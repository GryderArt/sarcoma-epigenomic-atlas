# Publish this repository to GitHub. Windows / PowerShell.
#
#   .\setup-github.ps1 -User <your-github-username>
#
# Installs the GitHub CLI if it is missing, signs you in (a browser window, once),
# creates the repository, pushes, and turns on GitHub Pages for the gap map.
#
# Your credential is stored by the GitHub CLI in Windows Credential Manager. It is
# never typed into a file, a script, or a chat window.

[CmdletBinding()]
param(
  [Parameter(Mandatory = $true)][string]$User,
  [string]$Repo    = "sarcoma-epigenomic-atlas",
  [ValidateSet("public", "private")][string]$Visibility = "public"
)

$ErrorActionPreference = "Stop"
function Step($n, $m) { Write-Host "`n[$n] $m" -ForegroundColor Cyan }
function Ok($m)       { Write-Host "    $m" -ForegroundColor Green }
function Warn($m)     { Write-Host "    $m" -ForegroundColor Yellow }

# ---------------------------------------------------------------- native commands
# Windows PowerShell 5.1 turns ANY stderr output from a native program into a terminating
# NativeCommandError when $ErrorActionPreference is "Stop". gh and git both report
# ordinary, non-error status on stderr: `gh auth status` says "You are not logged into any
# GitHub hosts" there, and `git push` writes its entire progress there. Under "Stop" that
# aborts the script at precisely the moments it exists to handle -- the first run of this
# script died on the sign-in check, before it could offer to sign you in.
#
# So every native call goes through one of these two helpers, which relax the preference
# for the duration of the call and report the exit code instead. PowerShell 7 does not
# behave this way, which is why this only ever bites on the shipped Windows console.

function Invoke-Native {
  # Run a native command with its output visible. Returns the command's OUTPUT; the exit
  # code lands in $NativeExit. Returning both down the pipeline would mean a caller doing
  # `| Select-Object -Last 1` picks up the exit code instead of the value it wanted --
  # which is how `gh api user --jq .login` first came back as "0".
  param([Parameter(Mandatory = $true)][scriptblock]$Cmd)
  $prev = $ErrorActionPreference
  $ErrorActionPreference = "Continue"
  try { & $Cmd; $script:NativeExit = $LASTEXITCODE } finally { $ErrorActionPreference = $prev }
}

function Test-Native {
  # Run a native command silently; return its exit code. For probes.
  param([Parameter(Mandatory = $true)][scriptblock]$Cmd)
  $prev = $ErrorActionPreference
  $ErrorActionPreference = "Continue"
  try { & $Cmd 2>&1 | Out-Null; return $LASTEXITCODE } finally { $ErrorActionPreference = $prev }
}

# ---------------------------------------------------------------- 0. Sanity
if (-not (Test-Path ".git")) {
  throw "Run this from inside the cloned repository (the folder containing README.md and pipeline\)."
}
if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
  throw "Git for Windows is not installed. Get it from https://git-scm.com/download/win and run this again."
}

# ---------------------------------------------------------------- 1. GitHub CLI
Step 1 "Checking for the GitHub CLI"
if (-not (Get-Command gh -ErrorAction SilentlyContinue)) {
  Warn "gh not found - installing via winget"
  Invoke-Native { winget install --id GitHub.cli --silent --accept-package-agreements --accept-source-agreements }
  $env:Path = [Environment]::GetEnvironmentVariable("Path", "Machine") + ";" +
              [Environment]::GetEnvironmentVariable("Path", "User")
  if (-not (Get-Command gh -ErrorAction SilentlyContinue)) {
    throw "gh installed but is not on PATH yet. Close this window, open a new PowerShell, and re-run."
  }
}
$ver = (Invoke-Native { gh --version } | Select-Object -First 1) -replace 'gh version ', ''
Ok "gh $ver"

# ---------------------------------------------------------------- 2. Sign in
Step 2 "Checking authentication"
if ((Test-Native { gh auth status }) -ne 0) {
  Warn "Not signed in. A browser window will open - approve it, then come back here."
  # Interactive: no redirection, so the device code and prompts reach you intact.
  Invoke-Native { gh auth login --hostname github.com --git-protocol https --web }
  if ($NativeExit -ne 0) { throw "Sign-in did not complete." }
}
$who = (Invoke-Native { gh api user --jq .login } | Select-Object -Last 1)
if (-not $who) { throw "Signed in, but could not read your account." }
$who = "$who".Trim()
Ok "signed in as $who"
if ($who -ne $User) {
  Warn "You passed -User $User but you are signed in as $who. Using $who."
  $User = $who
}
Invoke-Native { gh auth setup-git } | Out-Null

# git refuses to commit without an author, and a fresh Windows install has none. Derive
# one from the account just signed in, scoped to this repository so nothing global is
# touched. The noreply address is the one GitHub issues for exactly this purpose: commits
# attribute correctly on the web without publishing a real mailbox.
if (-not (Invoke-Native { git config --get user.email })) {
  $id = (Invoke-Native { gh api user --jq .id } | Select-Object -Last 1)
  git config user.email "$("$id".Trim())+$who@users.noreply.github.com"
  git config user.name  "$who"
  Warn "no git author was set - using $(git config --get user.email) for this repository"
}

# ---------------------------------------------------------------- 3. Placeholders
Step 3 "Filling documentation placeholders"
foreach ($f in @("README.md", "pipeline\README.tmpl.md", "CITATION.cff", "LICENSE-DATA")) {
  if (Test-Path $f) {
    $t = Get-Content $f -Raw
    $t = $t -replace '<user>', $User -replace '<you>', $User
    $t = $t -replace 'gryderart\.github\.io/sarcoma-epigenomic-atlas', "$($User.ToLower()).github.io/$Repo"
    Set-Content $f $t -NoNewline -Encoding UTF8
  }
}
Invoke-Native { git add -A } | Out-Null
if ((Test-Native { git diff --cached --quiet }) -ne 0) {
  # Piping a native call to Out-Null discards its output but not its outcome. The first
  # version of this did not read $NativeExit, so a commit that died on a missing git
  # author still printed "committed" and the run carried on with nothing committed.
  Invoke-Native { git commit -q -m "Point documentation at $User/$Repo" }
  if ($NativeExit -ne 0) { throw "Could not commit the placeholder changes - see above." }
  Ok "committed"
} else { Ok "nothing to change" }

# ---------------------------------------------------------------- 4. Create + push
Step 4 "Creating $User/$Repo ($Visibility) and pushing"
if ((Test-Native { gh repo view "$User/$Repo" }) -eq 0) {
  Warn "$User/$Repo already exists - pushing to it"
} else {
  Invoke-Native {
    gh repo create "$User/$Repo" "--$Visibility" `
      --description "A sample-level census of every public epigenomic experiment on sarcoma, all ages, built to find the holes."
  }
  if ($NativeExit -ne 0) { throw "Could not create $User/$Repo - see above." }
}

# Create and push are deliberately separate. `gh repo create --source=. --remote=origin`
# fails with "Unable to add remote origin" whenever an origin already exists -- and this
# repository is normally reached by cloning a bundle, which leaves origin pointing at the
# bundle file. Setting the remote ourselves works whether or not one is already there,
# and takes the same path whether the repository is new or not.
Test-Native { git remote remove origin } | Out-Null
Invoke-Native { git remote add origin "https://github.com/$User/$Repo.git" } | Out-Null
Invoke-Native { git push -u origin main }
if ($NativeExit -ne 0) { throw "Push failed. See the message above." }
Ok "pushed"

# ---------------------------------------------------------------- 5. Pages
Step 5 "Enabling GitHub Pages from main /docs"
$body = '{"source":{"branch":"main","path":"/docs"}}'
$code = Test-Native { $body | gh api -X POST "repos/$User/$Repo/pages" --input - }
if ($code -ne 0) {
  $code = Test-Native { $body | gh api -X PUT "repos/$User/$Repo/pages" --input - }
}
if ($code -eq 0) { Ok "Pages enabled" }
else { Warn "Could not enable Pages automatically. Settings -> Pages -> main, /docs" }

Write-Host "`nDone." -ForegroundColor Green
Write-Host "  Repository  https://github.com/$User/$Repo"
Write-Host "  Gap map     https://$($User.ToLower()).github.io/$Repo/   (a minute or two to build)"
Write-Host "`nWorth doing next: connect https://zenodo.org/account/settings/github/ and cut a"
Write-Host "v1.0.0 release, which mints a DOI you can cite in the white paper."
