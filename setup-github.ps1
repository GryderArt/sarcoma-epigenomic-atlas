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

if (-not (Test-Path ".git")) {
  throw "Run this from inside the cloned repository (the folder containing README.md and pipeline\)."
}

# ---------------------------------------------------------------- 1. GitHub CLI
Step 1 "Checking for the GitHub CLI"
if (-not (Get-Command gh -ErrorAction SilentlyContinue)) {
  Warn "gh not found - installing via winget"
  winget install --id GitHub.cli --silent --accept-package-agreements --accept-source-agreements
  $env:Path = [Environment]::GetEnvironmentVariable("Path", "Machine") + ";" +
              [Environment]::GetEnvironmentVariable("Path", "User")
  if (-not (Get-Command gh -ErrorAction SilentlyContinue)) {
    throw "gh installed but is not on PATH yet. Close this window, open a new PowerShell, and re-run."
  }
}
Ok "gh $((gh --version | Select-Object -First 1) -replace 'gh version ','')"

# ---------------------------------------------------------------- 2. Sign in
Step 2 "Checking authentication"
gh auth status 2>&1 | Out-Null
if ($LASTEXITCODE -ne 0) {
  Warn "Not signed in. A browser window will open - approve it, then come back here."
  gh auth login --hostname github.com --git-protocol https --web
  if ($LASTEXITCODE -ne 0) { throw "Sign-in did not complete." }
}
$who = (gh api user --jq .login 2>$null)
if (-not $who) { throw "Signed in, but could not read your account." }
Ok "signed in as $who"
if ($who -ne $User) {
  Warn "You passed -User $User but you are signed in as $who. Using $who."
  $User = $who
}

# ---------------------------------------------------------------- 3. Placeholders
Step 3 "Filling documentation placeholders"
foreach ($f in @("README.md", "pipeline\README.tmpl.md", "CITATION.cff", "LICENSE-DATA")) {
  if (Test-Path $f) {
    $t = Get-Content $f -Raw
    $t = $t -replace '<user>', $User -replace '<you>', $User
    $t = $t -replace 'gryderart\.github\.io/sarcoma-epigenomic-atlas', "$User.github.io/$Repo"
    Set-Content $f $t -NoNewline -Encoding UTF8
  }
}
git add -A
git diff --cached --quiet
if ($LASTEXITCODE -ne 0) { git commit -q -m "Point documentation at $User/$Repo"; Ok "committed" }
else { Ok "nothing to change" }

# ---------------------------------------------------------------- 4. Create + push
Step 4 "Creating $User/$Repo ($Visibility) and pushing"
$exists = $false
gh repo view "$User/$Repo" 2>&1 | Out-Null
if ($LASTEXITCODE -eq 0) { $exists = $true }

if ($exists) {
  Warn "$User/$Repo already exists - pushing to it"
  git remote remove origin 2>$null
  git remote add origin "https://github.com/$User/$Repo.git"
  git push -u origin main
} else {
  gh repo create "$User/$Repo" "--$Visibility" --source=. --remote=origin --push `
    --description "A sample-level census of every public epigenomic experiment on sarcoma, all ages, built to find the holes."
}
if ($LASTEXITCODE -ne 0) { throw "Push failed. See the message above." }
Ok "pushed"

# ---------------------------------------------------------------- 5. Pages
Step 5 "Enabling GitHub Pages from main /docs"
$body = '{"source":{"branch":"main","path":"/docs"}}'
$body | gh api -X POST "repos/$User/$Repo/pages" --input - 2>&1 | Out-Null
if ($LASTEXITCODE -ne 0) {
  $body | gh api -X PUT "repos/$User/$Repo/pages" --input - 2>&1 | Out-Null
}
if ($LASTEXITCODE -eq 0) { Ok "Pages enabled" }
else { Warn "Could not enable Pages automatically. Settings -> Pages -> main, /docs" }

Write-Host "`nDone." -ForegroundColor Green
Write-Host "  Repository  https://github.com/$User/$Repo"
Write-Host "  Gap map     https://$User.github.io/$Repo/   (a minute or two to build)"
Write-Host "`nWorth doing next: connect https://zenodo.org/account/settings/github/ and cut a"
Write-Host "v1.0.0 release, which mints a DOI you can cite in the white paper."
