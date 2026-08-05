#!/usr/bin/env bash
# Point this repository at your GitHub account and push it.
#   ./setup-github.sh <github-username> [repo-name]
set -euo pipefail

USER="${1:?usage: ./setup-github.sh <github-username> [repo-name]}"
REPO="${2:-sarcoma-epigenomic-atlas}"
URL="https://github.com/$USER/$REPO"

echo "Targeting $URL"

# fill the placeholders left in the docs
for f in README.md CITATION.cff LICENSE-DATA; do
  [ -f "$f" ] && sed -i.bak -e "s#<user>#$USER#g" -e "s#<you>#$USER#g" \
                            -e "s#gryderart.github.io/sarcoma-epigenomic-atlas#$USER.github.io/$REPO#g" \
                            "$f" && rm -f "$f.bak"
done
git add -A
git diff --cached --quiet || git commit -q -m "Point documentation at $USER/$REPO"

git remote remove origin 2>/dev/null || true
git remote add origin "$URL.git"

if command -v gh >/dev/null 2>&1 && gh auth status >/dev/null 2>&1; then
  echo "Creating the repository with the GitHub CLI..."
  gh repo create "$USER/$REPO" --public --source=. --remote=origin --push \
    --description "A sample-level census of every public epigenomic experiment on sarcoma, all ages, built to find the holes."
  echo "Enabling GitHub Pages from main /docs..."
  gh api -X POST "repos/$USER/$REPO/pages" -f "source[branch]=main" -f "source[path]=/docs" \
    >/dev/null 2>&1 || echo "  (enable Pages by hand: Settings -> Pages -> main, /docs)"
else
  echo
  echo "The GitHub CLI is not installed or not signed in."
  echo "Create an EMPTY repository at https://github.com/new named '$REPO'"
  echo "  - public, and do NOT add a README, .gitignore or licence"
  echo "Then run:"
  echo "    git push -u origin main"
fi

echo
echo "Done. Next:"
echo "  Settings -> Pages -> Source: main, folder /docs   (serves the gap map)"
echo "  $URL"
