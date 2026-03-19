# Branch Protection Recommendations

Recommended GitHub branch protection settings for the `main` branch of FSRM.

## Required Status Checks

The following CI jobs (from `.github/workflows/ci.yml`) should be configured as
**required status checks** that must pass before a pull request can be merged:

| Status Check Name | Required | Notes |
|-------------------|----------|-------|
| `Build` | ✅ Yes | Compilation must never silently break |
| `Unit Tests` | ✅ Yes | Fast-feedback tests for core components |
| `Functional Tests` | ✅ Yes | Feature-level correctness tests |
| `Physics Tests` | ✅ Yes | MMS verification tests |
| `Integration Tests` | ✅ Yes | End-to-end simulation tests |
| `CI Summary` | ✅ Yes | Aggregates all required check results |
| `Performance Tests` | ❌ No | Only runs on push to main (informational) |
| `Disabled Test Report` | ❌ No | Informational tracking only |
| `Static Analysis / cppcheck` | ❌ No | Informational initially; promote to required after warning cleanup |

## Pull Request Settings

| Setting | Recommended Value | Notes |
|---------|-------------------|-------|
| Require pull request reviews | Yes, 1 reviewer | Ensure code review before merge |
| Dismiss stale reviews on new push | Yes | Re-review after changes |
| Require review from code owners | Optional | Set up CODEOWNERS file if desired |
| Require branches to be up to date | Yes | Ensures CI runs on merged result |
| Require conversation resolution | Optional | Useful but not critical |
| Require signed commits | Optional | Nice-to-have for supply chain security |
| Require linear history | No | Merge commits are fine for this project |

## Administrator Settings

| Setting | Recommended Value | Notes |
|---------|-------------------|-------|
| Include administrators | No | Allow admin override for hotfixes |
| Allow force pushes | No | Protect commit history |
| Allow deletions | No | Prevent accidental branch deletion |

## How to Configure

1. Go to **Settings → Branches → Add branch protection rule**
2. Branch name pattern: `main`
3. Check **"Require a pull request before merging"**
4. Check **"Require status checks to pass before merging"**
5. Search for and add each required status check listed above
6. Check **"Require branches to be up to date before merging"**
7. Save changes

## Promoting Static Analysis to Required

Once the existing cppcheck warnings have been addressed, promote the static
analysis check to required:

1. Go to **Settings → Branches → main protection rule → Edit**
2. Under **"Require status checks to pass"**, add `Static Analysis / cppcheck`
3. Save changes

## Promoting Warning-Free Builds

To enforce warning-free builds in the future:

1. Clean up existing compiler warnings (check the CI build job summary for counts)
2. Add `-DCMAKE_CXX_FLAGS="-Werror"` to the cmake configure step in `ci.yml`
3. Add `Build` as a required status check (already recommended above)
