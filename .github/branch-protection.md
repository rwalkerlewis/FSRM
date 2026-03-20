# Branch Protection Recommendations for `main`

## Required Status Checks

The following CI jobs should be **required to pass** before merging to `main`:

| Job | Required | Rationale |
|-----|----------|-----------|
| **Build** | ✅ Yes | Code must compile |
| **Unit Tests** | ✅ Yes | Core component correctness |
| **Functional Tests** | ✅ Yes | Feature-level correctness |
| **Physics/MMS Tests** | ✅ Yes | PDE solver correctness |
| **Integration Tests** | ✅ Yes | End-to-end correctness |
| **Performance Tests** | ❌ No | Informational; only runs on main push |
| **Static Analysis** | ❌ No | Informational initially; promote later |

## Pull Request Rules

- **Require pull request reviews**: At least 1 approving review before merge.
- **Dismiss stale reviews**: When new commits are pushed.
- **Require branches to be up-to-date**: Merge only when the branch is current with `main`.
- **Require conversation resolution**: All review comments must be resolved.

## Additional Settings

- **Restrict force pushes**: No one can force-push to `main`.
- **Restrict deletions**: `main` cannot be deleted.
- **Require signed commits**: Optional but recommended.

## How to Configure (GitHub UI)

1. Go to **Settings → Branches → Branch protection rules**
2. Click **Add rule** for branch name pattern: `main`
3. Enable:
   - ☑ Require a pull request before merging
   - ☑ Require status checks to pass before merging
     - Search and add: `Build`, `Unit Tests`, `Functional Tests`, `Physics/MMS Tests`, `Integration Tests`
   - ☑ Require branches to be up to date before merging
   - ☑ Do not allow bypassing the above settings

## Future Enhancements

- **Enable `-Werror`**: Once warnings are cleaned up, add compiler-warning-free as a requirement.
- **Require Static Analysis**: Once cppcheck findings are addressed, promote to required.
- **Code coverage**: Add coverage reporting and set minimum thresholds.
