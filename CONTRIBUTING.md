# Contributing to FIRE

PRs are welcome! Please run the validation checks, the test suite, and the
formatters before submitting:

```bash
pixi run test-dry
pixi run test
pixi run test-verify
pixi run fmt
```

If your change touches manifest or reference handling, also run the
multi-sample test:

```bash
pixi run test-multi
```

## Conventional commits and releases

Releases are automated with
[release-please](https://github.com/googleapis/release-please): conventional
commits on `main` accumulate into a release PR that bumps the version in
`pixi.toml` and updates `CHANGELOG.md`; merging that PR tags the release.

PRs are squash-merged, so **the PR title must follow
[conventional commits](https://www.conventionalcommits.org)** (CI enforces
this):

- `feat: ...` — new feature; bumps the minor version (`0.2.x` → `0.3.0`).
  Note that FIRE output directories are versioned by major.minor (e.g.
  `-v0.2-`), so a minor bump changes output paths.
- `fix: ...` — bug fix; bumps the patch version (`0.2.0` → `0.2.1`) and does
  not change output paths.
- `feat!: ...` or `fix!: ...` — breaking change; also bumps the minor version
  while we are pre-1.0.
- `docs: ...`, `ci: ...`, `chore: ...`, `test: ...`, `refactor: ...` — no
  release.
