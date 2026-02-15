# Contributing

## Cutting a release

Preview what will happen:

```bash
cargo release patch
cargo dist plan
```

Do the release (merge to main first):

```bash
git checkout main
git merge <branch>
cargo release patch -x
```

`cargo release` bumps the version in `Cargo.toml`, commits, tags, and pushes.
The GitHub Actions release workflow then builds binaries and creates the GitHub Release automatically.
