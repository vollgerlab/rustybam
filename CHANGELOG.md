# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0](https://github.com/vollgerlab/rustybam/compare/v0.2.2...v0.3.0) - 2026-08-16

### Added

- add seq-content command and seq-stats improvements ([#28](https://github.com/vollgerlab/rustybam/pull/28))

### Other

- use the default changelog template
## [0.2.2](https://github.com/vollgerlab/rustybam/compare/v0.2.1...v0.2.2) - 2026-08-16

### Fixed

- move to rust-htslib 1.0 and bio 4.0, drop bio-io ([03f6e7e](https://github.com/vollgerlab/rustybam/commit/03f6e7eedb2f8307bf3d50744f780672a67064f9))
- bump rust-htslib to 0.45 for current hts-sys bindings ([807c402](https://github.com/vollgerlab/rustybam/commit/807c4020e2e77f08068a0e3ce1091ea8485b00b7))
- PAF parsing panics and inconsistent empty truncation ([49d8cd9](https://github.com/vollgerlab/rustybam/commit/49d8cd96552211c7556fd1aeec3253ff797cc972))

### Other

- link changelog entries to pull requests ([5536206](https://github.com/vollgerlab/rustybam/commit/5536206037c34671f3443bca5dbc1f5f8e8b1982))
- Merge pull request #26 from vollgerlab/fix/rust-htslib-045 ([4bb033c](https://github.com/vollgerlab/rustybam/commit/4bb033c9cbd95d666975d1b82ef50b44a25ff915))

## [0.2.1](https://github.com/vollgerlab/rustybam/compare/v0.2.0...v0.2.1) - 2026-08-15

### Fixed

- correct trailing-indel trimming, cs-tag inversion, and cigar parsing panics

### Other

- Merge pull request #22 from vollgerlab/fix/ultracode-2
- use the draft-release and dispatch pattern for cargo-dist
- add release-plz for automated releases
