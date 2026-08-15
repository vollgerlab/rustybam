# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.1](https://github.com/vollgerlab/rustybam/compare/v0.2.0...v0.2.1) - 2026-08-15

### Fixed

- PAF parsing panics and inconsistent empty truncation
- correct trailing-indel trimming, cs-tag inversion, and cigar parsing panics

### Other

- Merge pull request #22 from vollgerlab/fix/ultracode-2
- use the draft-release and dispatch pattern for cargo-dist
- add release-plz for automated releases
