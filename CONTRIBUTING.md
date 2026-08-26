# Contributing to spinDecon

Contributions are welcome. Please keep changes focused, reproducible, and free of generated build artifacts.

## Development workflow

1. Create or update the developer environment with `make`.
2. Make source changes in `applications/`, `native/`, `scripts/`, or other tracked source directories.
3. Run `make test` for the Python test suite.
4. Run `make` before submitting changes so the doctor checks pass.
5. For macOS packaging changes, also run `make app` and, when relevant, `make dmg`.

Do not commit `.venv`, `.build`, `dist`, generated native binaries/libraries, Python caches, DMGs, or application bundles. The repository `.gitignore` covers these outputs.

Please preserve third-party copyright and license notices when modifying or updating vendored dependencies.
