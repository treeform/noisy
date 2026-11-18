# Development workflow

- Use Nimby for dependency installation. Never use Nimble commands.
- Run the complete test suite with `nim r tests/tests.nim`.
- Keep benchmarks in `tests/bench_*.nim` and generators in `tools/gen_*.nim`.
- Commit completed work and push directly to `master`, as requested by the owner.
