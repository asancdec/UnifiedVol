# Contributing

Thank you for considering a contribution to UnifiedVol.

This project is a C++20 quantitative finance library. Contributions should
favor correctness, clear numerical behavior, focused tests, and small,
reviewable changes.

The project aims to follow front-office quantitative engineering standards:
numerical behavior should be explainable, testable, reproducible from public
sources, and performant enough for pricing and calibration workflows.

Performance matters for this project. For changes that may affect runtime,
allocation behavior, calibration speed, or pricing throughput, run the Linux
performance test workflow described in [BUILD.md](BUILD.md) when practical.

Feedback, review comments, and recommendations on existing code are very
welcome. If you see a simpler design, a more standard numerical method, a safer
API, or a better test strategy, please open an issue or mention it in the pull
request. This project is under active development and benefits from careful
outside review.

## Development Workflow

- Keep changes focused on one topic or issue.
- Prefer existing project patterns over introducing new abstractions.
- Add or update tests for behavioral changes.
- Keep public APIs small and explicit.
- Avoid unrelated formatting or cleanup in functional pull requests.

For all build, test, and performance-test commands, see [BUILD.md](BUILD.md).

Linux is the encouraged development and performance-testing environment for
this project. Builds may be extended to support other operating systems, but
platform-specific changes should preserve the Linux workflows documented in
[BUILD.md](BUILD.md).

## Testing Expectations

Before opening a pull request, run the relevant test suite described in
[BUILD.md](BUILD.md).

For most changes, unit tests are expected. Broader changes may also need
integration, regression, or performance tests.

The project aims to maintain at least 90% code coverage. Coverage is a guide,
not a substitute for meaningful numerical tests, but new code should generally
come with enough focused tests to preserve that target.

Good test targets include:

- validation failures and edge cases
- numerical invariants
- known values or regression fixtures
- API behavior for invalid input
- performance-sensitive paths when relevant

## Sources, Data, and Citations

All external data, formulas, algorithms, model descriptions, benchmark values,
and reference results contributed to this repository must come from public
sources that can be cited.

- Do not contribute proprietary, confidential, licensed, or client-owned market
  data.
- Do not contribute data scraped or copied from sources whose terms do not
  allow reuse.
- Cite papers, books, documentation, public datasets, and public web sources
  used to derive implementations, fixtures, tests, or examples.
- Add bibliographic references to `docs/citations.bib` when appropriate.
- Document dataset provenance in `docs/DATA.md` when adding or changing data
  files.
- Keep derived fixtures reproducible enough that reviewers can understand where
  the numbers came from.

## Code Style

Use the repository `.clang-format` file for formatting. Editors that use
clangd should pick it up automatically. The project uses an Allman brace style,
4-space indentation, a 90-column limit, and split multiline arguments rather
than tightly packed parameter lists.

General conventions:

- Use C++20.
- Keep headers and `.inl` files consistent with nearby code.
- Keep namespaces explicit and local to the component, for example
  `uv::core`, `uv::math`, `uv::models`, `uv::opt`, and `uv::io`.
- Prefer `std::span` for non-owning inputs.
- Prefer owning storage inside reusable model, surface, curve, optimizer, and
  interpolation objects.
- Use `Vector<T>` and `Complex<T>` from `Base/Types.hpp` where the surrounding
  code already uses the project aliases.
- Use project validation macros for user-facing precondition checks.
- Use `noexcept` only when a function is not expected to allocate, validate by
  throwing, call throwing APIs, or invoke arbitrary user code.
- Strongly prefer self-documenting code over comments. Comments are discouraged
  unless they explain a non-obvious numerical formula, source reference,
  invariant, performance trade-off, or external constraint that the code cannot
  make clear by naming and structure alone.

Architecture conventions:

- Public template declarations usually live in `.hpp` files and definitions in
  matching `Detail/*.inl` files included at the bottom of the header.
- Non-template implementation belongs in `src/`.
- Core data structures belong under `uv/Core`.
- Mathematical utilities belong under `uv/Math`.
- Model-specific code belongs under `uv/Models/<ModelName>`.
- Optimizer wrappers and calibration helpers belong under `uv/Optimization`.
- Tests mirror the library structure under `tests/Unit`, `tests/Integration`,
  `tests/Regression`, and `tests/Performance`.

Numerical code conventions:

- Validate dimensions, monotonicity, finiteness, and positivity/non-negativity
  at API boundaries.
- Keep tight pricing/calibration loops simple and allocation-free where
  practical.
- Prefer focused helper functions only when they reduce duplication without
  hiding important numerical formulas.
- Preserve existing fast paths when adding broader functionality.

## Pull Requests

Pull requests should include:

- a short description of the motivation
- a summary of the implementation
- tests added or updated
- the build/test commands run, referencing [BUILD.md](BUILD.md)
- any public sources, papers, datasets, or formulas used
- any known limitations or follow-up work

If the change affects numerical behavior, also include:

- what numerical results changed
- why the change is expected
- whether tolerances, fixtures, or regression values were updated
- any impact on calibration stability, arbitrage constraints, pricing accuracy,
  or performance

If the change affects public APIs, also include:

- the old behavior
- the new behavior
- whether the change is backward-compatible
- any migration notes for downstream users

## Issues

When opening an issue, please include:

- the observed behavior
- the expected behavior
- minimal reproduction steps when possible
- relevant input data, fixtures, or model parameters
- public sources or citations for any external data or formulas
- build/test context if the issue is environment-specific

## License

By contributing, you agree that your contribution will be licensed under the
same license as this repository.
