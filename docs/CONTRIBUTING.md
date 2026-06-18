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

## Guidelines

- Keep changes focused on one topic or issue.
- Prefer existing project patterns over introducing new abstractions.
- Add or update tests for behavioral changes.
- Avoid unrelated formatting or cleanup in functional pull requests.
- Use public sources for external formulas, data, algorithms, and benchmark
  values, and cite them in `citations.bib` when applicable.
- Do not contribute proprietary, confidential, licensed, or client-owned market
  data.

## Style

Use the repository `.clang-format` file.

## License

By contributing, you agree that your contribution will be licensed under the
same license as this repository.
