# Build Instructions

## Cloning the Repository

This project uses **git submodules**.

Clone the repository with submodules enabled:

```bash
git clone --recurse-submodules https://github.com/asancdec/UnifiedVol.git
```

If you already cloned the repository without submodules, initialize them with:

```bash
git submodule update --init --recursive
```

## Standard Build (Recommended)

```bash
rm -rf build/linux-gcc-release

cmake --preset linux-gcc-release
cmake --build --preset linux-gcc-release

./build/linux-gcc-release/unifiedvol_example
```

## Run All Tests

The test suite is split into **unit**, **integration**, **regression**, and
**performance** tests:

- **Unit:** correctness of one component.
- **Integration:** components working together.
- **Regression:** results match committed reference values.
- **Performance:** speed, allocation, and scalability guardrails.

Build the tests with:

```bash
cmake --build --preset linux-gcc-release-tests
```

Run the normal non-performance suite with CTest:

```bash
ctest --preset linux-gcc-release-nonperformance
```

Run all discovered tests, including performance guardrails:

```bash
ctest --preset linux-gcc-release-full
```

Run a specific test category with CTest presets:

```bash
ctest --preset linux-gcc-release-unit
ctest --preset linux-gcc-release-integration
ctest --preset linux-gcc-release-regression
cmake --build --preset linux-gcc-perf
ctest --preset linux-gcc-perf-performance
```

You can also run the test executables directly:

```bash
./build/linux-gcc-release/tests/unifiedvol_unit_tests
./build/linux-gcc-release/tests/unifiedvol_integration_tests
./build/linux-gcc-release/tests/unifiedvol_regression_tests
./build/linux-gcc-release/tests/unifiedvol_performance_tests
```

## Suggested CI Split

- Debug or release: `linux-gcc-release-nonperformance` for ordinary correctness checks.
- Perf build: `linux-gcc-perf-performance` when performance guardrails are requested.
- Optional matrix: GCC and Clang release builds.
- Golden files under `tests/Golden/` are manual reference data; normal test runs
  must not rewrite them.
