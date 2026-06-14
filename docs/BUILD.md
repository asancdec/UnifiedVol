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

## Standard Build

Configure and build the standard release target:

```bash
rm -rf build/linux-gcc-release

cmake --preset linux-gcc-release
cmake --build --preset linux-gcc-release
```

Run the example executable:

```bash
./build/linux-gcc-release/unifiedvol_example
```

## Tests

The test suite is split into **unit**, **integration**, **regression**, and
**performance** tests:

- **Unit tests:** validate one component in isolation.
- **Integration tests:** validate multiple components working together.
- **Regression tests:** compare results against committed reference values.
- **Performance tests:** check speed, allocation behavior, and scalability guardrails.

## Run Correctness Tests

Correctness tests use the standard release build.

Configure the release build:

```bash
cmake --preset linux-gcc-release
```

Build the correctness test targets:

```bash
cmake --build --preset linux-gcc-release-tests
```

Run all non-performance tests:

```bash
ctest --preset linux-gcc-release-nonperformance --output-on-failure
```

Run a specific correctness test category:

```bash
ctest --preset linux-gcc-release-unit --output-on-failure
ctest --preset linux-gcc-release-integration --output-on-failure
ctest --preset linux-gcc-release-regression --output-on-failure
```

## Run Performance Tests

Performance tests use a separate performance-oriented build. This avoids mixing
correctness testing with benchmark-style checks and keeps performance results
more stable.

Configure the performance build:

```bash
cmake --preset linux-gcc-perf
```

Build the performance test targets:

```bash
cmake --build --preset linux-gcc-perf-tests
```

Run the performance test suite:

```bash
ctest --preset linux-gcc-perf-performance --output-on-failure
```

## Run All Tests

To run the full test set, run the correctness suite first, then the performance
suite:

```bash
cmake --preset linux-gcc-release
cmake --build --preset linux-gcc-release-tests
ctest --preset linux-gcc-release-nonperformance --output-on-failure

cmake --preset linux-gcc-perf
cmake --build --preset linux-gcc-perf-tests
ctest --preset linux-gcc-perf-performance --output-on-failure
```

## Run Test Executables Directly

The test executables can also be run directly:

```bash
./build/linux-gcc-release/tests/unifiedvol_unit_tests
./build/linux-gcc-release/tests/unifiedvol_integration_tests
./build/linux-gcc-release/tests/unifiedvol_regression_tests
./build/linux-gcc-perf/tests/unifiedvol_performance_tests
```