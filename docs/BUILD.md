# Build Instructions

## Prerequisites

The documented presets target Linux with GCC, Ninja, CMake, and vcpkg.

- CMake >= 3.22
- Ninja
- GCC 13, available as `/usr/bin/g++-13`
- vcpkg checked out at `$HOME/dev/vcpkg`
- initialized git submodules
- ccache, optional; used automatically when installed

The vcpkg toolchain path is set in `CMakePresets.json`. If vcpkg is installed
elsewhere, update `CMAKE_TOOLCHAIN_FILE` in the presets before configuring.

When `ccache` is installed, CMake uses it automatically to speed up repeated
local builds. To disable this behavior for a build tree, configure with:

```bash
cmake --preset linux-gcc-release -DUNIFIEDVOL_ENABLE_CCACHE=OFF
```

## Quick Start

Configure, build, and run the standard release correctness suite:

```bash
cmake --preset linux-gcc-release
cmake --build --preset linux-gcc-release-tests
ctest --preset linux-gcc-release-nonperformance --output-on-failure
```

Build and run the example executable:

```bash
cmake --build --preset linux-gcc-release
./build/linux-gcc-release/unifiedvol_example
```

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

## Build Options

Common CMake options:

- `UNIFIEDVOL_BUILD_TESTS=ON/OFF`
- `UNIFIEDVOL_BUILD_EXAMPLE=ON/OFF`
- `UNIFIEDVOL_ENABLE_COVERAGE=ON/OFF`
- `UNIFIEDVOL_ENABLE_CCACHE=ON/OFF`

Example:

```bash
cmake --preset linux-gcc-release -DUNIFIEDVOL_BUILD_EXAMPLE=OFF
```

### Standard Build

Configure and build the standard release target:

```bash
cmake --preset linux-gcc-release
cmake --build --preset linux-gcc-release
```

Run the example executable:

```bash
./build/linux-gcc-release/unifiedvol_example
```

### Debug Build

```bash
cmake --preset linux-gcc-debug
cmake --build --preset linux-gcc-debug
```

## Tests

The test suite is split into **unit**, **integration**, **regression**, and
**performance** tests:

- **Unit tests:** validate one component in isolation.
- **Integration tests:** validate multiple components working together.
- **Regression tests:** compare results against committed reference values.
- **Performance tests:** check speed, allocation behavior, and scalability guardrails.

### Run Correctness Tests

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

Run the full release test suite, including performance-labelled tests in the
release build:

```bash
ctest --preset linux-gcc-release-full --output-on-failure
```

Run a specific correctness test category:

```bash
ctest --preset linux-gcc-release-unit --output-on-failure
ctest --preset linux-gcc-release-integration --output-on-failure
ctest --preset linux-gcc-release-regression --output-on-failure
```

### Run Performance Tests

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

Run performance tests on a quiet Linux machine when possible. The perf preset
uses `RelWithDebInfo`, `-O3`, debug symbols, and frame pointers for profiling.

### Run Coverage Tests

Coverage is GCC-only in this project. The coverage preset enables
`UNIFIEDVOL_ENABLE_COVERAGE`.

```bash
cmake --preset linux-gcc-coverage
cmake --build --preset linux-gcc-coverage-tests
ctest --preset linux-gcc-coverage-nonperformance --output-on-failure
```

Use local coverage-report tooling, such as `gcovr`, `gcov`, or `lcov`,
to generate reports from `build/linux-gcc-coverage`.

### Run All Tests

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

### Run Test Executables Directly

The test executables can also be run directly:

```bash
./build/linux-gcc-release/tests/unifiedvol_unit_tests
./build/linux-gcc-release/tests/unifiedvol_integration_tests
./build/linux-gcc-release/tests/unifiedvol_regression_tests
./build/linux-gcc-perf/tests/unifiedvol_performance_tests
```

## Static Analysis

Configure the release build and generate the compilation database used by
`clang-tidy`:

```bash
cmake --preset linux-gcc-release \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

cmake --build --preset linux-gcc-release -j "$(nproc)"
```

Run Clang-Tidy 18 with the repository `.clang-tidy` configuration:

```bash
run-clang-tidy-18 \
  -p build/linux-gcc-release \
  -config-file .clang-tidy \
  -header-filter "^$(pwd)/(uv|src|tests|examples)/.*" \
  -warnings-as-errors '*' \
  -j "$(nproc)" \
  "^$(pwd)/(uv|src|tests|examples)/.*\.(c|cc|cpp|cxx)$"
```

Only first-party translation units are analysed directly, excluding vendored
code under `external/`. First-party headers and `.inl` files are still checked
when included by those translation units. All enabled findings are treated as
errors.