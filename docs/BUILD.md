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

## Run All Unit Tests

```bash
ctest --test-dir build/linux-gcc-release --output-on-failure
```