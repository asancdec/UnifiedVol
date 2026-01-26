# Build Instructions

## Standard Build (Recommended)

```bash
rm -rf build

cmake -S . -B build -G "Ninja Multi-Config" \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_MANIFEST_MODE=ON \
  -DUNIFIEDVOL_BUILD_EXAMPLE=ON

cmake --build build --config Release

./build/Release/unifiedvol_example
```

---

## Profile-Guided Optimization (Optimized)

### Generate profile data (instrumented build)

```bash
rm -rf build-pgo-gen build-pgo-use pgo-data
mkdir -p pgo-data

cmake -S . -B build-pgo-gen -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_MANIFEST_MODE=ON \
  -DVCPKG_TARGET_TRIPLET=x64-linux \
  -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \
  -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=native -fprofile-generate=$(pwd)/pgo-data -DNDEBUG"

cmake --build build-pgo-gen -j
```

### Run representative workloads (repeat if possible)

```bash
./build-pgo-gen/unifiedvol_example
./build-pgo-gen/unifiedvol_example
```

### Use profile data (optimized build)

```bash
cmake -S . -B build-pgo-use -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_MANIFEST_MODE=ON \
  -DVCPKG_TARGET_TRIPLET=x64-linux \
  -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \
  -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=native -fprofile-use=$(pwd)/pgo-data -fprofile-correction -DNDEBUG"

cmake --build build-pgo-use -j

./build-pgo-use/unifiedvol
```