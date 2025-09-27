cmake_minimum_required(VERSION 3.22)
project(UnifiedVol VERSION 0.1.0 LANGUAGES CXX)

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Sources & library ---
file(GLOB_RECURSE UNIFIEDVOL_SOURCES CONFIGURE_DEPENDS src/*.cpp)
add_library(UnifiedVol STATIC ${UNIFIEDVOL_SOURCES})
target_include_directories(UnifiedVol PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

# --- NLopt dependency (prefer config; fallback to raw path) ---
# If configured with vcpkg toolchain, this usually succeeds:
find_package(NLopt CONFIG QUIET)
if(TARGET NLOPT::nlopt)
  target_link_libraries(UnifiedVol PUBLIC NLOPT::nlopt)
else()
  # Fallback: locate headers/lib by name (works with vcpkg classic layout too)
  find_path(NLOPT_INCLUDE_DIR nlopt.hpp REQUIRED)
  find_library(NLOPT_LIBRARY NAMES nlopt REQUIRED)
  target_include_directories(UnifiedVol PUBLIC "${NLOPT_INCLUDE_DIR}")
  target_link_libraries(UnifiedVol PUBLIC "${NLOPT_LIBRARY}")
endif()

# --- Example program (optional) ---
option(UNIFIEDVOL_BUILD_EXAMPLE "Build example program (examples/main.cpp)" ON)
if(UNIFIEDVOL_BUILD_EXAMPLE)
  add_executable(unifiedvol_example examples/main.cpp)
  target_link_libraries(unifiedvol_example PRIVATE UnifiedVol)
  target_include_directories(unifiedvol_example PRIVATE ${CMAKE_SOURCE_DIR}/include)

  # Copy data/ next to the exe so relative paths work
  add_custom_command(TARGET unifiedvol_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_SOURCE_DIR}/data"
            "$<TARGET_FILE_DIR:unifiedvol_example>/data"
  )
endif()