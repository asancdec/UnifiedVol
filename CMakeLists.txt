cmake_minimum_required(VERSION 3.22)
project(UnifiedVol VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

# --- UnifiedVol ---
file(GLOB_RECURSE UNIFIEDVOL_SOURCES CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/src/*.cpp
)

add_library(UnifiedVol STATIC ${UNIFIEDVOL_SOURCES})
set_target_properties(UnifiedVol PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(UnifiedVol PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/uv>
)

# --- Let's Be Rational (submodule) ---
set(LBR_DIR ${CMAKE_SOURCE_DIR}/external/lets_be_rational)

if(NOT EXISTS ${LBR_DIR})
  message(FATAL_ERROR
    "Missing submodule: external/lets_be_rational\n"
    "Run: git submodule update --init --recursive\n"
    "Or clone with: git clone --recurse-submodules https://github.com/asancdec/UnifiedVol.git"
  )
endif()

file(GLOB_RECURSE LBR_SOURCES CONFIGURE_DEPENDS
  ${LBR_DIR}/*.cpp
)

if(NOT LBR_SOURCES)
  message(FATAL_ERROR "No .cpp sources found under: ${LBR_DIR}")
endif()

add_library(lets_be_rational STATIC ${LBR_SOURCES})
set_target_properties(lets_be_rational PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(lets_be_rational PUBLIC
  ${LBR_DIR}
  ${LBR_DIR}/lets_be_rational
  ${LBR_DIR}/src
)

# --- Dependencies ---
find_package(NLopt CONFIG REQUIRED)
find_package(Ceres CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS math)

target_link_libraries(UnifiedVol PUBLIC
  lets_be_rational
  NLopt::nlopt
  Ceres::ceres
  Boost::boost
  Boost::math
)

# --- Example ---
option(UNIFIEDVOL_BUILD_EXAMPLE "Build example program (examples/main.cpp)" ON)
if(UNIFIEDVOL_BUILD_EXAMPLE)
  add_executable(unifiedvol_example ${CMAKE_SOURCE_DIR}/examples/main.cpp)
  target_link_libraries(unifiedvol_example PRIVATE UnifiedVol)

  add_custom_command(TARGET unifiedvol_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/data
            $<TARGET_FILE_DIR:unifiedvol_example>/data
  )
endif()
