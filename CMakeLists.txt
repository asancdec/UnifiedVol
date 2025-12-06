cmake_minimum_required(VERSION 3.22)
project(UnifiedVol VERSION 0.1.0 LANGUAGES CXX)

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- MSVC runtime: /MD (Release), /MDd (Debug) ---
if(MSVC)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

# --- Sources ---
file(GLOB_RECURSE UNIFIEDVOL_SOURCES CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/src/*.cpp
)

# --- Library ---
add_library(UnifiedVol STATIC ${UNIFIEDVOL_SOURCES})
set_target_properties(UnifiedVol PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(UnifiedVol PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

# --- Dependencies (via vcpkg, CONFIG mode) ---
find_package(NLopt CONFIG REQUIRED)       # provides NLopt::nlopt
find_package(Ceres CONFIG REQUIRED)       # provides Ceres::ceres 
find_package(Boost CONFIG REQUIRED COMPONENTS math)  # Boost::math, Boost::boost

target_link_libraries(UnifiedVol PUBLIC
  NLopt::nlopt
  Ceres::ceres
  Boost::boost
  Boost::math
)

# --- Example (optional) ---
option(UNIFIEDVOL_BUILD_EXAMPLE "Build example program (examples/main.cpp)" ON)
if(UNIFIEDVOL_BUILD_EXAMPLE)
  add_executable(unifiedvol_example ${CMAKE_SOURCE_DIR}/examples/main.cpp)
  target_link_libraries(unifiedvol_example PRIVATE UnifiedVol)
  target_include_directories(unifiedvol_example PRIVATE ${CMAKE_SOURCE_DIR}/include)

  # --- Copy the data folder next to the example binary ---
  add_custom_command(TARGET unifiedvol_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/data
            $<TARGET_FILE_DIR:unifiedvol_example>/data
  )
endif()

# --- Python bindings (pybind11) ---
option(UNIFIEDVOL_BUILD_PYTHON "Build Python bindings (pybind11)" ON)

if(UNIFIEDVOL_BUILD_PYTHON)
  find_package(pybind11 CONFIG REQUIRED)

  pybind11_add_module(interp
    ${CMAKE_SOURCE_DIR}/python/Bindings.cpp
  )

  target_link_libraries(interp PRIVATE UnifiedVol)

  target_include_directories(interp PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
  )

  # Put the compiled Python module next to your Python sources
  set_target_properties(interp PROPERTIES
    PREFIX ""  # avoid libinterp.so, just interp*.so
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/python"
  )
endif()