find_package(GTest CONFIG REQUIRED)

file(GLOB_RECURSE TEST_SOURCES
    CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_executable(unifiedvol_tests ${TEST_SOURCES})

target_link_libraries(unifiedvol_tests PRIVATE
  UnifiedVol
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(unifiedvol_tests)