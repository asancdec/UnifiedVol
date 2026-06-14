find_package(GTest CONFIG REQUIRED)

include(GoogleTest)

function(add_unifiedvol_test_executable target source_dir)
  file(GLOB_RECURSE TEST_SOURCES
      CONFIGURE_DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/${source_dir}/*.cpp
  )

  if(TEST_SOURCES)
    add_executable(${target} ${TEST_SOURCES})

    target_link_libraries(${target} PRIVATE
      UnifiedVol
      GTest::gtest_main
    )

    gtest_discover_tests(${target}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
  endif()
endfunction()

add_unifiedvol_test_executable(unifiedvol_unit_tests Unit)
add_unifiedvol_test_executable(unifiedvol_integration_tests Integration)
add_unifiedvol_test_executable(unifiedvol_regression_tests Regression)

add_custom_target(unifiedvol_tests
  DEPENDS
    unifiedvol_unit_tests
    unifiedvol_integration_tests
    unifiedvol_regression_tests
)
