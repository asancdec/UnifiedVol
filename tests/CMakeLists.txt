find_package(GTest CONFIG REQUIRED)

include(GoogleTest)

function(add_unifiedvol_test_executable target source_dir label)
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

    target_include_directories(${target} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}
    )

    gtest_discover_tests(${target}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      PROPERTIES LABELS ${label}
    )
  endif()
endfunction()

add_unifiedvol_test_executable(unifiedvol_unit_tests Unit unit)
add_unifiedvol_test_executable(unifiedvol_integration_tests Integration integration)
add_unifiedvol_test_executable(unifiedvol_regression_tests Regression regression)
add_unifiedvol_test_executable(unifiedvol_performance_tests Performance performance)

add_custom_target(unifiedvol_tests
  DEPENDS
    unifiedvol_unit_tests
    unifiedvol_integration_tests
    unifiedvol_regression_tests
    unifiedvol_performance_tests
)
