set(CMAKE_CXX_STANDARD 14)

add_executable(test_range_stp test_range_stp.cpp)
target_link_libraries(test_range_stp ${PROJECT_NAME})

add_executable(test_range_time test_range_time.cpp)
target_link_libraries(test_range_time ${PROJECT_NAME})

add_executable(test_build_range_index test_build_range_index.cpp)
target_link_libraries(test_build_range_index ${PROJECT_NAME})

add_executable(test_index_search_time test_index_search_time.cpp)
target_link_libraries(test_index_search_time ${PROJECT_NAME})

add_executable(test_index_size test_index_size.cpp)
target_link_libraries(test_index_size ${PROJECT_NAME})

add_executable(test_index_search_range test_index_search_range.cpp)
target_link_libraries(test_index_search_range ${PROJECT_NAME})