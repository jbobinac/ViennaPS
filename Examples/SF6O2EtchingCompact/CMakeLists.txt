cmake_minimum_required(VERSION 3.4)

project("SF6O2EtchingCompact")

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${VIENNAPS_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${VIENNAPS_LIBRARIES})

configure_file(parameters.hpp.in parameters.hpp)
set(dimension_extraction "DimensionExtraction")
add_executable(${dimension_extraction} ${dimension_extraction}.cpp)
target_include_directories(${dimension_extraction} PUBLIC ${VIENNAPS_INCLUDE_DIRS} ${CMAKE_BINARY_DIR})
target_link_libraries(${dimension_extraction} ${VIENNAPS_LIBRARIES})

add_dependencies(buildExamples ${PROJECT_NAME})

# configure_file(PyVtmHDimensions.csv ${CMAKE_CURRENT_BINARY_DIR}/PyVtmHDimensions.csv COPYONLY)
configure_file(dimensionsP&y35ptGrid.csv ${CMAKE_CURRENT_BINARY_DIR}/dimensionsP&y35ptGrid.csv COPYONLY)
#configure_file(dimensionsV&t.csv ${CMAKE_CURRENT_BINARY_DIR}/dimensionsV&t.csv COPYONLY)