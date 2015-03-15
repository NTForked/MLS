project(MLS)

cmake_minimum_required(VERSION 2.6.2)

find_package(PCL 1.3 REQUIRED common io visualization surface features search filters)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

include_directories( /opt/local/include/eigen3 )

add_executable(mls mls2.cpp)
target_link_libraries(mls ${PCL_LIBRARIES})
add_executable(clustering clustering.cpp)
target_link_libraries(clustering ${PCL_LIBRARIES})
add_executable(rbf rbf.cpp)
target_link_libraries(rbf ${PCL_LIBRARIES})