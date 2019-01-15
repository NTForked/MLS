project(MLS)

cmake_minimum_required(VERSION 2.6.2)

find_package(PCL 1.3 REQUIRED common io visualization surface features search filters)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable(mls mls2.cpp)
target_link_libraries(mls ${PCL_LIBRARIES})
add_executable(clustering clustering.cpp)
target_link_libraries(clustering ${PCL_LIBRARIES})
target_link_libraries(clustering Eigen3::Eigen)
add_executable(rbf rbf.cpp)
target_link_libraries(rbf ${PCL_LIBRARIES})
target_link_libraries(rbf Eigen3::Eigen)