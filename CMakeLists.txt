project(MLS)

cmake_minimum_required(VERSION 2.6.2)

find_package(CGAL QUIET COMPONENTS Core )
find_package(PCL 1.3 REQUIRED common io visualization surface features search filters)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories( /opt/local/include/eigen3 )
  
  create_single_source_cgal_program( "mls1.cpp")
  target_link_libraries(mls1 ${PCL_LIBRARIES})
  
  create_single_source_cgal_program( "mls2.cpp")
  target_link_libraries(mls2 ${PCL_LIBRARIES})

  create_single_source_cgal_program( "eigenexp.cpp")
  target_link_libraries(eigenexp ${PCL_LIBRARIES})
  

endif()
