project(Simplify)

cmake_minimum_required (VERSION 2.6)

find_package(CGAL QUIET COMPONENTS Core )
  include( ${CGAL_USE_FILE} )

  include_directories (BEFORE "../../include")

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIR})

add_executable(Simplify simplify.cpp)
target_link_libraries(Simplify ${GLUT_LIBRARY} ${OPENGL_LIBRARY})   
target_link_libraries(Simplify ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})                               