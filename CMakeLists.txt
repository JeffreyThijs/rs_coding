project(rs_coding)

include_directories(src/)

file(GLOB SOURCES "src/*.cpp")
 
add_executable(testStudent ${SOURCES})