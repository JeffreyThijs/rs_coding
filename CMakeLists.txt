cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 17)
project(rs_coding)

set (CMAKE_CXX_FLAGS "-Wall -pedantic -ansi")

find_package( Boost REQUIRED COMPONENTS program_options unit_test_framework)
include_directories(${Boost_INCLUDE_DIRS}
                    src/
                    galois/)

file(GLOB SOURCES 
     "src/rs_*.cpp"
     "src/helper_*.cpp"
     "galois/GaloisField.cpp"
     "galois/GaloisFieldElement.cpp"
     "galois/GaloisFieldPolynomial.cpp")
 
add_executable(rs "src/main.cpp" ${SOURCES})
target_link_libraries(rs ${Boost_PROGRAM_OPTIONS_LIBRARY})

# tests
file(GLOB TEST_SOURCES
     tests/*.cpp
)

add_executable (run_rs_tests ${TEST_SOURCES} ${SOURCES})
target_link_libraries (run_rs_tests ${Boost_LIBRARIES})