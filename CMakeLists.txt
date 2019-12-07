cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 17)
project(rs_coding)

find_package( Boost REQUIRED COMPONENTS program_options)

include_directories(${Boost_INCLUDE_DIRS}
                    src/
                    galois/)

file(GLOB SOURCES 
     "src/*.cpp"
     # "src/rs_*.cpp"
     # "tests/test_decoder.cpp"
     "galois/GaloisField.cpp"
     "galois/GaloisFieldElement.cpp"
     "galois/GaloisFieldPolynomial.cpp")
 
add_executable(rs ${SOURCES})
target_link_libraries( rs ${Boost_PROGRAM_OPTIONS_LIBRARY})