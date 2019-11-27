cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 17)
project(rs_coding)

include_directories(src/
                    galois/)

file(GLOB SOURCES 
     "src/*.cpp"
     "galois/GaloisField.cpp"
     "galois/GaloisFieldElement.cpp"
     "galois/GaloisFieldPolynomial.cpp")
 
add_executable(rs ${SOURCES})