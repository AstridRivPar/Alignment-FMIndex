cmake_minimum_required(VERSION 3.20)
project(FMIndex)

message(STATUS "[x] CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")



set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set (CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#find_package(Poco REQUIRED COMPONENTS Foundation Util)


add_subdirectory(src)
add_subdirectory(submodules)
add_subdirectory(libFMIndex)
