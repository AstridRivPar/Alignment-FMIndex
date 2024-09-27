cmake_minimum_required(VERSION 3.10)
project(bidirectionalsearch)




add_executable(Main main.cpp)
# add_executable(test runTests.cpp)

set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories( $ENV{HOME}/include ${CMAKE_CURRENT_SOURCE_DIR}/ReadingInput)

set(myincludes $ENV{HOME}/lib/libdivsufsort.a $ENV{HOME}/lib/libdivsufsort64.a $ENV{HOME}/lib/libsdsl.a )
target_link_libraries(Main PUBLIC ${myincludes})
# target_link_libraries(test PUBLIC ${myincludes})

add_custom_target(run_Main
    COMMAND Main
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "RUNNING approximate matching"
)
# add_custom_target(run_test
    
#     COMMAND test -O3 -std=c++11 ${CMAKE_CURRENT_SOURCE_DIR}/Resources/Synthetic/Models/E.txt ${CMAKE_CURRENT_SOURCE_DIR}/Resources/Unique/uniqueE.txt ${CMAKE_CURRENT_SOURCE_DIR}/OutputFiles/uniqueEStats_Numbers.txt 14
#     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#     COMMENT "RUNNING tests"
# )
# COMMAND server_app -i ../../tools/resources/${input_file}_to_ints_bwt_wm.txt -p 12345 -d 10
