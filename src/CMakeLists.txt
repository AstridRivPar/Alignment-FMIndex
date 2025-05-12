configure_file(${CMAKE_CURRENT_SOURCE_DIR}/FMIndex.properties FMIndex.properties COPYONLY)


add_executable(FMIndex
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

target_link_libraries(FMIndex PRIVATE submodules
                                    libFMIndex)
target_include_directories(FMIndex PUBLIC 
    ${PROJECT_BINARY_DIR})   

add_custom_target(run_FMIndex
    COMMAND FMIndex -f -p -i ${PROJECT_SOURCE_DIR}/Resources/testI.txt -q ${PROJECT_SOURCE_DIR}/Resources/testQ.txt -k 50 -o ${PROJECT_SOURCE_DIR}/OutputFiles/testO.txt

    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Running Main")

