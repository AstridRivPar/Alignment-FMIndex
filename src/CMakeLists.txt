configure_file(${CMAKE_CURRENT_SOURCE_DIR}/FMIndex.properties FMIndex.properties COPYONLY)

# find_package(glog REQUIRED)
# find_package(gflags REQUIRED)
add_executable(FMIndex
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
# add_definitions(-DC10_USE_GLOG)
# set(GFLAGS_INCLUDE_DIR /usr/lib/x86_64-linux-gnu)
target_link_libraries(FMIndex PRIVATE submodules
                                    libFMIndex)
target_include_directories(FMIndex PUBLIC 
    ${PROJECT_BINARY_DIR})   
# target_link_libraries(Main PUBLIC submodules glog::glog gflags libFMIndex)
add_custom_target(run_FMIndex
    COMMAND FMIndex -f -p -i ${PROJECT_SOURCE_DIR}/Resources/testI.txt -q ${PROJECT_SOURCE_DIR}/Resources/testQ.txt -k 50 -o ${PROJECT_SOURCE_DIR}/OutputFiles/testO.txt

    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Running Main")

# install(TARGETS FMIndex DESTINATION ${CMAKE_INSTALL_BINDIR})
# install(FILES ${CMAKE_CURRENT_BINARY_DIR}/FMIndex.properties DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT config)
