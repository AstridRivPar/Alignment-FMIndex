add_library(libFMIndex)

target_sources(libFMIndex
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Bi_FM_Index.cpp
    # ${CMAKE_CURRENT_SOURCE_DIR}/src/Helpers.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Reader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Matrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ED.cpp
)

target_include_directories(libFMIndex
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/include/
    ${BOOST_INCLUDE_DIRS}

) 

target_link_libraries(libFMIndex PUBLIC 
                        submodules
                        ${Boost_LIBRARIES}
                        Poco::Foundation
                        Poco::Util)
