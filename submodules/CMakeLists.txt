add_subdirectory(sdsl-lite)
add_library(submodules INTERFACE)
target_link_libraries(submodules INTERFACE sdsl-lite::sdsl-lite)