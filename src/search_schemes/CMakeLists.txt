add_subdirectory(generator)

add_library (oss expand.cpp isValid.cpp)
target_link_libraries (oss PUBLIC generator)
