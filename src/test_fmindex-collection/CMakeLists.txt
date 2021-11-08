cmake_minimum_required (VERSION 3.8)

project(test_fmindex-collection)

# fmindex-collectionunittests
add_executable(${PROJECT_NAME}
    checkOccTables.cpp
    checkBiFMIndex.cpp
    main.cpp
)
target_link_libraries(${PROJECT_NAME}
    Catch2::Catch2
    fmindex-collection
)
