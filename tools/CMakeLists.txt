message("PROJECT SOURCE DIR IS ${PROJECT_SOURCE_DIR}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

INCLUDE_DIRECTORIES(
	${CXXOPTS_INCLUDE}
	${ZLIB_INCLUDE}
	${BAMTOOLS_INCLUDE}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(verify
        verify.cpp
)

target_link_libraries(verify
	${CORE_LIB}
)

add_dependencies(verify ${MUTATION_VERIFIER_EXTERNAL_PROJECT})
