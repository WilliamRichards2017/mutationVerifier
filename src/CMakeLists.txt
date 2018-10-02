set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

message("ZLIB INCLUDE: ${ZLIB_LIBRARY}")

ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)
INCLUDE_DIRECTORIES(
	${ZLIB_INCLUDE}
	${BAMTOOLS_INCLUDE}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")


add_library(mutation_verifier_core STATIC
				   util.cpp	
)

SET(CORE_LIB mutation_verifier_core CACHE STRING "the core library")

target_link_libraries(${CORE_LIB}
	${BAMTOOLS_LIB}
	${BAMTOOLS_UTIL_LIB}
	${ZLIB_LIBRARY}
)

add_dependencies(${CORE_LIB} ${MUTATION_VERIFIER_EXTERNAL_PROJECT})
