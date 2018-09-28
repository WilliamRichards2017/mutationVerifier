CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(mutationVerifier)

SET(VERSION_MAJOR "0.8")
SET(VERSION_MINOR "0")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories("${PROJECT_SOURCE_DIR}/src")
add_subdirectory(src)
include_directories("${PROJECT_SOURCE_DIR}/externals")
add_subdirectory(externals)

#UNSET(CORE_LIBS CACHE)
#SET(CORE_LIBS CAHCE LIST)

INCLUDE_DIRECTORIES(
	${ZLIB_INCLUDE}
	${BAMTOOLS_INCLUDE}
)


add_executable(verify verify.cpp)
target_link_libraries(verify 
	mutation_verifier_core
)
		 
ADD_SUBDIRECTORY(gtest)
