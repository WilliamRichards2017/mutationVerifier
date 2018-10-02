CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(mutationVerifier)

SET(VERSION_MAJOR "0.8")
SET(VERSION_MINOR "0")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

UNSET(CORE_LIBS CACHE)
SET(CORE_LIBS CAHCE LIST)

include_directories("${PROJECT_SOURCE_DIR}/externals")
add_subdirectory(externals)

include_directories("${PROJECT_SOURCE_DIR}/src")
add_subdirectory(src)

include_directories("${PROJECT_SOURCE_DIR}/tools")
add_subdirectory(tools)

ADD_SUBDIRECTORY(gtest)
