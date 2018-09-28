# Setting up external library bamtools, we don't build it because we only need the include directories

SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(BAMTOOLS_PROJECT bamtools_project CACHE INTERNAL "bamtools project name")
SET(BAMTOOLS_DIR ${CMAKE_BINARY_DIR}/externals/bamtools CACHE INTERNAL "bamtools project directory")
SET(BAMTOOLS_LIB)
ExternalProject_Add(${BAMTOOLS_PROJECT}
	GIT_REPOSITORY https://github.com/dillonl/bamtools.git
	GIT_TAG 45129fa9035847f54baa2c1750d4d3f526793e82 #lock in the commit id so we don't this doesn't break in the future
	DEPENDS ${ZLIB_PROJECT}
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	PREFIX ${BAMTOOLS_DIR}
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
		-DZLIB_LIBRARY_PATH:PATH=${ZLIB_LIBRARY_PATH}
				-DZLIB_INCLUDE:PATH=${ZLIB_INCLUDE}
)

ExternalProject_Get_Property(${BAMTOOLS_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${BAMTOOLS_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${BAMTOOLS_PROJECT} BINARY_DIR)

SET(BAMTOOLS_LIB ${SOURCE_DIR}/lib/libbamtools.a CACHE INTERNAL "Bamtools Lib")
SET(BAMTOOLS_UTIL_LIB ${SOURCE_DIR}/lib/libbamtools-utils.a CACHE INTERNAL "Bamtools-util Lib")
SET(BAMTOOLS_INCLUDE ${SOURCE_DIR}/include CACHE INTERNAL "Bamtools Include")
