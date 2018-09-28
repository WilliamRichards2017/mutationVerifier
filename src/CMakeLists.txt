INCLUDE_DIRECTORIES(
	${ZLIB_INCLUDE}
	${BAMTOOLS_INCLUDE}
)


add_library(mutation_verifier_core STATIC
	util.cpp	
)

SET(CORE_LIB mutation_verifier_core)

target_link_libraries(${CORE_LIB}
	${ZLIB_LIBRARY}
	${BAMTOOLS_LIB}
	${BAMTOOLS_UTIL_LIB}

)

