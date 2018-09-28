set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(verify
        verify.cpp
)

target_link_libraries(verify
	${CORE_LIB}
)
