cmake_minimum_required(VERSION 3.8)
project(rodsim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(rodsim_cli
	src/main.cpp
	src/assembly.cpp
	src/stochastic.cpp
	src/event.cpp
	deps/pugixml/pugixml.cpp)

target_include_directories(rodsim_cli PUBLIC src/)
target_include_directories(rodsim_cli PUBLIC deps/pugixml/)

find_package(Threads REQUIRED)
target_link_libraries(rodsim_cli PRIVATE Threads::Threads)
