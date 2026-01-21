cmake_minimum_required(VERSION 3.16)

project(Mininec3_object_test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Mininec3_object_test main.cpp)

include(GNUInstallDirs)
install(TARGETS Mininec3_object_test
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
