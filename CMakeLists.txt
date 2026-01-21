cmake_minimum_required(VERSION 3.16)

project(Mininec3_object_test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Mininec3_object_test
    GeometryInput.hpp
    Connections.hpp
    GeometryOutput.hpp
    Connections.cpp
    GeometryInput.cpp
    GeometryOutput.cpp
    GeometryData.hpp
    SimulationState.hpp
    PrintCurrents.cpp
    PrintCurrents.hpp
    ImpedanceMatrix.cpp
    ImpedanceMatrix.hpp
    Beta_247_315.cpp
    Beta_247_315.hpp
    Sub_273_312.cpp
    Sub_273_312.hpp
    Gosub102.cpp
    Gosub102.hpp
    Gosub113.cpp
    Gosub113.hpp
    Gosub135.cpp
    Gosub135.hpp
    FrequencyInput.cpp
    FrequencyInput.hpp
    Gosub87.cpp
    Gosub87.hpp
    Kerneleval.cpp
    KernelEval.hpp
    EllipticIntegral.cpp
    EllipticIntegral.hpp
    SourceData.cpp
    SourceData.hpp
    SourceData.hpp
    ExcitationInput.cpp
    ExcitationInput.hpp
    MatSolve.cpp
    MatSolve.hpp
    Factorization.cpp
    Factorization.hpp
    main.cpp
    buildZ.hpp
    buildZ.cpp
)

include(GNUInstallDirs)
install(TARGETS Mininec3_object_test
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
