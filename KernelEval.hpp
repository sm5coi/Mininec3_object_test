#ifndef KERNELEVAL_HPP
#define KERNELEVAL_HPP

#include "SimulationState.hpp"

void KernelEval(
    SimulationState& S,
    double X2, double Y2, double Z2,
    double V1, double V2, double V3,
    double T,
    double I6u,
    double& T3,
    double& T4
    );

#endif
