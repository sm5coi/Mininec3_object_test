#include "Factorization.hpp"
#include "MatSolve.hpp"
#include "SourceData.hpp"
#include "buildZ.hpp"
#include <cmath>
#include <iostream>

void ImpedanceMatrix(SimulationState& S)
{
    const GeometryData g = S.geom;

    // init real & imag Z parts
    S.ZR.assign(g.N+1, std::vector<double>(g.N+1, 0.0));
    S.ZI.assign(g.N+1, std::vector<double>(g.N+1, 0.0));

    std::cout << "Entry of buildZ" << std::endl;
    buildZ( S);

    // AdditionOfLoads()

    std::cout << "Entry of Factorization" << std::endl;
    Factorization( S);

    std::cout << "Entry of MatSolve" << std::endl;
    MatSolve( S);

    std::cout << "Entry of SourceData" << std::endl;
    SourceData( S);

    std::cout << "Size of S.CurrX: " << S.CurrX.size() << "\n";


}
