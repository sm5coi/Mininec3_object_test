//#include "ImpedanceMatrix.hpp"
#include "Factorization.hpp"
#include "MatSolve.hpp"
#include "SourceData.hpp"
#include "buildZ.hpp"
#include <cmath>
#include <iostream>
//#include <iomanip>
//#include <complex>




void ImpedanceMatrix(SimulationState& S)
{
    const GeometryData g = S.geom;

    // init real & imag Z parts
    S.ZR.assign(g.N+1, std::vector<double>(g.N+1, 0.0));
    S.ZI.assign(g.N+1, std::vector<double>(g.N+1, 0.0));

    buildZ( S);

    // AdditionOfLoads()

    Factorization( S);

    MatSolve( S);

    SourceData( S);

    std::cout << "Size of S.CurrX: " << S.CurrX.size() << "\n";


}
