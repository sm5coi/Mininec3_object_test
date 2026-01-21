#include "SimulationState.hpp"
#include <cmath>

void Gosub113(
    SimulationState& S,
    double X1, double Y1, double Z1,
    double& X2, double& Y2, double& Z2,
    double& V1, double& V2, double& V3
    )
{
    const GeometryData g = S.geom;

    // ---- S(U) – S(M)
    // 112 REM ----- S(U)-S(M) GOES IN (X2,Y2,Z2)
    int I4 = static_cast<int>(std::floor(S.P2));  // 113

    if (static_cast<double>(I4) != S.P2)          // 114
    {
        int I5 = I4 + 1;                        // 115
        X2 = (g.Xa[I4] + g.Xa[I5]) / 2.0 - X1;      // 116
        Y2 = (g.Ya[I4] + g.Ya[I5]) / 2.0 - Y1;      // 117
        Z2 = S.K * (g.Za[I4] + g.Za[I5]) / 2.0 - Z1;// 118
    }
    else
    {
        X2 = g.Xa[(int)S.P2] - X1;                  // 120
        Y2 = g.Ya[(int)S.P2] - Y1;                  // 121
        Z2 = S.K * g.Za[(int)S.P2] - Z1;            // 122
    }

    // ---- S(V) – S(M)
    // 123 REM ----- S(V)-S(M) GOES IN (V1,V2,V3)
    I4 = static_cast<int>(std::floor(S.P3));      // 124

    if (static_cast<double>(I4) != S.P3)          // 125
    {
        int I5 = I4 + 1;                        // 126
        V1 = (g.Xa[I4] + g.Xa[I5]) / 2.0 - X1;      // 127
        V2 = (g.Ya[I4] + g.Ya[I5]) / 2.0 - Y1;      // 128
        V3 = S.K * (g.Za[I4] + g.Za[I5]) / 2.0 - Z1;// 129
    }
    else
    {
        V1 = g.Xa[(int)S.P3] - X1;                  // 131
        V2 = g.Ya[(int)S.P3] - Y1;                  // 132
        V3 = S.K * g.Za[(int)S.P3] - Z1;            // 133
    }
}
