#include "Gosub87.hpp"
#include "Gosub113.hpp"
#include "Gosub135.hpp"
#include <cmath>

void Gosub87( SimulationState& S, double& T1, double& T2)
{
    const GeometryData g = S.geom;

    // 84 REM ----- ENTRIES REQUIRED FOR IMPEDANCE MATRIX CALCULATION
    // 85 REM ----- S(M) GOES IN (X1,Y1,Z1) FOR SCALAR POTENTIAL
    // 86 REM ----- MOD FOR SMALL RADIUS TO WAVE LENGTH RATIO

    int FVS = 1;                // 87

    if (S.K >= 1)               // 88
    {
        if (g.A[S.P4] <= S.SRM) // 89
        {
            if ((S.P3 == S.P2 + 1.0) && (S.P1 == (S.P2 + S.P3) / 2.0))  // 90
            {
                T1 = 2.0 * std::log(g.Sa[S.P4] / g.A[S.P4]);            // 91
                T2 = -S.W * g.Sa[S.P4];                                 // 92
                return;                                                 // 93
            }
        }
    }

    int I4 = static_cast<int>(std::floor(S.P1));    // 94
    int I5 = I4 + 1;                                // 95

    double X1 = (g.Xa[I4] + g.Xa[I5]) / 2.0;        // 96
    double Y1 = (g.Ya[I4] + g.Ya[I5]) / 2.0;        // 97
    double Z1 = (g.Za[I4] + g.Za[I5]) / 2.0;        // 98

    double X2, Y2, Z2;
    double V1, V2, V3;

    Gosub113(
        S,
        X1, Y1, Z1,
        X2, Y2, Z2,
        V1, V2, V3
        );

    Gosub135(
        S,
        X2, Y2, Z2,
        V1, V2, V3,
        FVS,
        T1, T2
        );
}
