#include "Gosub102.hpp"
#include "Gosub113.hpp"
#include "Gosub135.hpp"
#include <cmath>

void Gosub102( SimulationState& S, double& T1, double& T2)
{
    const GeometryData g = S.geom;

    // 100 REM ----- S(M) GOES IN (X1,Y1,Z1) FOR VECTOR POTENTIAL
    // 101 REM ----- MOD FOR SMALL RADIUS TO WAVE LENGTH RATIO

    int FVS = 0;                                        // 102

    if (S.K >= 1)                                       // 103
    {
        if (g.A[S.P4] < S.SRM)                          // 104
        {
            if ((S.I == S.J) && (S.P3 - (S.P2 + 0.5)))  // 105
            {
                T1 = std::log(g.Sa[S.P4] / g.A[S.P4]);  // 106
                T2 = -S.W * g.Sa[S.P4] / 2.0;           // 107
                return;
            }
        }
    }

    double X1 = g.Xa[(int)S.P1];    // 109
    double Y1 = g.Ya[(int)S.P1];    // 110
    double Z1 = g.Za[(int)S.P1];    // 111

    double X2, Y2, Z2;
    double V1, V2, V3;


    // 112 REM ----- S(U)-S(M) GOES IN (X2,Y2,Z2)
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
