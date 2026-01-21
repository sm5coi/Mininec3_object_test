#include "Gosub135.hpp"
#include "KernelEval.hpp"
#include <cmath>

void Gosub135(
    SimulationState& S,
    double X2, double Y2, double Z2,
    double V1, double V2, double V3,
    int FVS,
    double& T1,
    double& T2)
{
    const GeometryData g = S.geom;

    // Gaussian quadrature coefficients Q
    static const double Q[15] = {0.0,
        0.288675135, 0.5, 0.430568156, 0.173927423,
        0.169990522, 0.326072577, 0.480144928, 0.050614268,
        0.398333239, 0.111190517, 0.262766205, 0.156853323,
        0.091717321, 0.181341892
    };

    // 134 REM ----- MAGNITUDE OF S(U) - S(M)
    double D0 = std::sqrt(X2*X2 + Y2*Y2 + Z2*Z2);   // 135, 137
    // 136 REM ----- MAGNITUDE OF S(V) - S(M)
    double D3 = std::sqrt(V1*V1 + V2*V2 + V3*V3);   // 138, 139

    // 140 REM ----- SQUARE OF WIRE RADIUS
    // 141 A2 = A(P4) * A(P4) In KernelEval()
    // 142 REM ----- MAGNITUDE OF S(V) - S(U)

    double S4 = (S.P3 - S.P2) * g.Sa[S.P4];        // 143

    // 144 REM ----- ORDER OF INTEGRATION
    // 145 REM ----- LTH ORDER GAUSSIAN QUADRATURE

    T1 = 0.0;                           // 146
    T2 = 0.0;                           // 147
    double I6u = 0.0;                   // 148
    double F2 = 1.0;                    // 149

    int L = 7;                          // 150

    double T = (D0 + D3) / g.Sa[S.P4];  // 151

    // criteria for exact kernel
    // 152 REM ----- CRITERIA FOR EXACT KERNEL
    if ((T > 1.1) || (S.CSd == 'N'))      // 153, 154
    {
        if (T > 6)   L = 3;             // 165
        if (T > 10)  L = 1;             // 166
    }
    else
    {
        // Check J2a coincidence conditions
        bool cond =
            (g.J2a[g.Wp[S.I]][1] == g.J2a[g.Wp[S.J]][1]) || // 155
            (g.J2a[g.Wp[S.I]][1] == g.J2a[g.Wp[S.J]][2]) || // 156
            (g.J2a[g.Wp[S.I]][2] == g.J2a[g.Wp[S.J]][1]) || // 157
            (g.J2a[g.Wp[S.I]][2] == g.J2a[g.Wp[S.J]][2]);   // 158

        if (cond)
        {
            if (g.A[S.P4] <= S.SRM)                         // 160
            {
                if (FVS == 1)                       // 161
                {
                    T1 = 2 * std::log(g.Sa[S.P4] / g.A[S.P4]);   // 91
                    T2 = -S.W * g.Sa[S.P4];                  // 92
                    return;                         // 93
                }
                else
                {
                    T1 = std::log(g.Sa[S.P4] / g.A[S.P4]);       // 106
                    T2 = -S.W * g.Sa[S.P4] / 2.0;            // 107
                    return;                         // 108
                }
            }

            F2 = 2 * (S.P3 - S.P2);                                     // 162
            I6u = (1 - std::log(S4 / F2 / 8.0 / g.A[S.P4])) / M_PI / g.A[S.P4]; // 163
        }
        else                        // 159
        {
            if (T > 6)   L = 3;     // 165
            if (T > 10)  L = 1;     // 166
        }
    }

    int I5 = L + L;                 // 167

    while (true)
    {
        double T3 = 0.0;            // 168
        double T4 = 0.0;            // 169

        double TT = (Q[L] + 0.5) / F2;          // 170
        KernelEval(S, X2, Y2, Z2, V1, V2, V3,   // 171
                   TT, I6u,
                   T3, T4);

        TT = (0.5 - Q[L]) / F2;                 // 172
        KernelEval(S, X2, Y2, Z2, V1, V2, V3,   // 173
                   TT, I6u,
                   T3, T4);
        L++;                // 174
        T1 += Q[L] * T3;    // 175
        T2 += Q[L] * T4;    // 176
        L++;                // 177

        if (L >= I5) break; // 178
    }

    T1 = S4 * (T1 + I6u);   // 179
    T2 = S4 * T2;           // 180
}
