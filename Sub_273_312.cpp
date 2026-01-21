#include "Sub_273_312.hpp"
#include "Gosub87.hpp"
#include <cmath>

// Helper to safely fetch Sa for J1/J2
static inline double SafeSa(const std::vector<double>& Sa, int idx)
{
    if (idx <= 0) return 0.0;
    if (idx >= (int)Sa.size()) return 0.0;
    return Sa[idx];
}

void Sub_273_312(
    SimulationState& S,
    GeometryData g,
    int J1, int J2,
    double T1, double T2,
    double& U1, double& U2)
{
    //    auto& Sa  = g.Sa;
    // auto& SRM = S.SRM;
    // auto& Xa  = g.Xa;
    // auto& Ya  = g.Ya;
    // auto& Za  = g.Za;
    // auto& W   = S.W;
    auto& A   = g.A;

    // 273 REM ----- COMPUTE PSI(M+1/2,N,N+1)
    S.P1 = S.P1 + 0.5;              // 274
    if (S.F8 == 2) S.P1 = S.P1 - 1; // 275
    S.P2 = S.P3;                    // 276
    S.P3 = S.P3 + 1.0;              // 277
    S.P4 = J2;                      // 278

    double T1_local = T1;
    double T2_local = T2;

    double U5 = 0.0;
    double U6 = 0.0;

    // ===== FIRST CALL =====
    if (S.F8 != 1)                              // 279
    {
        Gosub87(S, T1_local, T2_local);         // 283

        if (S.F8 < 2)                           // 284
        {
            U5 = T1_local;                      // 288
            U6 = T2_local;                      // 289
        }
        else
        {
            U1 = (2*T1_local - 4*U1*S.F5) / g.Sa[J1];   // 285
            U2 = (2*T2_local - 4*U2*S.F5) / g.Sa[J1];   // 286
            return;
        }
    }
    else
    {
        U5 = S.F5 * U1 + T1_local;      // 280
        U6 = S.F5 * U2 + T2_local;      // 281
    }

    // ===== SECOND CALL =====

    // 290 REM ----- COMPUTE PSI(M-1/2,N,N+1)
    S.P1 = S.P1 - 1;                        // 291

    Gosub87(S, T1_local, T2_local);         // 292

    U1 = (T1_local - U5) / g.Sa[J2];    //293
    U2 = (T2_local - U6) / g.Sa[J2];    // 294
    //}

    // ===== THIRD CALL =====
    // BASIC 296–299
    // 295 REM ----- COMPUTE PSI(M+1/2,N-1,N)
    S.P1 = S.P1 + 1;        // 296
    S.P3 = S.P2;            // 297
    S.P2 = S.P2 - 1.0;      // 298
    S.P4 = J1;              // 299

    Gosub87(S, T1_local, T2_local);             // 300

    double U3 = T1_local;                       // 301
    double U4 = T2_local;                       // 302

    // BASIC 304–309
    // 303 REM ----- COMPUTE PSI(M-1/2,N-1,N)
    if (S.F8 < 1)                               // 304
    {
        S.P1 = S.P1 - 1;                        // 308

        Gosub87(S, T1_local, T2_local);         // 309
    }
    else
    {
        T1_local = U5;                          // 305
        T2_local = U6;                          // 306
    }

    // 310 REM ----- GRADIENT OF SCALAR POTENTIAL CONTRIBUTION

    U1 = U1 + (U3 - T1_local) / g.Sa[J1];       // 311
    U2 = U2 + (U4 - T2_local) / g.Sa[J1];       // 312
}
