#pragma once
#include "MyStdLib.h"
#include "LinearAlg.h"
#include "cudacommons.h"
#include "Points.h"
typedef struct str_stats {
    double max[3];
    double min[3];
    double avgVel;
    double maxVol;
    double minVol;
    double volIni;
    double volMaxDelta;
    double maxDeltaVol;
    int contSave;
    double tol;
    double err;
} typ_stats;

__host__ double* statsGetMax();
__host__ double* statsGetMin();
__host__ double statsGetAvgVel();
__host__ double statsGetMaxVol();
__host__ double statsGetMinVol();
__host__ double statsGetVolIni();
__host__ double statsGetVolMaxDelta();
__host__ double statsGetMaxDeltaVol();
__host__ int statsGetContSave();
__host__ double statsGetTol();
__host__ double statsGetErr();
__host__ void statsSetTol(double);
__host__ void statsSetErr(double);


__host__ double statsGetVolIni();

__host__ void allocStats();
__host__ void deallocStats();
__host__ void iniStats();
__host__ void statsSetVolIni(double vol);
__host__ void statsIncContSave();
