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

double* statsGetMax();
double* statsGetMin();
double statsGetAvgVel();
double statsGetMaxVol();
double statsGetMinVol();
double statsGetVolIni();
double statsGetVolMaxDelta();
double statsGetMaxDeltaVol();
int statsGetContSave();
double statsGetTol();
double statsGetErr();
void statsSetTol(double);
void statsSetErr(double);


double statsGetVolIni();

__host__ void allocStats();
__host__ void deallocStats();
__host__ void iniStats();
__host__ void statsSetVolIni(double vol);
void statsIncContSave();
