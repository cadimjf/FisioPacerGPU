#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "MyStdLib.h"
#include "cudacommons.h"
using namespace std;
//stimulus
typedef struct str_sim {
    double iniTime;
    double period;
    double iniX;
    double endX;
    double iniY;
    double endY;
    double iniZ;
    double endZ;
} t_stim;


__host__ void readStimFile(string strFileSt);
__host__ double stimGetIniTime(int i);
__host__ double stimGetPeriod(int i);
__host__ double stimGetIniX(int i);
__host__ double stimGetEndX(int i);
__host__ double stimGetIniY(int i);
__host__ double stimGetEndY(int i);
__host__ double stimGetIniZ(int i);
__host__ double stimGetEndZ(int i);
__host__ int stimGetSize();

__host__ void stimDealloc();
__host__ __device__ int isStimulationTime(int pmRegion, double t, double dt);