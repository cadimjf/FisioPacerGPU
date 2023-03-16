#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "MyStdLib.h"
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


void readStimFile(string strFileSt);

double stimGetIniTime(int i);
double stimGetPeriod(int i);
double stimGetIniX(int i);
double stimGetEndX(int i);
double stimGetIniY(int i);
double stimGetEndY(int i);
double stimGetIniZ(int i);
double stimGetEndZ(int i);
int stimGetSize();

void stimDealloc();
int isStimulationTime(int pmRegion, double t, double dt);