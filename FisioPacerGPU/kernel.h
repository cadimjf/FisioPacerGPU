#pragma once
#include "cudacommons.h"
#include <vector>

#include <stdio.h>
#include "Constants.h"
#include "Geometry.h"
#include "DataStructures.h"
#include "Mechanics.h"
#include "MyStdLib.h"
#include "LinearAlg.h"
#include "ReadMesh.h"
//#include "Stopwatch.h"
#include "WriteOutput.h"
#include "NumericalMethods.h"
#include "Stats.h"
#include "Stimulus.h"

int simulate(typ_ca* CA, bool save);
void initializeCA(typ_ca* CA);
void computeStats(double);
int startCA(string paramFile, bool save);
void simulationStep(typ_ca* CA, double* forcesOnPts);
void allocCA(typ_ca* CA);
void deallocCA(typ_ca* CA);
