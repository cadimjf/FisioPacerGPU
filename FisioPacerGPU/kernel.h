#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

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

int simulate(typ_ca* CA, bool save);
void initializeCA(typ_ca* CA);
void stats(typ_ca* CA, double* forcesOnPts);
int startCA(string paramFile, bool save);

void allocCA(typ_ca* CA);
void deallocCA(typ_ca* CA);
