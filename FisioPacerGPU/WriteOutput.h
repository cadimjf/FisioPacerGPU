/* 
 * File:   WriteOutput.h
 * Author: ricardo
 *
 * Created on May 7, 2013, 5:33 PM
 */

#ifndef WRITEOUTPUT_H
#define	WRITEOUTPUT_H

#include <stdio.h>
#include "DataStructures.h"
#include "Constants.h"
#include "MyStdLib.h"
#include "CellularAutomaton.h"
#include "pressure.h"
using namespace std;
#include "WriteOutput.h"

void saveCSV(typ_ca *CA, string outFolder);
void saveVTK_V_cell( int t, typ_ca *CA, string outFolder);
void saveVTK_V_point( typ_ca *CA, string outFile);

void saveVTK_Simples( int t, typ_ca *CA, string outFolder);
void saveDebug(typ_ca *CA, string outFile, double *forcesOnPts);
void save_step(FILE *fileDt, typ_ca *CA, string strFolderOut, double *forcesOnPts);
#endif	/* WRITEOUTPUT_H */
