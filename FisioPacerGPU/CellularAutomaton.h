/* 
 * File:   CellularAutomaton.h
 * Author: ricardo
 *
 * Created on April 25, 2013, 9:26 AM
 */

#ifndef CELLULARAUTOMATON_H
#define	CELLULARAUTOMATON_H

#include "Constants.h"
#include "Geometry.h"
#include "DataStructures.h"
#include "MyStdLib.h"
#include "Stimulus.h"
//#include "pressure.h"
#include <iostream>
#include <stdio.h>
#include "omp.h"

using namespace std;



/**
 * 
 * @param el
 * @param i
 * @param iRegion
 * @return 
 */
double getActiveTension(int i, typ_ca *CA);
double getActiveTensionNormalized(int i, typ_ca *CA);
double getActiveTensionDiscrete( int i, typ_ca *CA);

double getV(int i ,typ_ca *CA);
double getVdiscret( int i,typ_ca *CA);

double getPropagationTime(int iElemNeighbor, int i, typ_ca *CA);
void cellActivation(int i, typ_ca *CA);
void incStates(int i, double dt, typ_ca *CA);
void CAStep_i(int i, typ_ca* CA);

int getNumThreads(int numThreads, int *threadsByIndividual);


void computeNewAPDElectroTonic(int i, typ_ca *CA);
void restartAPDElectroTonic(int i, typ_ca *CA);
double getAVGNeighborAPDElectroTonic(int i, typ_ca *CA);
double getIoutElectroTonic(int i ,typ_ca *CA);
double newAPDElectroTonic(double iout, double apd, double mapd);
#endif	/* CELLULARAUTOMATON_H */

