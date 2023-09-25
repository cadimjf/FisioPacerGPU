/* 
 * File:   Mechanics.h
 * Author: ricardo
 *
 * Created on April 25, 2013, 9:27 AM
 */

#ifndef MECHANICS_H
#define	MECHANICS_H
#include <omp.h>
#include "Geometry.h"
#include "DataStructures.h"
#include "CellularAutomaton.h"
//#include "MyStdLib.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <cstdlib>
#include <iostream>
#include <stdio.h>


#include "MyStdLib.h"
#include "pressure.h"

double getForce_HookesLawAxis(double k, double lengthT, double length0) ;
double getForce_Damping(double rVel[3], double kdamp, double axis[3]) ;
void getForce_HookesLawAngular(double force[3], double deltaA, double kAng, double axis[3]) ;
double getForce_VolPreserving(double volT, double vol0, double kvol);
int MecStep_i(typ_ca *CA, int i, double *forcesOnPts);
void getPassiveForce(int i, int iAxis, double force[3][2], double vel1[3], double vel2[3] , typ_ca *CA);
void getExternalForce(double force[3][2], int i, typ_ca *CA);
void getActiveForce( int i, double force[3][2], typ_ca *CA);

void computeForceOnElement(typ_ca *CA, double *forcesOnPts, int i);
void computeForceIntermediaria(typ_ca *CA, double *forcesOnPts);
void getForceOnSurface(typ_pressureface *face, typ_ca *CA, double *aForce);
void computePressurePoints(typ_ca *CA, double *forcesOnPts);

void applyPressureForcePoint(typ_ca *CA, double surfaceForces[3], double *forcesOnPts, int iPnt);
double getFaceAreaNormal(typ_pressureface *face, typ_ca *CA, double normal[3], double bary[3]);
//void getForceOnSurfaceVirtualHex(typ_ca *CA, int i, int iAxis, double *aForce,typ_face *face);
//void getPressIntPnt(typ_ca *CA, int i, double forceFiber[3][2],
//                    double forceSheet[3][2], double forceNSheet[3][2]);
#endif	/* MECHANICS_H */

