/* 
 * File:   Geometry.h
 * Author: Ricardo
 *
 * Created on 6 de Março de 2014, 19:11
 */

#ifndef GEOMETRY_H
#define	GEOMETRY_H
#include "DataStructures.h"
#include "CellularAutomaton.h"
#include "MyStdLib.h"
#include "ContinuumMech.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <iostream>

//#ifdef	__cplusplus
//extern "C" {
//#endif
void updateAlpha(int i, typ_ca *CA);
//
double areaTriangulo(double A[3], double B[3], double C[3]);
//
double getVolumeTetrahedron( typ_ca *CA, int i);

void getBarycenter( typ_ca *CA, int index);
//
 void iniGeometry(int i, typ_ca *CA);

void updateAxisDirectionOLD(typ_ca *CA, int i);
//
void findAllNormals(typ_ca *CA, int i, double n[4][3]);
//
void findVelocities(typ_ca *CA, int i, 
         double velsIter[6][3]);
//
void updateGeometry(
    int i, typ_ca *CA);

void debugpnts(typ_ca *CA, int);
    
void getInterceptPtsByInterpolation(typ_ca *CA, int i);


void getElementPtsInVector(int i, typ_ca *CA, 
        double pt1[3], double pt2[3], double pt3[3], double pt4[3]);
bool checksAreasGetEthaKsi(double ptI[3], double ptA[3], double ptB[3], 
        double ptC[3], double *ksi, double *etha, double *relEr);

void getProjection(double inA[3], double inB[3], double out[3], double *);

void findAxis(typ_ca *CA, int i, 
         double axis1[3], double axis2[3], double axis3[3], double interPts[3][6]);

void getPtsInVector(typ_point *pts, int i, double pt[3]);
    
void findFaceNormal(typ_ca *CA, int pt1, int pt2, int pt3, double n[3]);
void findIntersectionPlaneLine(int index, typ_ca *CA, 
        double reta[3][2], 
        double n[4][3],
        int iTipoFibra);
 void computeKs(typ_ca *CA, int iElem);
void fixNormal(double vPt1[3], double vPt2[3], double vPt3[3], double normal[3], double bary[3]);
void getFaceCenter(double vPt1[3], double vPt2[3], double vPt3[3], double bary[3]);
void iniVolumes(int iElem, typ_ca *CA);
void findDataByInterpol(double inData[3][4], double ck[4][6], double dataOut[3][6]);
 //#ifdef	__cplusplus
//}
//}
//#endif

#endif	/* GEOMETRY_H */

