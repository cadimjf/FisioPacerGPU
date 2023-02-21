/* 
 * File:   ContinuumMech.h
 * Author: ricardo
 *
 * Created on December 22, 2014, 5:16 PM
 */

#ifndef CONTINUUMMECH_H
#define	CONTINUUMMECH_H

#include "LinearAlg.h"
#include "MyStdLib.h"
#include <iostream>
#include <stdio.h>
 
void assemblyXMatrix(typ_point *points, double X[3][3], int i, typ_ca *CA);
void getDeformationTensor(double X0[3][3], double X1[3][3], double F[3][3]);
void computePolarDecompositition(double F[3][3], double R[3][3], double U[3][3]);
void getTensors(double X0[3][3], double X1[3][3], double F[3][3], 
         double R[3][3], double U[3][3] );
#endif	/* CONTINUUMMECH_H */

