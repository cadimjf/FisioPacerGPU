/* 
 * File:   LinearAlg.h
 * Author: ricardo
 *
 * Created on December 28, 2014, 9:07 PM
 */
#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <math.h>

#include "Constants.h"
#include "MyStdLib.h"
using namespace std;

#ifndef LINEARALG_H
#define	LINEARALG_H

#ifdef	__cplusplus
extern "C" {
#endif
//http://arxiv.org/abs/physics/0610206
//http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html

inline void HouseHolder(double A[3][3], double Q[3][3], double d[3], double e[2]);
int eigen(double A[3][3], double Q[3][3], double w[3]);
double det(double a[3][3]);
void printMatrix(double a[3][3]);
void computeSqrtMatrix(double U2[3][3], double U[3][3]);
void cross(double A[3], double B[3], double crossProd[3]);
double my_norm (double v[3]);
double dot(double[3], double[3]);
double normalizeVector(double v[3], int i);
void transpose(double M[3][3], double N[3][3]);
void invertMatrix(double M[3][3], double N[3][3]);

void invertMatrixEigenValue(double M[3][3], double N[3][3]);
void invertMatrixFormula(double M[3][3], double N[3][3], double detM);
void matrixMultiplication(double A[3][3], double B[3][3], double C[3][3]);
bool matrixCheck(double M[3][3], double N[3][3]);

#ifdef	__cplusplus
}
#endif

#endif	/* LINEARALG_H */

