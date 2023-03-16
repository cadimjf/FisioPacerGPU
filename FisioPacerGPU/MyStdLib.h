/* 
 * File:   MyStdLib.h
 * Author: ricardo
 *
 * Created on April 23, 2013, 1:33 PM
 */

#ifndef MYSTDLIB_H
#define	MYSTDLIB_H

//typedef double myNumType;

#include "DataStructures.h"
#include "Constants.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <iostream>


using namespace std;
double my_max(double a, double b);
double my_min(double a, double b);

int midpoint(int imin, int imax);
int binary_search(int A[], int key, int imin, int imax);

void printParameters(typ_ca *CA);
void setDefaultFolders();
void copyPoint(typ_point *from, typ_point *to);

#endif	/* MYSTDLIB_H */

