/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   pressure.h
 * Author: Ricardo
 *
 * Created on 24 de Janeiro de 2020, 16:03
 */
#include "DataStructures.h"
#include "cudacommons.h"

#ifndef PRESSURE_H
#define PRESSURE_H


 //pressure ap
typedef struct str_press {
    double vals[4];
    double times[4];
    int state;
    double time;
    //pressure faces
    int numFaces;
    double pressure;

}typ_press;

void readPressFile(string strFilePress, typ_ca* CA);
int pressureGetNumFaces();
double pressureGetPressure();
void pressureSetNumFaces(int n);
void pressureSetPressure(double p);



void gpucopyPressureCA();
void allocPressureCA();
void deallocPressureCA();
void iniPressureHost();
double getPressurePercent();
//void incPressureStates(double dt);
double getPressureDiscrete();
void pressureStep(double dt);
void readPressFile(string strFilePress, typ_ca* CA);
#endif /* PRESSURE_H */

