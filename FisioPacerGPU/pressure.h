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
void gpucopyPressureCA();

void allocPressureCA();
void deallocPressureCA();
void iniPressure();
double getPressurePercent();
void incPressureStates(double dt);
double getPressureDiscrete();

#endif /* PRESSURE_H */

