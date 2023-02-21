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

#ifndef PRESSURE_H
#define PRESSURE_H


void iniPressure(typ_ca *CA);
double getPressurePercent(typ_ca *CA);
void incPressureStates(typ_ca *CA);
double getPressureDiscrete(typ_ca *CA);

#endif /* PRESSURE_H */

