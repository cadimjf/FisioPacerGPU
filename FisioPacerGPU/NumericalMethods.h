/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   NumericalMethods.h
 * Author: Ricardo
 *
 * Created on 18 de Fevereiro de 2019, 13:41
 */

#ifndef NUMERICALMETHODS_H
#define NUMERICALMETHODS_H
//returns euler time step
#define EULER_STEP(y_old, right_hand_side, dt) y_old + dt*right_hand_side;
#define ACELERATION(force, mass) force/mass;
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Mechanics.h"
#include "DataStructures.h"
#include "Constants.h"
#include "MyStdLib.h"
void EulerMethod(typ_ca *CA, double *forcesOnPts);
void EulerOnPoint(typ_ca *CA, int k, double *forcesOnPts);
void VelocityVerletMethod(typ_ca *CA, double *forcesOnPts, double*);
void VelocityVerletMethod_pt1(typ_ca *CA,int k, double *forcesOnPts, double *, double *);
void VelocityVerletMethod_pt2(typ_ca *CA,
        int k, double *forcesOnPts, double *forcesOnPts_interm) ;
void getNewDt(double error, double tol, typ_ca *CA);

double getMaxTol(double err);
#endif /* NUMERICALMETHODS_H */

