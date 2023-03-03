/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "NumericalMethods.h"

/**
 * 
 * @param CA
 * @param forcesOnPts
 */
void EulerMethod(typ_ca *CA, double *forcesOnPts){
    for(int k=0;k<CA->params->pointsNum;k++){
        //Euler Method
        EulerOnPoint(CA, k, forcesOnPts);

    }
}

/**
 * 
 * @param CA
 * @param k
 * @param forcesOnPts
 */
void EulerOnPoint(typ_ca *CA, int k, double *forcesOnPts) 
{
    //
    if (!CA->pnts_old[k].xRestr) {
        double aX = ACELERATION(forcesOnPts[I2d(k,0,3)], CA->pnts_old[k].mass);
        CA->pnts_new[k].xV  = EULER_STEP(CA->pnts_old[k].xV, aX, CA->params->dt);
        CA->pnts_new[k].x  = EULER_STEP(CA->pnts_old[k].x, CA->pnts_old[k].xV, CA->params->dt);
    } else CA->pnts_new[k].x = CA->pnts_old[k].x + CA->pnts_old[k].xPreMov;
    //
    if (!CA->pnts_old[k].yRestr) {
        double aY = ACELERATION(forcesOnPts[I2d(k,1,3)], CA->pnts_old[k].mass);
        CA->pnts_new[k].yV  = EULER_STEP(CA->pnts_old[k].yV, aY, CA->params->dt);
        CA->pnts_new[k].y   = EULER_STEP(CA->pnts_old[k].y, CA->pnts_old[k].yV, CA->params->dt);
    } else CA->pnts_new[k].y = CA->pnts_old[k].y + CA->pnts_old[k].yPreMov;
    //
    if (!CA->pnts_old[k].zRestr) {
        double aZ = ACELERATION(forcesOnPts[I2d(k,2,3)], CA->pnts_old[k].mass);
        CA->pnts_new[k].zV  = EULER_STEP(CA->pnts_old[k].zV, aZ, CA->params->dt);
        CA->pnts_new[k].z   = EULER_STEP(CA->pnts_old[k].z, CA->pnts_old[k].zV, CA->params->dt);
    } else CA->pnts_new[k].z = CA->pnts_old[k].z + CA->pnts_old[k].zPreMov;

    
}
/**
 * 
 * https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
 * http://www.physics.udel.edu/~bnikolic/teaching/phys660/numerical_ode/node5.html
 * @param CA
 * @param forcesOnPts
 */
void VelocityVerletMethod(typ_ca *CA, double *forcesOnPts,  double *forcesOnPts_interm){
    double greatError = 0.0, greatTol=0.0;
    //computer the first step of velocity verlet method:
    //s(t+1) = s(t) + v(t)*h + 0.5*a(t)*h^2
    double error, tol;
    for(int k=0;k<CA->params->pointsNum;k++){
        VelocityVerletMethod_pt1(CA, k, forcesOnPts, &error, &tol);
        greatError  = my_max(error, greatError);
        greatTol    = my_max(tol,   greatTol);
        //copy velocities from p_old to points_interm, in order to find the proper forces on next step
        copyPoint(&(CA->pnts_old[k]), &(CA->pnts_intrm[k]));
        //update the points_interm x,y and z values, with the 1st step from velocity verlet
        CA->pnts_intrm[k].x = CA->pnts_new[k].x;
        CA->pnts_intrm[k].y = CA->pnts_new[k].y;
        CA->pnts_intrm[k].z = CA->pnts_new[k].z;
        CA->pnts_intrm[k].xV = (CA->pnts_new[k].x-CA->pnts_old[k].x)/CA->params->dt;
        CA->pnts_intrm[k].yV = (CA->pnts_new[k].y-CA->pnts_old[k].y)/CA->params->dt;
        CA->pnts_intrm[k].zV = (CA->pnts_new[k].z-CA->pnts_old[k].z)/CA->params->dt;        
        //initialize intermediate forces        
        forcesOnPts_interm[I2d(k,0,3)]=forcesOnPts_interm[I2d(k,1,3)]=forcesOnPts_interm[I2d(k,2,3)]=0.0f;
    }
    //swap intermediate step to old, in order to compute next step
    typ_point *aux1 = CA->pnts_old;
    CA->pnts_old = CA->pnts_intrm;
    CA->pnts_intrm = aux1;
    //now computes the forces again with s(t+1)
    computeForceIntermediaria(CA, forcesOnPts_interm);
    for(int k=0;k<CA->params->pointsNum;k++){
        VelocityVerletMethod_pt2(CA, k, forcesOnPts, forcesOnPts_interm);
        //forcesOnPts[I2d(k,0,3)]=forcesOnPts[I2d(k,1,3)]=forcesOnPts[I2d(k,2,3)]=0.0f;
    }
    //swap old to intermediate step, for the next iteration
    typ_point *aux2 = CA->pnts_intrm;
    CA->pnts_intrm = CA->pnts_old;
    CA->pnts_old   = aux2;
    getNewDt(greatError, greatTol, CA);
    CA->stats->tol = greatTol;
    CA->stats->err = greatError;
}
/**
 * 
 * @param error
 * @param tol
 * @param CA
 */
void getNewDt(double error, double tol, typ_ca *CA){
    double newDt = CA->params->dt;
    double deltaDT =0.0;
    if(fabs(error)<=0.00001){
        deltaDT=1.005;
    }else{
        deltaDT = sqrt(tol/error);
    }
    newDt = CA->params->dt*deltaDT;
    
    if(newDt>=CA->params->dtIni){
        CA->params->dt = newDt;
        CA->params->dt_sq = newDt*newDt;
    }else{
        CA->params->dt = CA->params->dtIni;
        CA->params->dt_sq = CA->params->dtIni*CA->params->dtIni;
    }
}
/**
 * 
 * @param CA
 * @param k
 * @param forcesOnPts
 */
void VelocityVerletMethod_pt1(typ_ca *CA,int k, double *forcesOnPts, double *error, double *tol) {
    double aErr[3]={0.0}, aTol[3]={0.0};
    if (!CA->pnts_old[k].xRestr) {
        double aX = ACELERATION(forcesOnPts[I2d(k,0,3)], CA->pnts_old[k].mass);
        CA->pnts_new[k].x = CA->pnts_old[k].x + CA->params->dt * CA->pnts_old[k].xV + 0.5 * aX*CA->params->dt_sq;
        aErr[0] = CA->pnts_new[k].x - EULER_STEP(CA->pnts_old[k].x, aX, CA->params->dt);
        aErr[0] = fabs(aErr[0]);
        aTol[0] = getMaxTol(fabs(CA->pnts_new[k].x));
    } else{
        CA->pnts_new[k].x = CA->pnts_old[k].x + CA->pnts_old[k].xPreMov;
    }
    //
    if (!CA->pnts_old[k].yRestr) {
        double aY = ACELERATION(forcesOnPts[I2d(k,1,3)], CA->pnts_old[k].mass);
        CA->pnts_new[k].y = CA->pnts_old[k].y + CA->params->dt * CA->pnts_old[k].yV + 0.5 * aY*CA->params->dt_sq;
        aErr[1] = CA->pnts_new[k].y - EULER_STEP(CA->pnts_old[k].y, aY, CA->params->dt);
        aErr[1] = fabs(aErr[1]);
        aTol[1] = getMaxTol(fabs(CA->pnts_new[k].y));
    } else{
        CA->pnts_new[k].y = CA->pnts_old[k].y + CA->pnts_old[k].yPreMov;
    } 
    //
    if (!CA->pnts_old[k].zRestr) {
        double aZ = ACELERATION(forcesOnPts[I2d(k,2,3)], CA->pnts_old[k].mass);
        CA->pnts_new[k].z = CA->pnts_old[k].z + CA->params->dt * CA->pnts_old[k].zV + 0.5 * aZ*CA->params->dt_sq;
        aErr[2] = CA->pnts_new[k].z - EULER_STEP(CA->pnts_old[k].z, aZ, CA->params->dt);
        aErr[2] = fabs(aErr[2]);
        aTol[2] = getMaxTol(fabs(CA->pnts_new[k].z));
    } else{
        CA->pnts_new[k].z = CA->pnts_old[k].z + CA->pnts_old[k].zPreMov;
    }
    *(error) = my_max(aErr[2], my_max(aErr[1], aErr[0]));
    *tol     = my_max(aTol[2], my_max(aTol[1], aTol[0]));
}
/**
 * 
 * @param val
 * @return 
 */
double getMaxTol(double val){
    return my_max(RELTOL_DT*fabs(val), ABSTOL_DT);
    
}
/**
 * 
 * @param CA
 * @param k
 * @param forcesOnPts
 * @param forcesOnPts_interm
 */
void VelocityVerletMethod_pt2(typ_ca *CA,
        int k, double *forcesOnPts, double *forcesOnPts_interm) {
    
    if (!CA->pnts_old[k].xRestr) {
        double aaX = ACELERATION(forcesOnPts[I2d(k,0,3)]+forcesOnPts_interm[I2d(k,0,3)], CA->pnts_old[k].mass);
        CA->pnts_new[k].xV = CA->pnts_old[k].xV + 0.5*(aaX)*CA->params->dt;
    } 
    if (!CA->pnts_old[k].yRestr) {
        double aaY = ACELERATION(forcesOnPts[I2d(k,1,3)]+forcesOnPts_interm[I2d(k,1,3)], CA->pnts_old[k].mass);
        CA->pnts_new[k].yV = CA->pnts_old[k].yV + 0.5*(aaY)*CA->params->dt;
    } 
    if (!CA->pnts_old[k].zRestr) {
        double aaZ = ACELERATION(forcesOnPts[I2d(k,2,3)]+forcesOnPts_interm[I2d(k,2,3)], CA->pnts_old[k].mass);
        CA->pnts_new[k].zV = CA->pnts_old[k].zV + 0.5*(aaZ)*CA->params->dt;
    }
}