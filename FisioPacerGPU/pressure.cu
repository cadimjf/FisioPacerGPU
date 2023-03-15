
#include "pressure.h"
#define P1 0
#define P2 1
#define P3 2
#define P4 3
typ_press* pressureCA;
typ_press* devicePressureCA;

__host__ void gpucopyPressureCA(){

    chker(cudaMemcpy(devicePressureCA, pressureCA, sizeof(typ_press), cudaMemcpyHostToDevice));

}


__host__ void allocPressureCA() {
    pressureCA = (typ_press*)malloc(sizeof(typ_press));
    if (pressureCA == NULL) {
        throw MyException("Allocation failure for stats structure.", __FILE__, __LINE__);
    }
    chker(cudaMalloc((void**)&devicePressureCA, sizeof(typ_press)));
}

__host__ void deallocPressureCA() {
    if (pressureCA != NULL) free(pressureCA);
    chker(cudaFree(devicePressureCA));
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
void iniPressure(){
    pressureCA->state = P1;
    pressureCA->time = 0.0;
    pressureCA->times[P1] = 0.05;//0.21;
    pressureCA->times[P2] = 0.01;//0.35;
    pressureCA->times[P3] = 0.20;//0.5;
    pressureCA->times[P4] = 0.72;//0.68
    pressureCA->vals[P1] = 0.01;//0.1;
    pressureCA->vals[P2] = 0.90;//1.0;
    pressureCA->vals[P3] = 1.00;//0.95;
    pressureCA->vals[P4] = 0.1; 
}

double interpolVals(double v1, double v2, double t1, double t2, double t){
    double m = (v2 - v1) / (t2 - t1);
    return v1 + m * (t - t1);
}

double interpolPress( int i1, int i2){
    double v2 = pressureCA->vals[i2];
    double t2 = pressureCA->times[i2];
    if(i1==-1)
        return interpolVals(0.0, v2, 0.0, t2, pressureCA->time);
    else{
        double v1 = pressureCA->vals[i1];
        double t1 = pressureCA->times[i1];
        return interpolVals(v1, v2, t1, t2, pressureCA->time);
    }
        
}

double getPressureDiscrete(){
    return pressureCA->vals[pressureCA->state];
}

double getPressurePercent(){
   
    switch(pressureCA->state){
        case P1:
            return interpolPress(-1, P1);
        case P2:
            return interpolPress(P1, P2);
        case P3:
            return interpolPress(P2, P3);
        case P4:
            return interpolPress(P3, P4);
    }
    return 0.0;
}

void incPressureStates(double dt)
{   
    pressureCA->time += dt;
    switch(pressureCA->state){
        case P1:
            pressureCA->state = (pressureCA->time >= pressureCA->times[P1])?P2:P1;
            break;
        case P2:
            pressureCA->state = (pressureCA->time >= pressureCA->times[P2])?P3:P2;
            break;
        case P3:
            pressureCA->state = (pressureCA->time >= pressureCA->times[P3])?P4:P3;
            break;
        case P4:
            if(pressureCA->time >= pressureCA->times[P4]){
                pressureCA->state = P1;
                pressureCA->time = 0.0;
            }
            break;
    }

}