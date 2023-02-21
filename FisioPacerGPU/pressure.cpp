
#include "pressure.h"
#define P1 0
#define P2 1
#define P3 2
#define P4 3
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
void iniPressure(typ_ca *CA){
    CA->pressureCA->state = P1;
    CA->pressureCA->time = 0.0;
    CA->pressureCA->times[P1] = 0.05;//0.21;
    CA->pressureCA->times[P2] = 0.01;//0.35;
    CA->pressureCA->times[P3] = 0.20;//0.5;
    CA->pressureCA->times[P4] = 0.72;//0.68
    CA->pressureCA->vals[P1] = 0.01;//0.1;
    CA->pressureCA->vals[P2] = 0.90;//1.0;
    CA->pressureCA->vals[P3] = 1.00;//0.95;
    CA->pressureCA->vals[P4] = 0.1; 
}

double interpolVals(double v1, double v2, double t1, double t2, double t){
    double m = (v2 - v1) / (t2 - t1);
    return v1 + m * (t - t1);
}

double interpolPress(typ_ca *CA, int i1, int i2){
    double v2 = CA->pressureCA->vals[i2];
    double t2 = CA->pressureCA->times[i2];
    if(i1==-1)
        return interpolVals(0.0, v2, 0.0, t2, CA->pressureCA->time);
    else{
        double v1 = CA->pressureCA->vals[i1];
        double t1 = CA->pressureCA->times[i1];
        return interpolVals(v1, v2, t1, t2, CA->pressureCA->time);
    }
        
}

double getPressureDiscrete(typ_ca *CA){
    return CA->pressureCA->vals[CA->pressureCA->state];
}

double getPressurePercent(typ_ca *CA){
   
    switch(CA->pressureCA->state){
        case P1:
            return interpolPress(CA, -1, P1);
        case P2:
            return interpolPress(CA, P1, P2);
        case P3:
            return interpolPress(CA, P2, P3);
        case P4:
            return interpolPress(CA, P3, P4);
    }
    return 0.0;
}

void incPressureStates(typ_ca *CA)
{
    
    typ_press *pr = CA->pressureCA;
    pr->time += CA->params->dt;
    switch(pr->state){
        case P1:
            pr->state = (pr->time >= pr->times[P1])?P2:P1;
            break;
        case P2:
            pr->state = (pr->time >= pr->times[P2])?P3:P2;
            break;
        case P3:
            pr->state = (pr->time >= pr->times[P3])?P4:P3;
            break;
        case P4:
            if(pr->time >= pr->times[P4]){
                pr->state = P1;
                pr->time = 0.0;
            }
            break;
    }

}