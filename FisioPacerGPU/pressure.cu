
#include "pressure.h"
#define P1 0
#define P2 1
#define P3 2
#define P4 3
typ_press* pressureHost;
typ_press* pressureDevice;

/*
* gpu copy
*/
__host__ void gpucopyPressureCA() {
    if (GPUMODE == 1)
    {
        chker(cudaMemcpy(pressureDevice, pressureHost, sizeof(typ_press), cudaMemcpyHostToDevice));
    }
}

/**
*
*/
__host__ void allocPressureCA()
{
    allocateHostVar(&pressureHost);
    allocateDeviceVar(&pressureDevice);
}
/**
*/
__host__ void deallocPressureCA() {//
    dealloc(pressureHost, pressureDevice);
}

/*
 */
__host__ void iniPressureHost()
{
    pressureHost->state = P1;
    pressureHost->time = 0.0;
    pressureHost->times[P1] = 0.05;//0.21;
    pressureHost->times[P2] = 0.01;//0.35;
    pressureHost->times[P3] = 0.20;//0.5;
    pressureHost->times[P4] = 0.72;//0.68
    pressureHost->vals[P1] = 0.01;//0.1;
    pressureHost->vals[P2] = 0.90;//1.0;
    pressureHost->vals[P3] = 1.00;//0.95;
    pressureHost->vals[P4] = 0.1;
}
/**
*/
__device__ __host__ double interpolVals(double v1, double v2, double t1, double t2, double t)
{
    double m = (v2 - v1) / (t2 - t1);
    return v1 + m * (t - t1);
}
/**
*/
__host__ double interpolPress(int i1, int i2)
{
    double v2 = pressureHost->vals[i2];
    double t2 = pressureHost->times[i2];
    if (i1 == -1)
        return interpolVals(0.0, v2, 0.0, t2, pressureHost->time);
    else {
        double v1 = pressureHost->vals[i1];
        double t1 = pressureHost->times[i1];
        return interpolVals(v1, v2, t1, t2, pressureHost->time);
    }
}
/*
*/
__host__ double getPressureDiscrete()
{
    return pressureHost->vals[pressureHost->state];
}
/**
*/
__host__ double getPressurePercent()
{
    switch (pressureHost->state) {
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
/*
*/
__device__ __host__ void incPressureStates(typ_press* pressure, double dt)
{
    pressure->time += dt;
    switch (pressure->state) {
    case P1:
        pressure->state = (pressure->time >= pressure->times[P1]) ? P2 : P1;
        break;
    case P2:
        pressure->state = (pressure->time >= pressure->times[P2]) ? P3 : P2;
        break;
    case P3:
        pressure->state = (pressure->time >= pressure->times[P3]) ? P4 : P3;
        break;
    case P4:
        if (pressure->time >= pressure->times[P4]) {
            pressure->state = P1;
            pressure->time = 0.0;
        }
        break;
    }
}

__global__ void incPressureStatesincGPU(typ_press* pressure, double dt) {
    incPressureStates(pressure, dt);
}
/*
* 
*/
void pressureStep(double dt) {
    if (GPUMODE == 0) {
        incPressureStates(pressureHost, dt);
    } else {
        incPressureStatesincGPU <<<1, 1 >>> (pressureDevice, dt);
        chker(cudaMemcpy(pressureHost, pressureDevice, sizeof(typ_press), cudaMemcpyDeviceToHost));
    }
}
