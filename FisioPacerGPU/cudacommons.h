#pragma once

#include <iostream>
#include <cstdio>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define GPUMODE 1
/*
 * EXECMODE:
 * 0 stands for HOST 
 * 1 stadands for device
*/
using namespace std;


#define chker(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true)
{
    if (code != cudaSuccess)
    {
        cout << "GPUassert:"<< cudaGetErrorString(code) << " in " << file << " : " <<line <<endl;
        if (abort) exit(code);
    }
}

/*
    Receives two pointers - one to host and other to device. 
    According to EXECMODE, return the proper pointer 
*/
template< typename T > T* getDataByMode(T* host, T* device)
{
    //HOST execution mode ==0
    if (GPUMODE == 0){
        return host;
    }else if (GPUMODE==1){//DEVICE
        return device;
    }else{
        return NULL;
    }
}
/**
* allocate device and host variables  
*/
template <typename T> void allocateHostVar(T** host) {
    *host = (T*)malloc(sizeof(T));
    if (*host == NULL) {
        throw MyException("Allocation failure for stats structure.", __FILE__, __LINE__);
    }
}

template <typename T> void allocateDeviceVar(T** device, int size=1) {
    if (GPUMODE == 1) {
        chker(cudaMalloc((void**)device, sizeof(T)*size));
    }    
}
/*
free memory
*/
template <typename T> void dealloc(T* host, T* device) {
    if (host != NULL) free(host);
    if (GPUMODE == 1) {
        chker(cudaFree(device));
    }
}

