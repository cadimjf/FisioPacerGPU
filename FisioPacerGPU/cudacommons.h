#pragma once

#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
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