
#include "Pressure.h"
#define P1 0
#define P2 1
#define P3 2
#define P4 3
typ_press* pressureHost;
__device__ typ_press* pressureDevice;
typ_press *gpuPntr;



__host__ __device__ typ_press* getPressureVar() {
    //is it device ou host
#if defined(__CUDA_ARCH__)
    // Device
    return pressureDevice;
#else
    // Host
    return pressureHost;
#endif    
}


int pressureGetNumFaces() { return getPressureVar()->numFaces; }
double pressureGetPressure(){ return getPressureVar()->pressure; }
void pressureSetNumFaces(int n)     { getPressureVar()->numFaces = n; }
void pressureSetPressure(double p)  { getPressureVar()->pressure=p;   }
/*
* gpu copy
*/
__host__ void gpucopyPressureCA() {
    if (GPUMODE == 1)
    {
        chker(cudaMemcpy(gpuPntr, pressureHost, sizeof(typ_press), cudaMemcpyHostToDevice));
    }
}

/**
*
*/
__host__ void allocPressureCA()
{
    cudaGetSymbolAddress((void**)&gpuPntr, pressureDevice);
    allocateHostVar(&pressureHost);
    allocateDeviceVar(&gpuPntr);
    

}
/**
*/
__host__ void deallocPressureCA() {//
    dealloc(pressureHost, pressureDevice);
    //if (aFaces != NULL) free(aFaces);

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
__device__ __host__ void incPressureStates(double dt)
{
    typ_press* pressure = getPressureVar();

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

__global__ void incPressureStatesincGPU( double dt) {
    incPressureStates( dt);
}
/*
* 
*/
void pressureStep(double dt) {
    if (GPUMODE == 0) {
        incPressureStates(dt);
    } else {
        incPressureStatesincGPU <<<1, 1 >>> (dt);
        //chker(cudaMemcpy(pressureHost, pressureDevice, sizeof(typ_press), cudaMemcpyDeviceToHost));
       // cout << pressureHost->time << endl;
    }
}



/**
 *
 * @param strFilePress
 * @param CA
 */
void readPressFile(string strFilePress, typ_ca* CA) {
    pressureSetNumFaces(0);
    //temp string to read lines
    string line;
    ifstream myfile2(strFilePress.c_str());

    if (myfile2.is_open()) {
        //first line contains the file's number of lines
        if (!getline(myfile2, line)) {
            stringstream ss;
            ss << "PRess File: failed to read line - " << strFilePress.c_str();
            string str = ss.str();
            throw MyException(str, __FILE__, __LINE__);
        }
        pressureSetNumFaces(atoi(line.c_str()));
        if (pressureGetNumFaces() <= 0) {
            stringstream ss;
            ss << "Invalid number of faces: " << pressureGetNumFaces();
            string str = ss.str();
            throw MyException(str, __FILE__, __LINE__);
            return;
        }
        getline(myfile2, line);
        pressureSetPressure(atof(line.c_str()));
        CA->params->aFaces = (typ_pressureface**)malloc(sizeof(typ_pressureface*) * pressureGetNumFaces());
        if (CA->params->aFaces == NULL) {
            stringstream ss;
            ss << "CA->params->aFaces: allocation: " << line;
            throw MyException(ss.str(), __FILE__, __LINE__);
        }
        int cont = 0;
        while (myfile2.good())
        {
            if (!getline(myfile2, line)) {
                stringstream ss;
                ss << "Press File: failed to read line: " << line;
                string str = ss.str();
                throw MyException(str, __FILE__, __LINE__);
                return;
            }
            stringstream strstream(line);
            string token;
            int countColumns = 0;
            typ_pressureface* face = (typ_pressureface*)malloc(sizeof(typ_pressureface));
            if (face == NULL) {
                stringstream ss;
                ss << "face File: allocation: " << line;
                throw MyException(ss.str(), __FILE__, __LINE__);
            }
            while (getline(strstream, token, ' ')) {
                switch (countColumns) {
                case 0:
                    face->pt1 = atoi(token.c_str());
                    CA->pnts_new[face->pt1].presFcs++;
                    CA->pnts_old[face->pt1].presFcs++;
                    break;
                case 1:
                    face->pt2 = atoi(token.c_str());
                    CA->pnts_new[face->pt2].presFcs++;
                    CA->pnts_old[face->pt2].presFcs++;
                    break;
                case 2:
                    face->pt3 = atoi(token.c_str());
                    CA->pnts_new[face->pt3].presFcs++;
                    CA->pnts_old[face->pt3].presFcs++;
                    break;
                default:
                    cout << "Invalid press file: " << countColumns << " columns!" << endl;
                }
                countColumns++;
            }

            CA->params->aFaces[cont] = face;
            if (cont >= pressureGetNumFaces()) {
                stringstream ss;
                ss << "There are more pres lines than allowed: " << pressureGetNumFaces();
                string str = ss.str();
                cout << str << endl;
                //throw MyException(str, __FILE__, __LINE__);
                break;
            }
            cont++;
        }
        myfile2.close();
    }
    else {
        cout << "Unable to open file " << strFilePress << endl << " running without pressures files" << endl;
        pressureSetNumFaces(0);
    }
}