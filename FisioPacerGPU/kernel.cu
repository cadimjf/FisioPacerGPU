#include "kernel.h"



typ_ca* deviceCA;
typ_param* deviceParams;
typ_t0_element* deviceIniElem;
double* deviceForcesOnPtsInterm;
double* deviceForcesOnPts;


__global__ void deviceIniForces(double* devForcesOnPtsInterm, double* devForcesOnPts, typ_param* devParams) {
    for (int k = 0; k < devParams->pointsNum; k++) {
        devForcesOnPts[I2d(k, 0, 3)] = devForcesOnPts[I2d(k, 1, 3)] = devForcesOnPts[I2d(k, 2, 3)] = 0.0f;
        devForcesOnPtsInterm[I2d(k, 0, 3)] = devForcesOnPtsInterm[I2d(k, 1, 3)] = devForcesOnPtsInterm[I2d(k, 2, 3)] = 0.0f;
    }
    
}


void deviceAlloc(typ_ca *CA) {

    chker(cudaMalloc((void**)&deviceCA,         sizeof(typ_ca)));
    
    chker(cudaMalloc((void**)&deviceParams,     sizeof(typ_param)));
    chker(cudaMalloc((void**)&deviceIniElem,        sizeof(typ_t0_element) * CA->params->elementsNum));
    chker(cudaMalloc((void**)&deviceForcesOnPts,sizeof(double) * CA->params->pointsNum * 3));
    chker(cudaMalloc((void**)&deviceForcesOnPtsInterm, sizeof(double) * CA->params->pointsNum * 3));
    /*
    typ_face** aFaces;

    int stimSize;
    t_stim** aStim; //
    int nRegions;
    t_par_ac** aParam;//
    //pressure faces
    int numFaces;
    typ_face** aFaces;//*/



    //gpuErrchk(cudaMalloc((void**)&CA_dev->params, sizeof(typ_param));
    //CA->params->aParam = (t_par_ac**)malloc(sizeof(t_par_ac*) * CA->params->nRegions);

    //gpuErrchk(cudaMalloc((void**)&CA_dev->params->aParam));
}

void deviceCopy(typ_ca *hCA) {
    chker(cudaMemcpy(deviceCA,         hCA,                 sizeof(typ_ca),     cudaMemcpyHostToDevice));
    gpucopyPressureCA();
    chker(cudaMemcpy(deviceParams,     hCA->params,         sizeof(typ_param), cudaMemcpyHostToDevice));
    chker(cudaMemcpy(deviceIniElem,    hCA->ini,            sizeof(typ_t0_element)*hCA->params->elementsNum, cudaMemcpyHostToDevice));
    for (int i = 0; i < pressureGetNumFaces(); i++) {
       // chker(cudaMemcpy(deviceAFaces[i], hCA->params->aFaces[i], sizeof(typ_face), cudaMemcpyHostToDevice));

    }
}


__global__ void teste(
    typ_ca* devCA)
{
    //int i = threadIdx.x;
    //devCA->time = deviceAFaces[0]->pt1;
//    devAFaces[0]-> = 666;;
}


void deviceDealloc( ) {
      
    /* TODO FIXME
    chker(cudaFree(deviceParams));
    chker(cudaFree(deviceIniElem));
    chker(cudaFree(deviceCA));
    chker(cudaFree(deviceForcesOnPts));
    chker(cudaFree(deviceForcesOnPtsInterm));
    */
}

/*
typedef struct str_cellularautomata{
    typ_t0_element   **ini;//
    typ_dt_element **t_old;//
    typ_dt_element **t_new;//
    lst_item **omega_b;  //
    lst_item **omega_a;//
    typ_param* params;//
    typ_point **pnts_new;//
    typ_point **pnts_old;//
    typ_point **pnts_intrm;//
}typ_ca;// automato cellular

*/


/**
 *
 * @param CA
 */
void initializeCA(typ_ca* CA)
{
    for (int i = 0; i < CA->params->elementsNum; i++)
    {
        iniVolumes(i, CA);
    }
    double mass = 0.0;
    //iterate over points to find their masses
    for (int i = 0; i < CA->params->pointsNum; i++) {
        mass = 0.0;
        lst_item* cur = CA->omega_b[i];
        while (cur != NULL) {
            mass += CA->ini[cur->value].volCel_ini;
            cur = cur->next;
        }
        CA->pnts_new[i].mass = CA->pnts_old[i].mass = ro_mass_dens * mass / 4.0;
    }
    for (int i = 0; i < CA->params->elementsNum; i++)
    {
        CA->t_old[i]->cellT = CA->t_new[i]->cellT = 0.0f;
        CA->t_old[i]->ECTNC_force_t_ini = CA->t_new[i]->ECTNC_force_t_ini = 0.0f;
        CA->t_old[i]->ECTNC_force_val_ini = CA->t_new[i]->ECTNC_force_val_ini = 0.0f;
        CA->t_old[i]->ECTNC_ap_t_ini = CA->t_new[i]->ECTNC_ap_t_ini = 0.0f;
        CA->t_old[i]->ECTNC_ap_val_ini = CA->t_new[i]->ECTNC_ap_val_ini = 0.0f;
        CA->t_old[i]->ECTNC_ap_t_end = CA->t_new[i]->ECTNC_ap_t_end = 0.0f;
        CA->t_old[i]->ECTNC_force_t_end = CA->t_new[i]->ECTNC_force_t_end = 0.0f;

        //initially, all cells are healthy
        CA->ini[i].cellCond = HEALTHY;
        //V state
        CA->t_new[i]->V_state = CA->t_old[i]->V_state = V0;
        //F state
        CA->t_new[i]->F_state = CA->t_old[i]->F_state = F0;
        iniGeometry(i, CA);
        //sets the pacemaker up
        for (int iS = 0; iS < stimGetSize(); iS++) {
            if (
                CA->t_old[i]->bary[0] > stimGetIniX(iS) && CA->t_old[i]->bary[0] < stimGetEndX(iS) &&
                CA->t_old[i]->bary[1] > stimGetIniY(iS) && CA->t_old[i]->bary[1] < stimGetEndY(iS) &&
                CA->t_old[i]->bary[2] > stimGetIniZ(iS) && CA->t_old[i]->bary[2] < stimGetEndZ(iS)
                ) {
                CA->ini[i].cellCond = paceMaker;
                CA->ini[i].pmRegion = iS;
            }
        }
        restartAPDElectroTonic(i, CA);
        CA->t_old[i]->APTime1 = CA->t_new[i]->APTime1;
        CA->t_old[i]->APTime2 = CA->t_new[i]->APTime2;
        CA->t_old[i]->APTime3 = CA->t_new[i]->APTime3;
        CA->t_old[i]->APTime4 = CA->t_new[i]->APTime4;

    }
    iniPressureHost();
}


/**
 *
 * @param CA
 * @param strFolderOut
 * @param finalTime
 * @param save
 * @param numThreads
 * @param threadsByIndividual
 * @return
 */
int simulate(typ_ca* CA, bool save) {
    //time units in SECONDS
    CA->time = 0.0;
    CA->timeSaving = 0.0;

    iniStats();

    int count = 0;
    typ_point* aux3;
    typ_dt_element** auxCA;
    char filename[255];
    sprintf(filename, "%sfisiopacer.txt", CA->params->outputFolder/*.c_str()*/);
    FILE* fileDt = fopen(filename, "w+");

    double* forcesOnPts_interm = (double*)malloc(sizeof(double) * CA->params->pointsNum * 3);
    if (forcesOnPts_interm == NULL) {
        throw MyException("Allocation failure for forcesOnPts_interm.", __FILE__, __LINE__);
    }
    double* forcesOnPts = (double*)malloc(sizeof(double) * CA->params->pointsNum * 3);
    if (forcesOnPts == NULL) {
        throw MyException("Allocation failure for forcesOnPts.", __FILE__, __LINE__);
    }
    for (int k = 0; k < CA->params->pointsNum; k++) {
        forcesOnPts[I2d(k, 0, 3)] = forcesOnPts[I2d(k, 1, 3)] = forcesOnPts[I2d(k, 2, 3)] = 0.0f;
        forcesOnPts_interm[I2d(k, 0, 3)] = forcesOnPts_interm[I2d(k, 1, 3)] = forcesOnPts_interm[I2d(k, 2, 3)] = 0.0f;
    }
   

    initializeCA(CA);
    double sumDt = 0.0;
    
    double vol = 0.0;
    for (int i = 0; i < CA->params->elementsNum; i++) {
         vol += CA->t_old[i]->volCel;
    }
    CA->volume = vol;
    statsSetVolIni(vol);
    int retValFinal = 0;

    

    computeStats(CA->volume);

    deviceAlloc(CA);
    deviceCopy(CA);
    if (GPUMODE == 1) {
        //zero the forces vector on gpu
        deviceIniForces <<<1, 1 >>> (deviceForcesOnPts, deviceForcesOnPtsInterm, deviceParams);

        teste <<< 1, 1 >>> (deviceCA);
        cudaMemcpy(CA, deviceCA, sizeof(typ_ca), cudaMemcpyDeviceToHost);
    }
    cout << "time: " << CA->time << endl;
    

    // exit(0);

    //cout<<"VOlume ini " << CA->volume<<endl;
    while (CA->time <= CA->params->simulationTime)
    {
        simulationStep(CA, forcesOnPts);
        if (CA->params->mecSim == 1) {
            // EulerMethod(CA, forcesOnPts);
            VelocityVerletMethod(CA, forcesOnPts, forcesOnPts_interm);
        }
        pressureStep(CA->params->dt);

        
           
        if (save) {
            save_step(fileDt, CA, CA->params->outputFolder, forcesOnPts);
        }
        for (int k = 0; k < CA->params->pointsNum; k++) {
            forcesOnPts[I2d(k, 0, 3)] = forcesOnPts[I2d(k, 1, 3)] = forcesOnPts[I2d(k, 2, 3)] = 0.0f;
        }
        if (GPUMODE == 1) {
            //zero the forces vector on gpu
            deviceIniForces <<<1, 1 >>> (deviceForcesOnPts, deviceForcesOnPtsInterm, deviceParams);
        }
        CA->time += CA->params->dt;
        auxCA = CA->t_old;
        CA->t_old = CA->t_new;
        CA->t_new = auxCA;
        //
        aux3 = CA->pnts_old;
        CA->pnts_old = CA->pnts_new;
        CA->pnts_new = aux3;
        //
        count++;
        sumDt += CA->params->dt;
    }// fim do while
    fclose(fileDt);
    //saveAPD(CA, strFolderOut, filename);
    if (forcesOnPts != NULL)
        free(forcesOnPts);
    forcesOnPts = NULL;
    if (forcesOnPts_interm != NULL)
        free(forcesOnPts_interm);
    forcesOnPts_interm = NULL;
    if (CA->params->printOutput == 1) {
        printf("final time: %.3f\n", CA->time);
        double iniVol = statsGetVolIni();

        printf("Volume: [%.3f %.3f] [%.3e %.3e] \n", statsGetMinVol() / iniVol * 100.0, statsGetMaxVol() / iniVol * 100.0, iniVol, CA->volume);
        printf("iterações: %d. Dt medio %g\n", count, sumDt / count);
    }
    deviceDealloc();
    return retValFinal;
}

/**
 *
 * @param nThreads
 * @param POINTS_OLD
 * @param forcesOnPts
 * @param time
 * @return
 */
void simulationStep(typ_ca* CA,
    double* forcesOnPts)
{
    double volT = 0.0;
    int contCA = 0;
    //#pragma omp parallel for schedule(static) num_threads(nThreads) reduction(+:volT)
    for (int i = 0; i < CA->params->elementsNum; i++) {
        if (CA->params->paSim == 1) {
            CAStep_i(i, CA);
            contCA++;
            computeNewAPDElectroTonic(i, CA);
        }
        computeForceOnElement(CA, forcesOnPts, i);
        //verifica se o volume da celula é menor que 1% do inicial - netste caso, mata o processo pq é sinal de erro.
        if (CA->t_new[i]->volCel < 0.001 * CA->ini[i].volCel_ini) {
            CA->params->mecSim = 0;
            CA->params->paSim = 0;
            if (CA->params->printOutput == 1) {
                cout << "Mata por volume pequeno[" << i << "]: " << CA->time << endl;
                cout << "Inicial: " << CA->ini[i].volCel_ini << " | Atual: " << CA->t_new[i]->volCel << endl;
                cout << CA->pnts_old[CA->ini[i].iPt1].x << " " <<
                    CA->pnts_old[CA->ini[i].iPt1].y << " " <<
                    CA->pnts_old[CA->ini[i].iPt1].z << endl;
                cout << CA->pnts_old[CA->ini[i].iPt2].x << " " <<
                    CA->pnts_old[CA->ini[i].iPt2].y << " " <<
                    CA->pnts_old[CA->ini[i].iPt2].z << endl;

                throw MyException("1percent volume.", __FILE__, __LINE__);

            }
        }
        volT += CA->t_new[i]->volCel;
    }//fim pragma
    //exit(0);
    if (CA->params->mecSim == 1) {
        computePressurePoints(CA, forcesOnPts);
    }
    CA->volume = volT;
}





/**
 *
 * @param numThreads
 * @param simulationTime
 * @param inputCarpFolder
 * @param outputFolder
 * @param paramFile
 * @param save
 * @param numSolutions
 * @param numRefPts
 * @param aRefIDPts
 * @param timeSolutions
 * @param aRefPts
 * @param threadsByIndividual
 * @return
 */
int startCA(string paramFile, bool save)
{
    try {


        //
        typ_ca* CA = (typ_ca*)malloc(sizeof(typ_ca));
        if (CA == NULL) {
            throw MyException("Allocation failure for CA.", __FILE__, __LINE__);
        }
        CA->params = (typ_param*)malloc(sizeof(typ_param));
        if (CA->params == NULL) {
            throw MyException("Allocation failure for parameter structure.", __FILE__, __LINE__);
        }
        readParameterFile(paramFile, CA);
        allocCA(CA);
        int retSim = simulate(CA, save);

        deallocCA(CA);
        return retSim;
    }
    catch (MyException& caught) {
        std::cout << caught.getMessage() << std::endl;
    }
    return 0;
}
/**
 *
 * @param CA
 */
void deallocCA(typ_ca* CA) {
    try {
        myArrayDeallocation<typ_dt_element>(CA->t_old, CA->params->elementsNum);
        myArrayDeallocation<typ_dt_element>(CA->t_new, CA->params->elementsNum);
        if (CA->ini         != NULL) free(CA->ini);
        if (CA->pnts_new    != NULL) free(CA->pnts_new);
        if (CA->pnts_old    != NULL) free(CA->pnts_old);
        if (CA->pnts_intrm  != NULL) free(CA->pnts_intrm);
        //
        freeAList(CA->omega_a, CA->params->elementsNum);
        freeAList(CA->omega_b, CA->params->pointsNum);
        stimDealloc();
        //
        //

        for (int i = 0; i < CA->params->nRegions; i++) {
            if (CA->params->aParam != NULL) {
                if (CA->params->aParam[i] != NULL) {
                    free(CA->params->aParam[i]);
                }
            }
        }
        if (CA->params->aParam != NULL) free(CA->params->aParam);
        if (CA->params != NULL) free(CA->params);
        deallocStats();
        deallocPressureCA();
        if (CA != NULL) free(CA);
    }
    catch (MyException& caught) {
        std::cout << caught.getMessage() << std::endl;
    }
}
/**
 *
 * @param strPtsFile
 * @param strEleFile
 * @param paramFile
 * @param strStmFile
 * @param strFibFile
 * @param strBoundFile
 * @param strPressFile
 * @return
 */
void allocCA(typ_ca* CA) {
    //input data
    string strPtsFile = CA->params->inputFolder;
    strPtsFile.append(".pts");
    //
    string strEleFile = CA->params->inputFolder;
    strEleFile.append(".elem");
    //
    string strFibFile = CA->params->inputFolder;
    strFibFile.append(".fib");
    //
    string strStmFile = CA->params->inputFolder;
    strStmFile.append(".stim");
    //
    string strBoundFile = CA->params->inputFolder;
    strBoundFile.append(".bound");
    //
    string strPressFile = CA->params->inputFolder;
    strPressFile.append(".press");

    try {
        //le o tamanho das malhas nos arquivos de entrada
        CA->params->pointsNum = readSize(strPtsFile);
        CA->params->elementsNum = readSize(strEleFile);
        
        
        CA->omega_a = iniAList(CA->params->elementsNum);
        if (CA->omega_a == NULL) {
            throw MyException("Allocation failure for CA->omega_a.", __FILE__, __LINE__);
        }
        CA->omega_b = iniAList(CA->params->pointsNum);
        if (CA->omega_b == NULL) {
            throw MyException("Allocation failure for CA->omega_b.", __FILE__, __LINE__);
        }
        //points
        CA->pnts_new = (typ_point*)malloc(CA->params->pointsNum * sizeof(typ_point));
        if (CA->pnts_new == NULL) {
            throw MyException("Allocation failure for CA->pnts_new.", __FILE__, __LINE__);
        }
        CA->pnts_old = (typ_point*)malloc(CA->params->pointsNum * sizeof(typ_point));
        if (CA->pnts_old == NULL) {
            throw MyException("Allocation failure for CA->pnts_old.", __FILE__, __LINE__);
        }
        CA->pnts_intrm = (typ_point*)malloc(CA->params->pointsNum * sizeof(typ_point));
        if (CA->pnts_intrm == NULL) {
            throw MyException("Allocation failure for CA->pnts_intrm.", __FILE__, __LINE__);
        }
        //elements
        CA->t_old = myArrayAllocation<typ_dt_element>(CA->params->elementsNum);
        if (CA->t_old == NULL) {
            throw MyException("Allocation failure for CA->t_old.", __FILE__, __LINE__);
        }
        CA->t_new = myArrayAllocation<typ_dt_element>(CA->params->elementsNum);
        if (CA->t_new == NULL) {
            throw MyException("Allocation failure for CA->t_new.", __FILE__, __LINE__);
        }
        CA->ini     = (typ_t0_element*)malloc(CA->params->elementsNum * sizeof(typ_t0_element));
        if (CA->ini == NULL) {
            throw MyException("Allocation failure for CA->ini.", __FILE__, __LINE__);
        }

        //

        allocStats();
        //opens the files and fill arrays
        allocPressureCA();
        openFile(CA, strPtsFile, strEleFile, strFibFile, strBoundFile, strPressFile, strStmFile);
        
    }
    catch (MyException& caught) {
        std::cout << caught.getMessage() << std::endl;
    }
}

/*
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);

__global__ void addKernel(int *c, const int *a, const int *b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

int main()
{
    const int arraySize = 5;
    const int a[arraySize] = { 1, 2, 3, 4, 5 };
    const int b[arraySize] = { 10, 20, 30, 40, 50 };
    int c[arraySize] = { 0 };

    // Add vectors in parallel.
    cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
        c[0], c[1], c[2], c[3], c[4]);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}
*/
