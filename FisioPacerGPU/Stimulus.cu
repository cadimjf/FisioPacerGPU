#include "Stimulus.h"

//stimulus regions
int stimSize;
t_stim* aStimHost;
__device__ t_stim* aStimDevice;
t_stim* gpuPntr;

/*
* returns the array of stimulus
* checks if the index is inside the valid range
* returns the data on gpu or cpu according to the function caller (is it host or device?)
*/
__host__ __device__ t_stim getItemStimulus(int i) {
        //is it device ou host
        #if defined(__CUDA_ARCH__)
            // Device
            return aStimDevice[i];
        #else
            // Host
            if (i<0 || i>stimSize) {
                stringstream ss;
                ss << "Invalid index for a stimulus: " << i << endl;
                ss << "Valid range: [0: " << stimSize << "]" << endl;
                string str = ss.str();
                throw MyException(str, __FILE__, __LINE__);
            }
            else {
                return aStimHost[i];
            }
        #endif    
}
/*
*/
__host__ __device__ int isStimulationTime(int pmRegion, double t, double dt) {
    t_stim s = getItemStimulus(pmRegion);
    return (((t >= s.iniTime) && (((t - s.iniTime) - (floor(((t - s.iniTime) / s.period)) * s.period)) <= dt))) ? 1 : 0;
}

__host__ double stimGetIniTime(int i) { return getItemStimulus(i).iniTime; }
__host__ double stimGetPeriod(int i) { return getItemStimulus(i).period; }
__host__ double stimGetIniX(int i) { return getItemStimulus(i).iniX; }
__host__ double stimGetEndX(int i) { return getItemStimulus(i).endX; }
__host__ double stimGetIniY(int i) { return getItemStimulus(i).iniY; }
__host__ double stimGetEndY(int i) { return getItemStimulus(i).endY; }
__host__ double stimGetIniZ(int i) { return getItemStimulus(i).iniZ; }
__host__ double stimGetEndZ(int i) { return getItemStimulus(i).endZ; }
__host__ int stimGetSize() { return stimSize; }

/*
*/
__host__ void stimDealloc() {
    if (aStimHost != NULL) free(aStimHost);
    if (GPUMODE == 1) {
        // TODO FIXME chker(cudaFree(gpuPntr));
    }
}


/**
 *
 * @param strFileSt
 * @param CA
 */
__host__ void readStimFile(string strFileSt) {
    //temp string to read lines
    stimSize = 0;
    string line;
    ifstream myfile2(strFileSt.c_str());

    if (myfile2.is_open()) {
        //first line contains the file's number of lines
        if (!getline(myfile2, line))  throw MyException("Stim File: failed to read line", __FILE__, __LINE__);
        stimSize = atoi(line.c_str());
        if (stimSize < 0) {
            stringstream ss;
            ss << "Invalid number of stimulus: " << stimSize << endl << strFileSt;
            string str = ss.str();
            throw MyException(str, __FILE__, __LINE__);
            return;
        }
        aStimHost = (t_stim*)malloc(sizeof(t_stim) * stimSize);
        if (!aStimHost) {
            stringstream ss;
            ss << "Stim allocation problem: " << endl << strFileSt;
            string str = ss.str();
            throw MyException(str, __FILE__, __LINE__);
        }
        int cont = 0;
        while (myfile2.good() && stimSize > 0)
        {
            if (!getline(myfile2, line)) {
                stringstream ss;
                ss << "Number of stimulus: " << stimSize << endl;
                ss << strFileSt << endl;
                ss << "Line number: " << cont << endl;
                ss << "Stim File: failed to read line: " << line;
                string str = ss.str();
                throw MyException(str, __FILE__, __LINE__);
                return;
            }
            stringstream strstream(line);
            string token;
            int countColumns = 0;
            t_stim stim;
            while (getline(strstream, token, ' ')) {
                switch (countColumns) {
                case 0:
                    stim.iniTime = atof(token.c_str());
                    break;
                case 1:
                    stim.period = atof(token.c_str());
                    break;
                case 2:
                    stim.iniX = atof(token.c_str());
                    break;
                case 3:
                    stim.endX = atof(token.c_str());
                    break;
                case 4:
                    stim.iniY = atof(token.c_str());
                    break;
                case 5:
                    stim.endY = atof(token.c_str());
                    break;
                case 6:
                    stim.iniZ = atof(token.c_str());
                    break;
                case 7:
                    stim.endZ = atof(token.c_str());
                    break;
                default:
                    cout << "Error: " << strFileSt.c_str() << " - Invalid stim file: "
                        << countColumns << " columns! " << line << endl;
                }
                countColumns++;
            }
            aStimHost[cont] = stim;
            
            if (cont >= stimSize) {
                stringstream ss;
                ss << "There are more stim lines than allowed: " << stimSize << endl << strFileSt;
                string str = ss.str();
                throw MyException(str, __FILE__, __LINE__);
                break;
            }
            cont++;
        }
        myfile2.close();
    }
    else {
        stringstream ss;
        ss << "Unable to open file " << strFileSt;
        string str = ss.str();
        throw MyException(str, __FILE__, __LINE__);
    }

    if(GPUMODE==1){
        cudaGetSymbolAddress((void**)&gpuPntr, aStimDevice);
        allocateDeviceVar(&gpuPntr, stimSize);
        chker(cudaMemcpy(gpuPntr, aStimHost, sizeof(t_stim) * stimSize, cudaMemcpyHostToDevice));
    }
    
}
