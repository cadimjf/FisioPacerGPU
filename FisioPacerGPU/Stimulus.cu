#include "Stimulus.h"



//stimulus regions
int stimSize;
t_stim* aStim;

t_stim getItemStimulus(int i) {
    if (i>= stimSize || i< 0) {
        stringstream ss;
        ss << "Invalid index for a stimulus: " << i << endl;
        ss << "Valid range: [0: " << stimSize<<"]" << endl;
        string str = ss.str();
        throw MyException(str, __FILE__, __LINE__);
    }
    return aStim[i];

}
/*
*/
int isStimulationTime(int pmRegion, double t, double dt) {
    t_stim s = getItemStimulus(pmRegion);
    return (((t >= s.iniTime) && (((t - s.iniTime) - (floor(((t - s.iniTime) / s.period)) * s.period)) <= dt))) ? 1 : 0;
}

double stimGetIniTime(int i) { return getItemStimulus(i).iniTime; }
double stimGetPeriod(int i) { return getItemStimulus(i).period; }
double stimGetIniX(int i) { return getItemStimulus(i).iniX; }
double stimGetEndX(int i) { return getItemStimulus(i).endX; }
double stimGetIniY(int i) { return getItemStimulus(i).iniY; }
double stimGetEndY(int i) { return getItemStimulus(i).endY; }
double stimGetIniZ(int i) { return getItemStimulus(i).iniZ; }
double stimGetEndZ(int i) { return getItemStimulus(i).endZ; }
int stimGetSize() { return stimSize; }

/*
*/
void stimDealloc() {
    if (aStim != NULL) free(aStim);
}


/**
 *
 * @param strFileSt
 * @param CA
 */
void readStimFile(string strFileSt) {
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
        aStim = (t_stim*)malloc(sizeof(t_stim) * stimSize);
        if (!aStim) {
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
            aStim[cont] = stim;
            
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
}
