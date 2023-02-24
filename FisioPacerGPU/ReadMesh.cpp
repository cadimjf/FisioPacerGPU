#include "ReadMesh.h"
//#include <bits/stdc++.h>
/**
 * 
 * @param iElement
 * @param CA
 */
void omegaBAdjacency(int iElement, typ_ca *CA){
    addItem(&(CA->omega_b[CA->ini[iElement]->iPt1]), iElement);
    addItem(&(CA->omega_b[CA->ini[iElement]->iPt2]), iElement);
    addItem(&(CA->omega_b[CA->ini[iElement]->iPt3]), iElement);
    addItem(&(CA->omega_b[CA->ini[iElement]->iPt4]), iElement);
}
/**
 * 
 * @param fileName
 * @param CA
 */
void fillFiberFile(string fileName, typ_ca *CA)
{
    ifstream myfile(fileName.c_str());
    //temp string to read lines
    string line;
    double aux_line;
    if (myfile.is_open())
    {
        //first line contains the file's number of lines
        if(!getline(myfile, line))  throw MyException("Fiber File: failed to read first line.", __FILE__, __LINE__);
        int fSize = atoi(line.c_str());
         if(fSize!=CA->params->elementsNum)
        {
            stringstream ss;
            ss <<"The number of Fibers elements does not match with fiber num: ";
            ss << fSize;
            ss << " / ";
            ss << CA->params->elementsNum;
            string str = ss.str();
            cout<<str;
            throw MyException(str, __FILE__, __LINE__);
        }
        //count the lines number
        int countLine = 0;

        //iterates over the points file
        while ( myfile.good() )
        {
            if(!getline(myfile, line))  break;
            stringstream strstream(line);
            string token;
            int countCoord=0;
            while (getline(strstream, token, ' ')) {
                aux_line = atof(token.c_str());
                if(countCoord==0){
                    CA->t_new[countLine]->fiberDir[0]  = CA->t_old[countLine]->fiberDir[0] = aux_line;
                }else if(countCoord==1){
                    CA->t_new[countLine]->fiberDir[1]  = CA->t_old[countLine]->fiberDir[1] = aux_line;
                }else if(countCoord==2){
                    CA->t_new[countLine]->fiberDir[2]  = CA->t_old[countLine]->fiberDir[2] = aux_line;
                }else if(countCoord==3){
                    CA->t_new[countLine]->sheetDir[0]  = CA->t_old[countLine]->sheetDir[0] = aux_line;
                }else if(countCoord==4){
                    CA->t_new[countLine]->sheetDir[1]  = CA->t_old[countLine]->sheetDir[1] = aux_line;
                }else if(countCoord==5){
                    CA->t_new[countLine]->sheetDir[2]  = CA->t_old[countLine]->sheetDir[2] = aux_line;
                }else if(countCoord==6){
                    CA->t_new[countLine]->nsheetDir[0]  = CA->t_old[countLine]->nsheetDir[0] = aux_line;
                }else if(countCoord==7){
                    CA->t_new[countLine]->nsheetDir[1]  = CA->t_old[countLine]->nsheetDir[1] = aux_line;
                }else if(countCoord==8){
                    CA->t_new[countLine]->nsheetDir[2]  = CA->t_old[countLine]->nsheetDir[2] = aux_line;
                }
                countCoord++;
            }
            countLine++;
            if(countLine>CA->params->elementsNum)
            {
//                cout<< "There are more lines than allowed: "<<p->elementsNum<<endl;
                break;
            }
        }
        myfile.close();
        if(countLine!=CA->params->elementsNum) {
            stringstream ss;
            ss <<"The number of Fibers elements does not match with elementNum: ";
            ss << countLine;
            ss << " / ";
            ss << CA->params->elementsNum;
            string str = ss.str();
            cout<<str;
//            throw MyException(str, __FILE__, __LINE__);
        }
    }
    else{
        stringstream ss;
        ss <<"Unable to open file: ";
        ss << fileName;
        string str = ss.str();
        throw MyException(str, __FILE__, __LINE__);
    }
    
}

/**
 *
 * @param fileName
 * @param points
 */
void pointsFile(string fileName, typ_ca *CA)
{
    ifstream myfile(fileName.c_str());
    //temp string to read lines
    string line;
    double aux_line;
    if (myfile.is_open())
    {
        //first line contains the file's number of lines
        if(!getline(myfile, line))  throw MyException("Points File: failed to read line", __FILE__, __LINE__);
        int fSize = atoi(line.c_str());
        if(CA->params->pointsNum!=fSize){
            stringstream ss;
            ss <<"Allocation error: points size on file: "<<fSize<<" Code constant: "<<CA->params->pointsNum<<endl;
            string str = ss.str();
            throw MyException(str, __FILE__, __LINE__);

        }
        //count the lines number
        int iLine = 0;
        //iterates over the points file
        while ( myfile.good() )
        {
            if(!getline(myfile, line))  throw MyException("Points File: failed to read line", __FILE__, __LINE__);
            stringstream strstream(line);
            string token;
            int countCoord=0;
            while (getline(strstream, token, ' ')) {
                aux_line = atof(token.c_str());
                if(countCoord==0)                    
                    CA->pnts_old[iLine]->x = CA->pnts_new[iLine]->x = aux_line;
                else if(countCoord==1)
                    CA->pnts_old[iLine]->y = CA->pnts_new[iLine]->y = aux_line;
                else if(countCoord==2)
                    CA->pnts_old[iLine]->z = CA->pnts_new[iLine]->z = aux_line;
                countCoord++;
            }
            CA->pnts_old[iLine]->mass    = CA->pnts_new[iLine]->mass      = 0.0f;

            CA->pnts_old[iLine]->xV      = CA->pnts_new[iLine]->xV        = 0.0f;
            CA->pnts_old[iLine]->yV      = CA->pnts_new[iLine]->yV        = 0.0f;
            CA->pnts_old[iLine]->zV      = CA->pnts_new[iLine]->zV        = 0.0f;

            CA->pnts_old[iLine]->xRestr  = CA->pnts_new[iLine]->xRestr    = false;
            CA->pnts_old[iLine]->yRestr  = CA->pnts_new[iLine]->yRestr    = false;
            CA->pnts_old[iLine]->zRestr  = CA->pnts_new[iLine]->zRestr    = false;

            CA->pnts_old[iLine]->xPreMov = CA->pnts_new[iLine]->xPreMov   = 0.0f;
            CA->pnts_old[iLine]->yPreMov = CA->pnts_new[iLine]->yPreMov   = 0.0f;
            CA->pnts_old[iLine]->zPreMov = CA->pnts_new[iLine]->zPreMov   = 0.0f;

            CA->pnts_new[iLine]->x0      = CA->pnts_old[iLine]->x0        = CA->pnts_old[iLine]->x;
            CA->pnts_new[iLine]->y0      = CA->pnts_old[iLine]->y0        = CA->pnts_old[iLine]->y;
            CA->pnts_new[iLine]->z0      = CA->pnts_old[iLine]->z0        = CA->pnts_old[iLine]->z;
            CA->pnts_new[iLine]->presFcs = CA->pnts_old[iLine]->presFcs   = 0;
            iLine++;
            if(iLine>=CA->params->pointsNum)
            {
//                cout<< "There are more lines than allowed: "<<p->pointsNum<<endl;
                break;
            }
        }
        myfile.close();
    }
    else{
        stringstream ss;
        ss <<"Unable to open file: "<< fileName;
        string str = ss.str();
        throw MyException(str, __FILE__, __LINE__);
    }
    
}
/**
 * 
 * @param fileName
 * @param CA
 */
void elementsFile(string fileName, typ_ca *CA)
{
    //temp string to read lines
    string line;
    ifstream myfile2(fileName.c_str());
    //
    if (myfile2.is_open())
    {
        //first line contains the file's number of lines
        if(!getline(myfile2, line))  throw MyException("Element File: failed to read line", __FILE__, __LINE__);
        int fSize = atoi(line.c_str());
        if(CA->params->elementsNum!=fSize){
            stringstream ss;
            ss <<"Allocation error: elements size on file: "<<fSize<<". Code constant: "<<CA->params->elementsNum<<endl;
            string str = ss.str();
            throw MyException(str, __FILE__, __LINE__);
        }
        int countEl = 0;
        while ( myfile2.good() )
        {

            if(!getline(myfile2, line))  throw MyException("Element File: failed to read line", __FILE__, __LINE__);
            stringstream strstream(line);
            string token;
            int countColumns=0;
            while (getline(strstream, token, ' ')) {
                switch(countColumns){
                    case 0:
                        //the first column is the element type:
                        if( token.compare("Tt")==0 ){
                            //element = new
                        }else
                            throw MyException("This file is not Tt", __FILE__, __LINE__);
                    break;
                    case 1:
                        CA->ini[countEl]->iPt1 = atoi(token.c_str());
                    break;
                    case 2:
                        CA->ini[countEl]->iPt2 = atoi(token.c_str());
                    break;
                    case 3:
                        CA->ini[countEl]->iPt3 = atoi(token.c_str());
                    break;
                    case 4:
                        CA->ini[countEl]->iPt4 = atoi(token.c_str());
                    break;
                    case 5://regiao
                        CA->ini[countEl]->iRegion = atoi(token.c_str())-1;
                    break;
                    default:
                        cout<<"Invalid position for an element"<<endl;
                }
                countColumns++;
            }
            
            //fills the neighboring array
            omegaBAdjacency(countEl, CA);
            countEl++;
            if(countEl>=CA->params->elementsNum)
            {
//                cout<< "There are more lines than allowed: "<<p->elementsNum<<endl;
                break;
            }
        }
        myfile2.close();
    }
    else{
        stringstream ss;
        ss <<"Unable to open file "<<fileName;
        string str = ss.str();
        throw MyException(str, __FILE__, __LINE__);

    }
    
}
/**
 * 
 * @param CA
 * @param strFilePts
 * @param strFileEle
 * @param strFileFib
 * @param boundFile
 * @param strPressFile
 */
void openFile(typ_ca *CA,
        string strFilePts, string strFileEle, string strFileFib, string boundFile,
        string strPressFile)
{

    using namespace std;
    //reads the points file
    pointsFile(strFilePts,CA);
    //reads the boundary conditions
    readBoundaryFile(boundFile,CA);
    //reads the elements file
    elementsFile(strFileEle,CA);
    OmegaANeighborhood(CA);
    fillFiberFile(strFileFib, CA);
    readPressFile(strPressFile, CA);
}
/**
 * 
 *
 * @param base
 * @param b1
 * @param b2
 * @param b3
 * @param b4
 * @return
 */
int hasAtLeastOnePtEqual(int base, int b1, int b2, int b3, int b4)
{
    if(base==b1 || base==b2 || base==b3 || base==b4)
        return 1;
    else
        return 0;
}

/**
 * 
 * @param el_old
 * @param p
 */
void OmegaANeighborhood(typ_ca *CA){
    for(int iElBase=0;  iElBase < CA->params->elementsNum;  iElBase++)
    {
        for(int iElBusca=0;  iElBusca < CA->params->elementsNum;  iElBusca++)
        {
            if(iElBusca!=iElBase)
            {
                int check=0;
                check += hasAtLeastOnePtEqual(CA->ini[iElBase]->iPt1, CA->ini[iElBusca]->iPt1,
                    CA->ini[iElBusca]->iPt2, CA->ini[iElBusca]->iPt3, CA->ini[iElBusca]->iPt4);
                check += hasAtLeastOnePtEqual(CA->ini[iElBase]->iPt2, CA->ini[iElBusca]->iPt1,
                    CA->ini[iElBusca]->iPt2, CA->ini[iElBusca]->iPt3, CA->ini[iElBusca]->iPt4);
                check += hasAtLeastOnePtEqual(CA->ini[iElBase]->iPt3, CA->ini[iElBusca]->iPt1,
                    CA->ini[iElBusca]->iPt2, CA->ini[iElBusca]->iPt3, CA->ini[iElBusca]->iPt4);
                check += hasAtLeastOnePtEqual(CA->ini[iElBase]->iPt4, CA->ini[iElBusca]->iPt1,
                    CA->ini[iElBusca]->iPt2, CA->ini[iElBusca]->iPt3, CA->ini[iElBusca]->iPt4);
                if(check>=1){
                    addItem(&(CA->omega_a[iElBase]), iElBusca);
                }
            }
        }

    }

}
/**
 * 
 * @param strFileSt
 * @param CA
 */
void readStimFile(string strFileSt, typ_ca * CA){
    //temp string to read lines
    CA->params->stimSize =0;
    string line;
    ifstream myfile2(strFileSt.c_str());

    if (myfile2.is_open()){
        //first line contains the file's number of lines
        if(!getline(myfile2, line))  throw MyException("Stim File: failed to read line", __FILE__, __LINE__);
        CA->params->stimSize = atoi(line.c_str());
        if(CA->params->stimSize<0){
            stringstream ss;
            ss <<"Invalid number of stimulus: "<<CA->params->stimSize<<endl<<strFileSt;
            string str = ss.str();
            throw MyException(str, __FILE__, __LINE__);
            return;
        }
        CA->params->aStim = (t_stim**)malloc(sizeof(t_stim*)*CA->params->stimSize);
        int cont = 0;
        while ( myfile2.good() && CA->params->stimSize>0)
        {
            if(!getline(myfile2, line)){
                stringstream ss;
                ss <<"Number of stimulus: "<<CA->params->stimSize<<endl;
                ss<<strFileSt<<endl;
                ss <<"Line number: "<<cont<<endl;
                ss <<"Stim File: failed to read line: "<<line;
                string str = ss.str();
                throw MyException(str, __FILE__, __LINE__);
                return;
            }
            stringstream strstream(line);
            string token;
            int countColumns=0;
            t_stim *stim = (t_stim*)malloc(sizeof(t_stim));
            while (getline(strstream, token, ' ')) {
                switch(countColumns){
                    case 0:
                        stim->iniTime = atof(token.c_str());
                    break;
                    case 1:
                        stim->period = atof(token.c_str());
                    break;
                    case 2:
                        stim->iniX = atof(token.c_str());
                    break;
                    case 3:
                        stim->endX = atof(token.c_str());
                    break;
                    case 4:
                        stim->iniY = atof(token.c_str());
                    break;
                    case 5:
                        stim->endY = atof(token.c_str());
                    break;
                    case 6:
                        stim->iniZ = atof(token.c_str());
                    break;
                    case 7:
                        stim->endZ = atof(token.c_str());
                    break;
                    default:
                        cout<<"Error: "<<strFileSt.c_str()<<" - Invalid stim file: "
                                <<countColumns<<" columns! "<<line<<endl;
                }
                countColumns++;
            }
            //fills the neighboring array
            CA->params->aStim[cont] = stim;
            if(cont>=CA->params->stimSize){
                stringstream ss;
                ss <<"There are more stim lines than allowed: "<<CA->params->stimSize<<endl<<strFileSt;
                string str = ss.str();
                throw MyException(str, __FILE__, __LINE__);
                break;
            }
            cont++;
        }
        myfile2.close();
    }else{
        stringstream ss;
        ss <<"Unable to open file "<<strFileSt;
        string str = ss.str();
        throw MyException(str, __FILE__, __LINE__);
    }
}

/**
 *
 * @param file
 * @return
 */
int readSize(string file){
    ifstream myfile(file.c_str());
    string line;
    if (myfile.is_open())
    {
        //first line contains the file's number of lines
        getline(myfile, line);
        myfile.close();
        return atoi(line.c_str());
    }else{
        stringstream ss;
        ss<<"File not found: "<<file<<endl;
        string str = ss.str();
        throw MyException(str, __FILE__, __LINE__);
    }
    return 0;
}
/**
 *
 * @param p
 */
void readParameterFile(string fName, typ_ca *CA){
    
    ifstream myfile(fName.c_str());
    string line;
    if (myfile.is_open())
    {
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read mecSim", __FILE__, __LINE__);
        CA->params->mecSim = atoi(line.c_str());
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read paSim", __FILE__, __LINE__);
        CA->params->paSim = atoi(line.c_str());
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read dt", __FILE__, __LINE__);
        CA->params->dt = atof(line.c_str());
        CA->params->dt_sq = CA->params->dt * CA->params->dt;
        CA->params->dtIni = CA->params->dt;
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read dtSave", __FILE__, __LINE__);
        CA->params->dtSave = atof(line.c_str());
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read Nregions", __FILE__, __LINE__);
        CA->params->nRegions = atoi(line.c_str());
        if(CA->params->nRegions<1)  throw MyException("Nregions has to be >=1", __FILE__, __LINE__);
        //aloca memoria para os parametros
        CA->params->aParam = (t_par_ac**)malloc(sizeof(t_par_ac*)*CA->params->nRegions);
        int contParam = 0;
        while ( myfile.good() )
        {
            if(!getline(myfile, line)){
                throw MyException("Stim File: failed to read line: ", __FILE__, __LINE__);
            }
            string token;
            int countColumns=0;
            stringstream strstream(line);
            t_par_ac *parAC = (t_par_ac*)malloc(sizeof(t_par_ac));
            while (getline(strstream, token, ' ')) {
                switch(countColumns){
                    case 0:
                        parAC->VV0 = atof(token.c_str());
                    break;
                    case 1:
                        parAC->VV1 = atof(token.c_str());
                    break;
                    case 2:
                        parAC->VV2 = atof(token.c_str());
                    break;
                    case 3:
                        parAC->VV3 = atof(token.c_str());
                    break;
                    case 4:
                        parAC->VV4 = atof(token.c_str());
                    break;
                    case 5:
                        parAC->VF0 = atof(token.c_str());
                    break;
                    case 6:
                        parAC->VF1 = atof(token.c_str());
                    break;
                    case 7:
                        parAC->VF2 = atof(token.c_str());
                    break;
                    case 8:
                        parAC->VF3 = atof(token.c_str());
                    break;
                    case 9:
                        parAC->VF4 = atof(token.c_str());
                    break;
                    case 10:
                        parAC->APTime1ini = atof(token.c_str());
                    break;
                    case 11:
                        parAC->APTime2ini = atof(token.c_str());
                    break;
                    case 12:
                        parAC->APTime3ini = atof(token.c_str());
                    break;
                    case 13:
                        parAC->APTime4ini = atof(token.c_str());
                    break;
                    case 14:
                        parAC->FTime1 = atof(token.c_str());
                    break;
                    case 15:
                        parAC->FTime2 = atof(token.c_str());
                    break;
                    case 16:
                        parAC->FTime3 = atof(token.c_str());
                    break;
                    case 17:
                        parAC->FTime4 = atof(token.c_str());
                    break;
                    case 18:
                        parAC->v[FIBER] = atof(token.c_str());
                    break;
                    case 19:
                        parAC->v[SHEET] = atof(token.c_str());
                    break;
                    case 20:
                        parAC->v[NORMAL] = atof(token.c_str());
                    break;
                    case 21:
                        parAC->forceMultiplier = atof(token.c_str());
                    break;
                    case 22:
                        parAC->EAxl[FIBER] = atof(token.c_str());
                    break;
                    case 23:
                        parAC->EAxl[SHEET] = atof(token.c_str());
                    break;
                    case 24:
                        parAC->EAxl[NORMAL] = atof(token.c_str());
                    break;
                    case 25:
                        parAC->EVol[FIBER] = atof(token.c_str());
                    break;
                    case 26:
                        parAC->EVol[SHEET] = atof(token.c_str());
                    break;
                    case 27:
                        parAC->EVol[NORMAL] = atof(token.c_str());
                    break;
                    case 28:
                        parAC->kDamp = atof(token.c_str());
                    break;
                    case 29:
                        parAC->EAng[FIBER] = atof(token.c_str());
                    break;
                    case 30:
                        parAC->EAng[SHEET] = atof(token.c_str());
                    break;
                    case 31:
                        parAC->EAng[NORMAL] = atof(token.c_str());
                    break;
                    case 32:
                        parAC->pEletroTonic = atof(token.c_str());
                        break;
                    default:
                        cout<<"Invalid param file: "<<countColumns<<" columns!"<<endl;
                }
                countColumns++;
            }
            //fills the neighboring array
            CA->params->aParam[contParam] = parAC;
            contParam++;
            //se ja leu os parametros que tinha q ler, sai do laço
            if(contParam>=CA->params->nRegions){
                break;
            }
        }
        //
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read gravityX", __FILE__, __LINE__);
        CA->params->gravityX = atof(line.c_str());
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read gravityY", __FILE__, __LINE__);
        CA->params->gravityY = atof(line.c_str());
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read gravityZ", __FILE__, __LINE__);
        CA->params->gravityZ = atof(line.c_str());
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read printOutput", __FILE__, __LINE__);
        CA->params->printOutput = atoi(line.c_str());
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read numThreads", __FILE__, __LINE__);
        CA->params->numThreads = atoi(line.c_str());
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read simulationTime", __FILE__, __LINE__);
        CA->params->simulationTime = atof(line.c_str());
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read inputFolder", __FILE__, __LINE__);
        strcpy(CA->params->inputFolder, line.c_str());
        if(!getline(myfile, line))  throw MyException("Param File: Failed to read outputFolder", __FILE__, __LINE__);
        strcpy(CA->params->outputFolder, line.c_str());
        myfile.close();
    }else{
        stringstream ss;
        ss<<"File not found: "<<fName<<endl;
        string str = ss.str();
        throw MyException(str, __FILE__, __LINE__);

    }
}

/**
 * 
 * @param boundFile
 * @param CA
 */
void readBoundaryFile(string boundFile, typ_ca *CA){
    int boundSize=0;
    //temp string to read lines
    string line;
    ifstream myfile2(boundFile.c_str());
    if (myfile2.is_open()){
        //first line contains the file's number of lines
        if(!getline(myfile2, line))  throw MyException("Boundary File: failed to read line", __FILE__, __LINE__);
        boundSize = atoi(line.c_str());
        int cont = 0;
        while ( myfile2.good() )
        {
            if(!getline(myfile2, line)){
                stringstream ss;
                ss<<boundFile<<endl;
                ss<<line<<endl;
                ss <<"boundary File: failed to read line: "<<"["<<cont<<"]: "<<line;
                string str = ss.str();
                throw MyException(str, __FILE__, __LINE__);
                return;
            }
            stringstream strstream(line);
            string token;
            int countColumns=0;
            int iPnt, iAxis;
            double preMov;
            while (getline(strstream, token, ' ')) {
                switch(countColumns){
                    case 0:
                        iPnt = atoi(token.c_str());
                    break;
                    case 1:
                        iAxis = atoi(token.c_str());
                    break;
                    case 2:
                        preMov = atof(token.c_str());
                    break;
                    default:
                        cout<<"Invalid boundary file: "<<countColumns<<" columns!"<<endl;
                }
                countColumns++;
            }
            if(iAxis==0){
                CA->pnts_new[iPnt]->xRestr = CA->pnts_old[iPnt]->xRestr  = true;
                CA->pnts_new[iPnt]->xPreMov= CA->pnts_old[iPnt]->xPreMov = preMov;
            }else if(iAxis==1){
                CA->pnts_new[iPnt]->yRestr = CA->pnts_old[iPnt]->yRestr  = true;
                CA->pnts_new[iPnt]->yPreMov= CA->pnts_old[iPnt]->yPreMov = preMov;
            }else if(iAxis==2){
                CA->pnts_new[iPnt]->zRestr = CA->pnts_old[iPnt]->zRestr  = true;
                CA->pnts_new[iPnt]->zPreMov= CA->pnts_old[iPnt]->zPreMov = preMov;
            }
            if(cont>=boundSize){
                stringstream ss;
                ss <<"There are more boundaray lines than allowed: "<<boundSize;
                string str = ss.str();
                //throw MyException(str, __FILE__, __LINE__);
                cout<<str<<endl;
                break;
            }
            cont++;
        }
        myfile2.close();
    }else{
        stringstream ss;
        ss <<"Unable to open file "<<boundFile;
        cout<<ss.str()<<". Running with no boundary conditions"<<endl;
    }
}


/**
 * 
 * @param strFilePress
 * @param CA
 */
void readPressFile(string strFilePress, typ_ca *CA){
    CA->params->numFaces=0;
    //temp string to read lines
    string line;
    ifstream myfile2(strFilePress.c_str());

    if (myfile2.is_open()){
        //first line contains the file's number of lines
        if(!getline(myfile2, line)){
            stringstream ss;
            ss <<"PRess File: failed to read line - "<<strFilePress.c_str();
            string str = ss.str();
            throw MyException(str, __FILE__, __LINE__);                
        }
        CA->params->numFaces = atoi(line.c_str());
        if(CA->params->numFaces<=0){
            stringstream ss;
            ss <<"Invalid number of faces: "<<CA->params->numFaces;
            string str = ss.str();
            throw MyException(str, __FILE__, __LINE__);
            return;
        }
        getline(myfile2, line);
        CA->params->pressure = atof(line.c_str());        
        CA->params->aFaces = (typ_face**)malloc(sizeof(typ_face*)*CA->params->numFaces);
        int cont = 0;
        while ( myfile2.good() )
        {
            if(!getline(myfile2, line)){
                stringstream ss;
                ss <<"Press File: failed to read line: "<<line;
                string str = ss.str();
                throw MyException(str, __FILE__, __LINE__);
                return;
            }
            stringstream strstream(line);
            string token;
            int countColumns=0;
            typ_face *face = (typ_face*)malloc(sizeof(typ_face));
            while (getline(strstream, token, ' ')) {
                switch(countColumns){
                    case 0:
                        face->pt1 = atoi(token.c_str());
                        CA->pnts_new[face->pt1]->presFcs++;
                        CA->pnts_old[face->pt1]->presFcs++;
                    break;
                    case 1:
                        face->pt2 = atoi(token.c_str());
                        CA->pnts_new[face->pt2]->presFcs++;
                        CA->pnts_old[face->pt2]->presFcs++;
                    break;
                    case 2:
                        face->pt3 = atoi(token.c_str());
                        CA->pnts_new[face->pt3]->presFcs++;
                        CA->pnts_old[face->pt3]->presFcs++;
                    break;
                    default:
                        cout<<"Invalid press file: "<<countColumns<<" columns!"<<endl;
                }
                countColumns++;
            }
                                
            CA->params->aFaces[cont] = face;
            if(cont>=CA->params->numFaces){
                stringstream ss;
                ss <<"There are more pres lines than allowed: "<<CA->params->numFaces;
                string str = ss.str();
                cout<<str<<endl;
                //throw MyException(str, __FILE__, __LINE__);
                break;
            }
            cont++;
        }
        myfile2.close();
    }else{
        cout <<"Unable to open file "<<strFilePress<<endl<<" running without pressures files"<<endl;
        CA->params->numFaces=0;
    }
}
