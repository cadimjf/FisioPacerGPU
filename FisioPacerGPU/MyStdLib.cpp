#include "MyStdLib.h"
double my_max(double a, double b){
    if(a>b)return a;
    else return b;
}
double my_min(double a, double b){
    if(a<b)return a;
    else return b;
}


/**
 * 
 * @param imin
 * @param imax
 * @return 
 */
int midpoint(int imin, int imax){
     return (int)(imin + ((imax - imin) / 2));
}

/**
 * 
 * @param A
 * @param key
 * @param imin
 * @param imax
 * @return 
 */
int binary_search(int A[], int key, int imin, int imax)
{
    // continue searching while [imin,imax] is not empty
    while (imax >= imin)
    {
        /* calculate the midpoint for roughly equal partition */
        int imid = midpoint(imin, imax);
        // determine which subarray to search
        if (A[imid] < key)
            // change min index to search upper subarray
            imin = imid + 1;
        else if (A[imid] > key)
            // change max index to search lower subarray
            imax = imid - 1;
        else
            // key found at index imid
            return imid;
    }
    // key not found
    return KEY_NOT_FOUND;
}

void printParameters(typ_ca *CA){
    cout<<CA->params->mecSim<<endl;
    cout<<CA->params->paSim<<endl;
    cout<<CA->params->dt<<endl;
    cout<<CA->params->dtSave<<endl;
    cout<<CA->params->nRegions<<endl;
    for(int i=0;i<CA->params->nRegions;i++){
        t_par_ac* ap = CA->params->aParam[i];
        cout<<ap->VV0<<" ";
        cout<<ap->VV1<<" ";
        cout<<ap->VV2<<" ";
        cout<<ap->VV3<<" ";
        cout<<ap->VV4<<" ";
        cout<<ap->VF0<<" ";
        cout<<ap->VF1<<" ";
        cout<<ap->VF2<<" ";
        cout<<ap->VF3<<" ";
        cout<<ap->VF4<<" ";
        cout<<ap->APTime1ini<<" ";
        cout<<ap->APTime2ini<<" ";
        cout<<ap->APTime3ini<<" ";
        cout<<ap->APTime4ini<<" ";
        cout<<ap->FTime1<<" ";
        cout<<ap->FTime2<<" ";
        cout<<ap->FTime3<<" ";
        cout<<ap->FTime4<<" ";
        cout<<ap->v[FIBER]<<" ";
        cout<<ap->v[SHEET]<<" ";
        cout<<ap->v[NORMAL]<<" ";
        cout<<ap->forceMultiplier<<" ";
        cout<<ap->EAxl[FIBER]<<" ";
        cout<<ap->EAxl[SHEET]<<" ";
        cout<<ap->EAxl[NORMAL]<<" ";
        cout<<ap->EVol[FIBER]<<" ";
        cout<<ap->EVol[SHEET]<<" ";
        cout<<ap->EVol[NORMAL]<<" ";
        cout<<ap->kDamp<<" ";
        cout<<ap->EAng[FIBER]<<" ";
        cout<<ap->EAng[SHEET]<<" ";
        cout<<ap->EAng[NORMAL]<<" ";
        cout<<ap->pEletroTonic<<endl;
    }
    cout<<CA->params->gravityX<<endl;
    cout<<CA->params->gravityY<<endl;
    cout<<CA->params->gravityZ<<endl;
    cout<<CA->params->printOutput<<endl;
    cout<<CA->params->numThreads<<endl;
    cout<<CA->params->simulationTime<<endl;
    cout<<CA->params->inputFolder<<endl;
    cout<<CA->params->outputFolder<<endl;
}

/**
 * 
 */
void setDefaultFolders(){
    //DEFAULT_FOLDERIN  = "/home/ricardo/Dropbox/Doutorado/data/teste2ele/";
    //DEFAULT_FOLDERIN  = "/home/ricardo/Dropbox/Doutorado/data/1tetraedro/";
    //DEFAULT_FOLDERIN  = "/home/ricardo/Dropbox/Doutorado/data/cubit/";
    string defaultFolderIn  = "/home/ricardo/Dropbox/Doutorado/Tese/Validacao/cubo6/FibraX";
    //DEFAULT_FOLDEROUT = "C:\\Users\\Ricardo\\output_malhaCA\\";
    string defaultFolderOut  ="/home/ricardo/outputFisioPacer/";

    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    time (&rawtime);
    timeinfo = localtime (&rawtime);

    strftime (buffer,24,"%Y_%b_%d_%H_%M_%S",timeinfo);

    defaultFolderOut.append(buffer);
    defaultFolderOut.append("/");
    
    cout<<"dir: " <<defaultFolderOut<<endl;
    
}



/**
 * 
 * @param from
 * @param to
 */
void copyPoint(typ_point *from, typ_point *to){
    to->x        = from->x;
    to->y        = from->y;
    to->z        = from->z;
    to->mass     = from->mass;
    to->xRestr   = from->xRestr;
    to->yRestr   = from->yRestr;
    to->zRestr   = from->zRestr;
    to->xV       = from->xV;
    to->yV       = from->yV;
    to->zV       = from->zV;
    
    
}
