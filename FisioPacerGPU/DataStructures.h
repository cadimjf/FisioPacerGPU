/* 
 * File:   Globals.h
 * Author: Ricardo
 *
 * Created on 27 de Setembro de 2013, 06:28
 */

#ifndef GLOBALS_H
#define	GLOBALS_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
using namespace std;
//Contains the points ids that define a tetrahedron face on the surface, where is applied external pressure
typedef struct sSurfacePressure {
    int pt1;
    int pt2;
    int pt3;
} typ_face;
//stimulus
typedef struct str_sim{
    double iniTime;
    double period;
    double iniX;
    double endX;
    double iniY;
    double endY;
    double iniZ;
    double endZ;
} t_stim;
//cellular automata parameters
typedef struct str_par_ac{
    //
    double VV0;
    double VV1;
    double VV2;
    double VV3;
    double VV4;
    //
    double VF0;
    double VF1;
    double VF2;
    double VF3;
    double VF4;
    //AP time subdivisions 
    double APTime1ini;
    double APTime2ini;
    double APTime3ini;
    double APTime4ini;
    //Force time subdivisions
    double FTime1;
    double FTime2;
    double FTime3;
    double FTime4;
    ///velocities m/s
    double v[3];

    //kPa
    double forceMultiplier;
    ///N/m or Kg/s^2
    double EAxl[3];
    // preserving volume coefficient
    double EVol[3];
    //    
    double kDamp;
    double EAng[3];
    //parameter to the eletrotonic interaction
    double pEletroTonic;
} t_par_ac;
/*   
*
 */
typedef struct sParam {
    int mecSim;
    int paSim;
    double dtIni;
    double dt;
    double dtSave;
    double dt_sq;
    int pointsNum;
    int elementsNum;
    double gravityX;
    double gravityY;
    double gravityZ;
    int printOutput;
    //stimulus regions
    int stimSize;
    t_stim **aStim; //
    int nRegions;
    t_par_ac **aParam;//
    //pressure faces
    int numFaces;
    typ_face **aFaces;//
    double pressure;
    int numThreads;
    double simulationTime;
    char inputFolder[260];
    char outputFolder[260];
} typ_param;


/*
 *
 */
typedef struct str_point{
    //positions x, y and z
    double x;
    double y;
    double z;
    double x0;
    double y0;
    double z0;
    //mass
    double mass;
    //is it is fixed or not
    bool xRestr;
    bool yRestr;
    bool zRestr;
    //prescribed movement
    double xPreMov;
    double yPreMov;
    double zPreMov;
    //velocity
    double xV;
    double yV;
    double zV;
    int presFcs;
}typ_point;


typedef struct str_time_element{
     //cell time 
    double cellT;    
    //electro-tonic times
    double ECTNC_force_t_ini;
    double ECTNC_force_val_ini;
    double ECTNC_ap_t_ini;
    double ECTNC_ap_val_ini;
    double ECTNC_ap_t_end;
    double ECTNC_force_t_end;
    //states of CA
    int V_state;
    int F_state;    
    //axis
    double fiberDir[3];
    double sheetDir[3];
    double nsheetDir[3];
    //barycenter of the element
    double bary[3];
    //cells volume
    double volCel;    
    double intPts[3][6];    
    double axisLengthT[3];
    //angulos no tempo t entre eixos
    double alphaT_12;
    double alphaT_13;
    double alphaT_23;
    
    double deltaAlpha_12;
    double deltaAlpha_13;
    double deltaAlpha_23;
    //areas do hexaedro virtual
    double areaHexFST;
    double areaHexFNT;
    double areaHexSNT;
    //AP time subdivisions 
    double APTime1;
    double APTime2;
    double APTime3;
    double APTime4;
}typ_dt_element;

typedef struct str_inielement{
    //cell condition: healthy, ischemic, dead, pacemaker
    int cellCond;
    int pmRegion;  
    double volCel_ini;
    double axisLength0[3];
    //direcoes iniciais das fibras
    double fiberDirIni[3];
    double sheetDirIni[3];
    double nsheetDirIni[3];
    //angulos iniciais entre eixos
    double alpha0_12;
    double alpha0_13;
    double alpha0_23;
    //indices dos pontos
    int iPt1;
    int iPt2;
    int iPt3;
    int iPt4;
    //regiao
    int iRegion;
    //indexes of areas containing the instersection points
    int hasPressure;
    int iFaces[6];
    int pressFaces[4];
    int ipressFaces[4];
    double ck[4][6];
    //barycenter of the element
    double bary_ini[3];
    //amortecimento calculado
    //E transformed into K
    double KAxl[3];
    double KVol[3];
    double KAng[3];
    double damp[3];
    double intercMass[6];  
    double areaHexFSTIni;
    double areaHexFNTIni;
    double areaHexSNTIni;
            
    
}typ_t0_element;

//neighbor list
struct lst_item{
    int value;
    struct lst_item* next;
};

//pressure ap
typedef struct str_press{
    double vals[4];
    double times[4];    
    int state;
    double time;
}typ_press;


typedef struct str_stats{
    double max[3];
    double min[3];
    double avgVel;
    double maxVol;
    double minVol;
    double volIni;
    double volMaxDelta;
    double maxDeltaVol;
    int contSave;
    double tol;
    double err;
} typ_stats;

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
    double volume;
    double time;
    double timeSaving;
    typ_stats *stats;//
    typ_press *pressureCA;
}typ_ca;// automato cellular
struct MyException : public std::exception
{
    std::string msg;
    std::string fileName;
    int lineNumber;
    MyException(std::string ss, std::string fn, int l){
        this->msg          = ss;
        this->fileName     = fn;
        this->lineNumber   = l;
    }
    ~MyException() {} // Updated
    const char* what() const throw() { 
        stringstream ss;
        ss<<"Error in file "<<this->fileName<< " at line "<<this->lineNumber<<":"<<endl;
        ss<<"\t"<<this->msg;
        string str = ss.str();
        return str.c_str(); 
    }
    std::string getMessage() { 
        stringstream ss;
        ss<<"Error in file "<<this->fileName<< " at line "<<this->lineNumber<<":"<<endl;
        ss<<"\t"<<this->msg;
        return ss.str();
        
    }
};
/**
 * 
 * @param length
 * @return 
 */
template< typename T > T* myAllocation( int length )
{
   T* al = (T*) malloc(length * sizeof(T));
   if(al==NULL) throw ("Allocation Failure", __FILE__, __LINE__);
   else  return al;
}

template< typename T > T** myArrayAllocation( int length )
{
   T** array = (T**) malloc(length*sizeof(T*));
   if(array==NULL){
       throw MyException("Allocation Failure", __FILE__, __LINE__);
   }else{
        for(int i=0;i<length;i++){
            array[i] =  (T*) malloc(sizeof(T));
        }
        return array;
   }
}

template< typename T > void myArrayDeallocation(T** array, int length )
{
    for(int i=0;i<length;i++){
        if(array[i]!=NULL) free(array[i]);
    }
    if(array!=NULL) free(array);
}


lst_item** iniAList(int size);
void freeAList(lst_item** list, int size);
void free_list(lst_item *head);

lst_item* createListItem(int value);
void addItem(lst_item** head, int value);
void printList(lst_item* list);

#endif	/* GLOBALS_H */
