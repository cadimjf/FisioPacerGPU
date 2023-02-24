/* 
 * File:   ReadMesh.h
 * Author: ricardo
 *
 * Created on April 23, 2013, 1:34 PM
 */


#ifndef READMESH_H
#define	READMESH_H
#include <exception>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Constants.h"
#include "MyStdLib.h"
#include "DataStructures.h"
#include "LinearAlg.h"
using namespace std;

void omegaBAdjacency(int iElement, typ_ca *CA);
void fillFiberFile(string fileName, typ_ca *CA);
void pointsFile(string fileName, typ_ca *CA);
void elementsFile(string fileName, typ_ca *CA);

void OmegaANeighborhood();
int readSize(string file);
void readPressFile(string strFilePress, typ_ca *CA);
void omegaBAdjacency(int iElement, typ_ca *CA);
void fillFiberFile(string fileName, typ_ca *CA);
void pointsFile(string fileName, typ_ca *CA);
void elementsFile(string fileName, typ_ca *CA);
void openFile(typ_ca *CA,
        string strFilePts, string strFileEle, string strFileFib, string boundFile,
        string strPressFile, string stimFile);
void OmegaANeighborhood(typ_ca *CA);
void readStimFile(string strFileSt, typ_ca * CA);
void readParameterFile(string fName, typ_ca *CA);
void readBoundaryFile(string boundFile, typ_ca *CA);
void readPressFile(string strFilePress, typ_ca *CA);
#endif	/* READMESH_H */

