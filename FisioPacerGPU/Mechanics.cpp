#include <stdio.h>
#include "Mechanics.h"
#include "MyStdLib.h"

/**
 *
 * @param k
 * @param lengthT
 * @param length0
 * @return
 */
double getForce_HookesLawAxis(double k, double lengthT, double length0) {
    double deltaL = lengthT - length0;
    return -k * (deltaL);
    //    double deltaL = lengthT/length0;
    //    return -k*log10(deltaL);
}

/**
 *
 * @param vel1
 * @param vel2
 * @param kdamp
 * @param axis
 * @param axisLenT
 * @return
 */
double getForce_Damping( double relativeVel[3], double kdamp, double axis[3]) {
    return -kdamp * dot(relativeVel, axis);
}

/**
 *
 * @param k
 * @param deltaA
 * @return
 */
void getForce_HookesLawAngular(double force[3], double deltaA, double kAng, double axis[3]) {
    double mod = -deltaA*kAng;
    force[0] = mod*axis[0];
    force[1] = mod*axis[1];
    force[2] = mod*axis[2]; 
}

/**
 *
 * @param volT
 * @param vol0
 * @param kvol
 */
double getForce_VolPreserving(double volT, double vol0, double kvol) {
    //return -kvol * (volT - vol0) / vol0;
    return -kvol * (volT - vol0);
}
/**
 * 
 * @param forceAxis
 * @param f
 * @param col
 */
void pressIntPoint(double forceAxis[3][2], double fSurf[3], int col){
    for(int k=0;k<3;k++){
        forceAxis[k][col]+=fSurf[k];
    }
}

/**
 *
 * @param POINTS_OLD
 * @param elements_old
 * @param i
 * @param p
 * @param forcesOnPts
 */
int MecStep_i(typ_ca *CA, int i, double *forcesOnPts)
{

    double forceSheet[3][2], forceNSheet[3][2], forceFiber[3][2];
    for (int ind = 0; ind < 3; ind++) {
        forceSheet[ind][0] = forceSheet[ind][1] = 0.0f;
        forceNSheet[ind][0] = forceNSheet[ind][1] = 0.0f;
        forceFiber[ind][0] = forceFiber[ind][1] = 0.0f;
    }
    //find the velocities of the intersection points
    double velsIter[6][3];
    findVelocities(CA, i, velsIter);
    //find the active force on the fiber
    if(CA->params->paSim==1){
        getActiveForce( i, forceFiber, CA);
    }
    //find the passive force on the fiber direction
    getPassiveForce(i, FIBER, forceFiber, velsIter[0], velsIter[1], CA);
    //find the passive force on the sheet direction
    getPassiveForce(i, SHEET, forceSheet, velsIter[2], velsIter[3], CA);
    //find the passive force on the Nsheet direction
    getPassiveForce(i, NORMAL, forceNSheet, velsIter[4], velsIter[5], CA);
    getExternalForce(forceFiber, i, CA);
    getExternalForce(forceSheet, i, CA);
    getExternalForce(forceNSheet, i, CA);
    
    //finds the force on the faces intercepted
    double forcePts[3][4];
    for (int j = 0; j < 3; j++)
        for (int k = 0; k < 4; k++)
            forcePts[j][k] = 0.0f;
    

    //distribui aforça da intersecao para os pontos do elemento
    for (int k = 0; k < 4; k++) {
        for (int j = 0; j < 3; j++) {
            //força na fibra interseçao 1
            forcePts[j][k] += CA->ini[i]->ck[k][0] * forceFiber[j][0];
            //força na fibra interseção 2
            forcePts[j][k] += CA->ini[i]->ck[k][1] * forceFiber[j][1];
            //força na sheet interseçao 3
            forcePts[j][k] += CA->ini[i]->ck[k][2] * forceSheet[j][0];
            //força na sheet interseção 4
            forcePts[j][k] += CA->ini[i]->ck[k][3] * forceSheet[j][1];
            //força na Nsheet interseçao 5
            forcePts[j][k] += CA->ini[i]->ck[k][4] * forceNSheet[j][0];
            //força na Nsheet interseção 6
            forcePts[j][k] += CA->ini[i]->ck[k][5] * forceNSheet[j][1];

        }
    }
    //adds the force computed above to the total force of the point.
    //This is necessary to take into account the force of all elements in a single point
    for (int jj = 0; jj < 3; jj++) {
        forcesOnPts[I2d(CA->ini[i]->iPt1, jj, 3)] += forcePts[jj][0];
        forcesOnPts[I2d(CA->ini[i]->iPt2, jj, 3)] += forcePts[jj][1];
        forcesOnPts[I2d(CA->ini[i]->iPt3, jj, 3)] += forcePts[jj][2];
        forcesOnPts[I2d(CA->ini[i]->iPt4, jj, 3)] += forcePts[jj][3];

        if (isnan(forcePts[jj][0]) || isnan(forcePts[jj][1]) || isnan(forcePts[jj][2]) || isnan(forcePts[jj][3]) ||
                isinf(forcePts[jj][0]) || isinf(forcePts[jj][1]) || isinf(forcePts[jj][2]) || isinf(forcePts[jj][3])) {
            if(CA->params->printOutput==1){
                printf("[%d] %d %f %f %f %f\n", i, jj, forcePts[jj][0], forcePts[jj][1], forcePts[jj][2], forcePts[jj][3]);
            }
            return -1;
        }
    }

    return 0;
}
/**
 * 
 * @param iElem
 * @param iAxis
 * @param CA
 * @param fAng1
 * @param fAnf2
 */
void getAngularForces(int iElem, int iAxis, typ_ca *CA, double fAng1[3], double fAng2[3], double axis[3]){
    double deltaAlpha12 = CA->t_old[iElem]->deltaAlpha_12;
    double deltaAlpha13 = CA->t_old[iElem]->deltaAlpha_13;
    double deltaAlpha23 = CA->t_old[iElem]->deltaAlpha_23;
    double k0 = CA->ini[iElem]->KAng[0];
    double k1 = CA->ini[iElem]->KAng[1];
    double k2 = CA->ini[iElem]->KAng[2];
    if(iAxis==FIBER){
        for(int i=0;i<3;i++){
            axis[i]  = CA->t_old[iElem]->fiberDir[i];
        }
        getForce_HookesLawAngular(fAng1, deltaAlpha12, k0, CA->t_old[iElem]->sheetDir);
        getForce_HookesLawAngular(fAng2, deltaAlpha13, k1, CA->t_old[iElem]->nsheetDir);
    }else if(iAxis==SHEET){
        for(int i=0;i<3;i++){
            axis[i]  = CA->t_old[iElem]->sheetDir[i];
        }
        getForce_HookesLawAngular(fAng1, deltaAlpha12, k0, CA->t_old[iElem]->fiberDir);
        getForce_HookesLawAngular(fAng2, deltaAlpha23, k2, CA->t_old[iElem]->nsheetDir);
    }if(iAxis==NORMAL){
        for(int i=0;i<3;i++){
            axis[i]  = CA->t_old[iElem]->nsheetDir[i];
        }
        getForce_HookesLawAngular(fAng1, deltaAlpha13, k1, CA->t_old[iElem]->fiberDir);
        getForce_HookesLawAngular(fAng2, deltaAlpha23, k2, CA->t_old[iElem]->sheetDir);
    }
}

/**
 * 
 * @param iElem
 * @param iAxis
 * @param force
 * @param vel1
 * @param vel2
 * @param CA
 */
void getPassiveForce(
        int iElem, int iAxis,
        double force[3][2], double vel1[3], double vel2[3], typ_ca *CA)
{
    double axis[3]={0.0};
    double fAng1[3]={0.0}, fAng2[3]={0.0};
    getAngularForces(iElem, iAxis, CA, fAng1, fAng2, axis);    
    //find the hooke's law resulting force
    double k=1.0;
    double FhookeAxis = getForce_HookesLawAxis(CA->ini[iElem]->KAxl[iAxis]*k, CA->t_old[iElem]->axisLengthT[iAxis], CA->ini[iElem]->axisLength0[iAxis]);
    //volume force
    double Fvol = getForce_VolPreserving(CA->t_old[iElem]->volCel, CA->ini[iElem]->volCel_ini, CA->ini[iElem]->KVol[iAxis]);
    double forceJ = 0.0;
    double modForca = FhookeAxis + Fvol;
    for (int j = 0; j < 3; j++) {        
        forceJ = axis[j]*(modForca) + fAng1[j] + fAng2[j];
        // - relativeVel[j]*ap->DAng[iAxis];
        /*area*(axis2[j]*(FhookeA1+FdampShear1) + axis3[j]*(FhookeA2+FdampShear2));*/
        force[j][0] += -(forceJ) - CA->ini[iElem]->damp[j]*vel1[j] ;
        force[j][1] +=  (forceJ) - CA->ini[iElem]->damp[j]*vel2[j] ;
    }
}

/**
 * @param force
 * @param i
 */
void getExternalForce(double force[3][2], int i, typ_ca *CA)
{
    // 9,80665 m/s2
    double mass = (ro_mass_dens * CA->ini[i]->volCel_ini);

    force[0][0] += CA->params->gravityX*mass;
    force[0][1] += CA->params->gravityX*mass;

    force[1][0] += CA->params->gravityY*mass;
    force[1][1] += CA->params->gravityY*mass;

    force[2][0] += CA->params->gravityZ*mass;
    force[2][1] += CA->params->gravityZ*mass;

}
/**
 * 
 * @param i
 * @param force
 * @param CA
 */
void getActiveForce( int i, double force[3][2], typ_ca *CA)
{
    double t = getActiveTension(i, CA);
    //double mod = t*CA->t_old[i]->areaHexSNT;
    double mod = t*CA->ini[i]->areaHexSNTIni;
    for (int j = 0; j < 3; j++) {
        force[j][0] += +mod * CA->t_old[i]->fiberDir[j];
        force[j][1] += -mod * CA->t_old[i]->fiberDir[j];
    }
}
/**
 *
 * @param nThreads
 * @param POINTS_OLD
 * @param forcesOnPts
 * @param time
 * @return
 */
void simulationStep( int nThreads,
                typ_ca* CA,
                double *forcesOnPts)
{
    double volT=0.0;
    int contCA=0;
    #pragma omp parallel for schedule(static) num_threads(nThreads) reduction(+:volT)
    for(int i=0;i<CA->params->elementsNum;i++){
        if(CA->params->paSim==1){
            CAStep_i(i, CA);
            contCA++;
            computeNewAPDElectroTonic(i, CA);
        }
        computeForceOnElement(CA, forcesOnPts, i);
        //verifica se o volume da celula é menor que 1% do inicial - netste caso, mata o processo pq é sinal de erro.
        if(CA->t_new[i]->volCel< 0.001 * CA->ini[i]->volCel_ini){
            CA->params->mecSim=0;
            CA->params->paSim=0;
            if(CA->params->printOutput==1){
                cout<<"Mata por volume pequeno["<<i<<"]: "<<CA->time<<endl;
                cout<<"Inicial: "<<CA->ini[i]->volCel_ini<<" | Atual: "<<CA->t_new[i]->volCel<<endl;
                cout<<  CA->pnts_old[CA->ini[i]->iPt1]->x<<" "<<
                        CA->pnts_old[CA->ini[i]->iPt1]->y<<" "<<
                        CA->pnts_old[CA->ini[i]->iPt1]->z<<endl;
                cout<<  CA->pnts_old[CA->ini[i]->iPt2]->x<<" "<<
                        CA->pnts_old[CA->ini[i]->iPt2]->y<<" "<<
                        CA->pnts_old[CA->ini[i]->iPt2]->z<<endl;

                throw MyException("1percent volume.", __FILE__, __LINE__);

            }
        }
        volT+= CA->t_new[i]->volCel;      
    }//fim pragma
    //exit(0);
    if(CA->params->mecSim==1){
       computePressurePoints(CA, forcesOnPts);
    }
    CA->volume=volT;
}
/**
 *
 * @param POINTS_OLD
 * @param forcesOnPts
 * @param i
 */
void computeForceOnElement(typ_ca *CA, double *forcesOnPts, int i)
{
    if(CA->params->mecSim==1){
        int retVal = MecStep_i(CA, i, forcesOnPts);
        if(retVal==-1){
            //The method didn't converge, so it changes the flags in order to no long compute the simulation.
            //This is done due the prohibition to return nor break loops inside omp clauses.
            CA->params->mecSim=0;
            CA->params->paSim=0;
        }else{
            //if the mechanical simulation is ok, update geometry
            updateGeometry(i, CA);
        }
    }
}

/**
 *
 * @param p
 * @param forcesOnPts
 */
void computeForceIntermediaria(
                typ_ca *CA,
                double *forcesOnPts){
    for(int i=0;i<CA->params->elementsNum;i++)
    {
        computeForceOnElement(CA, forcesOnPts, i);
    }
    if(CA->params->mecSim==1){
        computePressurePoints(CA, forcesOnPts);
    }
}


/*
* para cada area marcada no contorno
*      obtem a normal da area
*      obtem a area
*      verifica se a normal aponta para dentro ou para fora da malha
*      obtem a força = tensao * area
*      aplica em cada ponto
*      forcesOnPts[I2d(k,0,3)] += 1/3 * força
*
*/
void getForceOnSurface(typ_face *face, typ_ca *CA, double *aForce)
{
    double force=0.0;
    double normal[3]={0.0}, bary[3]={0.0};
    
    double area = getFaceAreaNormal(face, CA, normal, bary);
    force = -/*getPressurePercent(CA) * TODO FIXME retirar o comeentario para o modelo de pressao*/CA->params->pressure* area/3.0;
    //
    aForce[0] = force*normal[0];
    aForce[1] = force*normal[1];
    aForce[2] = force*normal[2];    
}

/**
 * 
 * @param face
 * @param CA
 * @param normal
 * @param bary
 * @return 
 */
double getFaceAreaNormal(typ_face *face, typ_ca *CA, double normal[3], double bary[3]){
    double vPt1[3], vPt2[3], vPt3[3];
    getPtsInVector(CA->pnts_old, face->pt1, vPt1);
    getPtsInVector(CA->pnts_old, face->pt2, vPt2);
    getPtsInVector(CA->pnts_old, face->pt3, vPt3);
    findFaceNormal(CA, face->pt1, face->pt2, face->pt3, normal);
    fixNormal(vPt1, vPt2, vPt3, normal, bary);
    //getFaceCenter(vPt1, vPt2, vPt3, bary);
    return areaTriangulo(vPt1, vPt2, vPt3);
}

/**
 *
 * @param POINTS_OLD
 * @param forcesOnPts
 */
void computePressurePoints(typ_ca *CA, double *forcesOnPts){
    typ_face *face;    
    double surfaceForces[3] = {0.0,0.0,0.0};
    double max=DBL_MIN, min=DBL_MAX, sum=0.0;
    for(int iF=0; iF < CA->params->numFaces; iF++){
        face = CA->params->aFaces[iF];
        getForceOnSurface(face, CA, surfaceForces);
        applyPressureForcePoint(CA, surfaceForces, forcesOnPts, face->pt1);
        applyPressureForcePoint(CA, surfaceForces, forcesOnPts, face->pt2);
        applyPressureForcePoint(CA, surfaceForces, forcesOnPts, face->pt3);
        double normF = my_norm(surfaceForces);
        if(normF>max) max=normF;
        if(normF<min) min=normF;
        sum+=normF;
    }   
    //exit(0);*
    //cout<<min<<"::"<<max<<"::"<<sum<<endl;
}
/**
 * 
 * @param CA
 * @param surfaceForces
 * @param forcesOnPts
 * @param iPnt
 */
void applyPressureForcePoint(typ_ca *CA,double surfaceForces[3], double *forcesOnPts, int iPnt){
    
    if (!CA->pnts_old[iPnt]->xRestr) {
        forcesOnPts[I2d(iPnt, 0, 3)] += surfaceForces[0];
    }
    //
    if (!CA->pnts_old[iPnt]->yRestr) {
        forcesOnPts[I2d(iPnt, 1, 3)] += surfaceForces[1];
    }
    //
    if (!CA->pnts_old[iPnt]->zRestr) {
        forcesOnPts[I2d(iPnt, 2, 3)] += surfaceForces[2];
    }

}


