#include <cfloat>
#include <float.h>

#include "CellularAutomaton.h"
/**
 *
 * @param el
 * @param i
 * @param iRegion
 * @return
 */
double getActiveTension(int i, typ_ca *CA)
{
    int iRegion = CA->ini[i]->iRegion;
    t_par_ac* aP = CA->params->aParam[iRegion];
    double f = getActiveTensionNormalized(i, CA);
    return f * aP->forceMultiplier;
}

double getActiveTensionNormalized(int i, typ_ca *CA){

    int iRegion = CA->ini[i]->iRegion;
    t_par_ac* aP = CA->params->aParam[iRegion];
    if(CA->ini[i]->cellCond==DEAD)return 0.0;

    if(CA->params->paSim==0 ) //Se o automato estÃƒÂ¡ desligado, retorna sempre o parametro do modulo da forca
       return 0.0;//aP->forceMultiplier;

    double f=0.0;
    double m;
    double tcell = CA->t_old[i]->cellT;
    switch(CA->t_old[i]->F_state){
        case F0:
            f = aP->VF0;
            break;
        case F1:
            f = aP->VF1;
            break;
        case F2:
            m = (aP->VF3 - aP->VF2) / (aP->FTime2 - aP->FTime1);
            f= aP->VF2 + m*(tcell-aP->FTime1);
            break;
        case F3:
            m = (aP->VF4 - aP->VF3) / (aP->FTime3 - aP->FTime2);
            f = aP->VF3 + m*(tcell - aP->FTime2);
            break;
        case F4:
            m = (aP->VF0 - aP->VF4) / (aP->FTime4 - aP->FTime3);
            f =  aP->VF4 + m*(tcell - aP->FTime3);
            break;
        default:
            f= 0.0;
    }
    return f ;
}

/**
 *
 * @param el
 * @param i
 * @return
 */
double getActiveTensionDiscrete( int i, typ_ca *CA)
{
    int iRegion = CA->ini[i]->iRegion;
    t_par_ac* aP = CA->params->aParam[iRegion];
    if(CA->ini[i]->cellCond==DEAD)return 0.0;
    if(CA->params->paSim==0 ) //Se o automato estÃƒÂ¡ desligado, retorna sempre o parametro do modulo da forca
       return aP->forceMultiplier;
    double f=0.0;

    //it is not in FR 1 state, behaves normally
    double tcell = CA->t_old[i]->cellT;
    switch(CA->t_old[i]->F_state){
        case F0:
            f = aP->VF0;
            break;
        case F1:
            f = aP->VF1;
            break;
        case F2:
            f = aP->VF2;
            break;
        case F3:
            f =  aP->VF3;
            break;
        case F4:
            f = aP->VF4;
            break;
        default:
            f= 0.0;
    }
    return f;
}
/**
 *
 * @param CA
 */
void initializeCA( typ_ca *CA )
{
    
    for(int i=0;i<CA->params->elementsNum;i++)
    {
        iniVolumes(i, CA);
    }
    double mass = 0.0;
    //iterate over points to find their masses
    for(int i=0;i<CA->params->pointsNum;i++){
        mass =0.0;
        lst_item *cur = CA->omega_b[i];
        while(cur != NULL){
            mass += CA->ini[cur->value]->volCel_ini;
            cur = cur->next;
        }
        CA->pnts_new[i]->mass = CA->pnts_old[i]->mass = ro_mass_dens*mass/4.0;   
    }
    for(int i=0;i<CA->params->elementsNum;i++)
    {
        CA->t_old[i]->cellT = CA->t_new[i]->cellT       = 0.0f;
        CA->t_old[i]->ECTNC_force_t_ini = CA->t_new[i]->ECTNC_force_t_ini       = 0.0f;
        CA->t_old[i]->ECTNC_force_val_ini  = CA->t_new[i]->ECTNC_force_val_ini        = 0.0f;
        CA->t_old[i]->ECTNC_ap_t_ini = CA->t_new[i]->ECTNC_ap_t_ini       = 0.0f;
        CA->t_old[i]->ECTNC_ap_val_ini  = CA->t_new[i]->ECTNC_ap_val_ini        = 0.0f;
        CA->t_old[i]->ECTNC_ap_t_end = CA->t_new[i]->ECTNC_ap_t_end       = 0.0f;
        CA->t_old[i]->ECTNC_force_t_end = CA->t_new[i]->ECTNC_force_t_end       = 0.0f;

        //initially, all cells are healthy
        CA->ini[i]->cellCond = HEALTHY;
        //V state
        CA->t_new[i]->V_state     = CA->t_old[i]->V_state   = V0;
        //F state
        CA->t_new[i]->F_state     = CA->t_old[i]->F_state   = F0;
        iniGeometry( i, CA);
        /*if(i==0){
            cout << "initializeCA"<<endl;
            cout<<CA->t_old[i]->axisLengthT[0]<<" "<<
                  CA->t_old[i]->axisLengthT[1]<<" "<<
                  CA->t_old[i]->axisLengthT[2]<<endl;
            cout<<"__________"<<endl;
        }*/
        //sets the pacemaker up
        for(int iS=0; iS<CA->params->stimSize; iS++){
            t_stim *s = CA->params->aStim[iS];
            if(
                    CA->t_old[i]->bary[0] > s->iniX && CA->t_old[i]->bary[0] < s->endX &&
                    CA->t_old[i]->bary[1] > s->iniY && CA->t_old[i]->bary[1] < s->endY &&
                    CA->t_old[i]->bary[2] > s->iniZ && CA->t_old[i]->bary[2] < s->endZ
            ){
                CA->ini[i]->cellCond = paceMaker;
                CA->ini[i]->pmRegion = iS;
            }
        }
        restartAPDElectroTonic(i, CA);
        CA->t_old[i]->APTime1 = CA->t_new[i]->APTime1;
        CA->t_old[i]->APTime2 = CA->t_new[i]->APTime2;
        CA->t_old[i]->APTime3 = CA->t_new[i]->APTime3;
        CA->t_old[i]->APTime4 = CA->t_new[i]->APTime4;
        
    }
    iniPressure(CA);
}

/**
 * return V in miliVolts
 */
double getV(int i ,typ_ca *CA)
{
    int iRegion = CA->ini[i]->iRegion;
    t_par_ac* ap = CA->params->aParam[iRegion];

    double m;
    switch(CA->t_old[i]->V_state){
        case V0:
            return ap->VV0;
        case V1:
            m = (ap->VV2 - ap->VV1)/(CA->t_old[i]->APTime1);
            return ap->VV1 + m*(CA->t_old[i]->cellT);
        case V2:
            m = (ap->VV3 - ap->VV2)/(CA->t_old[i]->APTime2 - CA->t_old[i]->APTime1);
            return ap->VV2 + m*(CA->t_old[i]->cellT - CA->t_old[i]->APTime1);
        case V3:
            m = (ap->VV4 - ap->VV3)/(CA->t_old[i]->APTime3 - CA->t_old[i]->APTime2);
            return ap->VV3 + m*(CA->t_old[i]->cellT - CA->t_old[i]->APTime2);
        case V4:
            m = (ap->VV0 - ap->VV4)/(CA->t_old[i]->APTime4 - CA->t_old[i]->APTime3);
            return ap->VV4 + m*(CA->t_old[i]->cellT - CA->t_old[i]->APTime3);
        default:
            return ap->VV0;
    }
}
/**
 *
 * @param cellsStates
 * @param i
 * @param p
 * @return
 */
double getVdiscret( int i,typ_ca *CA)
{
    int iRegion = CA->ini[i]->iRegion;
    t_par_ac* ap = CA->params->aParam[iRegion];
    switch(CA->t_old[i]->V_state){
        case V0:
            return ap->VV0;
        case V1:
            return ap->VV1;
        case V2:
            return ap->VV2;
        case V3:
            return ap->VV3;
        case V4:
            return ap->VV4;
        default:
            return ap->VV0;
    }
}
/**
 *
 * @param t
 * @param dt
 * @return
// */
int isStimulationTime(int i, typ_ca *CA){
    double t = CA->time;
    t_stim *s = CA->params->aStim[ CA->ini[i]->pmRegion];
    return (((t>=s->iniTime)&&(((t-s->iniTime)-(floor(((t-s->iniTime)/s->period))*s->period))<=CA->params->dt)))?1:0;

}
/**
 *
 * @param cellsPosT_old
 * @param iCurrent
 * @param iElemNeighbor
 * @return
 */
double getPropagationTime(int iElemNeighbor, int i, typ_ca *CA)
{
    t_par_ac* ap = CA->params->aParam[CA->ini[iElemNeighbor]->iRegion];

    double dir[3], dirN[3];
    //get the vector between the elements position
    dir[0] = CA->t_old[iElemNeighbor]->bary[0] - CA->t_old[i]->bary[0];
    dir[1] = CA->t_old[iElemNeighbor]->bary[1] - CA->t_old[i]->bary[1];
    dir[2] = CA->t_old[iElemNeighbor]->bary[2] - CA->t_old[i]->bary[2];
    //obtem a distancia entre o elemento e o vizinho
    double s = my_norm(dir);
    dirN[0] = dir[0]/s;
    dirN[1] = dir[1]/s;
    dirN[2] = dir[2]/s;
    //////////////
    double v1 = fabs(dot(dirN, CA->t_old[i]->fiberDir )) * ap->v[FIBER];
    double v2 = fabs(dot(dirN, CA->t_old[i]->sheetDir )) * ap->v[SHEET];
    double v3 = fabs(dot(dirN, CA->t_old[i]->nsheetDir)) * ap->v[NORMAL];
    return s/(v1+v2+v3);
}

/**
 *
 * @param i
 */
void cellActivation(int i, typ_ca *CA)
{
     //if the cell is dead, it doesnt increment any states
    if(CA->ini[i]->cellCond==DEAD)
        return;

    CA->t_new[i]->cellT    = 0.0f;
    CA->t_new[i]->V_state  = V1;
    CA->t_new[i]->F_state  = F1;
    CA->t_new[i]->ECTNC_force_t_ini    = 0.0f;
    CA->t_new[i]->ECTNC_force_val_ini     = 0.0f;
    CA->t_new[i]->ECTNC_ap_t_ini    = 0.0f;
    CA->t_new[i]->ECTNC_ap_val_ini     = 0.0f;
    CA->t_new[i]->ECTNC_ap_t_end    = 0.0f;
    CA->t_new[i]->ECTNC_force_t_end    = 0.0f;
}
/**
 *
 * @param cellsPosT_new
 * @param cellsPosT_old
 * @param cellsStates_new
 * @param cellsStates_old
 * @param dt
 */
void incStates(int i, double dt, typ_ca *CA)
{
    int iRegion = CA->ini[i]->iRegion;
    t_par_ac* ap = CA->params->aParam[iRegion];
    //if the cell is dead, it doesnt increment any states
    if(CA->ini[i]->cellCond==DEAD)
        return;
    if(CA->t_old[i]->V_state==V0 && CA->t_old[i]->F_state==F0){

    } else{
        CA->t_new[i]->cellT = CA->t_old[i]->cellT + dt;
    }
    switch(CA->t_old[i]->V_state){
        case V0:
            //keeps the same value
            CA->t_new[i]->V_state = V0;
            break;
        case V1:
            CA->t_new[i]->V_state = (CA->t_old[i]->cellT >= CA->t_old[i]->APTime1)?V2:V1;
            break;
        case V2:
            CA->t_new[i]->V_state = (CA->t_old[i]->cellT >= CA->t_old[i]->APTime2)?V3:V2;
            break;
        case V3:
            CA->t_new[i]->V_state = (CA->t_old[i]->cellT >= CA->t_old[i]->APTime3)?V4:V3;
            break;
        case V4:
            //CA->t_new[i]->V_state = (CA->t_old[i]->cellT >= CA->t_old[i]->APTime4)?V0:V4;
            
            if(CA->t_old[i]->cellT >= CA->t_old[i]->APTime4){
                CA->t_new[i]->V_state = V0;
                //restartAPDElectroTonic(i, CA);
            }else{
                CA->t_new[i]->V_state = V4;                
            }

            break;
        default:
            //keeps the same state
            CA->t_new[i]->V_state = CA->t_old[i]->V_state;
    }

    switch(CA->t_old[i]->F_state){
        case F0:
            CA->t_new[i]->F_state = F0;
            break;
        case F1:
            CA->t_new[i]->F_state = (CA->t_old[i]->cellT >= ap->FTime1)?F2:F1;
            break;
        case F2:
            CA->t_new[i]->F_state = (CA->t_old[i]->cellT >= ap->FTime2)?F3:F2;
            break;
        case F3:
            CA->t_new[i]->F_state = (CA->t_old[i]->cellT >= ap->FTime3)?F4:F3;
            break;
        case F4:
            if(CA->t_old[i]->cellT >= ap->FTime4){
                CA->t_new[i]->F_state = F0;
                CA->t_new[i]->cellT   = 0.0f;
                CA->t_new[i]->ECTNC_force_t_ini   = 0.0f;
                CA->t_new[i]->ECTNC_force_val_ini    = 0.0f;
                CA->t_new[i]->V_state = V0;
                CA->t_new[i]->ECTNC_ap_t_ini   = 0.0f;
                CA->t_new[i]->ECTNC_ap_val_ini    = 0.0f;
                CA->t_new[i]->ECTNC_ap_t_end   = 0.0f;
            }else{
                CA->t_new[i]->F_state = F4;
            }
            break;
        default:
            CA->t_new[i]->F_state = CA->t_old[i]->F_state;
    }
}
/**
 * 
 * @param i
 * @param CA
 */
void CAStep_i(int i, typ_ca *CA)
{
    int countStimulingNeighbors=0;

	lst_item *neighbor = CA->omega_a[i];
	//runs over the neighbors of each elements
	while(neighbor != NULL)
	{
            //0 pode ser estimulado
            //1 pode estimular os vizinhos
            //2 pode estimular os vizinhos
            //3 nao pode ser estimulado ou estimular
            //4 pode ser estimulado por 2 vizinhos ou mais mas nao pode estimular
            //se o elemento esta nos estados 0 e 4, onde o mesmo pode ser estimulado

            int neighbor_STATE_V = CA->t_old[neighbor->value]->V_state;
            //verifica se estimulo eh capaz de percorrer a distanica entre os elementos
            if( (neighbor_STATE_V==V1 || neighbor_STATE_V==V2) &&
                (CA->t_old[i]->V_state==V0 ||CA->t_old[i]->V_state==V4) )
            {
                double propTime		= getPropagationTime( neighbor->value, i, CA);
                countStimulingNeighbors	+= (CA->t_old[neighbor->value]->cellT>=propTime)?1:0;
            }
            neighbor = neighbor->next;
	}//end of visiting neighbors
    int activate=0;
    //verifies if it is time to activate cell
    if( (CA->ini[i]->cellCond==paceMaker) ){//if it is pacemaker and the stimulation time has come
        if( isStimulationTime(i, CA)==1 ){
            activate=1;
        }
    }
    if( (CA->t_old[i]->V_state==V0 && countStimulingNeighbors>=1) ||//if it is state V0 and it has at least 1 activated neighbour
        (CA->t_old[i]->V_state==V4 && countStimulingNeighbors>=2))//if it is state V4 and there are more than 2 neighbours
	{
            activate=1;
	}
    if(activate==1){
        cellActivation(i, CA);
    }else{
        incStates(i, CA->params->dt, CA);
    }
}
/**
 * 
 * @param numThreads
 * @param threadsByIndividual
 * @return 
 */
int getNumThreads(int numThreads, int *threadsByIndividual){
    return (threadsByIndividual!=NULL)?*threadsByIndividual:numThreads;
}
/**
 * 
 * @param CA
 * @param forcesOnPts
 */
void stats(typ_ca *CA, double *forcesOnPts){
    if(CA->volume>CA->stats->maxVol) CA->stats->maxVol=CA->volume;
    if(CA->volume<CA->stats->minVol) CA->stats->minVol=CA->volume;
    double deltaVol = fabs(CA->volume - CA->stats->volIni);
    if(deltaVol > CA->stats->maxDeltaVol){
        CA->stats->volMaxDelta = CA->volume;
        CA->stats->maxDeltaVol = deltaVol;
    }
    CA->stats->max[0]=CA->stats->max[1]=CA->stats->max[2]=-DBL_MAX;
    CA->stats->min[0]=CA->stats->min[1]=CA->stats->min[2]=+DBL_MAX;
    double sumVel=0.0;
    for(int k=0;k<CA->params->pointsNum;k++){
        //Zero the forces
        //Essa linha foi pra o método de verlet forcesOnPts[I2d(k,0,3)]=forcesOnPts[I2d(k,1,3)]=forcesOnPts[I2d(k,2,3)]=0.0f;
        //
        CA->stats->min[0] = min(CA->stats->min[0], CA->pnts_old[k]->x);
        CA->stats->min[1] = min(CA->stats->min[1], CA->pnts_old[k]->y);
        CA->stats->min[2] = min(CA->stats->min[2], CA->pnts_old[k]->z);
        CA->stats->max[0] = my_max(CA->stats->max[0], CA->pnts_old[k]->x);
        CA->stats->max[1] = my_max(CA->stats->max[1], CA->pnts_old[k]->y);
        CA->stats->max[2] = my_max(CA->stats->max[2], CA->pnts_old[k]->z);
        
        double v[3]={CA->pnts_old[k]->xV, CA->pnts_old[k]->yV, CA->pnts_old[k]->zV};
        sumVel += my_norm(v);
    }
    CA->stats->avgVel = sumVel / CA->params->pointsNum;
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
int simulate(  typ_ca *CA, bool save,int *threadsByIndividual){
    //time units in SECONDS
    CA->time  = 0.0;
    CA->timeSaving=0.0;
    CA->stats->contSave=0;

    CA->stats->minVol=DBL_MAX;
    CA->stats->maxVol=0.0;
    
    int count=0;
    typ_point **aux3;
    typ_dt_element **auxCA;
    char filename[255];
    sprintf(filename, "%sfisiopacer.txt", CA->params->outputFolder/*.c_str()*/);
    FILE *fileDt = fopen(filename, "w+");
        
    double *forcesOnPts_interm  = (double*)malloc(sizeof(double)*CA->params->pointsNum*3);
    if(forcesOnPts_interm==NULL) {
        throw MyException("Allocation failure for forcesOnPts_interm.", __FILE__, __LINE__);
    }
    double *forcesOnPts         = (double*)malloc(sizeof(double)*CA->params->pointsNum*3);
    if(forcesOnPts==NULL){
        throw MyException("Allocation failure for forcesOnPts.", __FILE__, __LINE__);
    }
    for(int k=0;k<CA->params->pointsNum;k++){
        forcesOnPts[I2d(k,0,3)]=forcesOnPts[I2d(k,1,3)]=forcesOnPts[I2d(k,2,3)]=0.0f;
        forcesOnPts_interm[I2d(k,0,3)]=forcesOnPts_interm[I2d(k,1,3)]=forcesOnPts_interm[I2d(k,2,3)]=0.0f;
    }
    initializeCA( CA );
    double sumDt=0.0; 
    CA->stats->volIni=0.0;
    for(int i=0;i<CA->params->elementsNum;i++){
        CA->stats->volIni+= CA->t_old[i]->volCel;
    }
    CA->volume = CA->stats->volIni;
    int retValFinal = 0;
    int nThreads =1;
    
    CA->stats->maxDeltaVol=0.0;
    CA->stats->volMaxDelta=0.0;
    CA->stats->avgVel=0.0;
    stats(CA, forcesOnPts);
    //cout<<"VOlume ini " << CA->volume<<endl;
    while(CA->time<=CA->params->simulationTime)
    {
        nThreads = getNumThreads(CA->params->numThreads, threadsByIndividual);
        simulationStep( nThreads, CA,  forcesOnPts);
        if(CA->params->mecSim==1){
           // EulerMethod(CA, forcesOnPts);
           VelocityVerletMethod(CA, forcesOnPts, forcesOnPts_interm);
        }
        incPressureStates(CA);
        if(save){
            save_step(fileDt, CA, CA->params->outputFolder, forcesOnPts);
        }
        for(int k=0;k<CA->params->pointsNum;k++){
            forcesOnPts[I2d(k,0,3)]=forcesOnPts[I2d(k,1,3)]=forcesOnPts[I2d(k,2,3)]=0.0f;
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
    if(forcesOnPts!=NULL)
        free(forcesOnPts);
    forcesOnPts = NULL;
    if(forcesOnPts_interm!=NULL)       
        free(forcesOnPts_interm);
    forcesOnPts_interm = NULL;
    if(CA->params->printOutput==1){
        printf("final time: %.3f\n",CA->time);
        printf("Volume: [%.3f %.3f] [%.3e %.3e] \n", CA->stats->minVol/CA->stats->volIni*100.0, CA->stats->maxVol/CA->stats->volIni*100.0, CA->stats->volIni , CA->volume);
        printf("iterações: %d. Dt medio %g\n", count, sumDt/count);
    }
    
    return retValFinal;
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
int startCA(string paramFile, bool save, int *threadsByIndividual)
{
    try {


        //
        typ_ca *CA = (typ_ca *) malloc(sizeof(typ_ca));
        if (CA == NULL) {
            throw MyException("Allocation failure for CA.", __FILE__, __LINE__);
        }
        readParameterFile(paramFile, CA);
        allocCA(CA);
        int retSim = simulate(CA, save, threadsByIndividual);

        deallocCA(CA);
        return retSim;
    }catch(MyException& caught){
        std::cout<<caught.getMessage()<<std::endl;
    }

}
/**
 * 
 * @param CA
 */
void deallocCA(typ_ca *CA){
    try{
        myArrayDeallocation<typ_dt_element>(CA->t_old, CA->params->elementsNum);
        myArrayDeallocation<typ_dt_element>(CA->t_new, CA->params->elementsNum);
        myArrayDeallocation<typ_t0_element>(CA->ini,   CA->params->elementsNum);
        //
        myArrayDeallocation<typ_point>(CA->pnts_old,   CA->params->pointsNum);
        myArrayDeallocation<typ_point>(CA->pnts_new,   CA->params->pointsNum);
        myArrayDeallocation<typ_point>(CA->pnts_intrm,   CA->params->pointsNum);
        //
        freeAList(CA->omega_a, CA->params->elementsNum);
        freeAList(CA->omega_b, CA->params->pointsNum);

        for(int i=0;i<CA->params->stimSize;i++){
            if(CA->params->aStim[i]!=NULL) free(CA->params->aStim[i]);
        }
        if(CA->params->aStim!=NULL) free(CA->params->aStim);
        //
        for(int i=0;i<CA->params->numFaces;i++){
            if(CA->params->aFaces[i]!=NULL) free(CA->params->aFaces[i]);
        }
        if(CA->params->aFaces!=NULL) free(CA->params->aFaces);
        //
        for(int i=0;i<CA->params->nRegions;i++){
            if(CA->params->aParam[i]!=NULL)
                free(CA->params->aParam[i]);
        }
        if(CA->params->aParam!=NULL) free(CA->params->aParam);                
        if(CA->params!=NULL) free(CA->params);
        if(CA->stats!=NULL) free(CA->stats);
        if(CA->pressureCA!=NULL) free(CA->pressureCA);
        if(CA!=NULL) free(CA);
    }catch(MyException& caught){
        std::cout<<caught.getMessage()<<std::endl;
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
void allocCA(typ_ca* CA){
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

    try{
        //le o tamanho das malhas nos arquivos de entrada
        CA->params->pointsNum    = readSize(strPtsFile);
        CA->params->elementsNum  = readSize(strEleFile);
        readStimFile(strStmFile, CA);
    	CA->omega_a = iniAList(CA->params->elementsNum);
        if(CA->omega_a==NULL) {
            throw MyException("Allocation failure for CA->omega_a.", __FILE__, __LINE__);
        }
        CA->omega_b = iniAList(CA->params->pointsNum);
        if(CA->omega_b==NULL){
            throw MyException("Allocation failure for CA->omega_b.", __FILE__, __LINE__);
        }
        //points
        CA->pnts_new  = myArrayAllocation<typ_point>(CA->params->pointsNum);
        if(CA->pnts_new==NULL){
            throw MyException("Allocation failure for CA->pnts_new.", __FILE__, __LINE__);
        }
        CA->pnts_old  = myArrayAllocation<typ_point>(CA->params->pointsNum);
        if(CA->pnts_old==NULL){
            throw MyException("Allocation failure for CA->pnts_old.", __FILE__, __LINE__);
        }
        CA->pnts_intrm= myArrayAllocation<typ_point>(CA->params->pointsNum);
        if(CA->pnts_intrm==NULL){
            throw MyException("Allocation failure for CA->pnts_intrm.", __FILE__, __LINE__);
        }
        //elements
        CA->t_old = myArrayAllocation<typ_dt_element>(CA->params->elementsNum);
        if(CA->t_old==NULL){
            throw MyException("Allocation failure for CA->t_old.", __FILE__, __LINE__);
        }
        CA->t_new = myArrayAllocation<typ_dt_element>(CA->params->elementsNum);
        if(CA->t_new==NULL){
            throw MyException("Allocation failure for CA->t_new.", __FILE__, __LINE__);
        }
        CA->ini   = myArrayAllocation<typ_t0_element>(CA->params->elementsNum);
        if(CA->ini==NULL){
            throw MyException("Allocation failure for CA->ini.", __FILE__, __LINE__);
        }
        //
               
        CA->stats = (typ_stats*)malloc(sizeof(typ_stats));
        if(CA->stats==NULL){
            throw MyException("Allocation failure for stats structure.", __FILE__, __LINE__);
        }
        //opens the files and fill arrays
	    openFile(CA, strPtsFile, strEleFile, strFibFile, strBoundFile, strPressFile);
        CA->pressureCA = (typ_press*)malloc(sizeof(typ_press));
        if(CA->pressureCA==NULL){
            throw MyException("Allocation failure for stats structure.", __FILE__, __LINE__);
        }
    }catch(MyException& caught){
        std::cout<<caught.getMessage()<<std::endl;
    }
}
/**
 * 
 * @param i
 * @param CA
 * @return 
 */
double getAVGNeighborAPDElectroTonic(int i, typ_ca *CA)
{
    lst_item *neighbor = CA->omega_a[i];
    int count=0;
    double sum=0.;
    //runs over the neighbors of each elements
    while(neighbor != NULL)
    {
        sum+=CA->t_old[neighbor->value]->APTime4;
        neighbor = neighbor->next;
        count++;
    }//end of visiting neighbors
    return sum/count;
}
/**
 * 
 * @param i
 * @param CA
 * @return 
 */
double getIoutElectroTonic(int i ,typ_ca *CA){
    int iRegion = CA->ini[i]->iRegion;
    t_par_ac* par = CA->params->aParam[iRegion];
    return par->pEletroTonic*getActiveTensionNormalized(i, CA);
}
/**
 * 
 * @param i
 * @param apd
 * @param mapd
 * @return 
 */
double newAPDElectroTonic(double iout, double apd, double mapd){
    return (1.0-iout)*apd + iout*mapd;
}
/**
 * 
 * @param i
 * @param CA
 * @return 
 */
void computeNewAPDElectroTonic(int i, typ_ca *CA){
    double iout = getIoutElectroTonic(i, CA);
    double mapd = getAVGNeighborAPDElectroTonic(i, CA);
    CA->t_new[i]->APTime4 = newAPDElectroTonic(iout, CA->t_old[i]->APTime4, mapd);
    double ratio = CA->t_new[i]->APTime4 / CA->t_old[i]->APTime4;
    CA->t_new[i]->APTime1 = ratio * CA->t_old[i]->APTime1;
    CA->t_new[i]->APTime2 = ratio * CA->t_old[i]->APTime2;
    CA->t_new[i]->APTime3 = ratio * CA->t_old[i]->APTime3;
    
}

/**
 * 
 * @param i
 * @param CA
 */
void restartAPDElectroTonic(int i, typ_ca *CA)
{
    int iRegion = CA->ini[i]->iRegion;
    t_par_ac* ap = CA->params->aParam[iRegion]; 
    CA->t_new[i]->APTime1  = ap->APTime1ini;
    CA->t_new[i]->APTime2  = ap->APTime2ini;
    CA->t_new[i]->APTime3  = ap->APTime3ini;
    CA->t_new[i]->APTime4  = ap->APTime4ini;
}