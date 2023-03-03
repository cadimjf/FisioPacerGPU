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
    int iRegion = CA->ini[i].iRegion;
    t_par_ac* aP = CA->params->aParam[iRegion];
    double f = getActiveTensionNormalized(i, CA);
    return f * aP->forceMultiplier;
}

double getActiveTensionNormalized(int i, typ_ca *CA){

    int iRegion = CA->ini[i].iRegion;
    t_par_ac* aP = CA->params->aParam[iRegion];
    if(CA->ini[i].cellCond==DEAD)return 0.0;

    if(CA->params->paSim==0 ) //Se o automato esta desligado, retorna sempre o parametro do modulo da forca
       return 0.0;//aP->forceMultiplier;

    double f=0.0;
    double m;
    double tcell = CA->t_old[i].cellT;
    switch(CA->t_old[i].F_state){
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
    int iRegion = CA->ini[i].iRegion;
    t_par_ac* aP = CA->params->aParam[iRegion];
    if(CA->ini[i].cellCond==DEAD)return 0.0;
    if(CA->params->paSim==0 ) //Se o automato estÃƒÂ¡ desligado, retorna sempre o parametro do modulo da forca
       return aP->forceMultiplier;
    double f=0.0;

    //it is not in FR 1 state, behaves normally
    double tcell = CA->t_old[i].cellT;
    switch(CA->t_old[i].F_state){
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
 * return V in miliVolts
 */
double getV(int i ,typ_ca *CA)
{
    int iRegion = CA->ini[i].iRegion;
    t_par_ac* ap = CA->params->aParam[iRegion];

    double m;
    switch(CA->t_old[i].V_state){
        case V0:
            return ap->VV0;
        case V1:
            m = (ap->VV2 - ap->VV1)/(CA->t_old[i].APTime1);
            return ap->VV1 + m*(CA->t_old[i].cellT);
        case V2:
            m = (ap->VV3 - ap->VV2)/(CA->t_old[i].APTime2 - CA->t_old[i].APTime1);
            return ap->VV2 + m*(CA->t_old[i].cellT - CA->t_old[i].APTime1);
        case V3:
            m = (ap->VV4 - ap->VV3)/(CA->t_old[i].APTime3 - CA->t_old[i].APTime2);
            return ap->VV3 + m*(CA->t_old[i].cellT - CA->t_old[i].APTime2);
        case V4:
            m = (ap->VV0 - ap->VV4)/(CA->t_old[i].APTime4 - CA->t_old[i].APTime3);
            return ap->VV4 + m*(CA->t_old[i].cellT - CA->t_old[i].APTime3);
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
    int iRegion = CA->ini[i].iRegion;
    t_par_ac* ap = CA->params->aParam[iRegion];
    switch(CA->t_old[i].V_state){
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
    t_stim *s = CA->params->aStim[ CA->ini[i].pmRegion];
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
    t_par_ac* ap = CA->params->aParam[CA->ini[iElemNeighbor].iRegion];

    double dir[3] = { 0., 0., 0. }, dirN[3] = {0., 0., 0.};
    //get the vector between the elements position
    dir[0] = CA->t_old[iElemNeighbor].bary[0] - CA->t_old[i].bary[0];
    dir[1] = CA->t_old[iElemNeighbor].bary[1] - CA->t_old[i].bary[1];
    dir[2] = CA->t_old[iElemNeighbor].bary[2] - CA->t_old[i].bary[2];
    //obtem a distancia entre o elemento e o vizinho
    double s = my_norm(dir);
    dirN[0] = dir[0]/s;
    dirN[1] = dir[1]/s;
    dirN[2] = dir[2]/s;
    //////////////
    double v1 = fabs(dot(dirN, CA->t_old[i].fiberDir )) * ap->v[FIBER];
    double v2 = fabs(dot(dirN, CA->t_old[i].sheetDir )) * ap->v[SHEET];
    double v3 = fabs(dot(dirN, CA->t_old[i].nsheetDir)) * ap->v[NORMAL];
    return s/(v1+v2+v3);
}

/**
 *
 * @param i
 */
void cellActivation(int i, typ_ca *CA)
{
     //if the cell is dead, it doesnt increment any states
    if(CA->ini[i].cellCond==DEAD)
        return;

    CA->t_new[i].cellT    = 0.0f;
    CA->t_new[i].V_state  = V1;
    CA->t_new[i].F_state  = F1;
    CA->t_new[i].ECTNC_force_t_ini    = 0.0f;
    CA->t_new[i].ECTNC_force_val_ini     = 0.0f;
    CA->t_new[i].ECTNC_ap_t_ini    = 0.0f;
    CA->t_new[i].ECTNC_ap_val_ini     = 0.0f;
    CA->t_new[i].ECTNC_ap_t_end    = 0.0f;
    CA->t_new[i].ECTNC_force_t_end    = 0.0f;
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
    int iRegion = CA->ini[i].iRegion;
    t_par_ac* ap = CA->params->aParam[iRegion];
    //if the cell is dead, it doesnt increment any states
    if(CA->ini[i].cellCond==DEAD)
        return;
    if(CA->t_old[i].V_state==V0 && CA->t_old[i].F_state==F0){

    } else{
        CA->t_new[i].cellT = CA->t_old[i].cellT + dt;
    }
    switch(CA->t_old[i].V_state){
        case V0:
            //keeps the same value
            CA->t_new[i].V_state = V0;
            break;
        case V1:
            CA->t_new[i].V_state = (CA->t_old[i].cellT >= CA->t_old[i].APTime1)?V2:V1;
            break;
        case V2:
            CA->t_new[i].V_state = (CA->t_old[i].cellT >= CA->t_old[i].APTime2)?V3:V2;
            break;
        case V3:
            CA->t_new[i].V_state = (CA->t_old[i].cellT >= CA->t_old[i].APTime3)?V4:V3;
            break;
        case V4:
            //CA->t_new[i].V_state = (CA->t_old[i].cellT >= CA->t_old[i].APTime4)?V0:V4;
            
            if(CA->t_old[i].cellT >= CA->t_old[i].APTime4){
                CA->t_new[i].V_state = V0;
                //restartAPDElectroTonic(i, CA);
            }else{
                CA->t_new[i].V_state = V4;                
            }

            break;
        default:
            //keeps the same state
            CA->t_new[i].V_state = CA->t_old[i].V_state;
    }

    switch(CA->t_old[i].F_state){
        case F0:
            CA->t_new[i].F_state = F0;
            break;
        case F1:
            CA->t_new[i].F_state = (CA->t_old[i].cellT >= ap->FTime1)?F2:F1;
            break;
        case F2:
            CA->t_new[i].F_state = (CA->t_old[i].cellT >= ap->FTime2)?F3:F2;
            break;
        case F3:
            CA->t_new[i].F_state = (CA->t_old[i].cellT >= ap->FTime3)?F4:F3;
            break;
        case F4:
            if(CA->t_old[i].cellT >= ap->FTime4){
                CA->t_new[i].F_state = F0;
                CA->t_new[i].cellT   = 0.0f;
                CA->t_new[i].ECTNC_force_t_ini   = 0.0f;
                CA->t_new[i].ECTNC_force_val_ini    = 0.0f;
                CA->t_new[i].V_state = V0;
                CA->t_new[i].ECTNC_ap_t_ini   = 0.0f;
                CA->t_new[i].ECTNC_ap_val_ini    = 0.0f;
                CA->t_new[i].ECTNC_ap_t_end   = 0.0f;
            }else{
                CA->t_new[i].F_state = F4;
            }
            break;
        default:
            CA->t_new[i].F_state = CA->t_old[i].F_state;
    }
}
/**
 * 
 * @param i
 * @param CA
 */

void CAStep_i(int i, typ_ca* CA)
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

            int neighbor_STATE_V = CA->t_old[neighbor->value].V_state;
            //verifica se estimulo eh capaz de percorrer a distanica entre os elementos
            if( (neighbor_STATE_V==V1 || neighbor_STATE_V==V2) &&
                (CA->t_old[i].V_state==V0 ||CA->t_old[i].V_state==V4) )
            {
                double propTime		= getPropagationTime( neighbor->value, i, CA);
                countStimulingNeighbors	+= (CA->t_old[neighbor->value].cellT>=propTime)?1:0;
            }
            neighbor = neighbor->next;
	}//end of visiting neighbors
    int activate=0;
    //verifies if it is time to activate cell
    if( (CA->ini[i].cellCond==paceMaker) ){//if it is pacemaker and the stimulation time has come
        if( isStimulationTime(i, CA)==1 ){
            activate=1;
        }
    }
    if( (CA->t_old[i].V_state==V0 && countStimulingNeighbors>=1) ||//if it is state V0 and it has at least 1 activated neighbour
        (CA->t_old[i].V_state==V4 && countStimulingNeighbors>=2))//if it is state V4 and there are more than 2 neighbours
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
        sum+=CA->t_old[neighbor->value].APTime4;
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
    int iRegion = CA->ini[i].iRegion;
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
    CA->t_new[i].APTime4 = newAPDElectroTonic(iout, CA->t_old[i].APTime4, mapd);
    double ratio = CA->t_new[i].APTime4 / CA->t_old[i].APTime4;
    CA->t_new[i].APTime1 = ratio * CA->t_old[i].APTime1;
    CA->t_new[i].APTime2 = ratio * CA->t_old[i].APTime2;
    CA->t_new[i].APTime3 = ratio * CA->t_old[i].APTime3;
    
}

/**
 * 
 * @param i
 * @param CA
 */
void restartAPDElectroTonic(int i, typ_ca *CA)
{
    int iRegion = CA->ini[i].iRegion;
    t_par_ac *ap = CA->params->aParam[iRegion]; 
    CA->t_new[i].APTime1  = ap->APTime1ini;
    CA->t_new[i].APTime2  = ap->APTime2ini;
    CA->t_new[i].APTime3  = ap->APTime3ini;
    CA->t_new[i].APTime4  = ap->APTime4ini;
}