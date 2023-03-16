#include "Stats.h"

typ_stats* stats;

double* statsGetMax() { return stats->max; }
double* statsGetMin() { return stats->min; }
double statsGetAvgVel() { return stats->avgVel; }
double statsGetMaxVol() { return stats->maxVol; }
double statsGetMinVol() { return stats->minVol; }
double statsGetVolIni() { return stats->volIni; }
double statsGetVolMaxDelta() { return stats->volMaxDelta; }
double statsGetMaxDeltaVol() { return stats->maxDeltaVol; }
int statsGetContSave() { return stats->contSave; }
double statsGetTol() { return stats->tol; }
double statsGetErr() { return stats->err; }

void statsSetTol(double tol) {
    stats->tol = tol;
}
void statsSetErr(double err) {
    stats->err = err;
}

/**
*
*/
__host__ void allocStats()
{
    allocateHostVar(&stats);
}
/**
*/
__host__ void deallocStats() {//
    if (stats != NULL) free(stats); 
}
/*
*/
__host__ void iniStats() {
    stats->contSave = 0;
    stats->minVol = DBL_MAX;
    stats->maxVol = 0.0;
    stats->maxDeltaVol = 0.0;
    stats->volMaxDelta = 0.0;
    stats->avgVel = 0.0;
    stats->volIni = 0.0;

}
/*
*/
__host__ void statsSetVolIni(double vol) {
    stats->volIni = vol;
}
/**
 *
 * @param CA
 */
void computeStats(double volume) {
    if (volume > stats->maxVol) stats->maxVol = volume;
    if (volume < stats->minVol) stats->minVol = volume;
    double deltaVol = fabs(volume - stats->volIni);
    if (deltaVol > stats->maxDeltaVol) {
        stats->volMaxDelta = volume;
        stats->maxDeltaVol = deltaVol;
    }
    pointsStats(stats->max, stats->min, &(stats->avgVel));
}


void statsIncContSave() {
    stats->contSave++;
}