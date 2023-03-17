#include "Stats.h"

typ_stats* stats;

__host__ double* statsGetMax() { return stats->max; }
__host__ double* statsGetMin() { return stats->min; }
__host__ double statsGetAvgVel() { return stats->avgVel; }
__host__ double statsGetMaxVol() { return stats->maxVol; }
__host__ double statsGetMinVol() { return stats->minVol; }
__host__ double statsGetVolIni() { return stats->volIni; }
__host__ double statsGetVolMaxDelta() { return stats->volMaxDelta; }
__host__ double statsGetMaxDeltaVol() { return stats->maxDeltaVol; }
__host__ int statsGetContSave() { return stats->contSave; }
__host__ double statsGetTol() { return stats->tol; }
__host__ double statsGetErr() { return stats->err; }

__host__ void statsSetTol(double tol) {
    stats->tol = tol;
}
__host__ void statsSetErr(double err) {
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
__host__ void computeStats(double volume) {
    if (volume > stats->maxVol) stats->maxVol = volume;
    if (volume < stats->minVol) stats->minVol = volume;
    double deltaVol = fabs(volume - stats->volIni);
    if (deltaVol > stats->maxDeltaVol) {
        stats->volMaxDelta = volume;
        stats->maxDeltaVol = deltaVol;
    }
    pointsStats(stats->max, stats->min, &(stats->avgVel));
}

/*
*/
__host__ void statsIncContSave() {
    stats->contSave++;
}