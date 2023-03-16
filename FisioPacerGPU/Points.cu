#include "Points.h"
/*
*/
__host__ void pointsStats(double max[3], double min[3], double* avgVel) {

    max[0] = max[1] = max[2] = -DBL_MAX;
    min[0] = min[1] = min[2] = +DBL_MAX;
    double sumVel = 0.0;
    /* TODO FIXME for (int k = 0; k < CA->params->pointsNum; k++) {
        //Zero the forces
        //Essa linha foi pra o método de verlet forcesOnPts[I2d(k,0,3)]=forcesOnPts[I2d(k,1,3)]=forcesOnPts[I2d(k,2,3)]=0.0f;
        //
        min[0] = my_min(min[0], CA->pnts_old[k].x);
        min[1] = my_min(min[1], CA->pnts_old[k].y);
        min[2] = my_min(min[2], CA->pnts_old[k].z);
        max[0] = my_max(max[0], CA->pnts_old[k].x);
        max[1] = my_max(max[1], CA->pnts_old[k].y);
        max[2] = my_max(max[2], CA->pnts_old[k].z);

        double v[3] = { CA->pnts_old[k].xV, CA->pnts_old[k].yV, CA->pnts_old[k].zV };
        sumVel += my_norm(v);
    }
    *avgVel = sumVel / CA->params->pointsNum;
    */
}