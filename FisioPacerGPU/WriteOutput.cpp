#include "WriteOutput.h"


/**
 * 
 * @param fileDt
 * @param CA
 * @param forcesOnPts
 * @param strFolderOut
 */
void save_step(FILE *fileDt, typ_ca *CA, string strFolderOut, double *forcesOnPts){
    double diff =  CA->time - CA->timeSaving;
    if(diff>=0.0)
    {//checks if diff is too close to zero
        fprintf(fileDt, "%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", 
            CA->stats->contSave , 
            CA->params->dt, CA->stats->avgVel, CA->stats->tol, CA->stats->err, 
            getActiveTensionNormalized(0, CA),
            getActiveTensionDiscrete( 0, CA),
            getV(1, CA), getV(0, CA), 
            getV(2, CA), getVdiscret(0, CA), 
            getPressurePercent(CA), getPressureDiscrete(CA),
            CA->stats->min[0], CA->stats->min[1], CA->stats->min[2], 
            CA->stats->max[0], CA->stats->max[1], CA->stats->max[2],
            CA->volume);
        //saveDebug(CA, strFolderOut, forcesOnPts);
        saveVTK_V_point(CA, strFolderOut);
        /*saveCSV(contSave, CA, strFolderOut);
        saveVTK_V_cell(contSave, POINTS_OLD, strFolderOut);*/
        CA->timeSaving += CA->params->dtSave;
        CA->stats->contSave++;
    }
}


void saveCSV( typ_ca *CA, string outFolder)
{
    char filename[255];
    sprintf(filename, "%sfile.%05d.csv", outFolder.c_str(), CA->stats->contSave);
    FILE *file = fopen(filename, "w+");
    if(file==NULL){
        stringstream ss;
        ss<<"Failed to open file: "<<filename;
        throw MyException(ss.str(), __FILE__, __LINE__);
    }
    //fprintf(file, "%d\n",CA->params->pointsNum);
//    fprintf(file, "x y z\n");
    for(int i=0;i<CA->params->pointsNum;i++){
        fprintf(file, "%f %f %f\n", CA->pnts_old[i]->x, CA->pnts_old[i]->y, CA->pnts_old[i]->z);
    }
    if(fclose(file)!=0)throw MyException("Failed to close file.", __FILE__, __LINE__);
}

void saveVTK_V_cell(
        int t,
        typ_ca *CA,
        string outFolder
        )
{
    char filename[255];
    sprintf(filename, "%sfile.%d.vtu", outFolder.c_str(), t);


    FILE *file = fopen(filename, "w+");

    fprintf(file, "<?xml version=\"1.0\"?>\n");
    fprintf(file, "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >\n");
    fprintf(file, "<UnstructuredGrid>\n");
    fprintf(file, "<Piece  NumberOfPoints=\"%d\"  NumberOfCells=\"%d\">\n", CA->params->pointsNum, CA->params->elementsNum);
    fprintf(file, "<Points>\n");
    fprintf(file, "<DataArray  type=\"Float32\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->pointsNum;i++){
        fprintf(file, "%g ",CA->pnts_old[i]->x );
        fprintf(file, "%g ",CA->pnts_old[i]->y );
        fprintf(file, "%g ",CA->pnts_old[i]->z );

    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "</Points>\n");
    fprintf(file, "<Cells>\n");
    fprintf(file, "<DataArray  type=\"UInt32\"  Name=\"connectivity\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%d ", CA->ini[i]->iPt1);
        fprintf(file, "%d ", CA->ini[i]->iPt2);
        fprintf(file, "%d ", CA->ini[i]->iPt3);
        fprintf(file, "%d ", CA->ini[i]->iPt4);
        fprintf(file, " ");
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "<DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"ascii\">\n");
    int countOffset=elementsCols;
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%d ",countOffset);
        countOffset +=elementsCols;
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "10 ");
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "</Cells>\n");

    fprintf(file, "<CellData>\n");
    fprintf(file, "<DataArray type=\"Float32\" Name=\"v\" format=\"ascii\">\n");
    float v;
    for(int i=0;i<CA->params->elementsNum;i++){
        v= getV(i, CA);
        fprintf(file, "%.7f ", v);
    }
    fprintf(file, "</DataArray>\n");
    
    
    /*fprintf(file, "<DataArray type=\"Float32\" Name=\"f\" format=\"ascii\">\n");
    float f;
    for(int i=0;i<p->elementsNum;i++){
        if(p->forceMultiplier!=0.0)
            f= getActiveTension(el_old, i, p, 0.0)/p->forceMultiplier;
        else f=0.0;
        fprintf(file, "%g ", f);
    }
    fprintf(file, "</DataArray>\n");
    */
    /*fprintf(file, "<DataArray type=\"Float32\" Name=\"F_State\" format=\"ascii\">\n");
    for(int i=0;i<p->elementsNum;i++){
        fprintf(file, "%d ", el_old[i]->V_state);
    }
    fprintf(file, "</DataArray>\n");*/

    /*fprintf(file, "<DataArray type=\"Float32\" Name=\"vol\" format=\"ascii\">\n");
    for(int i=0;i<p->elementsNum;i++){
        fprintf(file, "%f ", elGeo_old[i]->volCel);
    }
    fprintf(file, "</DataArray>\n");*/

    fprintf(file, "</CellData>\n");
    fprintf(file, "</Piece>\n");
    fprintf(file, "</UnstructuredGrid>\n");
    fprintf(file, "</VTKFile>\n");

    fclose(file);

}
void saveNormals(
        typ_ca *CA, string outFile)
{
    char filename[255];
    sprintf(filename, "%snormal_%d.vtk", outFile.c_str(), CA->stats->contSave);
    FILE *file = fopen(filename, "w+");
    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "vtk output\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(file, "POINTS %d double\n", CA->params->numFaces);
    typ_face *face;    
    for(int iF=0; iF < CA->params->numFaces; iF++){
        face = CA->params->aFaces[iF];
        double normal[3]={0.0}, bary[3]={0.0};
        getFaceAreaNormal(face, CA, normal,bary);
        fprintf(file, "%f %f %f\n", bary[0], bary[1], bary[2]);
    }   
    fprintf(file, "CELLS %d %d\n", CA->params->numFaces,CA->params->numFaces*2);
    for(int iF=0; iF < CA->params->numFaces; iF++){
        fprintf(file, "1 %d\n",iF);
    }  
    fprintf(file, "CELL_TYPES %d\n", CA->params->numFaces);
    for(int iF=0; iF < CA->params->numFaces; iF++){
        fprintf(file, "1\n");
    }  
    fprintf(file, "POINT_DATA %d\n", CA->params->numFaces);
    fprintf(file, "NORMALS n double\n");
    for(int iF=0; iF < CA->params->numFaces; iF++){
        face = CA->params->aFaces[iF];
        double normal[3]={0.0}, bary[3]={0.0};
        getFaceAreaNormal(face, CA, normal,bary);
        fprintf(file, "%f %f %f\n", normal[0], normal[1], normal[2]);
    }

    

    fclose(file);
}


void saveVTK_V_point(
        typ_ca *CA,
        string outFile)
{
    char filename[255];
    sprintf(filename, "%sfile_%d.vtu", outFile.c_str(), CA->stats->contSave);

    FILE *file = fopen(filename, "w+");

    fprintf(file, "<?xml version=\"1.0\"?>\n");
    fprintf(file, "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >\n");
    fprintf(file, "<UnstructuredGrid>\n");
    fprintf(file, "<Piece  NumberOfPoints=\"%d\"  NumberOfCells=\"%d\">\n", CA->params->pointsNum, CA->params->elementsNum);
    fprintf(file, "<Points>\n");
    fprintf(file, "<DataArray  type=\"Float32\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");
    
    for(int i=0;i<CA->params->pointsNum;i++){
        fprintf(file, "%f %f %f ",
                CA->pnts_old[i]->x/*/10.0*//**1000.0*/,
                CA->pnts_old[i]->y/*/10.0*//**1000.0*/,
                CA->pnts_old[i]->z/*/10.0*//**1000.0*/);

    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "</Points>\n");
    fprintf(file, "<Cells>\n");
    fprintf(file, "<DataArray  type=\"UInt32\"  Name=\"connectivity\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%d ",CA->ini[i]->iPt1);
        fprintf(file, "%d ",CA->ini[i]->iPt2);
        fprintf(file, "%d ",CA->ini[i]->iPt3);
        fprintf(file, "%d ",CA->ini[i]->iPt4);
        fprintf(file, " ");
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "<DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"ascii\">\n");
    int countOffset=elementsCols;
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%d ",countOffset);
        countOffset +=elementsCols;
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "10 ");
    }
    fprintf(file, "</DataArray>\n");
    
    fprintf(file, "</Cells>\n");
    fprintf(file, "<PointData>\n");
   
    fprintf(file, "</PointData>\n");

    fprintf(file, "<CellData>\n");

    fprintf(file, "<DataArray type=\"Float32\" Name=\"V\" format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%.2f ", getV(i,CA));
    }
    fprintf(file, "</DataArray>\n");
    
    fprintf(file, "<DataArray type=\"Float32\" Name=\"APD\" format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%.2f ", CA->t_old[i]->APTime4);
    }
    fprintf(file, "</DataArray>\n");
    
    fprintf(file, "</CellData>\n");
    fprintf(file, "</Piece>\n");

    fprintf(file, "</UnstructuredGrid>\n");
    fprintf(file, "</VTKFile>\n");

    fclose(file);
}

void saveDebug(
        typ_ca *CA,
        string outFile, double *forcesOnPts)
{
    saveNormals(CA, outFile);
            
    char filename[255];
    sprintf(filename, "%sfile_%d.vtu", outFile.c_str(), CA->stats->contSave);

    FILE *file = fopen(filename, "w+");

    fprintf(file, "<?xml version=\"1.0\"?>\n");
    fprintf(file, "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >\n");
    fprintf(file, "<UnstructuredGrid>\n");
    fprintf(file, "<Piece  NumberOfPoints=\"%d\"  NumberOfCells=\"%d\">\n", CA->params->pointsNum, CA->params->elementsNum);
    fprintf(file, "<Points>\n");
    fprintf(file, "<DataArray  type=\"Float32\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");
    
    for(int i=0;i<CA->params->pointsNum;i++){
        fprintf(file, "%f %f %f ",
                CA->pnts_old[i]->x/*/10.0*//**1000.0*/,
                CA->pnts_old[i]->y/*/10.0*//**1000.0*/,
                CA->pnts_old[i]->z/*/10.0*//**1000.0*/);

    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "</Points>\n");
    fprintf(file, "<Cells>\n");
    fprintf(file, "<DataArray  type=\"UInt32\"  Name=\"connectivity\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%d ",CA->ini[i]->iPt1);
        fprintf(file, "%d ",CA->ini[i]->iPt2);
        fprintf(file, "%d ",CA->ini[i]->iPt3);
        fprintf(file, "%d ",CA->ini[i]->iPt4);
        fprintf(file, " ");
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "<DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"ascii\">\n");
    int countOffset=elementsCols;
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%d ",countOffset);
        countOffset +=elementsCols;
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "10 ");
    }
    fprintf(file, "</DataArray>\n");
    
    fprintf(file, "</Cells>\n");
    fprintf(file, "<PointData>\n");
   
    fprintf(file, "<DataArray Name=\"displacement\" type=\"Float32\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->pointsNum;i++){
        double displacement[3];
        displacement[0] = CA->pnts_old[i]->x - CA->pnts_old[i]->x0;
        displacement[1] = CA->pnts_old[i]->y - CA->pnts_old[i]->y0;
        displacement[2] = CA->pnts_old[i]->z - CA->pnts_old[i]->z0;
        //fprintf(file, "%e ", my_norm(aForce));
        fprintf(file, "%e %e %e ", displacement[0], displacement[1], displacement[2]);
    }
    fprintf(file, "</DataArray>\n");
    
    fprintf(file, "<DataArray Name=\"force\" type=\"Float32\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->pointsNum;i++){
        double aForce[3];
        aForce[0] = forcesOnPts[I2d(i,0,3)];
        aForce[1] = forcesOnPts[I2d(i,1,3)];
        aForce[2] = forcesOnPts[I2d(i,2,3)];
        //fprintf(file, "%e ", my_norm(aForce));
        fprintf(file, "%e %e %e ", aForce[0], aForce[1], aForce[2]);
    }
    fprintf(file, "</DataArray>\n");
    
    fprintf(file, "<DataArray Name=\"vel\" type=\"Float32\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->pointsNum;i++){
        fprintf(file, "%.7f %.7f %.7f ", CA->pnts_old[i]->xV, CA->pnts_old[i]->yV, CA->pnts_old[i]->zV);
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "<DataArray Name=\"restr\" type=\"Float32\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->pointsNum;i++){
        fprintf(file, "%d %d %d ", CA->pnts_old[i]->xRestr, CA->pnts_old[i]->yRestr, CA->pnts_old[i]->zRestr);
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "</PointData>\n");

    fprintf(file, "<CellData>\n");
            
    fprintf(file, "<DataArray type=\"Float32\" Name=\"f\" format=\"ascii\">\n");
    float f;
    for(int i=0;i<CA->params->elementsNum;i++){
        f= getActiveTension(i, CA);
        fprintf(file, "%g ", f);
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "<DataArray type=\"Float32\" Name=\"v_cell\" format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%d ", CA->t_old[i]->V_state);
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "<DataArray type=\"Float32\" Name=\"V\" format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%.2f ", getV(i,CA));
    }
    fprintf(file, "</DataArray>\n");
    
    //Fibra
    fprintf(file, "<DataArray type=\"Float32\" Name=\"Fiber\" format=\"ascii\" NumberOfComponents=\"3\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%g %g %g ", CA->t_old[i]->fiberDir[0], CA->t_old[i]->fiberDir[1],CA->t_old[i]->fiberDir[2]);
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "<DataArray type=\"Float32\" Name=\"sheet\" format=\"ascii\" NumberOfComponents=\"3\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%g %g %g ", CA->t_old[i]->sheetDir[0], CA->t_old[i]->sheetDir[1], CA->t_old[i]->sheetDir[2]);
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "<DataArray type=\"Float32\" Name=\"nsheet\" format=\"ascii\" NumberOfComponents=\"3\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%g %g %g ", CA->t_old[i]->nsheetDir[0], CA->t_old[i]->nsheetDir[1], CA->t_old[i]->nsheetDir[2]);
    }
    fprintf(file, "</DataArray>\n");
    
    fprintf(file, "<DataArray type=\"Float32\" Name=\"KAng\" format=\"ascii\" NumberOfComponents=\"3\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%g %g %g ", CA->ini[i]->KAng[0], CA->ini[i]->KAng[1],CA->ini[i]->KAng[2]);
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "<DataArray type=\"Float32\" Name=\"KAxl\" format=\"ascii\" NumberOfComponents=\"3\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%g %g %g ", CA->ini[i]->KAxl[0], CA->ini[i]->KAxl[1],CA->ini[i]->KAxl[2]);
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "<DataArray type=\"Float32\" Name=\"vol\" format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%g ", CA->t_old[i]->volCel);
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "<DataArray type=\"Float32\" Name=\"alpha12\" format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%g ", CA->t_old[i]->alphaT_12);
    }
    fprintf(file, "</DataArray>\n");
    //
    fprintf(file, "<DataArray type=\"Float32\" Name=\"alpha13\" format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%g ", CA->t_old[i]->alphaT_13);
    }
    fprintf(file, "</DataArray>\n");
    //
    fprintf(file, "<DataArray type=\"Float32\" Name=\"alpha23\" format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%g ", CA->t_old[i]->alphaT_23);
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "<DataArray type=\"Float32\" Name=\"axisLength\" format=\"ascii\"  NumberOfComponents=\"3\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%g %g %g ",  CA->t_old[i]->axisLengthT[0], CA->t_old[i]->axisLengthT[1], CA->t_old[i]->axisLengthT[2]);
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "</CellData>\n");
    fprintf(file, "</Piece>\n");

    fprintf(file, "</UnstructuredGrid>\n");
    fprintf(file, "</VTKFile>\n");

    fclose(file);
}


void saveVTK_Simples(
        int t,
        typ_ca *CA,
        string outFolder)
{
    char filename[255];
    sprintf(filename, "%sfile_%d.vtu", outFolder.c_str(), t);

    FILE *file = fopen(filename, "w+");

    fprintf(file, "<?xml version=\"1.0\"?>\n");
    fprintf(file, "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >\n");
    fprintf(file, "<UnstructuredGrid>\n");
    fprintf(file, "<Piece  NumberOfPoints=\"%d\"  NumberOfCells=\"%d\">\n", CA->params->pointsNum, CA->params->elementsNum);
    fprintf(file, "<Points>\n");
    fprintf(file, "<DataArray  type=\"Float32\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->pointsNum;i++){
        fprintf(file, "%g %g %g ",CA->pnts_old[i]->x, CA->pnts_old[i]->y, CA->pnts_old[i]->z);

    }
    fprintf(file, "</DataArray>\n");
    fprintf(file, "</Points>\n");
    fprintf(file, "<Cells>\n");
    fprintf(file, "<DataArray  type=\"UInt32\"  Name=\"connectivity\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%d ",CA->ini[i]->iPt1);
        fprintf(file, "%d ",CA->ini[i]->iPt2);
        fprintf(file, "%d ",CA->ini[i]->iPt3);
        fprintf(file, "%d ",CA->ini[i]->iPt4);

        fprintf(file, " ");
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "<DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"ascii\">\n");
    int countOffset=elementsCols;
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "%d ",countOffset);
        countOffset +=elementsCols;
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->elementsNum;i++){
        fprintf(file, "10 ");
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "</Cells>\n");
    fprintf(file, "<PointData>\n");
    fprintf(file, "<DataArray type=\"Float32\" Name=\"v\" format=\"ascii\">\n");
    float V,vol, numerador, denominador;
    int adjacentElement;
    float maior;
    for(int i=0;i<CA->params->pointsNum;i++){
        maior = -99999.0;
        numerador=0.0;
        denominador=0.0;
        lst_item *cur = CA->omega_b[i];
        while(cur != NULL){
            adjacentElement = cur->value;
            V = getV(adjacentElement, CA);
            vol = CA->t_old[adjacentElement]->volCel;
            numerador   += V*vol;
            denominador += vol;
            if(V>maior) maior=V;
            cur = cur->next;
        }
        if(denominador!=0.0){
            //printf("%d\n",i);
            //fprintf(file, "%g ", numerador/denominador);
            fprintf(file, "%g ", maior);
        }else {
            fprintf(file, "%g ", 0.0);
        }
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "<DataArray Name=\"restriction\" type=\"UInt32\"  NumberOfComponents=\"1\"  format=\"ascii\">\n");
    for(int i=0;i<CA->params->pointsNum;i++){
        fprintf(file, "%d ",CA->pnts_old[i]->zRestr);
//        fprintf(file, "%d ",points[i]->isRestrictedY);
//        fprintf(file, "%d ",points[i]->isRestrictedZ);

    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "</PointData>\n");

    fprintf(file, "</Piece>\n");

    fprintf(file, "</UnstructuredGrid>\n");
    fprintf(file, "</VTKFile>\n");

    fclose(file);
}
