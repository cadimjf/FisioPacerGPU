#include <float.h>

#include "Geometry.h"

/**
 * 
 * @param etha
 * @param ksi
 * @return 
 */
double getN1(double etha, double ksi)
{
    return 1.0f -ksi-etha;
}
/**
 * 
 * @param etha
 * @param ksi
 * @return 
 */
double getN2(double etha, double ksi)
{
    return ksi;
}
/**
 * 
 * @param etha
 * @param ksi
 * @return 
 */
double getN3(double etha, double ksi)
{
    return etha;
}
/**
 * 
 * @param CA
 * @param index
 * @return 
 */
 double getVolumeTetrahedron(typ_ca *CA, int index)
{
    //assembly a matrix (a-d, b-d, c-d)
    double aux[3][3] = { { 0., 0., 0. } , {0., 0., 0.} , {0., 0., 0.} };
    int pt1, pt2, pt3, pt4;
    pt4=CA->ini[index]->iPt4;
    //
    pt1=CA->ini[index]->iPt1;
    aux[0][0]= CA->pnts_old[pt1].x - CA->pnts_old[pt4].x;
    aux[1][0]= CA->pnts_old[pt1].y - CA->pnts_old[pt4].y;
    aux[2][0] = CA->pnts_old[pt1].z - CA->pnts_old[pt4].z;
    //
    pt2=CA->ini[index]->iPt2;
    aux[0][1]= CA->pnts_old[pt2].x - CA->pnts_old[pt4].x;
    aux[1][1]= CA->pnts_old[pt2].y - CA->pnts_old[pt4].y;
    aux[2][1]= CA->pnts_old[pt2].z - CA->pnts_old[pt4].z;
    //
    pt3=CA->ini[index]->iPt3;
    aux[0][2]= CA->pnts_old[pt3].x - CA->pnts_old[pt4].x;
    aux[1][2]= CA->pnts_old[pt3].y - CA->pnts_old[pt4].y;
    aux[2][2]= CA->pnts_old[pt3].z - CA->pnts_old[pt4].z;
    
    return fabs(det(aux))/6.0f;
 }
    

/**
 * 
 * @param points
 * @param index
 */
void getBarycenter(typ_ca *CA, int index)
{
    CA->t_new[index]->bary[0]=0.0f;
    CA->t_new[index]->bary[1]=0.0f;
    CA->t_new[index]->bary[2]=0.0f;
    int pt1, pt2, pt3, pt4;
    pt1=CA->ini[index]->iPt1;
    pt2=CA->ini[index]->iPt2;
    pt3=CA->ini[index]->iPt3;
    pt4=CA->ini[index]->iPt4;
    
    //
    CA->t_new[index]->bary[0] = CA->pnts_old[pt1].x + CA->pnts_old[pt2].x +
           CA->pnts_old[pt3].x + CA->pnts_old[pt4].x;
    //
    CA->t_new[index]->bary[1] = CA->pnts_old[pt1].y + CA->pnts_old[pt2].y +
           CA->pnts_old[pt3].y + CA->pnts_old[pt4].y;
    //
    CA->t_new[index]->bary[2] = CA->pnts_old[pt1].z + CA->pnts_old[pt2].z +
           CA->pnts_old[pt3].z + CA->pnts_old[pt4].z;

    CA->t_new[index]->bary[0] = CA->t_new[index]->bary[0]/elementsCols;
    CA->t_new[index]->bary[1] = CA->t_new[index]->bary[1]/elementsCols;
    CA->t_new[index]->bary[2] = CA->t_new[index]->bary[2]/elementsCols;
}

/*
# obtem a area do triangulo definido pelos pontos A, B, C
*/
double areaTriangulo(double A[3], double B[3], double C[3])
{
    double aux1[3] = { 0., 0., 0. }; //A-B
    double aux2[3] = { 0., 0., 0. }; //B-C
    double aux3[3] = { 0., 0., 0. }; //cross product
    for(int i=0;i<3;i++){
        aux1[i] = A[i] - B[i];
        aux2[i] = B[i] - C[i];
    }
    cross(aux1,aux2,aux3);
    return 0.5 * my_norm( aux3 );
}


/**
 * Finds the line parallel to axis that pass through the barycenter:
 * (X-Xb)=tL, so 
 * X=tL+Xb, where t is a parameter,  L is the direction normalized and Xb is the barycenter.
 * the tL is stored at column 0 of retaParalel: [x][0] and the second term, at retaPAralel[x][1]
 * @param bary
 * @param dir
 * @param retaParal
 */
void findLineAxisBary(double bary[3], double dir[3], double retaParal[3][2]){
    for(int i=0;i<3;i++){
        retaParal[i][0] = dir[i];
        retaParal[i][1] = bary[i];
    }
    
}
/**
 * Computes the normal the plane of a face defined by three points 
 * n=(X1-X2) x (X2-X3)
 * and then n is normalized.
 * @param point_old
 * @param pt1
 * @param pt2
 * @param pt3
 * @param n
 */
void findFaceNormal(typ_ca *CA, int pt1, int pt2, int pt3, double n[3]){
    double aux1[3]={0., 0., 0.}, aux2[3] = { 0., 0., 0. };
    aux1[0]=CA->pnts_old[pt1].x - CA->pnts_old[pt2].x;
    aux1[1]=CA->pnts_old[pt1].y - CA->pnts_old[pt2].y;
    aux1[2] = CA->pnts_old[pt1]. z - CA->pnts_old[pt2].z;
    //
    aux2[0]=CA->pnts_old[pt2].x - CA->pnts_old[pt3].x;
    aux2[1]=CA->pnts_old[pt2].y - CA->pnts_old[pt3].y;
    aux2[2]=CA->pnts_old[pt2].z - CA->pnts_old[pt3].z;
    //finds the cross product n=aux1 x aux2
    cross(aux1, aux2, n);
 //normalizes the vector
    double norma = normalizeVector(n, -2);
    if(norma<TOLERANCE){
        stringstream ss;
        ss<<" PAU AQUIIIIIIIIIIIII"<<endl;
        ss<<"NORMAZERO: "<<norma<<endl;
        ss<<". Pontos: ";
        ss<<pt1<<" "<<pt2<<" "<<pt3<<endl;
        ss<<"vetores: "<<endl;
        ss<<aux1[0]<<" "<<aux1[1]<<" "<<aux1[2]<<endl;
        ss<<aux2[0]<<" "<<aux2[1]<<" "<<aux2[2]<<endl;
        ss<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
        ss<<"- -=-=-=-=-=-"<<endl;
        ss<<CA->pnts_old[pt1].x<<" "<<CA->pnts_old[pt1].y<<" "<<CA->pnts_old[pt1].z<<endl;
        ss<<CA->pnts_old[pt2].x<<" "<<CA->pnts_old[pt2].y<<" "<<CA->pnts_old[pt2].z<<endl;
        ss<<CA->pnts_old[pt3].x<<" "<<CA->pnts_old[pt3].y<<" "<<CA->pnts_old[pt3].z<<endl;
        
        string str = ss.str();
        if(CA->params->printOutput==1){
            throw MyException(str, __FILE__, __LINE__);
        }
    }
}

/**
 * t is given by
 * t=((X0-Xb).N ) / (L.N)
 * @param eqLine
 * @param randomPtFace
 * @param nFace
 * @return 
 */
double findT(double eqLine[3][2], double randomPtFace[3], double nFace[3], double *checkDot){
    double aux[3] = { 0., 0., 0. };
    //X0-Xb (Xb is at eqLine[:][1])
    aux[0] = randomPtFace[0] - eqLine[0][1];
    aux[1] = randomPtFace[1] - eqLine[1][1];
    aux[2] = randomPtFace[2] - eqLine[2][1];
    //
    double aux2[3] = { 0., 0., 0. };
    aux2[0] = eqLine[0][0];
    aux2[1] = eqLine[1][0];
    aux2[2] = eqLine[2][0];
    //t=((X0-Xb).N ) / (L.N)
    *checkDot = dot(aux2, nFace);
    return dot(aux, nFace)/ *checkDot;
}

/**
 * 
 * @param reta
 * @param t
 * @param pt
 */
void findIntersectionPoint(double reta[3][2], double t, double pt[3])
{
    pt[0]= reta[0][0]*t + reta[0][1];
    pt[1]= reta[1][0]*t + reta[1][1];
    pt[2]= reta[2][0]*t + reta[2][1];
}

/**
 * 
 * @param sj12
 * @param sj13
 * @param s123
 * @param etha
 * @param ksi
 */
void findInterpolationFunction(double sj12, double sj13, double s123, 
        double *etha, double *ksi){
    *etha = sj12/ s123;
    *ksi  = sj13/ s123;
}
/**
 * @param ptI
 * @param ptA
 * @param ptB
 * @param ptC
 * @param reta
 * @param t
 * @param ksi
 * @param etha
 * @return 
 */
bool checksAreasGetEthaKsi(double ptI[3], double ptA[3], double ptB[3], 
        double ptC[3], double *ksi, double *etha, double *relEr)
{
    
    double sABC = areaTriangulo(ptA, ptB, ptC);
    double sIAB = areaTriangulo(ptI, ptA, ptB);
    double sIAC = areaTriangulo(ptI, ptA, ptC);
    double sIBC = areaTriangulo(ptI, ptB, ptC);
    
    double absEr = fabs(sABC -(sIAB + sIAC + sIBC));
    double myRelEr = fabs(absEr/sABC);
    findInterpolationFunction(sIAB, sIAC, sABC, etha, ksi);    
    if(myRelEr<=RELTOL && absEr<=ABSTOL) 
    {
        return true;
    }else{
        *relEr = myRelEr;
        return false;
    }
}
/**
 * 
 * @param reta
 * @param ptA
 * @param ptB
 * @param ptC
 * @param n
 * @param ptI
 * @return 
 */
int checkPtBelongsToFace(double reta[3][2], double ptA[3], double ptB[3], double ptC[3], 
        double n[3], double ptI[3], double *ksi, double *etha, double *relEr)
{
    double t;
    //checkDot shall contain the value of the dot project: L.N, if it is zero, there's not intersection.
    //If it is not zero, the line crosses the plane, but another check is necessary to verify if the point is in the face
    double checkDot;
    
    //face 1
    t = findT(reta, ptA, n, &checkDot);
    //computes the intersection point
    findIntersectionPoint(reta, t, ptI);
    if(fabs(checkDot)>=TOLERANCE ){
        if(checksAreasGetEthaKsi(ptI, ptA, ptB, ptC, ksi, etha, relEr))
            return 1;
    } 
    return 0;
}
/**
 * 
 * @param ptA
 * @param ptB
 * @return 
 */
bool isPtsEqual(double ptA[3], double ptB[3])
{
    if(fabs(ptA[0]-ptB[0])>=TOLERANCE) return false;
    if(fabs(ptA[1]-ptB[1])>=TOLERANCE) return false;
    if(fabs(ptA[2]-ptB[2])>=TOLERANCE) return false;
    return true;
}
/**
 * 
 * @param ptIValid1
 * @param ptIValid2
 * @param ptI
 * @param check
 * @return 
 */
bool addPoint(double ptIValid1[3], double ptIValid2[3], double ptI[3], int *check)
{
    if(*check==0){
        ptIValid1[0] = ptI[0];
        ptIValid1[1] = ptI[1];
        ptIValid1[2] = ptI[2];
        (*check)++;
        return true; 
    }else if(*check==1){
        if(!isPtsEqual(ptI, ptIValid1))
        {
            ptIValid2[0] = ptI[0];
            ptIValid2[1] = ptI[1];
            ptIValid2[2] = ptI[2];
            (*check)++;
            return true;
        }
    }
    return false;
}

void getPtsInVector(typ_point *pts, int i, double pt[3]){
    pt[0]=pts[i].x;
    pt[1]=pts[i].y;
    pt[2]=pts[i].z;    
}
/**
 * 
 * @param el
 * @param pts
 * @param pt1
 * @param pt2
 * @param pt3
 * @param pt4
 */
void getElementPtsInVector(int i, typ_ca *CA, 
        double pt1[3], double pt2[3], double pt3[3], double pt4[3])
{
    getPtsInVector(CA->pnts_old, CA->ini[i]->iPt1, pt1);
    getPtsInVector(CA->pnts_old, CA->ini[i]->iPt2, pt2);
    getPtsInVector(CA->pnts_old, CA->ini[i]->iPt3, pt3);
    getPtsInVector(CA->pnts_old, CA->ini[i]->iPt4, pt4);
    
}
/**
 * @param el
 * @param POINTS_OLD
 * @param reta
 * @param n1
 * @param n2
 * @param n3
 * @param n4
 * @param iTipoFibra
 */
void findIntersectionPlaneLine(int index, typ_ca *CA, 
        double reta[3][2], 
        double n[4][3],
        int iTipoFibra)
{
    
    double pt1[3];
    double pt2[3];
    double pt3[3];
    double pt4[3];
    getElementPtsInVector(index, CA, pt1, pt2, pt3, pt4);
    int check=0;
    double ptI1[3], ptI2[3], ptI3[3], ptI4[3];
    double ptIValid1[3] = { 0., 0., 0. }, ptIValid2[3] = { 0., 0., 0. };
    ptIValid1[0]=ptIValid2[0]=0.0;
    ptIValid1[1]=ptIValid2[1]=0.0;
    ptIValid1[2]=ptIValid2[2]=0.0;
    
    //returning interpolation function
    double ksi, etha, relError=0.0;
    //index for insert the interpolation function on the proper 
    int col;
    int smallestErrorID=-1;
    double smallestError=DBL_MAX;
    //face 1: 123
    if(checkPtBelongsToFace(reta, pt1, pt2, pt3, n[0], ptI1, &ksi, &etha, &relError)==1)
    {
        if(addPoint(ptIValid1, ptIValid2, ptI1, &check))
        {
            col = iTipoFibra+check-1;
            CA->ini[index]->iFaces[col]= 0;
            CA->ini[index]->ck[0][col] = getN1(etha, ksi);
            CA->ini[index]->ck[1][col] = getN2(etha, ksi);
            CA->ini[index]->ck[2][col] = getN3(etha, ksi);
            CA->ini[index]->ck[3][col] = 0.0f;
        }
    }else{
        //o ponto não pertece à face. Guarda o índice 
        if(relError<smallestError){
            smallestError   = relError;
            smallestErrorID = 0;
        }
    }
    //face 2: 124
    if(checkPtBelongsToFace(reta, pt1, pt2, pt4, n[1], ptI2, &ksi, &etha, &relError)==1)
    {
        if(addPoint(ptIValid1, ptIValid2, ptI2, &check))
        {
            col = iTipoFibra+check-1;
            CA->ini[index]->iFaces[col]= 1;
            CA->ini[index]->ck[0][col] = getN1(etha, ksi);
            CA->ini[index]->ck[1][col] = getN2(etha, ksi);
            CA->ini[index]->ck[2][col] = 0.0f;
            CA->ini[index]->ck[3][col] = getN3(etha, ksi);
        }
    }else{
        //o ponto não pertece à face. Guarda o índice 
        if(relError<smallestError){
            smallestError   = relError;
            smallestErrorID = 1;
        }
    }
    //face 3: 134
    if(checkPtBelongsToFace(reta, pt1, pt3, pt4, n[2], ptI3, &ksi, &etha, &relError)==1)
    {
        if(addPoint(ptIValid1, ptIValid2, ptI3, &check))
        {
            col = iTipoFibra+check-1;
            CA->ini[index]->iFaces[col]= 2;
            CA->ini[index]->ck[0][col] = getN1(etha, ksi);
            CA->ini[index]->ck[1][col] = 0.0f;
            CA->ini[index]->ck[2][col] = getN2(etha, ksi);
            CA->ini[index]->ck[3][col] = getN3(etha, ksi);
        }
    }else{
        //o ponto não pertece à face. Guarda o índice 
        if(relError<smallestError){
            smallestError   = relError;
            smallestErrorID = 2;
        }
    }
    //face 4: 234
    if(checkPtBelongsToFace(reta, pt2, pt3, pt4, n[3], ptI4, &ksi, &etha, &relError)==1)
    {
        if(addPoint(ptIValid1, ptIValid2, ptI4, &check))
        {
            col = iTipoFibra+check-1;
            CA->ini[index]->iFaces[col]= 3;
            CA->ini[index]->ck[0][col] = 0.0f;
            CA->ini[index]->ck[1][col] = getN1(etha, ksi);
            CA->ini[index]->ck[2][col] = getN2(etha, ksi);
            CA->ini[index]->ck[3][col] = getN3(etha, ksi);
        }
    }else{
        //o ponto não pertece à face. Guarda o índice 
        if(relError<smallestError){
            smallestError   = relError;
            smallestErrorID = 3;
        }
    }
    //verifica se não há exatamente dois pontos encontrados
    if(check==1){//se há um ponto, considera como outra interseção o ponto com menor erro relativo
        if(smallestErrorID==0){
            checkPtBelongsToFace(reta, pt1, pt2, pt3, n[0], ptI1, &ksi, &etha, &relError);
            if(addPoint(ptIValid1, ptIValid2, ptI1, &check))
            {
                col = iTipoFibra+check-1;
                CA->ini[index]->iFaces[col]= 0;
                CA->ini[index]->ck[0][col] = getN1(etha, ksi);
                CA->ini[index]->ck[1][col] = getN2(etha, ksi);
                CA->ini[index]->ck[2][col] = getN3(etha, ksi);
                CA->ini[index]->ck[3][col] = 0.0f;
            }
        }else if(smallestErrorID==1){
            checkPtBelongsToFace(reta, pt1, pt2, pt4, n[1], ptI2, &ksi, &etha, &relError);
            if(addPoint(ptIValid1, ptIValid2, ptI2, &check))
            {
                col = iTipoFibra+check-1;
                CA->ini[index]->iFaces[col]= 1;
                CA->ini[index]->ck[0][col] = getN1(etha, ksi);
                CA->ini[index]->ck[1][col] = getN2(etha, ksi);
                CA->ini[index]->ck[2][col] = 0.0f;
                CA->ini[index]->ck[3][col] = getN3(etha, ksi);
            }
        }else if(smallestErrorID==2){
            checkPtBelongsToFace(reta, pt1, pt3, pt4, n[2], ptI3, &ksi, &etha, &relError);
            if(addPoint(ptIValid1, ptIValid2, ptI3, &check))
            {
                col = iTipoFibra+check-1;
                CA->ini[index]->iFaces[col]= 2;
                CA->ini[index]->ck[0][col] = getN1(etha, ksi);
                CA->ini[index]->ck[1][col] = 0.0f;
                CA->ini[index]->ck[2][col] = getN2(etha, ksi);
                CA->ini[index]->ck[3][col] = getN3(etha, ksi);
            }
        }else if(smallestErrorID==3){
            checkPtBelongsToFace(reta, pt2, pt3, pt4, n[3], ptI4, &ksi, &etha, &relError);
            if(addPoint(ptIValid1, ptIValid2, ptI4, &check))
            {
                col = iTipoFibra+check-1;
                CA->ini[index]->iFaces[col]= 3;
                CA->ini[index]->ck[0][col] = 0.0f;
                CA->ini[index]->ck[1][col] = getN1(etha, ksi);
                CA->ini[index]->ck[2][col] = getN2(etha, ksi);
                CA->ini[index]->ck[3][col] = getN3(etha, ksi);
            }
        }
    }
    if(check!=2){
        stringstream ss;
        ss << "PAUUUUUU CHECK: "<<check;
        ss << "\npt1\tpt2\tpt3\tpt4\tptI1\tptI2\n";
        for(int i=0;i<3;i++)
        {
            ss<<pt1[i]<<"\t"<<pt2[i]<<"\t"<<pt3[i]<<"\t"<<pt4[i]<<"\t"; 
            ss<<ptIValid1[i]<<"\t"<<ptIValid2[i]<<endl;
        }
        string str = ss.str();
        if(CA->params->printOutput==1){
            throw MyException(str, __FILE__, __LINE__);
        }
    }
}
/**
 * 
 * @param POINTS_OLD
 * @param i
 * @param normals
 */
void findAllNormals(typ_ca *CA, int i, double normals[4][3])
{
    
    //finds the normal of the 4 faces
    //face 1 points 123   
    findFaceNormal(CA, CA->ini[i]->iPt1, CA->ini[i]->iPt2, CA->ini[i]->iPt3, normals[0]);
    //face 2 points 124
    findFaceNormal(CA, CA->ini[i]->iPt1, CA->ini[i]->iPt2, CA->ini[i]->iPt4, normals[1]);
    //face 3 points 134
    findFaceNormal(CA, CA->ini[i]->iPt1, CA->ini[i]->iPt3, CA->ini[i]->iPt4, normals[2]);
    //face 4 points 234
    findFaceNormal(CA, CA->ini[i]->iPt2, CA->ini[i]->iPt3, CA->ini[i]->iPt4, normals[3]);
}   

/**
 * 
 * @param i
 * @param j
 * @param k
 */
void swapInterpolationCols(int i, int j, int k, typ_ca *CA)
{
    int aux1=CA->ini[i]->iFaces[k];
    CA->ini[i]->iFaces[k] = CA->ini[i]->iFaces[j];
    CA->ini[i]->iFaces[j]=aux1;
    //
    for(int ii=0;ii<4;ii++){
        double aux=CA->ini[i]->ck[ii][k];
        CA->ini[i]->ck[ii][k] = CA->ini[i]->ck[ii][j];
        CA->ini[i]->ck[ii][j]=aux;
    }    
}
/**
 * 
 * @param oldAx
 * @param newAx
 * @return 
 */
 int checkAxis(double oldAx[3], double newAx[3]){
    double absEr, relEr;
    for(int j=0;j<3;j++){        
        absEr   = fabs(oldAx[j]-newAx[j]);
        relEr   = absEr/fabs(oldAx[j]);
        if(absEr>ABSTOL && relEr>RELTOL && fabs(oldAx[j])>TOLERANCE 
                &&fabs(oldAx[j])/oldAx[j]!=fabs(newAx[j])/newAx[j] )
            return 1;
    }
    return 0;

 }
 /**
  * 
  * @param 
  * @param face
  * @return 
  */
 int checkFaceByPnt(int pt, typ_face *f)
 {
    if(f->pt1 == pt || f->pt2 == pt || f->pt3 == pt){
        return 1;
    } 
    return 0;
 }
 //Verifica se a face de pressao pertence a face 1: 123
 int checkFace1(int i, typ_ca *CA, typ_face *face){
     if(checkFaceByPnt(CA->ini[i]->iPt1, face)==0)
         return 0;
     if(checkFaceByPnt(CA->ini[i]->iPt2, face)==0)
         return 0;
     if(checkFaceByPnt(CA->ini[i]->iPt3, face)==0)
         return 0;
    return 1;
 }
 //Verifica se a face de pressao pertence a face 3: 124
 int checkFace2(int i, typ_ca *CA, typ_face *face){
     if(checkFaceByPnt(CA->ini[i]->iPt1, face)==0)
         return 0;
     if(checkFaceByPnt(CA->ini[i]->iPt2, face)==0)
         return 0;
     if(checkFaceByPnt(CA->ini[i]->iPt4, face)==0)
         return 0;
    return 1;
 }
  //Verifica se a face de pressao pertence a face 3: 134
 int checkFace3(int i, typ_ca *CA, typ_face *face){
     if(checkFaceByPnt(CA->ini[i]->iPt1, face)==0)
         return 0;
     if(checkFaceByPnt(CA->ini[i]->iPt3, face)==0)
         return 0;
     if(checkFaceByPnt(CA->ini[i]->iPt4, face)==0)
         return 0;
    return 1;
 }
 //Verifica se a face de pressao pertence a face 3: 234
 int checkFace4(int i, typ_ca *CA, typ_face *face){
     if(checkFaceByPnt(CA->ini[i]->iPt2, face)==0)
         return 0;
     if(checkFaceByPnt(CA->ini[i]->iPt3, face)==0)
         return 0;
     if(checkFaceByPnt(CA->ini[i]->iPt4, face)==0)
         return 0;
    return 1;
 }
/**
 * 
 * @param i
 * @param CA
 */
void interceptionPntPressure(int i, typ_ca *CA){
    for(int iF=0; iF < CA->params->numFaces; iF++){
        typ_face *face = CA->params->aFaces[iF];
        CA->ini[i]->pressFaces[0] = checkFace1(i, CA, face);
        CA->ini[i]->pressFaces[1] = checkFace2(i, CA, face);
        CA->ini[i]->pressFaces[2] = checkFace3(i, CA, face);
        CA->ini[i]->pressFaces[3] = checkFace4(i, CA, face);
        int soma=0;
        for(int j=0;j<4;j++){
            soma+=CA->ini[i]->pressFaces[j];
            if(CA->ini[i]->pressFaces[j]==1){
               CA->ini[i]->ipressFaces[j]=iF;
               CA->ini[i]->hasPressure=1;
            }else{
               CA->ini[i]->ipressFaces[j]=-1;
            }
        }
        if(CA->ini[i]->hasPressure==1) break;        
    }
}
void iniVolumes(int iElem, typ_ca *CA)
{
    //computes the volume
    CA->ini[iElem]->volCel_ini = CA->t_new[iElem]->volCel  = 
        CA->t_old[iElem]->volCel = getVolumeTetrahedron(CA, iElem);

}
 
void iniGeometry(int i, typ_ca *CA)
{   
    //Finds the barycenter
    getBarycenter(CA, i);
    CA->ini[i]->bary_ini[0] = CA->t_old[i]->bary[0] = CA->t_new[i]->bary[0];
    CA->ini[i]->bary_ini[1] = CA->t_old[i]->bary[1] = CA->t_new[i]->bary[1];
    CA->ini[i]->bary_ini[2] = CA->t_old[i]->bary[2] = CA->t_new[i]->bary[2];
    //Tt 0 1 2 3 1
    /*cout<<"Tt "<<CA->ini[i]->iPt1<<" "<<CA->ini[i]->iPt2;
    cout<<" "<<CA->ini[i]->iPt3<<" "<<CA->ini[i]->iPt4<<" ";
    if(CA->ini[i]->bary_ini[0]<3.0){
        cout<<"2";
    }else if(CA->ini[i]->bary_ini[0]>7.0){
        cout<<"3";
    }else{
        cout<<"1";
    }
    cout<<endl;*/
    //Finds the line parallel to an axis that pass by the barycenter
    double retaParalFibra[3][2], retaParalSheet[3][2], retaParalNSheet[3][2];
    findLineAxisBary(CA->t_old[i]->bary, CA->t_old[i]->fiberDir,  retaParalFibra);
    findLineAxisBary(CA->t_old[i]->bary, CA->t_old[i]->sheetDir,  retaParalSheet);
    findLineAxisBary(CA->t_old[i]->bary, CA->t_old[i]->nsheetDir, retaParalNSheet);
    double n[4][3] = { 0. };
    for(int ii=0;ii<4;ii++){
        n[ii][0]=n[ii][1]=n[ii][2]=0.0;
    }
    findAllNormals(CA, i, n);
    findIntersectionPlaneLine(i, CA, retaParalFibra,  n, 0);
    findIntersectionPlaneLine(i, CA, retaParalSheet,  n, 2);
    findIntersectionPlaneLine(i, CA, retaParalNSheet, n, 4);
    //find the intersection points and update the fiber directions
    updateAxisDirectionOLD(CA, i);
    // se alguma direação foi invertida. Por exemplo, a fibra era pra ser [1 0 0] e a interpolação encontrou [-1 0 0]
    //para o modelo mecanico nao faz diferença. apenas interfere na visualização das fibras no paraview
    int checkSwap=0;
    if( checkAxis(CA->t_old[i]->fiberDir, CA->t_new[i]->fiberDir) ==1)
    {
        swapInterpolationCols(i, 0, 1, CA);
        checkSwap++;
    }
    if( checkAxis(CA->t_old[i]->sheetDir, CA->t_new[i]->sheetDir) ==1)
    {
        swapInterpolationCols(i, 2, 3, CA);
        checkSwap++;
    }
    if( checkAxis(CA->t_old[i]->nsheetDir, CA->t_new[i]->nsheetDir) ==1 )
    {
        swapInterpolationCols(i, 4, 5, CA);
        checkSwap++;
    }
    //se houve troca, recomputa os eixos
    if(checkSwap!=0)    updateAxisDirectionOLD(CA, i);
    //
    CA->ini[i]->axisLength0[FIBER] = CA->t_old[i]->axisLengthT[FIBER] = CA->t_new[i]->axisLengthT[FIBER];
    CA->t_old[i]->fiberDir[0]    = CA->t_new[i]->fiberDir[0];
    CA->t_old[i]->fiberDir[1]    = CA->t_new[i]->fiberDir[1];
    CA->t_old[i]->fiberDir[2]    = CA->t_new[i]->fiberDir[2];    
    //
    CA->ini[i]->axisLength0[SHEET] = CA->t_old[i]->axisLengthT[SHEET] = CA->t_new[i]->axisLengthT[SHEET];
    CA->t_old[i]->sheetDir[0]    = CA->t_new[i]->sheetDir[0];
    CA->t_old[i]->sheetDir[1]    = CA->t_new[i]->sheetDir[1];
    CA->t_old[i]->sheetDir[2]    = CA->t_new[i]->sheetDir[2];
    //
    CA->ini[i]->axisLength0[NORMAL] = CA->t_old[i]->axisLengthT[NORMAL] = CA->t_new[i]->axisLengthT[NORMAL];
    CA->t_old[i]->nsheetDir[0]    = CA->t_new[i]->nsheetDir[0];
    CA->t_old[i]->nsheetDir[1]    = CA->t_new[i]->nsheetDir[1];
    CA->t_old[i]->nsheetDir[2]    = CA->t_new[i]->nsheetDir[2];
    /*if(i==0){
        cout<<  CA->ini[i]->axisLength0[FIBER]<<" "
            <<  CA->ini[i]->axisLength0[SHEET]<<" "
            << CA->ini[i]->axisLength0[NORMAL]<<endl;    
        debugpnts(CA, i);
    }*/
    CA->ini[i]->areaHexFSTIni = CA->ini[i]->axisLength0[FIBER]*CA->ini[i]->axisLength0[SHEET];
    CA->ini[i]->areaHexFNTIni = CA->ini[i]->axisLength0[FIBER]*CA->ini[i]->axisLength0[NORMAL];
    CA->ini[i]->areaHexSNTIni = CA->ini[i]->axisLength0[SHEET]*CA->ini[i]->axisLength0[NORMAL];
    
    CA->t_old[i]->areaHexFST = CA->t_new[i]->areaHexFST = CA->ini[i]->areaHexFSTIni;
    CA->t_old[i]->areaHexFNT = CA->t_new[i]->areaHexFNT = CA->ini[i]->areaHexFNTIni;
    CA->t_old[i]->areaHexSNT = CA->t_new[i]->areaHexSNT = CA->ini[i]->areaHexSNTIni;
    
    for(int j=0;j<3;j++){
        CA->ini[i]->fiberDirIni[j]   = CA->t_new[i]->fiberDir[j];
        CA->ini[i]->sheetDirIni[j]   = CA->t_new[i]->sheetDir[j];
        CA->ini[i]->nsheetDirIni[j]  = CA->t_new[i]->nsheetDir[j];
    }
    //
   /* CA->t_old[i]->alphaT_12 = CA->ini[i]->alpha0_12 = 0.0;
    CA->t_old[i]->alphaT_13 = CA->ini[i]->alpha0_13 = 0.0;
    CA->t_old[i]->alphaT_23 = CA->ini[i]->alpha0_23 = 0.0;*/
    updateAlpha( i, CA);
    CA->t_old[i]->alphaT_12 = CA->ini[i]->alpha0_12 = CA->t_new[i]->alphaT_12;
    CA->t_old[i]->alphaT_13 = CA->ini[i]->alpha0_13 = CA->t_new[i]->alphaT_13;
    CA->t_old[i]->alphaT_23 = CA->ini[i]->alpha0_23 = CA->t_new[i]->alphaT_23;
    
    CA->t_old[i]->deltaAlpha_12 = CA->t_new[i]->deltaAlpha_12=0.0;
    CA->t_old[i]->deltaAlpha_13 = CA->t_new[i]->deltaAlpha_13=0.0;
    CA->t_old[i]->deltaAlpha_23 = CA->t_new[i]->deltaAlpha_23=0.0;
    
    //
    getInterceptPtsByInterpolation(CA, i);
    for(int jj=0;jj<6;jj++){
        if(CA->ini[i]->iFaces[jj]<0 || CA->ini[i]->iFaces[jj]>4)
        {
            stringstream ss;
            ss<<"Area index ["<<CA->ini[i]->iFaces[jj]<<"] does not exist."<<endl;
            string str = ss.str();
            if(CA->params->printOutput==1) {
                throw MyException(str, __FILE__, __LINE__);
            }
        }
        CA->t_old[i]->intPts[0][jj]= CA->t_new[i]->intPts[0][jj];
        CA->t_old[i]->intPts[1][jj]= CA->t_new[i]->intPts[1][jj];
        CA->t_old[i]->intPts[2][jj]= CA->t_new[i]->intPts[2][jj];
    }
    computeKs(CA, i);
    CA->ini[i]->hasPressure=0;
    interceptionPntPressure(i,CA);
}
void getIntercMasses(typ_ca *CA, int iElem){
    double pntMass[4] = { 0., 0., 0. ,0.};
    
    int iPt1 = CA->ini[iElem]->iPt1;
    int iPt2 = CA->ini[iElem]->iPt2;
    int iPt3 = CA->ini[iElem]->iPt3;
    int iPt4 = CA->ini[iElem]->iPt4;
    pntMass[0] = CA->pnts_old[iPt1].mass;
    pntMass[1] = CA->pnts_old[iPt2].mass;
    pntMass[2] = CA->pnts_old[iPt3].mass;
    pntMass[3] = CA->pnts_old[iPt4].mass;
    
    double mass=0.0;
    for (int d = 0 ; d < 6 ; d++ ){
        for (int k = 0 ; k < 4 ; k++ ){                
            mass = mass + pntMass[k]*CA->ini[iElem]->ck[k][d];
        }
        CA->ini[iElem]->intercMass[d] = mass;
        mass = 0.0;
    }
}
/**
 * 
 * @param CA
 * @param iElem
 */
void computeKs(typ_ca *CA, int iElem){
    getIntercMasses(CA, iElem);
    double area_L0[3] = { 0., 0., 0. }, area[3] = { 0., 0., 0. };
    typ_t0_element *ini = CA->ini[iElem];
    //calcula o damp de um jeito novo para tentar amortecer o cisalhamento
    //fibra
    area_L0[0] = ini->areaHexSNTIni / ini->axisLength0[0];
    area[0]    = ini->areaHexSNTIni;
    //sheet
    area_L0[1] = ini->areaHexFNTIni / ini->axisLength0[1];
    area[1]    = ini->areaHexFNTIni;
    //normal
    area_L0[2] = ini->areaHexFSTIni / ini->axisLength0[2];
    area[2]    = ini->areaHexFSTIni;
    
    t_par_ac* ap = CA->params->aParam[ini->iRegion];
    for(int i =0;i<3;i++){
        ini->KAxl[i] = area_L0[i]*ap->EAxl[i];
        ini->KVol[i] = area_L0[i]*ap->EVol[i];
        ini->KAng[i] = area[i]*ap->EAng[i];
        ini->damp[i] = area[i]*ap->kDamp;
    }
}
  
 /**
  * 
  * @param inData
  * @param ck
  * @param dataOut
  */
 void findDataByInterpol(double inData[3][4], double ck[4][6], double dataOut[3][6])
 {
    double sum=0.0;
    for (int c = 0 ; c < 3 ; c++ )
    {
        for (int d = 0 ; d < 6 ; d++ )
        {
            for (int k = 0 ; k < 4 ; k++ )
            {
                sum = sum + inData[c][k]*ck[k][d];
            }
            dataOut[c][d] = sum;
            sum = 0.0;
        }
    }
}
 void debugpnts(typ_ca *CA, int i){
    int iPt1=CA->ini[i]->iPt1;
    int iPt2=CA->ini[i]->iPt2;
    int iPt3=CA->ini[i]->iPt3;
    int iPt4=CA->ini[i]->iPt4;
    cout<<"ponto: "<<iPt1<<" ";
    cout<<CA->pnts_old[iPt1].x<<" ";
    cout<<CA->pnts_old[iPt1].y<<" ";
    cout<<CA->pnts_old[iPt1].z<<endl;
    cout<<"ponto: "<<iPt2<<" ";
    cout<<CA->pnts_old[iPt2].x<<" ";
    cout<<CA->pnts_old[iPt2].y<<" ";
    cout<<CA->pnts_old[iPt2].z<<endl;
    cout<<"ponto: "<<iPt3<< " ";
    cout<<CA->pnts_old[iPt3].x<<" ";
    cout<<CA->pnts_old[iPt3].y<<" ";
    cout<<CA->pnts_old[iPt3].z<<endl;
    cout<<"ponto: "<<iPt4<< " ";
    cout<<CA->pnts_old[iPt4].x<<" ";
    cout<<CA->pnts_old[iPt4].y<<" ";
    cout<<CA->pnts_old[iPt4].z<<endl;
    cout<<"============="<<endl;
    for(int j=0; j<3; j++)
    {
        printf("%.e :: %e = %e - %e\n",  
               CA->t_new[i]->sheetDir[j] , 
               CA->t_new[i]->intPts[j][3]-CA->t_new[i]->intPts[j][2], 
               CA->t_new[i]->intPts[j][3], 
               CA->t_new[i]->intPts[j][2]);
    }
    cout<<"============="<<endl;
    for (int d = 0 ; d < 6 ; d++ )
    {
        for (int k = 0 ; k < 4 ; k++ )
        {
            printf("%.e ", CA->ini[i]->ck[k][d]);
        }
        printf("\n");
    }
    printf("ifaces: \n");
    for(int jj=0;jj<6;jj++){
        printf("%d ", CA->ini[i]->iFaces[jj]);
    }
    printf("\n intPts: \n");
    for(int jj=0;jj<6;jj++){
        printf("%e ",  CA->t_new[i]->intPts[0][jj]);
        printf("%e ",  CA->t_new[i]->intPts[1][jj]);
        printf("%e\n", CA->t_new[i]->intPts[2][jj]);
    }
    cout<<"============="<<endl;
    //getchar();  
 }
 /**
  * 
  * @param POINTS_OLD
  * @param elements_old
  * @param i
  * @param interPts
  */
 void getInterceptPtsByInterpolation(typ_ca *CA, int i)
 {
    double points[3][4] = { 0. };
    int iPt1=CA->ini[i]->iPt1;
    int iPt2=CA->ini[i]->iPt2;
    int iPt3=CA->ini[i]->iPt3;
    int iPt4=CA->ini[i]->iPt4;
    points[0][0]    = CA->pnts_old[iPt1].x;
    points[1][0]    = CA->pnts_old[iPt1].y;
    points[2][0]    = CA->pnts_old[iPt1].z;
    points[0][1]    = CA->pnts_old[iPt2].x;
    points[1][1]    = CA->pnts_old[iPt2].y;
    points[2][1]    = CA->pnts_old[iPt2].z;
    points[0][2]    = CA->pnts_old[iPt3].x;
    points[1][2]    = CA->pnts_old[iPt3].y;
    points[2][2]    = CA->pnts_old[iPt3].z;
    points[0][3]    = CA->pnts_old[iPt4].x;
    points[1][3]    = CA->pnts_old[iPt4].y;
    points[2][3]    = CA->pnts_old[iPt4].z;
    //    
    findDataByInterpol(points, CA->ini[i]->ck, CA->t_new[i]->intPts);
 }
 /**
  * 
  * @param POINTS_OLD
  * @param elements_old
  * @param i
  * @param axis1
  * @param axis2
  * @param axis3
  */
void findAxis(typ_ca *CA, int i)
{
    //find the intercepetion points
    getInterceptPtsByInterpolation(CA, i);
    //now find the directions
    for(int j=0; j<3; j++)
    {
        CA->t_new[i]->fiberDir[j] = CA->t_new[i]->intPts[j][1] - CA->t_new[i]->intPts[j][0];
        CA->t_new[i]->sheetDir[j] = CA->t_new[i]->intPts[j][3] - CA->t_new[i]->intPts[j][2];
        CA->t_new[i]->nsheetDir[j]= CA->t_new[i]->intPts[j][5] - CA->t_new[i]->intPts[j][4];
    }
    
 }

 
 void findVelocities(typ_ca *CA,  int iElem, 
         double velsIter[6][3])
{
    double velsIterTEMP[3][6] = { 0. };
    double velPts[3][4] = { 0. };
    
    int iPt1=CA->ini[iElem]->iPt1;
    int iPt2=CA->ini[iElem]->iPt2;
    int iPt3=CA->ini[iElem]->iPt3;
    int iPt4=CA->ini[iElem]->iPt4;
    
    velPts[0][0]    = CA->pnts_old[iPt1].xV;
    velPts[1][0]    = CA->pnts_old[iPt1].yV;
    velPts[2][0]    = CA->pnts_old[iPt1].zV;
    velPts[0][1]    = CA->pnts_old[iPt2].xV;
    velPts[1][1]    = CA->pnts_old[iPt2].yV;
    velPts[2][1]    = CA->pnts_old[iPt2].zV;
    velPts[0][2]    = CA->pnts_old[iPt3].xV;
    velPts[1][2]    = CA->pnts_old[iPt3].yV;
    velPts[2][2]    = CA->pnts_old[iPt3].zV;
    velPts[0][3]    = CA->pnts_old[iPt4].xV;
    velPts[1][3]    = CA->pnts_old[iPt4].yV;
    velPts[2][3]    = CA->pnts_old[iPt4].zV;
    if(velsIterTEMP==NULL)
        throw MyException("velsIterTEMP", __FILE__, __LINE__);
    if(velPts==NULL)
        throw MyException("velPts", __FILE__, __LINE__);
    if(CA->ini[iElem]->ck==NULL)
        throw MyException("CA->ini[iElem]->ck", __FILE__, __LINE__);
    //////// TODO FIXME nao sei porque este codigo estava replicado em vez de chamar a funcao
    findDataByInterpol(velPts, CA->ini[iElem]->ck, velsIterTEMP);
    /*double sum = 0.0;
    for (int c = 0 ; c < 3 ; c++ )
    {
        for (int d = 0 ; d < 6 ; d++ )
        {
            for (int k = 0 ; k < 4 ; k++ )
            {
                sum = sum + velPts[c][k]*CA->ini[iElem]->ck[k][d];
            }
            velsIterTEMP[c][d] = sum;
            sum = 0.0;
        }
    }*/
    //tranpose
    for(int ii=0;ii<3;ii++){
        for(int jj=0;jj<6;jj++){
            velsIter[jj][ii] = velsIterTEMP[ii][jj];
        }
    }

 }

 /**
 * 
 * @param POINTS_OLD
 * @param elements_old
 * @param elements_new
 */
void updateAxisDirectionOLD(typ_ca *CA, int i)
{
    findAxis(CA, i);
    CA->t_new[i]->axisLengthT[FIBER] = normalizeVector(CA->t_new[i]->fiberDir, i);
    CA->t_new[i]->axisLengthT[SHEET] = normalizeVector(CA->t_new[i]->sheetDir, i);
    CA->t_new[i]->axisLengthT[NORMAL]= normalizeVector(CA->t_new[i]->nsheetDir, i);
}

/**
 * 
 * @param el_new
 * @param el_old
 * @param i
 */
void updateAlpha(int i, typ_ca *CA){    
    
    CA->t_new[i]->alphaT_12 = dot(CA->t_old[i]->fiberDir, CA->t_old[i]->sheetDir);
    CA->t_new[i]->alphaT_13 = dot(CA->t_old[i]->fiberDir, CA->t_old[i]->nsheetDir);
    CA->t_new[i]->alphaT_23 = dot(CA->t_old[i]->sheetDir, CA->t_old[i]->nsheetDir);
    
    CA->t_new[i]->deltaAlpha_12 = CA->t_old[i]->alphaT_12 - CA->ini[i]->alpha0_12;
    CA->t_new[i]->deltaAlpha_13 = CA->t_old[i]->alphaT_13 - CA->ini[i]->alpha0_13;
    CA->t_new[i]->deltaAlpha_23 = CA->t_old[i]->alphaT_23 - CA->ini[i]->alpha0_23;
    
}

        
/**
 * 
 * @param elements_old
 * @param elements_new
 * @param i
 * @param POINTS_OLD
 * @param POINTS_NEW
 */
void updateGeometry(
    int i, typ_ca *CA)
{
    CA->t_new[i]->volCel = getVolumeTetrahedron( CA, i);    
    getBarycenter(CA, i);
    updateAxisDirectionOLD(CA, i);
    updateAlpha(i, CA);
    
    //atualiza areas
    CA->t_new[i]->areaHexFST = CA->t_old[i]->axisLengthT[FIBER]*CA->t_old[i]->axisLengthT[SHEET];
    CA->t_new[i]->areaHexFNT = CA->t_old[i]->axisLengthT[FIBER]*CA->t_old[i]->axisLengthT[NORMAL];
    CA->t_new[i]->areaHexSNT = CA->t_old[i]->axisLengthT[SHEET]*CA->t_old[i]->axisLengthT[NORMAL];
}
/**
 * 
 * @param inA
 * @param inB
 * @param out
 */
void getProjection(double inA[3], double inB[3], double out[3], double *dotAB)
{
    double bHat[3] = { 0., 0., 0. };
    bHat[0]=inB[0]; bHat[1]=inB[1]; bHat[2]=inB[2];
    normalizeVector(bHat, -1);
    *dotAB = dot(inA, bHat);
    out[0]=(*dotAB)*bHat[0];
    out[1]=(*dotAB)*bHat[1];
    out[2]=(*dotAB)*bHat[2];
}

/**
 * Esta função corrige a direção da normal. 
 * Faz com que a direção da normal seja para dentro do ventriculo, cujo ponto central é 
 * arbitrariamente determinado em [0 0 5].
 * Encontra um vetor entre o ponto central do triângulo e o ponto de referencia no centro do ventriculo.
 * Então encontra o produto escalar entre este vetor e a normal. Se for positivo, inverte a normal.
 * @param vPt1
 * @param vPt2
 * @param vPt3
 * @param normal
 */
void fixNormal(double vPt1[3], double vPt2[3], double vPt3[3], double normal[3], double bary[3]){
    
    double vRef[3]={0.0};
    //ponto de referencia central
    double refPnt[3] = {0.0, 0.0, 0.0};
    for(int i=0; i<3;i++){
        //media dos pontos
        bary[i]       = (vPt1[i]+vPt2[i]+vPt3[i])/3.0;
        //acha o vetor
        vRef[i] = bary[i] - refPnt[i];
    }
    //acha o dot product    
    double dotP = dot(vRef, normal);
    //se o produto escalar for positivo, inverte o sentido da normal
    if(dotP>0){
        normal[0] = -normal[0];
        normal[1] = -normal[1];
        normal[2] = -normal[2];
    }
}

void getFaceCenter(double vPt1[3], double vPt2[3], double vPt3[3], double bary[3]){
    for(int i=0; i<3;i++){
        //media dos pontos
        bary[i]       = (vPt1[i]+vPt2[i]+vPt3[i])/3.0;
    }    
}