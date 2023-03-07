#include "ContinuumMech.h"

 /**
  * Assembly the matrices defining the tetrahedron position by three vectors:
  * v1 = p2-p1
  * v2 = p3-p1
  * v3 = p4-p1
  * @param points
  * @param elements
  * @param X
  * @param i
  */  
 void assemblyXMatrix(typ_point *points, double X[3][3], int i, typ_ca *CA)
 {
    int iPt1=CA->ini[i].iPt1;
    int iPt2=CA->ini[i].iPt2;
    int iPt3=CA->ini[i].iPt3;
    int iPt4=CA->ini[i].iPt4;
    //first column
    X[0][0] = points[iPt2].x - points[iPt1].x;
    X[1][0] = points[iPt2].y - points[iPt1].y;
    X[2][0] = points[iPt2].z - points[iPt1].z;
    //Second column
    X[0][1] = points[iPt3].x - points[iPt1].x;
    X[1][1] = points[iPt3].y - points[iPt1].y;
    X[2][1] = points[iPt3].z - points[iPt1].z;
    //Third column
    X[0][2] = points[iPt4].x - points[iPt1].x;
    X[1][2] = points[iPt4].y - points[iPt1].y;
    X[2][2] = points[iPt4].z - points[iPt1].z;
 }
 /**
  * Finds the deformation tensor F with the initial and final configurations X0 and X1:
  * X1=F*X0,
  * F=X1*X0^-1
  * @param X0
  * @param X1
  * @param F
  * @return 
  */
 void getDeformationTensor(double X0[3][3], double X1[3][3], double F[3][3])
 {
    double X0Inv[3][3];
    invertMatrix(X0, X0Inv);
    matrixMultiplication(X1, X0Inv, F);
    
   
 }
 
 
/**
 * Polar decomposition:
 * Finds rotation and stretch tensors (R and) U from F.
 * @param F
 * @param R
 * @param U
 */
 void computePolarDecompositition(double F[3][3], double R[3][3], double U[3][3])
 {
    //U^2 = F'*F
    double FT[3][3], U2[3][3], UI[3][3];
    transpose(F, FT);
    matrixMultiplication(FT, F, U2);
    //U=sqrt(U^2)
    computeSqrtMatrix(U2, U);
    //invert U
    invertMatrix(U, UI);
    //R=FU^-1
    matrixMultiplication(F, UI, R);
 }
 
 /**
  * 
  * @param X0
  * @param X1
  * @param F
  * @param R
  * @param U
  */
 void getTensors(double X0[3][3], double X1[3][3], double F[3][3], 
         double R[3][3], double U[3][3] ){
    
    getDeformationTensor(X0, X1, F);
    computePolarDecompositition(F, R, U);
}
