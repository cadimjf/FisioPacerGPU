#include <string>

#include "LinearAlg.h"
//http://arxiv.org/abs/physics/0610206
//http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html

// Macros
#define SQR(x)      ((x)*(x))                        // x^2 


// ----------------------------------------------------------------------------
inline void HouseHolder(double A[3][3], double Q[3][3], double d[3], double e[2])
// ----------------------------------------------------------------------------
// Reduces a symmetric 3x3 matrix to tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// ---------------------------------------------------------------------------
{
  const int n = 3;
  double u[n], q[n];
  double omega, f;
  double K, h, g;
  
  // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
  for (int i=0; i < n; i++)
  {
    Q[i][i] = 1.0;
    for (int j=0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }
#endif

  // Bring first row and column to the desired form 
  h = SQR(A[0][1]) + SQR(A[0][2]);
  if (A[0][1] > 0)
    g = -sqrt(h);
  else
    g = sqrt(h);
  e[0] = g;
  f    = g * A[0][1];
  u[1] = A[0][1] - g;
  u[2] = A[0][2];
  
  omega = h - f;
  if (omega > 0.0)
  {
    omega = 1.0 / omega;
    K     = 0.0;
    for (int i=1; i < n; i++)
    {
      f    = A[1][i] * u[1] + A[i][2] * u[2];
      q[i] = omega * f;                  // p
      K   += u[i] * f;                   // u* A u
    }
    K *= 0.5 * SQR(omega);

    for (int i=1; i < n; i++)
      q[i] = q[i] - K * u[i];
    
    d[0] = A[0][0];
    d[1] = A[1][1] - 2.0*q[1]*u[1];
    d[2] = A[2][2] - 2.0*q[2]*u[2];
    
    // Store inverse Householder transformation in Q
#ifndef EVALS_ONLY
    for (int j=1; j < n; j++)
    {
      f = omega * u[j];
      for (int i=1; i < n; i++)
        Q[i][j] = Q[i][j] - f*u[i];
    }
#endif

    // Calculate updated A[1][2] and store it in e[1]
    e[1] = A[1][2] - q[1]*u[2] - u[1]*q[2];
  }
  else
  {
    for (int i=0; i < n; i++)
      d[i] = A[i][i];
    e[1] = A[1][2];
  }
}



// Macros
#define SQR(x)      ((x)*(x))                        // x^2 


// ----------------------------------------------------------------------------
int eigen(double A[3][3], double Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the QL algorithm with implicit shifts, preceded by a
// Householder reduction to tridiagonal form.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Dependencies:
//   dsytrd3()
// ----------------------------------------------------------------------------
{
  const int n = 3;
  double e[3];                   // The third element is used only as temporary workspace
  double g, r, p, f, b, s, c, t; // Intermediate storage
  int nIter;
  int m;

  // Transform A to real tridiagonal form by the Householder method
  HouseHolder(A, Q, w, e);
  
  // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
  // with the QL method
  //
  // Loop over all off-diagonal elements
  for (int l=0; l < n-1; l++)
  {
    nIter = 0;
    while (1)
    {
      // Check for convergence and exit iteration loop if off-diagonal
      // element e(l) is zero
      for (m=l; m <= n-2; m++)
      {
        g = fabs(w[m])+fabs(w[m+1]);
        if (fabs(e[m]) + g == g)
          break;
      }
      if (m == l)
        break;
      
      if (nIter++ >= 30)
        return -1;

      // Calculate g = d_m - k
      g = (w[l+1] - w[l]) / (e[l] + e[l]);
      r = sqrt(SQR(g) + 1.0);
      if (g > 0)
        g = w[m] - w[l] + e[l]/(g + r);
      else
        g = w[m] - w[l] + e[l]/(g - r);

      s = c = 1.0;
      p = 0.0;
      for (int i=m-1; i >= l; i--)
      {
        f = s * e[i];
        b = c * e[i];
        if (fabs(f) > fabs(g))
        {
          c      = g / f;
          r      = sqrt(SQR(c) + 1.0);
          e[i+1] = f * r;
          c     *= (s = 1.0/r);
        }
        else
        {
          s      = f / g;
          r      = sqrt(SQR(s) + 1.0);
          e[i+1] = g * r;
          s     *= (c = 1.0/r);
        }
        
        g = w[i+1] - p;
        r = (w[i] - g)*s + 2.0*c*b;
        p = s * r;
        w[i+1] = g + p;
        g = c*r - b;

        // Form eigenvectors
#ifndef EVALS_ONLY
        for (int k=0; k < n; k++)
        {
          t = Q[k][i+1];
          Q[k][i+1] = s*Q[k][i] + c*t;
          Q[k][i]   = c*Q[k][i] - s*t;
        }
#endif 
      }
      w[l] -= p;
      e[l]  = g;
      e[m]  = 0.0;
    }
  }

  return 0;
}

/**
 * 
 * @param a
 * @return 
 */
/*
double det(double a[3][3]){
    printf("DET;::::       %f %f %f\n",
             a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2]), 
            -a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2]),
            +a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]));
    
    
    return (a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
           -a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2])
           +a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]));
}*/
double det(double matrix[3][3]){
    double aux1 =  matrix[0][0]*matrix[1][1]*matrix[2][2] +
                matrix[0][1]*matrix[1][2]*matrix[2][0] +
                matrix[0][2]*matrix[1][0]*matrix[2][1];
    double aux2 =  matrix[0][2]*matrix[1][1]*matrix[2][0] +
                matrix[2][1]*matrix[1][2]*matrix[0][0] +
                matrix[2][2]*matrix[1][0]*matrix[0][1];
    return aux1-aux2;
}


double my_norm (double v[3]){
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/*
 * cross product
 */
 void cross(double A[3], double B[3], double crossProd[3]){
    crossProd[0] = A[1]*B[2] - A[2]*B[1];
    crossProd[1] = A[2]*B[0] - A[0]*B[2];
    crossProd[2] = A[0]*B[1] - A[1]*B[0];
}
/**
 * 
 * @param A
 * @param B
 * @return 
 */
 double dot(double A[3], double B[3]){
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}
 
 /**
 * Normalize a vector
 */
double normalizeVector(double v[3], int index){
    double n = my_norm(v);
    if(fabs(n)>TOLERANCE){
        v[0]    = v[0]/n;
        v[1]    = v[1]/n;
        v[2]    = v[2]/n;
    }else{
       // printf("Norm is too low[%d]: %e / %e\n",index,  n, TOLERANCE);
    }
    return n;
}
/**
 * 
 * @param a
 */
void printMatrix(double a[3][3]){
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++)
            printf("%e ",a[i][j]);
        printf("\n");
    }
    printf("---------------------------\n");

}

 /**
  * 
  * @param U2
  * @param U
  */
 void computeSqrtMatrix(double U2[3][3], double U[3][3])
 {
    double eigenVectors[3][3];
    double eigenValues[3];    
    eigen(U2, eigenVectors, eigenValues);
    double diagEigVal[3][3], invEigVec[3][3], aux[3][3];
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            if(i==j)
                diagEigVal[i][j]=sqrt(eigenValues[j]);
            else
                diagEigVal[i][j]=0.0;
        }
    }
    //Q*D^0.5*Q^-1
    invertMatrix(eigenVectors, invEigVec);
    matrixMultiplication(eigenVectors, diagEigVal, aux);
    matrixMultiplication(aux, invEigVec, U);
 }
 /**
  * 
  * @param M
  * @param N
  */
void invertMatrixEigenValue(double M[3][3], double N[3][3]){
    double eigenVectors[3][3];
    double eigenValues[3];

    eigen(M, eigenVectors, eigenValues);
    double diagEigVal[3][3], invEigVec[3][3], aux[3][3];
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++)
        {
            if(i==j)
                diagEigVal[i][j]=1.0/eigenValues[j];
            else
                diagEigVal[i][j]=0.0;
        }
    }
    //Q*D^-1*Q^-1
    invertMatrix(eigenVectors, invEigVec);
    matrixMultiplication(eigenVectors, diagEigVal, aux);
    matrixMultiplication(aux, invEigVec, N);
}
/**
 * 
 * @param M
 * @param N
 * @param detM
 */
void invertMatrixFormula(double M[3][3], double N[3][3], double detM){
    //First Row
    N[0][0] = (M[1][1] * M[2][2] - M[1][2] * M[2][1]) / detM;
    N[0][1] = (M[0][2] * M[2][1] - M[0][1] * M[2][2]) / detM;
    N[0][2] = (M[0][1] * M[1][2] - M[0][2] * M[1][1]) / detM;
    //Second Row
    N[1][0] = (M[1][2] * M[2][0] - M[1][0] * M[2][2]) / detM;
    N[1][1] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) / detM;
    N[1][2] = (M[0][2] * M[1][0] - M[0][0] * M[1][2]) / detM;
    //Third Row
    N[2][0] = (M[1][0] * M[2][1] - M[1][1] * M[2][0]) / detM;
    N[2][1] = (M[0][1] * M[2][0] - M[0][0] * M[2][1]) / detM;
    N[2][2] = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) / detM;
}
/**
 * N=M^-1
 * http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
 * @param M
 * @param N
 */
void invertMatrix(double M[3][3], double N[3][3]){
    double detM = det(M);
    //determinant too much close to zero. There's no way to invert the matrix
    if(fabs(detM)<TOLERANCE){
        stringstream ss;
        ss<<"Trying to invert a matrix with 0 determinant: det="<<detM<<";"<<endl;
        ss<<"Tolerance: "<<RELTOL<<endl;
        ss<<M[0][0]<<" "<<M[0][1]<<" "<<M[0][2]<<endl;
        ss<<M[1][0]<<" "<<M[1][1]<<" "<<M[1][2]<<endl;
        ss<<M[2][0]<<" "<<M[2][1]<<" "<<M[2][2]<<endl;
        string str = ss.str();
        throw MyException(str, __FILE__, __LINE__);
    }
    double tol=1.0e-4;
    //very small determinant, compute the inverse via eigen decomposition
    if(fabs(detM)>=TOLERANCE && fabs(detM)<tol){
        invertMatrixEigenValue(M, N);
    //determinant is big enough to compute the inverse via direct formula
    }else if(fabs(detM)>=tol){
        invertMatrixFormula(M, N, detM);       
    }
    //Check if N is really M^-1, by multiplying MN and verifying if it the Identity matrix.
    /*double testeI[3][3];
    matrixMultiplication(M, N, testeI);
    if(!matrixCheck(testeI, NULL)){
        stringstream ss;
        ss<<"Inversion Failed!!"<<endl;
        ss<<"Det: "<<detM<<endl;
        ss<<"M*M^-1!=I "<<endl;
        ss<<"teste: "<<endl;
        ss<<testeI[0][0]<<" "<<testeI[0][1]<<" "<<testeI[0][2]<<endl;
        ss<<testeI[1][0]<<" "<<testeI[1][1]<<" "<<testeI[1][2]<<endl;
        ss<<testeI[2][0]<<" "<<testeI[2][1]<<" "<<testeI[2][2]<<endl;
        ss<<"M: "<<endl;
        ss<<M[0][0]<<" "<<M[0][1]<<" "<<M[0][2]<<endl;
        ss<<M[1][0]<<" "<<M[1][1]<<" "<<M[1][2]<<endl;
        ss<<M[2][0]<<" "<<M[2][1]<<" "<<M[2][2]<<endl;
        ss<<"N: "<<endl;
        ss<<N[0][0]<<" "<<N[0][1]<<" "<<N[0][2]<<endl;
        ss<<N[1][0]<<" "<<N[1][1]<<" "<<N[1][2]<<endl;
        ss<<N[2][0]<<" "<<N[2][1]<<" "<<N[2][2]<<endl;
        string str = ss.str();
        throw MyException(str, __FILE__, __LINE__);
    }*/
}

/**
 * Matrix multiplication
 * C=A*B
 * @param A
 * @param B
 * @param C
 */
void matrixMultiplication(double A[3][3], double B[3][3], double C[3][3])
{
    double sum=0.0;
    for (int i = 0 ; i < 3 ; i++ )
    {
        for (int j = 0 ; j < 3 ; j++ )
        {
            for (int k = 0 ; k < 3 ; k++ )
            {
                sum = sum + A[i][k]*B[k][j];
            }
            C[i][j] = sum;
            sum = 0.0;
        }
    }
}
/**
 * Transpose matrix M
 * N = M'
 * @param M
 * @param N
 */
void transpose(double M[3][3], double N[3][3]){
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            N[j][i]=M[i][j];
}


/**
 * 
 * @param M
 * @param N
 * @return 
 */
bool matrixCheck(double M[3][3], double N[3][3]){
    
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            //if N it is null, compares to identity matrix
            if(N==NULL){
                if(i==j){//diagonal
                    if(fabs(1.0-M[i][j])>ABSTOL) 
                        return false;
                }else{  //nondiagonal
                    if(fabs(M[i][j])>ABSTOL)  
                        return false;
                }
            }else{ //if N it not null, compares M and N
                if(fabs(M[i][j])<TOLERANCE){
                    if(fabs((M[i][j]-N[i][j]))>ABSTOL) 
                        return false;
                }else{
                    if(fabs((M[i][j]-N[i][j])/M[i][j])>RELTOL && fabs((M[i][j]-N[i][j]))>ABSTOL) 
                        return false;
                }
            }
        }
    }
    return true;
}
