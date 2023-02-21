/* 
 * File:   Constants.h
 * Author: ricardo
 *
 * Created on April 24, 2013, 3:43 PM
 */

#ifndef CONSTANTS_H
#define	CONSTANTS_H

#define I2d(i, j, cols) i*cols+j
#define RELTOL_DT 0.001
#define ABSTOL_DT 0.0001

//#define RELTOL_DT 0.001
//#define ABSTOL_DT 0.01


#define KEY_NOT_FOUND -1
#define TOLERANCE 1e-12
#define LIN_ALG_TOL 2.0e-5
#define RELTOLAREA 1e-7
#define RELTOL 1e-5 //5e-7
#define ABSTOL 1e-7 //5e-7
#define elementsCols 4 //nao alterar, usado para encoontrar o baricentro
//cell conditions
#define HEALTHY 0
#define ISCHEMIC 1
#define DEAD	3
//cell pace maker
#define paceMaker 4
//cell V states
#define V0 0
#define V1 1
#define V2 2
#define V3 3
#define V4 4

//cell F states
#define F0 0
#define F1 1
#define F2 2
#define F3 3
#define F4 4

//mass density (kg/m^3) 
//#define ro 1000.0
//mass density (kg/cm^3) 
#define ro_mass_dens 0.001f //0.001g/mm^3
//mass density (kg/cm^3) 
// ESTAVA com este arqui #define ro 1.0e-15

#define FIBER 0
#define SHEET 1
#define NORMAL 2

#define THREADSPERBLOCK 512

//#define VEL_STOP 0.001 
//#define VEL_STOP 0.01//exp1
#define VEL_STOP 0.05 //exp 2 e 3
#define EMPTY_FLAG -1
#endif	/* CONSTANTS_H */
