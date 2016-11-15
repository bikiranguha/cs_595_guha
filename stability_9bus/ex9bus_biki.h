#ifndef PF_H
#define PF_H

#include <petscsnes.h>
#include <petscdmnetwork.h>

# define MAXLINE 1000
#define REF_BUS 3
#define PV_BUS 2
#define PQ_BUS 1
#define ISOLATED_BUS 4
#define NGEN_AT_BUS_MAX 15
#define NLOAD_AT_BUS_MAX 1

/* 2. Bus data */
/* 11 columns */
struct _p_VERTEXDATA{
  PetscInt      bus_i; /* Integer bus number .. used by some formats like Matpower */
  //char	 	i[20]; /* Bus Number */
  //char 		name[20]; /* Bus Name */
  PetscScalar 	basekV; /* Bus Base kV */
  //PetscInt 	ide; /* Bus type code */
  //PetscScalar 	gl; /* Active component of shunt admittance to ground */
  //PetscScalar 	bl; /* Reactive component of shunt admittance to ground */
  //PetscInt 	area; /* Area number */
  //PetscInt 	zone; /* Zone number */
  PetscScalar 	vr; /* Real part of bus voltage; in pu */
  PetscScalar 	vi; /* Imaginary part of bus voltage ; in pu */
 // PetscInt 	owner; /* Owner number */
  //PetscInt	internal_i; /* Internal Bus Number */
  PetscInt      ngen; /* Number of generators incident at this bus */
  PetscInt      gidx[NGEN_AT_BUS_MAX]; /* list of inndices for accessing the generator data in GEN structure */
  PetscInt      nload;
  PetscInt      lidx[NLOAD_AT_BUS_MAX];
  PetscScalar   yd[2];//yd[0]: imaginary part of yii; yd[1]: real part of yii
  
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_VERTEXDATA *VERTEXDATA;

/* 3. Load data */
/* 12 columns */
struct _p_LOAD{
  PetscInt      bus_i; /* Bus number */
 // char 		i[20]; /* Bus Number or extended bus name*/
 // char 		id[20]; /* Load identifier, in case of multiple loads. 1 by default */
 // PetscInt 	status; /* Load status */
 // PetscInt 	area; /* Area to which load is assigned */
 // PetscInt 	zone; /* Zone to which load is assigned */
 // PetscScalar 	pl; /* Active power component of constant MVA load */
 // PetscScalar 	ql; /* Reactive power component of constant MVA load */
//  PetscScalar 	ip; /* Active power component of constant current load: MW pu V */
//  PetscScalar 	iq; /* Reactive power component of constant current load: Mvar pu V */
//  PetscScalar 	yp; /* Active power component of constant admittance load: MW pu V */
 // PetscScalar 	yq; /* Reactive power component of constant admittance load: Mvar pu V */
//  PetscInt 	owner; /* Owner number */
//  PetscInt	internal_i; /* Internal Bus Number */
 // PetscScalar   scale_load;
 
 PetscScalar PD0;
 PetscScalar QD0;
 PetscInt    ld_nsegsp;
 PetscScalar ld_alphap;
 PetscScalar ld_betap;
 PetscInt    ld_nsegsq;
 PetscScalar ld_alphaq;
 PetscScalar ld_betaq; 
 
 /* ld_nsegsp,ld_nsegsq - Number of individual load models for real and reactive power loads
    ld_alphap,ld_alphap - Percentage contribution (weights) or loads
    P_D0                - Real power load
    Q_D0                - Reactive power load
    V_m(t)              - Voltage magnitude at time t
    V_m0                - Voltage magnitude at t = 0
    ld_betap, ld_betaq  - exponents describing the load model for real and reactive part */
	
	
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_LOAD *LOAD;

/* 4. Generator data */
/* 20+ columns */
/******* 20, USING ONLY 1 OWNER's WORTH OF DATA. COME BACK TO THIS LATER, if necessary ******/
struct _p_GEN{
  PetscInt      bus_i;
  char 		i[20]; /* Bus Number or extended bus name*/
 // char 		id[20]; /* Generator identifier, in case of multiple generators at same bus. 1 by default */
  PetscScalar 	pg; /* Generator active power output */
  PetscScalar 	qg; /* Generator reactive power output */
  //PetscScalar 	qt; /* Maximum reactive power output: Mvar */
 // PetscScalar 	qb; /* Minimum reactive power output: Mvar */
 // PetscScalar 	vs; /* Regulated voltage setpoint: pu */
  PetscInt 	ireg; /* Remote bus number/identifier */
  PetscScalar 	mbase; /* MVA base of the machine */
  
 // PetscScalar 	zr; /* Complex machine impedance ZSOURCE in pu on mbase */
  //PetscScalar 	zx; /* ----------------------"------------------------- */
  //PetscScalar 	rt; /* Step-up transformer impedance XTRAN in pu on mbase */
 // PetscScalar 	xt; /* -----------------------"-------------------------- */
 // PetscScalar 	gtap; /* Step-up transformer turns ratio */
  //PetscInt 	status; /* Machine status */
 // PetscScalar 	rmpct; /* Mvar % required to hold voltage at remote bus */
//  PetscScalar 	pt; /* Gen max active power output: MW */
 // PetscScalar 	pb; /* Gen min active power output: MW */
 // PetscInt 	o1; /* Owner number */
 // PetscScalar 	f1; /* Fraction of ownership */
 // PetscInt	internal_i; /* Internal Bus Number */
  //PetscScalar   scale_gen;
 /* Generator real and reactive powers (found via loadflow) */
/* Generator constants */
 PetscScalar H  ;   /* Inertia constant */ 
 PetscScalar Rs; /* Stator Resistance */
 PetscScalar Xd;  /* d-axis reactance */
 PetscScalar Xdp; /* d-axis transient reactance */
 PetscScalar Xq ; /* q-axis reactance Xq(1) set to 0.4360, value given in text 0.0969 */
 PetscScalar Xqp; /* q-axis transient reactance */
 PetscScalar Td0p ; /* d-axis open circuit time constant */
 PetscScalar Tq0p ; /* q-axis open circuit time constant */
PetscScalar M; /* M = 2*H/w_s */
PetscScalar D; /* D = 0.1*M */

PetscScalar TM; /* Mechanical Torque */
/* Exciter system constants */
PetscScalar KA;  /* Voltage regulartor gain constant */
PetscScalar TA;     /* Voltage regulator time constant */
PetscScalar KE;     /* Exciter gain constant */
PetscScalar TE; /* Exciter time constant */
PetscScalar KF;  /* Feedback stabilizer gain constant */
PetscScalar TF;    /* Feedback stabilizer time constant */
PetscScalar k1;
PetscScalar k2;  /* k1 and k2 for calculating the saturation function SE = k1*exp(k2*Efd) */

PetscScalar Vref; 
  
  
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_GEN *GEN;

/* 17+ columns */
struct _p_EDGEDATA{
  PetscInt      fbus;
  PetscInt      tbus;
 // char 		i[20]; /* Bus Number or extended bus name*/
 // char 		j[20]; /* Bus Number or extended bus name*/
 // char 		ckt[20]; /* Circuit identifier. 1 by default */
  //PetscScalar 	r; /* Branch resistance: pu */
  //PetscScalar 	x; /* Branch reactance: pu */
  //PetscScalar 	b; /* Branch charging susceptance: pu */
  //PetscScalar 	rateA; /* rate A in MVA */
 // PetscScalar 	rateB; /* rate B in MVA */
 // PetscScalar 	rateC; /* rate C in MVA */
 // PetscScalar   tapratio;
 // PetscScalar   phaseshift;
  //PetscScalar 	gi; /* Complex admittance at 'i' end: pu */
  //PetscScalar 	bi; /* Complex admittance at 'i' end: pu */
  //PetscScalar 	gj; /* Complex admittance at 'j' end: pu */
  //PetscScalar 	bj; /* Complex admittance at 'j' end: pu */
  //PetscInt 	status; /* Service status */
  //PetscScalar 	length; /* Line length */
  //PetscInt 	o1; /* Owner number */
  //PetscScalar 	f1; /* Fraction of ownership */
  //PetscInt	internal_i; /* Internal From Bus Number */
  //PetscInt	internal_j; /* Internal To Bus Number */
  //PetscScalar   yff[2],yft[2],ytf[2],ytt[2]; /* [G,B] */
  PetscScalar   yft[2];//yft[0]: imaginary part of yij;yft[1]: real part of yij 
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_EDGEDATA *EDGEDATA;


typedef struct{
  PetscScalar sbase; /* System base MVA */
  PetscInt    nbus,ngen,nbranch,nload; /* # of buses,gens,branches, and loads (includes elements which are
                                          out of service */
  VERTEXDATA  bus;
  LOAD        load;
  GEN         gen;
  EDGEDATA    branch;
}SYSDATA PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

extern PetscErrorCode ReadData(SYSDATA*, Num /* Make sure to complete this */);


 
// typedef struct{
  // PetscScalar sbase; /* System base MVA */
  // PetscInt    nbus,ngen,nbranch,nload; /* # of buses,gens,branches, and loads (includes elements which are
                                          // out of service */
  // VERTEXDATA  bus;
  // LOAD        load;
  // GEN         gen;
  // EDGEDATA    branch;
// }PFDATA PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

// extern PetscErrorCode PFReadMatPowerData(PFDATA*,char*);

#endif
