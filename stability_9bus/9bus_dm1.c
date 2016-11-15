static char help[] = "Power grid stability analysis of WECC 9 bus system.\n\
This example is based on the 9-bus (node) example given in the book Power\n\
Systems Dynamics and Stability (Chapter 7) by P. Sauer and M. A. Pai.\n\
The power grid in this example consists of 9 buses (nodes), 3 generators,\n\
3 loads, and 9 transmission lines. The network equations are written\n\
in current balance form using rectangular coordiantes.\n\n";

/* T
   Concepts: DMNetwork
   Concepts: PETSc TS solver
*/

//#include "9bus.h"
#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>

#define freq 60
#define w_s (2*PETSC_PI*freq)
#include <petscdmnetwork.h>
#include <petscksp.h>

typedef struct {
  PetscInt id; /* Bus Number or extended bus name*/
  PetscInt gbus; /* Bus Number or extended bus name*/
  //char 		id[20]; /* Generator identifier, in case of multiple generators at same bus. 1 by default */
  PetscScalar 	mbase; /* MVA base of the machine */
  PetscScalar PG; /* Generator active power output */
  PetscScalar QG; /* Generator reactive power output */
/* Generator constants */
PetscScalar H;   /* Inertia constant */
PetscScalar Rs; /* Stator Resistance */
PetscScalar Xd;  /* d-axis reactance */
PetscScalar Xdp; /* d-axis transient reactance */
PetscScalar Xq; /* q-axis reactance Xq(1) set to 0.4360, value given in text 0.0969 */
PetscScalar Xqp; /* q-axis transient reactance */
PetscScalar Td0p; /* d-axis open circuit time constant */
PetscScalar Tq0p; /* q-axis open circuit time constant */
PetscScalar M; /* M = 2*H/w_s */
PetscScalar D; /* D = 0.1*M */
PetscScalar TM; /* Mechanical Torque */

/* Exciter system constants */
PetscScalar KA ;  /* Voltage regulartor gain constant */
PetscScalar TA;     /* Voltage regulator time constant */
PetscScalar KE;     /* Exciter gain constant */
PetscScalar TE; /* Exciter time constant */
PetscScalar KF;  /* Feedback stabilizer gain constant */
PetscScalar TF;    /* Feedback stabilizer time constant */
PetscScalar k1;
PetscScalar k2;  /* k1 and k2 for calculating the saturation function SE = k1*exp(k2*Efd) */

PetscScalar Vref;
 
} Gen;

typedef struct {
  PetscInt    id; /* node id */
  PetscInt    nofgen;/* Number of generators at the bus*/
  PetscInt    nofload;/*  Number of load at the bus*/
 // PetscInt gbus; /* Bus Number or extended bus name*/
 // PetscScalar inj; /* current injection (A) */
  //PetscBool   gr; /* grounded ? */
  PetscScalar 	yff[2]; /* yff[0]= imaginary part of admittance, yff[1]=real part of admittance*/
  //PetscScalar 	vr; /* real part of Bus voltage; in pu */
  //PetscScalar 	vi; /* imaginary part of Bus voltage phase angle */
  PetscScalar   vr;
  PetscScalar   vi;
} Bus;

typedef struct {
  PetscInt    id; /* bus id */
  PetscInt lbus; /* Bus Number or extended bus name*/
  PetscScalar PD0;
  PetscScalar QD0;
  PetscInt    ld_nsegsp;
  //PetscScalar *ld_alphap;//ld_alphap=[1,0,0], an array, not a value, so use *ld_alphap;
  PetscScalar ld_alphap[3];//ld_alphap=[1,0,0], an array, not a value, so use *ld_alphap;
  //define: PetscScalar ld_alphap[20], then strcpy(load[i].ld_alphap, ld_alphap); 
  //define: PetscScalar *ld_alphap, then load[i].ld_alphap=ld_alphap.
  PetscScalar ld_betap[3];
  PetscInt    ld_nsegsq;
  PetscScalar ld_alphaq[3];
  PetscScalar ld_betaq[3];
} Load;

typedef struct {
  PetscInt    id; /* node id */
 // PetscScalar inj; /* current injection (A) */
  //PetscBool   gr; /* grounded ? */
  PetscScalar 	yft[2]; /* yft[0]= imaginary part of admittance, yft[1]=real part of admittance*/
  //PetscScalar 	vr; /* real part of Bus voltage; in pu */
  //PetscScalar 	vi; /* imaginary part of Bus voltage phase angle */
} Branch;

typedef struct {
  PetscReal   tfaulton,tfaultoff; /* Fault on and off times */
  PetscInt    faultbus; /* Fault bus */
  PetscScalar Rfault;
  PetscReal   t0,tmax;//t0: initial time,final time
  //PetscInt    neqs_gen,neqs_net,neqs_pgrid; // neqs_gen: No. of gen equations, and so on....
  Mat         Sol; /* Matrix to save solution at each time step */
  PetscInt    stepnum;
  PetscBool   alg_flg;
  PetscReal   t;
  IS          is_diff; /* indices for differential equations */
  IS          is_alg; /* indices for algebraic equations */
  PetscBool   setisdiff; /* TS computes truncation error based only on the differential variables */
  PetscScalar    ybusfault[18];
  
} Userctx;

#undef __FUNCT__
#define __FUNCT__ "read_data"
PetscErrorCode read_data(PetscInt ngen, PetscInt nload, PetscInt nbus, PetscInt nbranch, Mat Ybus,Vec V0,Gen **pgen,Load **pload,Bus **pbus, Branch **pbranch, PetscInt **pedgelist)
{
  PetscErrorCode    ierr;
  PetscInt          i,row[1],col[2];
  //PetscInt          j;
  Bus              *bus;
  Branch            *branch;
  Gen               *gen;
  Load              *load;
  //PetscInt          edgelist[18]={0,3,1,6,2,8,3,4,3,5,4,6,5,8,6,7,7,8};//ask
  PetscInt          *edgelist;
  PetscScalar    *varr;
  PetscScalar M[3], D[3];

  PetscInt nofgen[9] = {1,1,1,0,0,0,0,0,0}; /* Buses at which generators are incident */
  PetscInt nofload[9] = {0,0,0,0,1,1,0,1,0}; /* Buses at which loads are incident */
 
  /*10 parameters*//* Generator real and reactive powers (found via loadflow) */
  PetscInt gbus[3] = {0,1,2}; /* Buses at which generators are incident */
  PetscScalar PG[3] = {0.716786142395021,1.630000000000000,0.850000000000000};
  PetscScalar QG[3] = {0.270702180178785,0.066120127797275,-0.108402221791588};
  /* Generator constants */
  PetscScalar H[3]    = {23.64,6.4,3.01};   /* Inertia constant */
  PetscScalar Rs[3]   = {0.0,0.0,0.0}; /* Stator Resistance */
  PetscScalar Xd[3]   = {0.146,0.8958,1.3125};  /* d-axis reactance */
  PetscScalar Xdp[3]  = {0.0608,0.1198,0.1813}; /* d-axis transient reactance */
  PetscScalar Xq[3]   = {0.4360,0.8645,1.2578}; /* q-axis reactance Xq(1) set to 0.4360, value given in text 0.0969 */
  PetscScalar Xqp[3]  = {0.0969,0.1969,0.25}; /* q-axis transient reactance */
  PetscScalar Td0p[3] = {8.96,6.0,5.89}; /* d-axis open circuit time constant */
  PetscScalar Tq0p[3] = {0.31,0.535,0.6}; /* q-axis open circuit time constant */

  /* Exciter system constants (8 parameters)*/
  PetscScalar KA[3] = {20.0,20.0,20.0};  /* Voltage regulartor gain constant */
  PetscScalar TA[3] = {0.2,0.2,0.2};     /* Voltage regulator time constant */
  PetscScalar KE[3] = {1.0,1.0,1.0};     /* Exciter gain constant */
  PetscScalar TE[3] = {0.314,0.314,0.314}; /* Exciter time constant */
  PetscScalar KF[3] = {0.063,0.063,0.063};  /* Feedback stabilizer gain constant */
  PetscScalar TF[3] = {0.35,0.35,0.35};    /* Feedback stabilizer time constant */
  PetscScalar k1[3] = {0.0039,0.0039,0.0039};
  PetscScalar k2[3] = {1.555,1.555,1.555};  /* k1 and k2 for calculating the saturation function SE = k1*exp(k2*Efd) */
  
  /* Load constants(10 parameter)
  We use a composite load model that describes the load and reactive powers at each time instant as follows
  P(t) = \sum\limits_{i=0}^ld_nsegsp \ld_alphap_i*P_D0(\frac{V_m(t)}{V_m0})^\ld_betap_i
  Q(t) = \sum\limits_{i=0}^ld_nsegsq \ld_alphaq_i*Q_D0(\frac{V_m(t)}{V_m0})^\ld_betaq_i
  where
    id                  - index of the load
    lbus                - Buses at which loads are incident 
    ld_nsegsp,ld_nsegsq - Number of individual load models for real and reactive power loads
    ld_alphap,ld_alphap - Percentage contribution (weights) or loads
    P_D0                - Real power load
    Q_D0                - Reactive power load
    V_m(t)              - Voltage magnitude at time t
    V_m0                - Voltage magnitude at t = 0
    ld_betap, ld_betaq  - exponents describing the load model for real and reactive part

    Note: All loads have the same characteristic currently.
  */
   PetscInt lbus[3] = {4,5,7};
   PetscScalar PD0[3] = {1.25,0.9,1.0};
   PetscScalar QD0[3] = {0.5,0.3,0.35};
   PetscInt    ld_nsegsp[3] = {3,3,3};
   const PetscScalar ld_alphap[3] = {1,0,0};
   PetscScalar ld_betap[3]  = {2,1,0};
   PetscInt    ld_nsegsq[3] = {3,3,3};
   PetscScalar ld_alphaq[3] = {1,0,0};
   PetscScalar ld_betaq[3]  = {2,1,0};
   
   
  M[0] = 2*H[0]/w_s; M[1] = 2*H[1]/w_s; M[2] = 2*H[2]/w_s;
  D[0] = 0.1*M[0]; D[1] = 0.1*M[1]; D[2] = 0.1*M[2];
   
   
  
  PetscFunctionBeginUser;
  ierr = PetscCalloc4(nbus,&bus,ngen,&gen,nload,&load,nbranch,&branch);CHKERRQ(ierr);
  //Allocates 2 cleared (zeroed) arrays of memory both aligned to PETSC_MEMALIGN
  // PetscCalloc2(size_t m1,type **r1,size_t m2,type **r2)
  // Input:
  // m1- number of elements to allocate in 1st chunk (may be zero)
  // m2- number of elements to allocate in 2nd chunk (may be zero)
  // Output:
  // r1- memory allocated in first chunk
  // r2- memory allocated in second chunk
  //ask: what's "memory aligned to PETSC_MEMALIGN"?
  // without the command, out of memory range
  
   ierr = VecGetArray(V0,&varr);CHKERRQ(ierr);
  
 //read bus data
 for (i = 0; i < nbus; i++) {
    bus[i].id  = i;
    bus[i].nofgen=nofgen[i];
    bus[i].nofload=nofload[i];
	bus[i].vr=varr[2*i];
	bus[i].vi=varr[2*i+1];
    row[0]=2*i;
    col[0]=2*i;
    col[1]=2*i+1;
    ierr=MatGetValues(Ybus,1,row,2,col,bus[i].yff);CHKERRQ(ierr);//real and imaginary part of admittance from Ybus into yff
  }
  
  //read generator data
 for (i = 0; i < ngen; i++) {
    gen[i].id  = i;
    gen[i].gbus= gbus[i];
	gen[i].PG = PG[i];
	gen[i].QG = QG[i];
    gen[i].H= H[i];
    gen[i].Rs= Rs[i];
    gen[i].Xd= Xd[i];
    gen[i].Xdp= Xdp[i];
    gen[i].Xq= Xq[i];
    gen[i].Xqp= Xqp[i]; 
    gen[i].Td0p= Td0p[i];
    gen[i].Tq0p= Tq0p[i];
    gen[i].M = M[i];
    gen[i].D = D[i];	
    
    //exciter system
    gen[i].KA= KA[i];
    gen[i].TA= TA[i];
    gen[i].KE=KE[i];
    gen[i].TE=TE[i];
    gen[i].KF=KF[i];
    gen[i].TF=TF[i];
    gen[i].k1=k1[i];
    gen[i].k2=k2[i];  
  }
  
   for (i = 0; i < nload; i++) {
    load[i].id  = i;
    load[i].lbus  =lbus[i];
    load[i].PD0=PD0[i];
    load[i].QD0=QD0[i];
      load[i].ld_nsegsp=ld_nsegsp[i];
	  
	  load[i].ld_alphap[0]=ld_alphap[0];
	  load[i].ld_alphap[1]=ld_alphap[1];
	  load[i].ld_alphap[2]=ld_alphap[2];
	  
	   load[i].ld_alphaq[0]=ld_alphaq[0];
	  load[i].ld_alphaq[1]=ld_alphaq[1];
	  load[i].ld_alphaq[2]=ld_alphaq[2];
	  
	  load[i].ld_betap[0]=ld_betap[0];
	  load[i].ld_betap[1]=ld_betap[1];
      load[i].ld_betap[2]=ld_betap[2];
      load[i].ld_nsegsq=ld_nsegsq[i];

	  
	  load[i].ld_betaq[0]=ld_betaq[0];
	  load[i].ld_betaq[1]=ld_betaq[1];
      load[i].ld_betaq[2]=ld_betaq[2];
   
  }

 ierr = PetscCalloc1(2*nbranch,&edgelist);CHKERRQ(ierr);

 //ask
  for (i = 0; i < nbranch; i++) {
    switch (i) {
      case 0:
        edgelist[2*i]     = 0;
        edgelist[2*i + 1] = 3;
        break;
      case 1:
        edgelist[2*i]     = 1;
        edgelist[2*i + 1] = 6;
        break;
      case 2:
        edgelist[2*i]     = 2;
        edgelist[2*i + 1] = 8;
        break;
      case 3:
        edgelist[2*i]     = 3;
        edgelist[2*i + 1] = 4;
        break;
      case 4:
        edgelist[2*i]     = 3;
        edgelist[2*i + 1] = 5;
        break;
      case 5:
        edgelist[2*i]     = 4;
        edgelist[2*i + 1] = 6;
        break;
      case 6:
        edgelist[2*i]     = 5;
        edgelist[2*i + 1] = 8;
        break;
      case 7:
        edgelist[2*i]     = 6;
        edgelist[2*i + 1] = 7;
        break;
      case 8:
        edgelist[2*i]     = 7;
        edgelist[2*i + 1] = 8;
        break;
      default:
        break;
    }
  }
 
 
 //read branch data
  for (i = 0; i < nbranch; i++) {
    branch[i].id  = i;
    row[0]=edgelist[2*i]*2;
    col[0]=edgelist[2*i+1]*2;
    col[1]=edgelist[2*i+1]*2+1;
    ierr=MatGetValues(Ybus,1,row,2,col,branch[i].yft);CHKERRQ(ierr);//imaginary part of admittance
  }
  
  
  // *pngen     =ngen;
  // *pnload    =nload;
  // *pnbus     =nbus;
  // *pnbranch  =nbranch;
  *pgen     =gen;
  *pload    =load;
  *pbus     = bus;
  *pbranch  = branch;
  *pedgelist= edgelist;
  
  ierr = VecRestoreArray(V0,&varr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "SetInitialGuess"
PetscErrorCode SetInitialGuess(DM networkdm, Vec X, Vec V0)
{
  PetscErrorCode ierr;
  Bus            *bus;
  Gen            *gen;
  PetscInt       v, vStart, vEnd, offset;
  PetscBool      ghostvtex;
  Vec            localX;
  PetscScalar    *xarr;
  PetscInt       key,numComps,j,offsetd;
  DMNetworkComponentGenericDataType *arr;
  //PetscInt       i ;
  PetscInt		 idx=0;
  PetscScalar    Vr,Vi,IGr,IGi,Vm,Vm2;
  PetscScalar    Eqp,Edp,delta;
  PetscScalar    Efd,RF,VR; /* Exciter variables */
  PetscScalar    Id,Iq;  /* Generator dq axis currents */
  PetscScalar    theta,Vd,Vq,SE;

  PetscFunctionBegin;
  // M[0] = 2*H[0]/w_s; M[1] = 2*H[1]/w_s; M[2] = 2*H[2]/w_s;
  // D[0] = 0.1*M[0]; D[1] = 0.1*M[1]; D[2] = 0.1*M[2];
  
  ierr = DMNetworkGetVertexRange(networkdm,&vStart, &vEnd);CHKERRQ(ierr);

  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);

  ierr = VecSet(X,0.0);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArray(localX,&xarr);CHKERRQ(ierr);
  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);
  
   for (v = vStart; v < vEnd; v++) {
     ierr = DMNetworkIsGhostVertex(networkdm,v,&ghostvtex);CHKERRQ(ierr);
     if (ghostvtex) continue;
   
     ierr = DMNetworkGetVariableOffset(networkdm,v,&offset);CHKERRQ(ierr);
     ierr = DMNetworkGetNumComponents(networkdm,v,&numComps);CHKERRQ(ierr);
     for (j=0; j < numComps; j++) {
       ierr = DMNetworkGetComponentTypeOffset(networkdm,v,j,&key,&offsetd);CHKERRQ(ierr);
       if (key == 1) {
	 bus = (Bus*)(arr+offsetd);
	 xarr[offset] = bus->vr;
	 xarr[offset+1] = bus->vi;//->
	 Vr = bus->vr;
	 Vi= bus->vi;
         } else if(key == 2) {
	gen = (Gen*)(arr+offsetd);
    Vm  = PetscSqrtScalar(Vr*Vr + Vi*Vi); Vm2 = Vm*Vm;
    IGr = (Vr*gen->PG + Vi*gen->QG)/Vm2; // Real part of gen current
    IGi = (Vi*gen->PG - Vr*gen->QG)/Vm2; // Imaginary part of gen current

    
	/* Machine angle */
	delta = atan2(Vi+gen->Xq*IGr,Vr-gen->Xq*IGi);  // Link for atan2: https://en.wikipedia.org/wiki/Atan2 

    theta = PETSC_PI/2.0 - delta;

    Id = IGr*PetscCosScalar(theta) - IGi*PetscSinScalar(theta); /* d-axis stator current */
    Iq = IGr*PetscSinScalar(theta) + IGi*PetscCosScalar(theta); /* q-axis stator current */

    Vd = Vr*PetscCosScalar(theta) - Vi*PetscSinScalar(theta);
    Vq = Vr*PetscSinScalar(theta) + Vi*PetscCosScalar(theta);

    Edp = Vd + gen->Rs*Id - gen->Xqp*Iq; /* d-axis transient EMF */
    Eqp = Vq + gen->Rs*Iq + gen->Xdp*Id; /* q-axis transient EMF */

    gen->TM = gen->PG;
	idx=offset+2;
	xarr[idx]   = Eqp;
    xarr[idx+1] = Edp;
    xarr[idx+2] = delta;
    xarr[idx+3] = w_s;

    idx = idx + 4;

    xarr[idx]   = Id;
    xarr[idx+1] = Iq;

    idx = idx + 2;

    /* Exciter */
    Efd = Eqp + (gen->Xd - gen->Xdp)*Id;
    SE  = gen->k1*PetscExpScalar(gen->k2*Efd);
    VR  =  gen->KE*Efd + SE;
    RF  =  gen->KF*Efd/gen->TF;

    xarr[idx]   = Efd;
    xarr[idx+1] = RF;
    xarr[idx+2] = VR;

    gen->Vref = Vm + (VR/gen->KA);
	
	
      }
     }
   }
  ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
 }

 /* Converts from machine frame (dq) to network (phase a real,imag) reference frame */
#undef __FUNCT__
#define __FUNCT__ "dq2ri"
PetscErrorCode dq2ri(PetscScalar Fd,PetscScalar Fq,PetscScalar delta,PetscScalar *Fr, PetscScalar *Fi)
{
  PetscFunctionBegin;
  *Fr =  Fd*PetscSinScalar(delta) + Fq*PetscCosScalar(delta);
  *Fi = -Fd*PetscCosScalar(delta) + Fq*PetscSinScalar(delta);
  PetscFunctionReturn(0);
}

/* Converts from network frame ([phase a real,imag) to machine (dq) reference frame */
#undef __FUNCT__
#define __FUNCT__ "ri2dq"
PetscErrorCode ri2dq(PetscScalar Fr,PetscScalar Fi,PetscScalar delta,PetscScalar *Fd, PetscScalar *Fq)
{
  PetscFunctionBegin;
  *Fd =  Fr*PetscSinScalar(delta) - Fi*PetscCosScalar(delta);
  *Fq =  Fr*PetscCosScalar(delta) + Fi*PetscSinScalar(delta);
  PetscFunctionReturn(0);
}

 
 
 
 #undef __FUNCT__
 #define __FUNCT__ "FormIFunction"
  PetscErrorCode FormIFunction(TS ts, PetscReal t,Vec X,Vec Xdot,Vec F, Userctx *user){
	
  PetscErrorCode ierr;
  //Wash           wash=(Wash)ctx;
  DM             networkdm;
  Vec            localX,localXdot,localF;
  //const PetscInt *cone;
  PetscInt       vfrom,vto,offsetfrom,offsetto;
  //PetscInt     type,varoffset,voffset;
  PetscInt       v,vStart,vEnd,e;
  //PetscInt       eStart,eEnd/* pipeoffset */;
 // PetscBool      ghost;
  PetscScalar    *farr;
  //PetscScalar  *vf,*juncx,*juncf;
  //Pipe           pipe;
  //PipeField      *pipex,*pipexdot,*pipef;
  // DMDALocalInfo  info;
 // Junction       junction;
  MPI_Comm       comm;
  PetscMPIInt    rank,size;
  const PetscScalar *xarr,*xdotarr;
  DMNetworkComponentGenericDataType *arr;
  PetscScalar    Vd,Vq,SE;
  //Userctx        *user=(Userctx*)ctx;
    
	
  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)ts,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr); 
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr); 
 
  ierr = VecSet(F,0.0);CHKERRQ(ierr);
  
  ierr = TSGetDM(ts,&networkdm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localF);CHKERRQ(ierr); 
  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localXdot);CHKERRQ(ierr);
  ierr = VecSet(localF,0.0);CHKERRQ(ierr);

  /* update ghost values of locaX and locaXdot */
  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  
  ierr = DMGlobalToLocalBegin(networkdm,Xdot,INSERT_VALUES,localXdot);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,Xdot,INSERT_VALUES,localXdot);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localXdot,&xdotarr);CHKERRQ(ierr);
  ierr = VecGetArray(localF,&farr);CHKERRQ(ierr);
  
  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr); 
  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr); 
  
  for (v=vStart; v < vEnd; v++) { // Vertices contain all the info about each component (bus, generator, branch and load)
    PetscInt    i,j,offsetd, offset, key;
    PetscScalar Vr, Vi;
    Bus         *bus;
    Gen         *gen;
    Load        *load;
    PetscBool   ghostvtex;
    PetscInt    numComps;
	PetscScalar Yffr, Yffi;
	PetscScalar   Vm, Vm2,Vm0;
	PetscScalar  Vr0, Vi0;
	PetscScalar  PD,QD;

    ierr = DMNetworkIsGhostVertex(networkdm,v,&ghostvtex);CHKERRQ(ierr);  // Ghost vertices may be buses belonging to other processors whose info is needed by the local processor.
	
    ierr = DMNetworkGetNumComponents(networkdm,v,&numComps);CHKERRQ(ierr); // Buses, branches, gen and loads are the components here.
    ierr = DMNetworkGetVariableOffset(networkdm,v,&offset);CHKERRQ(ierr); // Get the offset for accessing the variable associated with the given vertex/edge from the local vector.
	
	
	//if (ghostvtex) continue;
	
	
    for (j = 0; j < numComps; j++) {
      ierr = DMNetworkGetComponentTypeOffset(networkdm,v,j,&key,&offsetd);CHKERRQ(ierr);
      if (key == 1) {  // bus structure
	  
	  PetscInt       nconnedges;
	  const PetscInt *connedges;

	bus = (Bus*)(arr+offsetd);  // Q: How does this directly go to the bus required?
	
	if (!ghostvtex) { // Calculate current due to diagonal elements in Ybus
	  Vr = xarr[offset];
	  Vi= xarr[offset+1];
	  
	  
	  
	  Yffr = bus->yff[1];
	  Yffi = bus->yff[0];
	  
	  
	  
	  
	  // if (user->alg_flg == PETSC_TRUE){
		  
		  // Yffr += user->ybusfault[bus->id*2+1];
		  // Yffi += user->ybusfault[bus->id*2];
	  // }
	  
	  Vr0 = bus->vr;
	  Vi0 = bus->vi;
	  
	  farr[offset] +=  Yffi*Vr + Yffr*Vi; // imaginary current due to diagonal elements
	  farr[offset+1] += Yffr*Vr - Yffi*Vi; // real current due to diagonal elements
	 	}
		
	ierr = DMNetworkGetSupportingEdges(networkdm,v,&nconnedges,&connedges);CHKERRQ(ierr);  // nconnedges: no. of branches connected to the bus, connedges: list of branches	
		
		
		for (i=0; i < nconnedges; i++) { // loop for connected branches to the specific bus
	  Branch       *branch;
	  PetscInt       keye;
          PetscScalar    Yfti, Yftr;
          const PetscInt *cone;
		  PetscScalar  Vfr, Vfi, Vtr, Vti;

	  e = connedges[i];  // The current branch
	  ierr = DMNetworkGetComponentTypeOffset(networkdm,e,0,&keye,&offsetd);CHKERRQ(ierr);
	  branch = (Branch*)(arr+offsetd);
	  
	  Yfti = branch->yft[0];  
	  Yftr = branch->yft[1];
	  

	  ierr = DMNetworkGetConnectedNodes(networkdm,e,&cone);CHKERRQ(ierr);
	  vfrom = cone[0];  // From bus (vertex)
	  vto   = cone[1];   // To bus (vertex)

	  ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
	  ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

	  
	  // From bus and to bus real and imaginary voltages
	  Vfr     = xarr[offsetfrom];
	  Vfi     = xarr[offsetfrom+1];
	  Vtr	  = xarr[offsetto];
	  Vti     = xarr[offsetto+1];
	

	  if (vfrom == v) {
	    farr[offsetfrom]   += Yftr*Vti + Yfti*Vtr;
	    farr[offsetfrom+1] += Yftr*Vtr - Yfti*Vti;
	  } else {
	    farr[offsetto]   += Yftr*Vfi + Yfti*Vfr;
	    farr[offsetto+1] += Yftr*Vfr - Yfti*Vfi;
	  }
	  
	}

   }
   
   else if (key == 2){
	   
	 if (!ghostvtex) {
	  gen = (Gen*)(arr+offsetd);
	  
  PetscScalar    Eqp,Edp,delta,w; /* Generator variables */
  PetscScalar    Efd,RF,VR; /* Exciter variables */
  PetscScalar    Id,Iq;  /* Generator dq axis currents */
  //PetscScalar    Vd,Vq,SE;
  PetscScalar    IGr,IGi;
  PetscScalar    Zdq_inv[4],det;
  PetscInt       idx;
  PetscScalar    Xd, Xdp, Td0p, Xq, Xqp, Tq0p, TM, D, M, Rs ; // Generator parameters
  PetscScalar    k1,k2,KE,TE,TF,KA,KF,Vref, TA; // Generator parameters

	 
    idx = offset + 2;	 

	 /* Generator state variables */
	Eqp   = xarr[idx];
    Edp   = xarr[idx+1];
    delta = xarr[idx+2];
    w     = xarr[idx+3];
    Id    = xarr[idx+4];
    Iq    = xarr[idx+5];
    Efd   = xarr[idx+6];
    RF    = xarr[idx+7];
    VR    = xarr[idx+8];
	
	/* Generator parameters */
	Xd = gen->Xd;
	Xdp = gen->Xdp;
	Td0p = gen->Td0p;
	Xq = gen->Xq;
    Xqp = gen->Xqp;
    Tq0p = gen->Tq0p;
	TM = gen->TM;
	D = gen->D;
	M = gen->M;
	Rs = gen->Rs;
    k1 = gen->k1;
	k2 = gen->k2;
	KE = gen->KE;
	TE = gen->TE;
	TF = gen->TF;
	KA = gen->KA;
	KF = gen->KF;
	Vref = gen->Vref;
	TA = gen->TA;

    /* Generator differential equations */
    farr[idx]   = (Eqp + (Xd - Xdp)*Id - Efd)/Td0p + xdotarr[idx];
    farr[idx+1] = (Edp - (Xq - Xqp)*Iq)/Tq0p  + xdotarr[idx+1];
    farr[idx+2] = -w + w_s + xdotarr[idx+2];
    farr[idx+3] = (-TM + Edp*Id + Eqp*Iq + (Xqp - Xdp)*Id*Iq + D*(w - w_s))/M  + xdotarr[idx+3];

    Vr = xarr[offset]; /* Real part of generator terminal voltage */
    Vi = xarr[offset+1]; /* Imaginary part of the generator terminal voltage */

    ierr = ri2dq(Vr,Vi,delta,&Vd,&Vq);CHKERRQ(ierr);
    /* Algebraic equations for stator currents */
    det = Rs*Rs + Xdp*Xqp;

    Zdq_inv[0] = Rs/det;
    Zdq_inv[1] = Xqp/det;
    Zdq_inv[2] = -Xdp/det;
    Zdq_inv[3] = Rs/det;

    farr[idx+4] = Zdq_inv[0]*(-Edp + Vd) + Zdq_inv[1]*(-Eqp + Vq) + Id;
    farr[idx+5] = Zdq_inv[2]*(-Edp + Vd) + Zdq_inv[3]*(-Eqp + Vq) + Iq;

    /* Add generator current injection to network */
    ierr = dq2ri(Id,Iq,delta,&IGr,&IGi);CHKERRQ(ierr);

    farr[offset]   -= IGi;
    farr[offset+1] -= IGr;

    Vm = PetscSqrtScalar(Vd*Vd + Vq*Vq);

    SE = k1*PetscExpScalar(k2*Efd);

    /* Exciter differential equations */
    farr[idx+6] = (KE*Efd + SE - VR)/TE + xdotarr[idx+6];
    farr[idx+7] = (RF - KF*Efd/TF)/TF + xdotarr[idx+7];
    farr[idx+8] = (VR - KA*RF + KA*KF*Efd/TF - KA*(Vref - Vm))/TA + xdotarr[idx+8];

	}  
	      
   }
   
   else if (key ==3){
	
	if (!ghostvtex) {
	  
	  PetscInt k;
	  PetscInt    ld_nsegsp;
	 PetscScalar *ld_alphap;//ld_alphap=[1,0,0], an array, not a value, so use *ld_alphap;
  PetscScalar *ld_betap;
  PetscInt    ld_nsegsq;
  PetscScalar *ld_alphaq;
  PetscScalar *ld_betaq;
  PetscScalar PD0, QD0, IDr,IDi;
  
	  
	  
	  load = (Load*)(arr+offsetd);
	  
	  
	  
	  /* Load Parameters */
	  
	  ld_nsegsp = load->ld_nsegsp;
	  ld_alphap = load->ld_alphap;
	  ld_betap = load->ld_betap;
	  ld_nsegsq = load->ld_nsegsq;
	  ld_alphaq = load->ld_alphaq;
	  ld_betaq = load->ld_betaq;
	  PD0 = load->PD0;
	  QD0 = load->QD0;
	  
	  
	Vr = xarr[offset]; /* Real part of generator terminal voltage */
    Vi = xarr[offset+1]; /* Imaginary part of the generator terminal voltage */
    Vm  = PetscSqrtScalar(Vr*Vr + Vi*Vi); Vm2 = Vm*Vm;
    Vm0 = PetscSqrtScalar(Vr0*Vr0 + Vi0*Vi0);
    PD  = QD = 0.0;
    for (k=0; k < ld_nsegsp; k++) PD += ld_alphap[k]*PD0*PetscPowScalar((Vm/Vm0),ld_betap[k]);
    for (k=0; k < ld_nsegsq; k++) QD += ld_alphaq[k]*QD0*PetscPowScalar((Vm/Vm0),ld_betaq[k]);

    /* Load currents */
    IDr = (PD*Vr + QD*Vi)/Vm2;
    IDi = (-QD*Vr + PD*Vi)/Vm2;

    farr[offset]   += IDi;
    farr[offset+1] += IDr;
	}
   }
	  
	  
	  
	}
	// farr[offset]  = 0;
    // farr[offset+1] = 0;
	
  }
  

  
  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(localXdot,&xdotarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(localF,&farr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localXdot);CHKERRQ(ierr);

  // VecView(localF,PETSC_VIEWER_STDOUT_WORLD);
  
  
  ierr = DMLocalToGlobalBegin(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localF);CHKERRQ(ierr);
  PetscFunctionReturn(0);
  
  
}




/* This function is used for solving the algebraic system only during fault on and
   off times. It computes the entire F and then zeros out the part corresponding to
   differential equations
 F = [0;g(y)];
*/
 #undef __FUNCT__
#define __FUNCT__ "AlgFunction"
//PetscErrorCode AlgFunction (SNES snes, Vec X, Vec F, Userctx *user)
PetscErrorCode AlgFunction (SNES snes, Vec X, Vec F, void *ctx) // the last argument should be 'void' to be compatible with SNESSetFunction
//PetscErrorCode AlgFunction (SNES snes, Vec X, Vec F)
{
	
	
  PetscErrorCode ierr;
 
  DM             networkdm;
  Vec            localX,localF;
  //const PetscInt *cone;
  PetscInt       vfrom,vto,offsetfrom,offsetto;
  //PetscInt       type,varoffset,voffset;
  PetscInt       v,vStart,vEnd,e;
  //PetscInt       eStart,eEnd/* pipeoffset */;
  //PetscBool      ghost;
  PetscScalar    *farr;
  //PetscScalar    *vf,*juncx,*juncf;
  Userctx        *user=(Userctx*)ctx;
  
 
  //DMDALocalInfo  info;
 
  MPI_Comm       comm;
  PetscMPIInt    rank,size;
  const PetscScalar *xarr;
  DMNetworkComponentGenericDataType *arr;
    
	
  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject) snes,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr); 
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr); 
 
  ierr = VecSet(F,0.0);CHKERRQ(ierr);
  
  ierr = SNESGetDM(snes,&networkdm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localF);CHKERRQ(ierr); 
  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);
  ierr = VecSet(localF,0.0);CHKERRQ(ierr);

  /* update ghost values of locaX and locaXdot */
  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);


  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecGetArray(localF,&farr);CHKERRQ(ierr);
  
  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr); 
  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr); 
  
  for (v=vStart; v < vEnd; v++) { // Vertices contain all the info about each component (bus, generator, branch and load)
    PetscInt    i,j,offsetd, offset, key;
    PetscScalar Vr, Vi;
    Bus         *bus;
    Gen         *gen;
    Load        *load;
    PetscBool   ghostvtex;
    PetscInt    numComps;
	PetscScalar Yffr, Yffi;
	PetscScalar   Vm, Vm2,Vm0;
	PetscScalar  Vr0, Vi0;
	PetscScalar  PD,QD;

    ierr = DMNetworkIsGhostVertex(networkdm,v,&ghostvtex);CHKERRQ(ierr);  // Ghost vertices may be buses belonging to other processors whose info is needed by the local processor.
	
    ierr = DMNetworkGetNumComponents(networkdm,v,&numComps);CHKERRQ(ierr); // Buses, branches, gen and loads are the components here.
    ierr = DMNetworkGetVariableOffset(networkdm,v,&offset);CHKERRQ(ierr); // Get the offset for accessing the variable associated with the given vertex/edge from the local vector.
	
	

	
	
    for (j = 0; j < numComps; j++) {
      ierr = DMNetworkGetComponentTypeOffset(networkdm,v,j,&key,&offsetd);CHKERRQ(ierr);
      if (key == 1) {  // bus structure
	  
	  PetscInt       nconnedges;
	  const PetscInt *connedges;

	bus = (Bus*)(arr+offsetd);  // Q: How does this directly go to the bus required?
	
	if (!ghostvtex) { // Calculate current due to diagonal elements in Ybus
	  Vr = xarr[offset];
	  Vi= xarr[offset+1];
	  Yffr = bus->yff[1]+ user->ybusfault[bus->id*2+1];
	  Yffi = bus->yff[0]+ user->ybusfault[bus->id*2];
	  //Yffr = bus->yff[1];
	  //Yffi = bus->yff[0];
	  Vr0 = bus->vr;   // have used these values in the load model
	  Vi0 = bus->vi;
	  
	  farr[offset] +=  Yffi*Vr + Yffr*Vi; // imaginary current due to diagonal elements
	  farr[offset+1] += Yffr*Vr - Yffi*Vi; // real current due to diagonal elements
	 	}
		
	ierr = DMNetworkGetSupportingEdges(networkdm,v,&nconnedges,&connedges);CHKERRQ(ierr);  // nconnedges: no. of branches connected to the bus, connedges: list of branches	
		
		
		for (i=0; i < nconnedges; i++) { // loop for connected branches to the specific bus
	  Branch       *branch;
	  PetscInt       keye;
          PetscScalar    Yfti, Yftr;
          const PetscInt *cone;
		  PetscScalar  Vfr, Vfi, Vtr, Vti;

	  e = connedges[i];  // The current branch
	  ierr = DMNetworkGetComponentTypeOffset(networkdm,e,0,&keye,&offsetd);CHKERRQ(ierr);
	  branch = (Branch*)(arr+offsetd);
	  
	  Yfti = branch->yft[0];  
	  Yftr = branch->yft[1];
	  

	  ierr = DMNetworkGetConnectedNodes(networkdm,e,&cone);CHKERRQ(ierr);
	  vfrom = cone[0];  // From bus (vertex)
	  vto   = cone[1];   // To bus (vertex)

	  ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
	  ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

	  
	  //From bus and to bus real and imaginary voltages
	  Vfr     = xarr[offsetfrom];
	  Vfi     = xarr[offsetfrom+1];
	  Vtr	  = xarr[offsetto];
	  Vti     = xarr[offsetto+1];
	

	  if (vfrom == v) {
	    farr[offsetfrom]   += Yftr*Vti + Yfti*Vtr;
	    farr[offsetfrom+1] += Yftr*Vtr - Yfti*Vti;
	  } else {
	    farr[offsetto]   += Yftr*Vfi + Yfti*Vfr;
	    farr[offsetto+1] += Yftr*Vfr - Yfti*Vfi;
	  }
	  
	}

   }
   
   else if (key == 2){
	   
	 if (!ghostvtex) {
	  gen = (Gen*)(arr+offsetd);
	  
  PetscScalar    Eqp,Edp,delta,w; /* Generator variables */
  PetscScalar    Efd,RF,VR; /* Exciter variables */
  PetscScalar    Id,Iq;  /* Generator dq axis currents */
  PetscScalar    Vd,Vq,SE;
  PetscScalar    IGr,IGi;
  PetscScalar    Zdq_inv[4],det;
  PetscInt       idx;
  PetscScalar    Xd, Xdp, Td0p, Xq, Xqp, Tq0p, TM, D, M, Rs ; // Generator parameters
  PetscScalar    k1,k2,KE,TE,TF,KA,KF,Vref, TA; // Generator parameters

	 
    idx = offset + 2;	 

	 /* Generator state variables */
	Eqp   = xarr[idx];
    Edp   = xarr[idx+1];
    delta = xarr[idx+2];
    w     = xarr[idx+3];
    Id    = xarr[idx+4];
    Iq    = xarr[idx+5];
    Efd   = xarr[idx+6];
    RF    = xarr[idx+7];
    VR    = xarr[idx+8];
	
	/* Generator parameters */
	Xd = gen->Xd;
	Xdp = gen->Xdp;
	Td0p = gen->Td0p;
	Xq = gen->Xq;
    Xqp = gen->Xqp;
    Tq0p = gen->Tq0p;
	TM = gen->TM;
	D = gen->D;
	M = gen->M;
	Rs = gen->Rs;
    k1 = gen->k1;
	k2 = gen->k2;
	KE = gen->KE;
	TE = gen->TE;
	TF = gen->TF;
	KA = gen->KA;
	KF = gen->KF;
	Vref = gen->Vref;
	TA = gen->TA;

    /* Set generator differential equation residual functions to zero */
    farr[idx]   = 0 ;
    farr[idx+1] = 0 ;
    farr[idx+2] = 0 ;
    farr[idx+3] = 0;
	
	
	

    Vr = xarr[offset]; /* Real part of generator terminal voltage */
    Vi = xarr[offset+1]; /* Imaginary part of the generator terminal voltage */

    ierr = ri2dq(Vr,Vi,delta,&Vd,&Vq);CHKERRQ(ierr);
    /* Algebraic equations for stator currents */
    det = Rs*Rs + Xdp*Xqp;

    Zdq_inv[0] = Rs/det;
    Zdq_inv[1] = Xqp/det;
    Zdq_inv[2] = -Xdp/det;
    Zdq_inv[3] = Rs/det;

    farr[idx+4] = Zdq_inv[0]*(-Edp + Vd) + Zdq_inv[1]*(-Eqp + Vq) + Id;
    farr[idx+5] = Zdq_inv[2]*(-Edp + Vd) + Zdq_inv[3]*(-Eqp + Vq) + Iq;

    /* Add generator current injection to network */
    ierr = dq2ri(Id,Iq,delta,&IGr,&IGi);CHKERRQ(ierr);

    farr[offset]   -= IGi;
    farr[offset+1] -= IGr;

    Vm = PetscSqrtScalar(Vd*Vd + Vq*Vq);

    SE = k1*PetscExpScalar(k2*Efd);

   	/* Set exciter differential equation residual functions equal to zero*/
    farr[idx+6] = 0;
    farr[idx+7] = 0;
    farr[idx+8] = 0;

	}  
	      
   }
   
   else if (key ==3){
	
	if (!ghostvtex) {
	  
	  PetscInt k;
	  PetscInt    ld_nsegsp;
	 PetscScalar *ld_alphap;//ld_alphap=[1,0,0], an array, not a value, so use *ld_alphap;
  PetscScalar *ld_betap;
  PetscInt    ld_nsegsq;
  PetscScalar *ld_alphaq;
  PetscScalar *ld_betaq;
  PetscScalar PD0, QD0, IDr,IDi;
  
	  
	  
	  load = (Load*)(arr+offsetd);
	  
	  
	  
	  /* Load Parameters */
	  
	  ld_nsegsp = load->ld_nsegsp;
	  ld_alphap = load->ld_alphap;
	  ld_betap = load->ld_betap;
	  ld_nsegsq = load->ld_nsegsq;
	  ld_alphaq = load->ld_alphaq;
	  ld_betaq = load->ld_betaq;
	  PD0 = load->PD0;
	  QD0 = load->QD0;
	  
	  
	Vr = xarr[offset]; /* Real part of generator terminal voltage */
    Vi = xarr[offset+1]; /* Imaginary part of the generator terminal voltage */
    Vm  = PetscSqrtScalar(Vr*Vr + Vi*Vi); Vm2 = Vm*Vm;
    Vm0 = PetscSqrtScalar(Vr0*Vr0 + Vi0*Vi0);
    PD  = QD = 0.0;
    for (k=0; k < ld_nsegsp; k++) PD += ld_alphap[k]*PD0*PetscPowScalar((Vm/Vm0),ld_betap[k]);
    for (k=0; k < ld_nsegsq; k++) QD += ld_alphaq[k]*QD0*PetscPowScalar((Vm/Vm0),ld_betaq[k]);

    /* Load currents */
    IDr = (PD*Vr + QD*Vi)/Vm2;
    IDi = (-QD*Vr + PD*Vi)/Vm2;

    farr[offset]   += IDi;
    farr[offset+1] += IDr;
	}
   }
	  
	  
	  
	}

	
  }
  

  
  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(localF,&farr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);

  //VecView(localF,PETSC_VIEWER_STDOUT_WORLD);
  
  
  ierr = DMLocalToGlobalBegin(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localF);CHKERRQ(ierr);
  PetscFunctionReturn(0); 
  
  }
  
  



 

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char ** argv)
{
  PetscErrorCode ierr;
  // PetscInt       numEdges=9,numVertices=9;
  // int            *edges = NULL;
  PetscInt       i,j,*edgelist= NULL;;  
  // UserCtx        User;
     PetscMPIInt    size,rank;
     PetscInt       neqs_gen=0,neqs_net=0,neqs_pgrid=0,ngen=0,nbus=0,nbranch=0,nload=0;//
  // PetscInt       eStart, eEnd, vStart, vEnd,j,neqs_gen,neqs_net,neqs_pgrid;//
  Vec            X,F,F_alg,Xdot,V0;//V0: initial real and imaginary voltage of all buses
  // Mat            J,A;
 //SNES           snes;
   TS                ts;
  SNES           snes_alg;
  PetscViewer    Xview,Ybusview;
  // PetscInt       i,idx,*idx2,row_loc,col_loc;
  // PetscScalar    *x,*mat,val,*amat;
  Vec            vatol;
   Mat         Ybus; /* Network admittance matrix */
   Bus        *bus;
   Branch     *branch;
   Gen        *gen;
   Load       *load;
   //PetscScalar    sv[2];
   //PetscInt       row[1],col[2];
   //PetscScalar val[2];
   DM             networkdm;
   PetscInt       componentkey[4];
   PetscLogStage  stage1;
   PetscInt       eStart, eEnd, vStart, vEnd;
   PetscInt       genj,loadj;
   PetscInt m=0,n=0;
   //PetscReal ftime = 2.0; // Final time to iterate to
   Userctx       user;
   //TSConvergedReason reason;
   //PetscBool alg_flg;
   //PetscScalar    ybusfault[18];


  
  ierr = PetscInitialize(&argc,&argv,"petscoptions",help);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  


 /* Read initial voltage vector and Ybus */
 if (!rank) {
    ngen=3;
    nbus=9;  
    nbranch=9;
    nload=3;
    neqs_gen   = 9*ngen; /* # eqs. for generator subsystem */
    neqs_net   = 2*nbus; /* # eqs. for network subsystem   */
    neqs_pgrid = neqs_gen + neqs_net;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"X.bin",FILE_MODE_READ,&Xview);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"Ybus.bin",FILE_MODE_READ,&Ybusview);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_SELF,&V0);CHKERRQ(ierr);
  ierr = VecSetSizes(V0,PETSC_DECIDE,neqs_net);CHKERRQ(ierr);
  ierr = VecLoad(V0,Xview);CHKERRQ(ierr);
  //ierr=VecView(V0,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_SELF,&Ybus);CHKERRQ(ierr);
  ierr = MatSetSizes(Ybus,PETSC_DECIDE,PETSC_DECIDE,neqs_net,neqs_net);CHKERRQ(ierr);
  ierr = MatGetLocalSize(Ybus,&m,&n);CHKERRQ(ierr);
  
  ierr = MatSetType(Ybus,MATBAIJ);CHKERRQ(ierr);
  /*  ierr = MatSetBlockSize(Ybus,2);CHKERRQ(ierr); */
  ierr = MatLoad(Ybus,Ybusview);CHKERRQ(ierr);
  
  //MatView(Ybus, PETSC_VIEWER_STDOUT_WORLD);
  //read data
  ierr = read_data(ngen, nload, nbus, nbranch, Ybus, V0, &gen, &load, &bus, &branch, &edgelist);CHKERRQ(ierr);
  
  /* Destroy unnecessary stuff */
  ierr = PetscViewerDestroy(&Xview);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&Ybusview);CHKERRQ(ierr);
 // ierr = MatDestroy(&Ybus);CHKERRQ(ierr);
 
 
 }
  
  ierr = DMNetworkCreate(PETSC_COMM_WORLD,&networkdm);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"branchstruct",sizeof(Branch),&componentkey[0]);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"busstruct",sizeof(Bus),&componentkey[1]);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"genstruct",sizeof(Gen),&componentkey[2]);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"loadstruct",sizeof(Load),&componentkey[3]);CHKERRQ(ierr);
  
  ierr = MPI_Barrier(PETSC_COMM_WORLD);CHKERRQ(ierr);
  
  ierr = PetscLogStageRegister("Create network",&stage1);CHKERRQ(ierr);

  PetscLogStagePush(stage1);
  
  //ierr = PetscLogStageRegister("Create network",&stage1);CHKERRQ(ierr);
  //PetscLogStagePush(stage1);
  // /* Set number of nodes/edges */
  ierr = DMNetworkSetSizes(networkdm,nbus,nbranch,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);

  /* Add edge connectivity */
  ierr = DMNetworkSetEdgeList(networkdm,edgelist);CHKERRQ(ierr);
  /* Set up the network layout */
   ierr = DMNetworkLayoutSetUp(networkdm);CHKERRQ(ierr);
  
  /* We don't use these data structures anymore since they have been copied to networkdm */
  if (!rank) {
       ierr = PetscFree(edgelist);CHKERRQ(ierr);
  }
  
  // /* Add network components: physical parameters of nodes and branches*/
  if (!rank) {
    ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
    genj=0; loadj=0;
    for (i = eStart; i < eEnd; i++) {
      ierr = DMNetworkAddComponent(networkdm,i,componentkey[0],&branch[i-eStart]);CHKERRQ(ierr);
    }

    ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
    for (i = vStart; i < vEnd; i++) {
      ierr = DMNetworkAddComponent(networkdm,i,componentkey[1],&bus[i-vStart]);CHKERRQ(ierr);
      /* Add number of variables */
      ierr = DMNetworkAddNumVariables(networkdm,i,2);CHKERRQ(ierr);
      if (bus[i-vStart].nofgen) {
      for (j = 0; j < bus[i-vStart].nofgen; j++) {
	      ierr = DMNetworkAddComponent(networkdm,i,componentkey[2],&gen[genj++]);CHKERRQ(ierr);
          ierr = DMNetworkAddNumVariables(networkdm,i,9);CHKERRQ(ierr);
      }
    }
    if (bus[i-vStart].nofload) {
      for (j=0; j < bus[i-vStart].nofload; j++) {
	      ierr = DMNetworkAddComponent(networkdm,i,componentkey[3],&load[loadj++]);CHKERRQ(ierr);
      }
    }
    }
  }
  
  ierr = DMSetUp(networkdm);CHKERRQ(ierr);
  
  if (!rank) {
    ierr = PetscFree(bus);CHKERRQ(ierr);
    ierr = PetscFree(gen);CHKERRQ(ierr);
    ierr = PetscFree(branch);CHKERRQ(ierr);
    ierr = PetscFree(load);CHKERRQ(ierr);
    //ierr = PetscFree(pfdata);CHKERRQ(ierr);
    //I have already created DM, so the above are not useful any longer.
  }

  // for parallel options
  if (size > 1) {
    DM distnetworkdm;
    /* Network partitioning and distribution of data */
    ierr = DMNetworkDistribute(networkdm,0,&distnetworkdm);CHKERRQ(ierr);
    ierr = DMDestroy(&networkdm);CHKERRQ(ierr);
    networkdm = distnetworkdm;
  }
  
  PetscLogStagePop();
  
  
  ierr = DMCreateGlobalVector(networkdm,&X);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(networkdm,&Xdot);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&F);CHKERRQ(ierr);
  ierr = SetInitialGuess(networkdm, X, V0); CHKERRQ(ierr);
  //VecView(X,PETSC_VIEWER_STDOUT_WORLD);
  
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"Transient stability fault options","");CHKERRQ(ierr);
  {
    user.tfaulton  = 1.0;
    user.tfaultoff = 1.2;
    user.Rfault    = 0.0001;
    user.setisdiff = PETSC_FALSE;
    user.faultbus  = 8;
    ierr           = PetscOptionsReal("-tfaulton","","",user.tfaulton,&user.tfaulton,NULL);CHKERRQ(ierr);
    ierr           = PetscOptionsReal("-tfaultoff","","",user.tfaultoff,&user.tfaultoff,NULL);CHKERRQ(ierr);
    ierr           = PetscOptionsInt("-faultbus","","",user.faultbus,&user.faultbus,NULL);CHKERRQ(ierr);
    user.t0        = 0.0;
    user.tmax      = 5.0;
    ierr           = PetscOptionsReal("-t0","","",user.t0,&user.t0,NULL);CHKERRQ(ierr);
    ierr           = PetscOptionsReal("-tmax","","",user.tmax,&user.tmax,NULL);CHKERRQ(ierr);
    ierr           = PetscOptionsBool("-setisdiff","","",user.setisdiff,&user.setisdiff,NULL);CHKERRQ(ierr);
	
	for (i = 0; i < 18; i++) {
    user.ybusfault[i]=0;
    }
    user.ybusfault[user.faultbus*2+1]=1/user.Rfault;
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  

  // for (i = 0; i < 18; i++) {
    // ybusfault[i]=0;
    // }
    // ybusfault[user.faultbus*2+1]=1/user.Rfault;

  
  /* Setup TS solver                                           */
  /*--------------------------------------------------------*/
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
   ierr = TSSetDM(ts,(DM)networkdm);CHKERRQ(ierr);
  //ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  //ierr = TSSetEquationType(ts,TS_EQ_DAE_IMPLICIT_INDEX1);CHKERRQ(ierr);
  //ierr = TSARKIMEXSetFullyImplicit(ts,PETSC_TRUE);CHKERRQ(ierr);
  //ierr = TSSetSolution(ts,x);CHKERRQ(ierr);
  
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetIFunction(ts,NULL, (TSIFunction) FormIFunction,&user);CHKERRQ(ierr);
   //VecView(F,PETSC_VIEWER_STDOUT_WORLD);
  
  
  ierr = TSSetDuration(ts,1000,user.tfaulton);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,0.0,0.01);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  //ierr = TSSetPostStep(ts,SaveSolution);CHKERRQ(ierr);// do the save solution

  // if(user.setisdiff) {
    // const PetscInt *idx;
    // PetscScalar *vatoli;
    // PetscInt k;
 //Create vector of absolute tolerances and set the algebraic part to infinity
    // ierr = VecDuplicate(X,&vatol);CHKERRQ(ierr);
    // ierr = VecSet(X,100000.0);CHKERRQ(ierr);
    // ierr = VecGetArray(vatol,&vatoli);CHKERRQ(ierr);
    // ierr = ISGetIndices(user.is_diff,&idx);
    // for(k=0; k < 7*ngen; k++) vatoli[idx[k]] = 1e-2;
    // ierr = VecRestoreArray(vatol,&vatoli);CHKERRQ(ierr);
  // }
  
  user.alg_flg = PETSC_FALSE;
  /* Prefault period */
   ierr = TSSolve(ts,X);CHKERRQ(ierr);
   //VecView(X,PETSC_VIEWER_STDOUT_WORLD);
  // ierr = TSGetConvergedReason(ts,&reason);CHKERRQ(ierr);
  
  
  /* Create the nonlinear solver for solving the algebraic system */
  /* Note that although the algebraic system needs to be solved only for
     Idq and V, we reuse the entire system including xgen. The xgen
     variables are held constant by setting their residuals to 0 and
     putting a 1 on the Jacobian diagonal for xgen rows
  */

 
  
  
  ierr = VecDuplicate(X,&F_alg);CHKERRQ(ierr);
   ierr = TSGetSNES(ts,&snes_alg);CHKERRQ(ierr);
  //ierr = SNESCreate(PETSC_COMM_WORLD,&snes_alg);CHKERRQ(ierr);
   ierr = SNESSetFunction(snes_alg,F_alg,AlgFunction,&user);CHKERRQ(ierr);
 //ierr = MatZeroEntries(J);CHKERRQ(ierr);
 //ierr = SNESSetJacobian(snes_alg,J,J,AlgJacobian,&user);CHKERRQ(ierr);
  ierr = SNESSetOptionsPrefix(snes_alg,"alg_");CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes_alg);CHKERRQ(ierr);
  
  
  
 /* Apply disturbance - resistive fault at user.faultbus */
  /* This is done by adding shunt conductance to the diagonal location
     in the Ybus matrix */


  user.alg_flg = PETSC_TRUE;
  /* Solve the algebraic equations */
  ierr = SNESSolve(snes_alg,NULL,X);CHKERRQ(ierr);
  VecView(X,PETSC_VIEWER_STDOUT_WORLD);

  
   ierr = VecDestroy(&F_alg);CHKERRQ(ierr);
   ierr = VecDestroy(&X);CHKERRQ(ierr);
   ierr = VecDestroy(&Xdot);CHKERRQ(ierr);
   ierr = VecDestroy(&F);CHKERRQ(ierr);
   
   if (rank == 0){
   ierr = MatDestroy(&Ybus);CHKERRQ(ierr);
   ierr = VecDestroy(&V0);CHKERRQ(ierr);
   }
  // ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = DMDestroy(&networkdm);CHKERRQ(ierr);
  ierr=TSDestroy(&ts); CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
 

 }
