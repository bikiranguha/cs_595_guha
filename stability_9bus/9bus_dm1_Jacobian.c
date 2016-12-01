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
  PetscInt    neqs_pgrid;
  
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
  //MPI_Comm       comm;
  //PetscMPIInt    rank,size;
  const PetscScalar *xarr,*xdotarr;
  DMNetworkComponentGenericDataType *arr;
  PetscScalar    Vd,Vq,SE;
  //Userctx        *user=(Userctx*)ctx;
    
	
  PetscFunctionBegin;
  // ierr = PetscObjectGetComm((PetscObject)ts,&comm);CHKERRQ(ierr);
  // ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr); 
  // ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr); 
 
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
	  
	  
	  
	  
	  if (user->alg_flg == PETSC_TRUE){
		  
		  Yffr += user->ybusfault[bus->id*2+1];
		  Yffi += user->ybusfault[bus->id*2];
	  }
	  
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
 
  // MPI_Comm       comm;
  // PetscMPIInt    rank,size;
  const PetscScalar *xarr;
  DMNetworkComponentGenericDataType *arr;
    
	
  PetscFunctionBegin;
  // ierr = PetscObjectGetComm((PetscObject) snes,&comm);CHKERRQ(ierr);
  // ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr); 
  // ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr); 
 
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
	  
	   Yffr = bus->yff[1];
	   Yffi = bus->yff[0];
	  
	  if (user->alg_flg == PETSC_TRUE){
	  
	  Yffr += user->ybusfault[bus->id*2+1];
	  Yffi += user->ybusfault[bus->id*2];	  
	  }
	  
	  
	  
	  
	
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
 

/*
   J = [-df_dx, -df_dy
        dg_dx, dg_dy]
*/
#undef __FUNCT__
#define __FUNCT__ "ResidualJacobian"
PetscErrorCode ResidualJacobian(SNES snes,Vec X,Mat J,Mat B,void *ctx)
{
  PetscErrorCode ierr;
  Userctx        *user=(Userctx*)ctx;
  //Vec            Xgen,Xnet;
  //PetscScalar    *xgen,*xnet;
  PetscInt       i,idx=0;
  PetscScalar    Vr,Vi,Vm,Vm2;
  PetscScalar    Eqp,Edp,delta; /* Generator variables */
  PetscScalar    Efd;
  PetscScalar    Id,Iq;  /* Generator dq axis currents */
  PetscScalar    Vd,Vq,SE;
  PetscScalar    val[10];
  PetscInt       row[2],col[10];
  PetscInt       net_start=user->neqs_gen;
  PetscScalar    Zdq_inv[4],det;
  PetscScalar    dVd_dVr,dVd_dVi,dVq_dVr,dVq_dVi,dVd_ddelta,dVq_ddelta;
  PetscScalar    dIGr_ddelta,dIGi_ddelta,dIGr_dId,dIGr_dIq,dIGi_dId,dIGi_dIq;
  PetscScalar    dSE_dEfd;
  PetscScalar    dVm_dVd,dVm_dVq,dVm_dVr,dVm_dVi;
  PetscInt          ncols;
  const PetscInt    *cols;
  const PetscScalar *yvals;
  PetscInt          k;
  PetscScalar PD,QD,Vm0,*v0,Vm4;
  PetscScalar dPD_dVr,dPD_dVi,dQD_dVr,dQD_dVi;
  PetscScalar dIDr_dVr,dIDr_dVi,dIDi_dVr,dIDi_dVi;
   DM             networkdm;
   Vec            localX; 
    PetscInt       vfrom,vto,offsetfrom,offsetto;
  PetscInt       v,vStart,vEnd,e;
  const PetscScalar *xarr;
  DMNetworkComponentGenericDataType *arr;

    
  PetscFunctionBegin;
  //ierr  = MatZeroEntries(B);CHKERRQ(ierr);
  //ierr  = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  //ierr  = DMCompositeScatter(user->dmpgrid,X,Xgen,Xnet);CHKERRQ(ierr);

  //ierr = VecGetArray(Xgen,&xgen);CHKERRQ(ierr);
  //ierr = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);

  ierr = MatZeroEntries(J); CHKERRQ(ierr);
  //Zeros all entries of a matrix

  ierr = SNESGetDM(snes,&networkdm);CHKERRQ(ierr);
  //(input, output);(snes, *dm)
  //Gets the DM that may be used by some preconditioners
  
  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);
  // Gets a Seq PETSc vector that may be used with the DMXXX routines. 
  //(input, output)
  //(dm,*g)

  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  //update ghost values by communicating with other processes

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);
  //(x, **a)
  //Get read-only pointer to contiguous array containing this processor's portion of the vector data.

  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);
  //Returns the component data array
  //(dm, **componentdataarray)
  //（input,output）
  
  
    for (v=vStart; v < vEnd; v++) {
    PetscInt    i,j,offsetd,key;
    PetscInt    offset,goffset;
    PetscScalar Vm;
    PetscScalar Sbase=User->Sbase;
    VERTEXDATA  bus;
    PetscBool   ghostvtex;
    PetscInt    numComps;

    ierr = DMNetworkIsGhostVertex(networkdm,v,&ghostvtex);CHKERRQ(ierr);
    
    ierr = DMNetworkGetNumComponents(networkdm,v,&numComps);CHKERRQ(ierr);
    
    for (j = 0; j < numComps; j++) {
      ierr = DMNetworkGetVariableOffset(networkdm,v,&offset);CHKERRQ(ierr);
      ierr = DMNetworkGetVariableGlobalOffset(networkdm,v,&goffset);CHKERRQ(ierr);
      //Get the global offset for the variable associated with the given vertex/edge from the global vector
      //(input,input,output)
      //(,the vertex point,the offset )
      
      ierr = DMNetworkGetComponentTypeOffset(networkdm,v,j,&key,&offsetd);CHKERRQ(ierr);

// key 1 is incomplete
      if (key == 1) {
        PetscInt       nconnedges;
	const PetscInt *connedges;

	bus = (VERTEXDATA)(arr+offsetd);
    //? why values for reference bus?
	if (!ghostvtex) {
	  /* Handle reference bus constrained dofs */
	  if (bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) {
	    row[0] = goffset; row[1] = goffset+1;
	    col[0] = goffset; col[1] = goffset+1; col[2] = goffset; col[3] = goffset+1;
	    values[0] = 1.0; values[1] = 0.0; values[2] = 0.0; values[3] = 1.0;
        //We know voltage angle and magnitude for reference bus.
        //residual function for reference bus:
        //thetam=thea, that is f1=theta-thetam=0;
        //vm=v, that is f2=v-vm=0;
        //take derivative of f1 with respect to theta=1, take derivative of f1 with respect to v=0;
        //take derivative of f2 with respect to theta=0, take derivative of f2 with respect to v=1;
        
	    ierr = MatSetValues(J,2,row,2,col,values,ADD_VALUES);CHKERRQ(ierr);
        //Inserts or adds a block of values into a matrix.
        //(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv)
        // mat- the matrix
        // m, idxm- the number of rows and their global indices
       // n, idxn- the number of columns and their global indices
       // v- a logically two-dimensional array of values
       //addv- either ADD_VALUES or INSERT_VALUES, 
       //row[0]: real power, P
       //row[1]: reactive power, Q
       //MatView(J,PETSC_VIEWER_STDOUT_WORLD);
	    break;
	  }
	  
	  Vm = xarr[offset+1];//?
	  
	  /* Shunt injections */
          row[0] = goffset; row[1] = goffset+1;
          col[0] = goffset; col[1] = goffset+1;
          values[0] = values[1] = values[2] = values[3] = 0.0;
          if (bus->ide != PV_BUS) {
            values[1] = 2.0*Vm*bus->gl/Sbase;
            values[3] = -2.0*Vm*bus->bl/Sbase;
          }
          ierr = MatSetValues(J,2,row,2,col,values,ADD_VALUES);CHKERRQ(ierr);
	}

	ierr = DMNetworkGetSupportingEdges(networkdm,v,&nconnedges,&connedges);CHKERRQ(ierr);
    //Return the supporting edges for this vertex point
    //(DM dm,PetscInt vertex,PetscInt *nedges,const PetscInt *edges[])
    // Input:
    // dm- The DMNetwork object
    // vertex- the vertex point
    // Output:
    // nedges- number of edges connected to this vertex point
    // edges- List of edge points
    
    // nconnedges- number of edges connected to this vertex point
	for (i=0; i < nconnedges; i++) {
	  EDGEDATA       branch;
	  VERTEXDATA     busf,bust;
	  PetscInt       offsetfd,offsettd,keyf,keyt;
          PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
          const PetscInt *cone;
          PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;

	  e = connedges[i];
	  ierr = DMNetworkGetComponentTypeOffset(networkdm,e,0,&key,&offsetd);CHKERRQ(ierr);
      //Gets the type along with the offset for indexing the component value from the component data array
      //(DM dm,PetscInt p, PetscInt compnum, PetscInt *compkey, PetscInt *offset)
      //Input Parameters
      // dm- The DMNetwork object
      // p- vertex/edge point
      // compnum- component number
      // Output Parameters
      // compkey- the key obtained when registering the component
     // offset- offset into the component data array associated with the vertex/edge point

	  branch = (EDGEDATA)(arr+offsetd);
	  if (!branch->status) continue;
	  
	  Gff = branch->yff[0];
	  Bff = branch->yff[1];
	  Gft = branch->yft[0];
	  Bft = branch->yft[1];
	  Gtf = branch->ytf[0];
	  Btf = branch->ytf[1];
	  Gtt = branch->ytt[0];
	  Btt = branch->ytt[1];

	  ierr = DMNetworkGetConnectedNodes(networkdm,e,&cone);CHKERRQ(ierr);
      //Return the connected vertices for this edge point
      //(DM dm,PetscInt e,const PetscInt *cone[])
      // Input:
      // dm- The DMNetwork object
      // e- the edge point
      // Output:
      // cone -vertices connected to this edge 
      
	  vfrom = cone[0];
	  vto   = cone[1];

	  ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
      //Get the offset for accessing the variable associated with the given vertex/edge from the local vector.
      //(DM dm,PetscInt p,PetscInt *offset)
      // Input:
      // p- the vertex point
      // Output:
      // offset -the offset 
      
	  ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);
	  ierr = DMNetworkGetVariableGlobalOffset(networkdm,vfrom,&goffsetfrom);CHKERRQ(ierr);
      // In networkdm, first part is edge, second part is vertex
      //Get the global offset for the variable associated with the given vertex/edge from the global vector.
      //(DM dm,PetscInt p,PetscInt *offsetg)
      //Input:
      // p- the vertex point
      // Output:
      // offsetg -the offset 
      
	  ierr = DMNetworkGetVariableGlobalOffset(networkdm,vto,&goffsetto);CHKERRQ(ierr);

	  if (goffsetto < 0) goffsetto = -goffsetto - 1;

	  thetaf = xarr[offsetfrom];
	  Vmf     = xarr[offsetfrom+1];
      //why not thetat=xarr[offsetfrom+2]? Is it possible offsetfrom+1=offsetto, Vmf and thetat are repeated?
	  thetat = xarr[offsetto];
	  Vmt     = xarr[offsetto+1];
	  thetaft = thetaf - thetat;
	  thetatf = thetat - thetaf;

	  ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&keyf,&offsetfd);CHKERRQ(ierr);
      //Gets the type along with the offset for indexing the component value from the component data array
      //(DM dm,PetscInt p, PetscInt compnum, PetscInt *compkey, PetscInt *offset)
      // Input:
      // p- vertex point
      // compnum- component number
      // Output:
      // compkey- the key obtained when registering the component
      // offset- offset into the component data array associated with the vertex/edge point

	  ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&keyt,&offsettd);CHKERRQ(ierr);
	  busf = (VERTEXDATA)(arr+offsetfd);
	  bust = (VERTEXDATA)(arr+offsettd);
      // one of vfrom, vto is v
	  if (vfrom == v) {
	    if (busf->ide != REF_BUS) {
	      /*    farr[offsetfrom]   += Gff*Vmf*Vmf + Vmf*Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));  */
          // row[0] is P
          //col[0,1,2,3]=thetaf, vf, thetav, vt;
	      row[0]  = goffsetfrom;
	      col[0]  = goffsetfrom; col[1] = goffsetfrom+1; col[2] = goffsetto; col[3] = goffsetto+1;
	      values[0] =  Vmf*Vmt*(Gft*-PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft)); /* df_dthetaf */    
	      values[1] =  2.0*Gff*Vmf + Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft)); /* df_dVmf */
	      values[2] =  Vmf*Vmt*(Gft*PetscSinScalar(thetaft) + Bft*-PetscCosScalar(thetaft)); /* df_dthetat */
	      values[3] =  Vmf*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft)); /* df_dVmt */
	      
	      ierr = MatSetValues(J,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
	    }
	    if (busf->ide != PV_BUS && busf->ide != REF_BUS) {
          // row[0] is Q
          //col[0,1,2,3]=thetaf, vf, thetav, vt;
	      row[0] = goffsetfrom+1;
	      col[0]  = goffsetfrom; col[1] = goffsetfrom+1; col[2] = goffsetto; col[3] = goffsetto+1;
	      /*    farr[offsetfrom+1] += -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft)); */
	      values[0] =  Vmf*Vmt*(Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
	      values[1] =  -2.0*Bff*Vmf + Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
	      values[2] =  Vmf*Vmt*(-Bft*PetscSinScalar(thetaft) + Gft*-PetscCosScalar(thetaft));
	      values[3] =  Vmf*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
	      
	      ierr = MatSetValues(J,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
	    }
	  } else {
        // vto == v; change Bft, Gft, thetaft,Vmf to Btf, Gtf, thetatf,Vmt
	    if (bust->ide != REF_BUS) {
          // row[0] is P
          //col[0,1,2,3]=thetat, Vmt, thetaf, Vmf;
	      row[0] = goffsetto;
	      col[0] = goffsetto; col[1] = goffsetto+1; col[2] = goffsetfrom; col[3] = goffsetfrom+1;
	      /*    farr[offsetto]   += Gtt*Vmt*Vmt + Vmt*Vmf*(Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf)); */
	      values[0] =  Vmt*Vmf*(Gtf*-PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetaft)); /* df_dthetat */
	      values[1] =  2.0*Gtt*Vmt + Vmf*(Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf)); /* df_dVmt */
	      values[2] =  Vmt*Vmf*(Gtf*PetscSinScalar(thetatf) + Btf*-PetscCosScalar(thetatf)); /* df_dthetaf */
	      values[3] =  Vmt*(Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf)); /* df_dVmf */
	      
	      ierr = MatSetValues(J,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
	    }
	    if (bust->ide != PV_BUS && bust->ide != REF_BUS) {
          // row[0] is Q
          //col[0,1,2,3]=thetat, Vmt, thetaf, Vmf;
	      row[0] = goffsetto+1;
	      col[0] = goffsetto; col[1] = goffsetto+1; col[2] = goffsetfrom; col[3] = goffsetfrom+1;
	      /*    farr[offsetto+1] += -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf)); */
	      values[0] =  Vmt*Vmf*(Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
	      values[1] =  -2.0*Btt*Vmt + Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
	      values[2] =  Vmt*Vmf*(-Btf*PetscSinScalar(thetatf) + Gtf*-PetscCosScalar(thetatf));
	      values[3] =  Vmt*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
	      
	      ierr = MatSetValues(J,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
	    }
	  }
	}
    
    // for PV bus, row[0] is Q
	if (!ghostvtex && bus->ide == PV_BUS) {
	  row[0] = goffset+1; col[0] = goffset+1;
      //rpw[0] is for Q, so row[0]=goffset+1;
      //col[0] is for v, so col[0]=gooffset+1;
      //first row is P, second row is Q; first column is theta, seconde column is v.
      //if row[0] for P, row[0]=goffset;
      //if col[0] for theta, col[0]=goffset;      
	  values[0]  = 1.0;
	  ierr = MatSetValues(J,1,row,1,col,values,ADD_VALUES);CHKERRQ(ierr);
      //PV bus, voltage is fixed
      //function: vm=v;
      //residual function: f=v-vm=0;
      //take derivative of f with respect to v=1;
	}
      }
    }
  }  else if (key == 2){
	   
	 if (!ghostvtex) {
	  gen = (Gen*)(arr+offsetd);
	  
  PetscScalar    Eqp,Edp,delta,w; /* Generator variables */
  PetscScalar    Efd,RF,VR; /* Exciter variables */
  PetscScalar    Id,Iq;  /* Generator dq axis currents */
  //PetscScalar    Vd,Vq,SE;
  PetscScalar    IGr,IGi;
  PetscScalar    Zdq_inv[4],det;
  PetscInt       idx, gidx;
  PetscScalar    Xd, Xdp, Td0p, Xq, Xqp, Tq0p, TM, D, M, Rs ; // Generator parameters
  PetscScalar    k1,k2,KE,TE,TF,KA,KF,Vref, TA; // Generator parameters
  PetscScalar    vr,vi;
	 
    idx = offset + 2;  // local indices
    gidx = goffset + 2;	// global indices

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

    
    /*    fgen[idx]   = (Eqp + (Xd[i] - Xdp[i])*Id - Efd)/Td0p[i]; */
    row[0] = gidx;
    col[0] = gidx;           col[1] = gidx+4;          col[2] = gidx+6;
    val[0] = 1/ Td0p; val[1] = (Xd - Xdp)/ Td0p; val[2] = -1/Td0p;

    ierr = MatSetValues(J,1,row,3,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*    fgen[idx+1] = (Edp - (Xq[i] - Xqp[i])*Iq)/Tq0p[i]; */
    row[0] = gidx + 1;
    col[0] = gidx + 1;       col[1] = gidx+5;
    val[0] = 1/Tq0p; val[1] = -(Xq - Xqp)/Tq0p;
    ierr   = MatSetValues(J,1,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*    fgen[idx+2] = - w + w_s; */
    row[0] = gidx + 2;
    col[0] = gidx + 2; col[1] = gidx + 3;
    val[0] = 0;       val[1] = -1;
    ierr   = MatSetValues(J,1,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*    fgen[gidx+3] = (-TM[i] + Edp*Id + Eqp*Iq + (Xqp[i] - Xdp[i])*Id*Iq + D[i]*(w - w_s))/M[i]; */
    row[0] = gidx + 3;
    col[0] = gidx; col[1] = gidx + 1; col[2] = gidx + 3;       col[3] = gidx + 4;                  col[4] = gidx + 5;
    val[0] = Iq/M;  val[1] = Id/M;      val[2] = D/M; val[3] = (Edp + (Xqp-Xdp)*Iq)/M; val[4] = (Eqp + (Xqp - Xdp)*Id)/M;
    ierr   = MatSetValues(J,1,row,5,col,val,INSERT_VALUES);CHKERRQ(ierr);
    
    
    Vr   = xarr[offset]; /* Real part of generator terminal voltage */
    Vi   = xarr[offset+1]; /* Imaginary part of the generator terminal voltage */
    ierr = ri2dq(Vr,Vi,delta,&Vd,&Vq);CHKERRQ(ierr);

    det = Rs*Rs + Xdp*Xqp;

    Zdq_inv[0] = Rs/det;
    Zdq_inv[1] = Xqp/det;
    Zdq_inv[2] = -Xdp/det;
    Zdq_inv[3] = Rs/det;

    dVd_dVr    = PetscSinScalar(delta); dVd_dVi = -PetscCosScalar(delta);
    dVq_dVr    = PetscCosScalar(delta); dVq_dVi = PetscSinScalar(delta);
    dVd_ddelta = Vr*PetscCosScalar(delta) + Vi*PetscSinScalar(delta);
    dVq_ddelta = -Vr*PetscSinScalar(delta) + Vi*PetscCosScalar(delta);


    /*    fgen[idx+4] = Zdq_inv[0]*(-Edp + Vd) + Zdq_inv[1]*(-Eqp + Vq) + Id; */
    row[0] = gidx+4;
    col[0] = gidx;         col[1] = gidx+1;        col[2] = gidx + 2;
    val[0] = -Zdq_inv[1]; val[1] = -Zdq_inv[0];  val[2] = Zdq_inv[0]*dVd_ddelta + Zdq_inv[1]*dVq_ddelta;
    col[3] = gidx + 4; col[4] = net_start+2*gbus[i];                     col[5] = net_start + 2*gbus[i]+1;  // col[4] and col[5] are with respect to Vr and Vi
    val[3] = 1;       val[4] = Zdq_inv[0]*dVd_dVr + Zdq_inv[1]*dVq_dVr; val[5] = Zdq_inv[0]*dVd_dVi + Zdq_inv[1]*dVq_dVi;
    ierr   = MatSetValues(J,1,row,6,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*  fgen[idx+5] = Zdq_inv[2]*(-Edp + Vd) + Zdq_inv[3]*(-Eqp + Vq) + Iq; */
    row[0] = gidx+5;
    col[0] = gidx;         col[1] = gidx+1;        col[2] = gidx + 2;
    val[0] = -Zdq_inv[3]; val[1] = -Zdq_inv[2];  val[2] = Zdq_inv[2]*dVd_ddelta + Zdq_inv[3]*dVq_ddelta;
    col[3] = gidx + 5; col[4] = net_start+2*gbus[i];                     col[5] = net_start + 2*gbus[i]+1;
    val[3] = 1;       val[4] = Zdq_inv[2]*dVd_dVr + Zdq_inv[3]*dVq_dVr; val[5] = Zdq_inv[2]*dVd_dVi + Zdq_inv[3]*dVq_dVi;
    ierr   = MatSetValues(J,1,row,6,col,val,INSERT_VALUES);CHKERRQ(ierr);
	
	
	dIGr_ddelta = Id*PetscCosScalar(delta) - Iq*PetscSinScalar(delta);
    dIGi_ddelta = Id*PetscSinScalar(delta) + Iq*PetscCosScalar(delta);
    dIGr_dId    = PetscSinScalar(delta);  dIGr_dIq = PetscCosScalar(delta);
    dIGi_dId    = -PetscCosScalar(delta); dIGi_dIq = PetscSinScalar(delta);
	
	
	/* fnet[goffset]   -= IGi; */
    row[0] = goffset;
    col[0] = gidx+2;        col[1] = gidx + 4;   col[2] = gidx + 5;
    val[0] = -dIGi_ddelta; val[1] = -dIGi_dId; val[2] = -dIGi_dIq;
    ierr = MatSetValues(J,1,row,3,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /* fnet[goffset+1]   -= IGr; */
    row[0] = goffset+1;
    col[0] = gidx+2;        col[1] = gidx + 4;   col[2] = gidx + 5;
    val[0] = -dIGr_ddelta; val[1] = -dIGr_dId; val[2] = -dIGr_dIq;
    ierr   = MatSetValues(J,1,row,3,col,val,INSERT_VALUES);CHKERRQ(ierr);

    



    Vm = PetscSqrtScalar(Vd*Vd + Vq*Vq);

   /*    fgen[idx+6] = (KE*Efd + SE - VR)/TE; */
    /*    SE  = k1*PetscExpScalar(k2*Efd); */

    dSE_dEfd = k1*k2*PetscExpScalar(k2*Efd);

    row[0] = gidx + 6;
    col[0] = gidx + 6;                     col[1] = gidx + 8;
    val[0] = (KE + dSE_dEfd)/TE;  val[1] = -1/TE;
    ierr   = MatSetValues(J,1,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /* Exciter differential equations */

    /*    fgen[idx+7] = (RF - KF*Efd/TF)/TF; */
    row[0] = gidx + 7;
    col[0] = gidx + 6;       col[1] = gidx + 7;
    val[0] = (-KF/TF)/TF;  val[1] = 1/TF;
    ierr   = MatSetValues(J,1,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*    fgen[idx+8] = (VR - KA*RF + KA*KF*Efd/TF - KA*(Vref - Vm))/TA; */
    /* Vm = (Vd^2 + Vq^2)^0.5; */

    dVm_dVd    = Vd/Vm; dVm_dVq = Vq/Vm;
    dVm_dVr    = dVm_dVd*dVd_dVr + dVm_dVq*dVq_dVr;
    dVm_dVi    = dVm_dVd*dVd_dVi + dVm_dVq*dVq_dVi;
    row[0]     = gidx + 8;
    col[0]     = gidx + 6;           col[1] = gidx + 7; col[2] = gidx + 8;
    val[0]     = (KA*KF/TF)/TA; val[1] = -KA/TA;  val[2] = 1/TA;
    col[3]     = goffset; col[4] = goffset+1;
    val[3]     = KA*dVm_dVr/TA;         val[4] = KA*dVm_dVi/TA;
    ierr       = MatSetValues(J,1,row,5,col,val,INSERT_VALUES);CHKERRQ(ierr);
 
	
	//stop here

   
	}  
	      
   }  else if (key ==3){
	
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
  
  
  /* Generator subsystem */
  for (i=0; i < ngen; i++) {
    Eqp   = xgen[idx];
    Edp   = xgen[idx+1];
    delta = xgen[idx+2];
    Id    = xgen[idx+4];
    Iq    = xgen[idx+5];
    Efd   = xgen[idx+6];

    

    Vr   = xnet[2*gbus[i]]; /* Real part of generator terminal voltage */
    Vi   = xnet[2*gbus[i]+1]; /* Imaginary part of the generator terminal voltage */
    ierr = ri2dq(Vr,Vi,delta,&Vd,&Vq);CHKERRQ(ierr);

    det = Rs[i]*Rs[i] + Xdp[i]*Xqp[i];

    Zdq_inv[0] = Rs[i]/det;
    Zdq_inv[1] = Xqp[i]/det;
    Zdq_inv[2] = -Xdp[i]/det;
    Zdq_inv[3] = Rs[i]/det;

    dVd_dVr    = PetscSinScalar(delta); dVd_dVi = -PetscCosScalar(delta);
    dVq_dVr    = PetscCosScalar(delta); dVq_dVi = PetscSinScalar(delta);
    dVd_ddelta = Vr*PetscCosScalar(delta) + Vi*PetscSinScalar(delta);
    dVq_ddelta = -Vr*PetscSinScalar(delta) + Vi*PetscCosScalar(delta);

    /*    fgen[idx+4] = Zdq_inv[0]*(-Edp + Vd) + Zdq_inv[1]*(-Eqp + Vq) + Id; */
    row[0] = idx+4;
    col[0] = idx;         col[1] = idx+1;        col[2] = idx + 2;
    val[0] = -Zdq_inv[1]; val[1] = -Zdq_inv[0];  val[2] = Zdq_inv[0]*dVd_ddelta + Zdq_inv[1]*dVq_ddelta;
    col[3] = idx + 4; col[4] = net_start+2*gbus[i];                     col[5] = net_start + 2*gbus[i]+1;  // col[4] and col[5] are with respect to Vr and Vi
    val[3] = 1;       val[4] = Zdq_inv[0]*dVd_dVr + Zdq_inv[1]*dVq_dVr; val[5] = Zdq_inv[0]*dVd_dVi + Zdq_inv[1]*dVq_dVi;
    ierr   = MatSetValues(J,1,row,6,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*  fgen[idx+5] = Zdq_inv[2]*(-Edp + Vd) + Zdq_inv[3]*(-Eqp + Vq) + Iq; */
    row[0] = idx+5;
    col[0] = idx;         col[1] = idx+1;        col[2] = idx + 2;
    val[0] = -Zdq_inv[3]; val[1] = -Zdq_inv[2];  val[2] = Zdq_inv[2]*dVd_ddelta + Zdq_inv[3]*dVq_ddelta;
    col[3] = idx + 5; col[4] = net_start+2*gbus[i];                     col[5] = net_start + 2*gbus[i]+1;
    val[3] = 1;       val[4] = Zdq_inv[2]*dVd_dVr + Zdq_inv[3]*dVq_dVr; val[5] = Zdq_inv[2]*dVd_dVi + Zdq_inv[3]*dVq_dVi;
    ierr   = MatSetValues(J,1,row,6,col,val,INSERT_VALUES);CHKERRQ(ierr);

    dIGr_ddelta = Id*PetscCosScalar(delta) - Iq*PetscSinScalar(delta);
    dIGi_ddelta = Id*PetscSinScalar(delta) + Iq*PetscCosScalar(delta);
    dIGr_dId    = PetscSinScalar(delta);  dIGr_dIq = PetscCosScalar(delta);
    dIGi_dId    = -PetscCosScalar(delta); dIGi_dIq = PetscSinScalar(delta);

    /* fnet[2*gbus[i]]   -= IGi; */
    row[0] = net_start + 2*gbus[i];
    col[0] = idx+2;        col[1] = idx + 4;   col[2] = idx + 5;
    val[0] = -dIGi_ddelta; val[1] = -dIGi_dId; val[2] = -dIGi_dIq;
    ierr = MatSetValues(J,1,row,3,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /* fnet[2*gbus[i]+1]   -= IGr; */
    row[0] = net_start + 2*gbus[i]+1;
    col[0] = idx+2;        col[1] = idx + 4;   col[2] = idx + 5;
    val[0] = -dIGr_ddelta; val[1] = -dIGr_dId; val[2] = -dIGr_dIq;
    ierr   = MatSetValues(J,1,row,3,col,val,INSERT_VALUES);CHKERRQ(ierr);

    Vm = PetscSqrtScalar(Vd*Vd + Vq*Vq);

    /*    fgen[idx+6] = (KE[i]*Efd + SE - VR)/TE[i]; */
    /*    SE  = k1[i]*PetscExpScalar(k2[i]*Efd); */

    dSE_dEfd = k1[i]*k2[i]*PetscExpScalar(k2[i]*Efd);

    row[0] = idx + 6;
    col[0] = idx + 6;                     col[1] = idx + 8;
    val[0] = (KE[i] + dSE_dEfd)/TE[i];  val[1] = -1/TE[i];
    ierr   = MatSetValues(J,1,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /* Exciter differential equations */

    /*    fgen[idx+7] = (RF - KF[i]*Efd/TF[i])/TF[i]; */
    row[0] = idx + 7;
    col[0] = idx + 6;       col[1] = idx + 7;
    val[0] = (-KF[i]/TF[i])/TF[i];  val[1] = 1/TF[i];
    ierr   = MatSetValues(J,1,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*    fgen[idx+8] = (VR - KA[i]*RF + KA[i]*KF[i]*Efd/TF[i] - KA[i]*(Vref[i] - Vm))/TA[i]; */
    /* Vm = (Vd^2 + Vq^2)^0.5; */

    dVm_dVd    = Vd/Vm; dVm_dVq = Vq/Vm;
    dVm_dVr    = dVm_dVd*dVd_dVr + dVm_dVq*dVq_dVr;
    dVm_dVi    = dVm_dVd*dVd_dVi + dVm_dVq*dVq_dVi;
    row[0]     = idx + 8;
    col[0]     = idx + 6;           col[1] = idx + 7; col[2] = idx + 8;
    val[0]     = (KA[i]*KF[i]/TF[i])/TA[i]; val[1] = -KA[i]/TA[i];  val[2] = 1/TA[i];
    col[3]     = net_start + 2*gbus[i]; col[4] = net_start + 2*gbus[i]+1;
    val[3]     = KA[i]*dVm_dVr/TA[i];         val[4] = KA[i]*dVm_dVi/TA[i];
    ierr       = MatSetValues(J,1,row,5,col,val,INSERT_VALUES);CHKERRQ(ierr);
    idx        = idx + 9;
  }

  for (i=0; i<nbus; i++) {// allocating Ybus indices in J for branch currents
    ierr   = MatGetRow(user->Ybus,2*i,&ncols,&cols,&yvals);CHKERRQ(ierr);
    row[0] = net_start + 2*i;
    for (k=0; k<ncols; k++) {// adding non zero columns for the branch currents elements in Ir (real current balance equations)
      col[k] = net_start + cols[k];
      val[k] = yvals[k];
    }
    ierr = MatSetValues(J,1,row,ncols,col,val,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatRestoreRow(user->Ybus,2*i,&ncols,&cols,&yvals);CHKERRQ(ierr);

    ierr   = MatGetRow(user->Ybus,2*i+1,&ncols,&cols,&yvals);CHKERRQ(ierr);
    row[0] = net_start + 2*i+1;
    for (k=0; k<ncols; k++) {// adding non zero columns for the branch currents elements in Ii (imaginary current balance equations)
      col[k] = net_start + cols[k];
      val[k] = yvals[k];
    }
    ierr = MatSetValues(J,1,row,ncols,col,val,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatRestoreRow(user->Ybus,2*i+1,&ncols,&cols,&yvals);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(J,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);

  ierr = VecGetArray(user->V0,&v0);CHKERRQ(ierr);
  for (i=0; i < nload; i++) {
    Vr      = xnet[2*lbus[i]]; /* Real part of load bus voltage */
    Vi      = xnet[2*lbus[i]+1]; /* Imaginary part of the load bus voltage */
    Vm      = PetscSqrtScalar(Vr*Vr + Vi*Vi); Vm2 = Vm*Vm; Vm4 = Vm2*Vm2;
    Vm0     = PetscSqrtScalar(v0[2*lbus[i]]*v0[2*lbus[i]] + v0[2*lbus[i]+1]*v0[2*lbus[i]+1]);
    PD      = QD = 0.0;
    dPD_dVr = dPD_dVi = dQD_dVr = dQD_dVi = 0.0;
    for (k=0; k < ld_nsegsp[i]; k++) {  // 
      PD      += ld_alphap[k]*PD0[i]*PetscPowScalar((Vm/Vm0),ld_betap[k]);
      dPD_dVr += ld_alphap[k]*ld_betap[k]*PD0[i]*PetscPowScalar((1/Vm0),ld_betap[k])*Vr*PetscPowScalar(Vm,(ld_betap[k]-2));
      dPD_dVi += ld_alphap[k]*ld_betap[k]*PD0[i]*PetscPowScalar((1/Vm0),ld_betap[k])*Vi*PetscPowScalar(Vm,(ld_betap[k]-2));
    }
    for (k=0; k < ld_nsegsq[i]; k++) {
      QD      += ld_alphaq[k]*QD0[i]*PetscPowScalar((Vm/Vm0),ld_betaq[k]);
      dQD_dVr += ld_alphaq[k]*ld_betaq[k]*QD0[i]*PetscPowScalar((1/Vm0),ld_betaq[k])*Vr*PetscPowScalar(Vm,(ld_betaq[k]-2));
      dQD_dVi += ld_alphaq[k]*ld_betaq[k]*QD0[i]*PetscPowScalar((1/Vm0),ld_betaq[k])*Vi*PetscPowScalar(Vm,(ld_betaq[k]-2));
    }

    /*    IDr = (PD*Vr + QD*Vi)/Vm2; */
    /*    IDi = (-QD*Vr + PD*Vi)/Vm2; */

    dIDr_dVr = (dPD_dVr*Vr + dQD_dVr*Vi + PD)/Vm2 - ((PD*Vr + QD*Vi)*2*Vr)/Vm4;
    dIDr_dVi = (dPD_dVi*Vr + dQD_dVi*Vi + QD)/Vm2 - ((PD*Vr + QD*Vi)*2*Vi)/Vm4;

    dIDi_dVr = (-dQD_dVr*Vr + dPD_dVr*Vi - QD)/Vm2 - ((-QD*Vr + PD*Vi)*2*Vr)/Vm4;
    dIDi_dVi = (-dQD_dVi*Vr + dPD_dVi*Vi + PD)/Vm2 - ((-QD*Vr + PD*Vi)*2*Vi)/Vm4;


    /*    fnet[2*lbus[i]]   += IDi; */
    row[0] = net_start + 2*lbus[i];
    col[0] = net_start + 2*lbus[i];  col[1] = net_start + 2*lbus[i]+1;
    val[0] = dIDi_dVr;               val[1] = dIDi_dVi;
    ierr   = MatSetValues(J,1,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
    /*    fnet[2*lbus[i]+1] += IDr; */
    row[0] = net_start + 2*lbus[i]+1;
    col[0] = net_start + 2*lbus[i];  col[1] = net_start + 2*lbus[i]+1;
    val[0] = dIDr_dVr;               val[1] = dIDr_dVi;
    ierr   = MatSetValues(J,1,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(user->V0,&v0);CHKERRQ(ierr);

  ierr = VecRestoreArray(Xgen,&xgen);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xnet,&xnet);CHKERRQ(ierr);

  ierr = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   J = [I, 0
        dg_dx, dg_dy]     Q: Why is df_dx = 1 for diagonal elements?
*/
#undef __FUNCT__
#define __FUNCT__ "AlgJacobian"
PetscErrorCode AlgJacobian(SNES snes,Vec X,Mat A,Mat B,void *ctx)
{
  PetscErrorCode ierr;
  Userctx        *user=(Userctx*)ctx;

  PetscFunctionBegin;
  ierr = ResidualJacobian(snes,X,A,B,ctx);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatZeroRowsIS(A,user->is_diff,1.0,NULL,NULL);CHKERRQ(ierr); // zeroes all the rows related to differential equations and put '1' for diagonal elements
  PetscFunctionReturn(0);
}
 
  
/*
   J = [a*I-df_dx, -df_dy
        dg_dx, dg_dy]                Q: Dont understand this function?
*/

#undef __FUNCT__
#define __FUNCT__ "IJacobian"
PetscErrorCode IJacobian(TS ts,PetscReal t,Vec X,Vec Xdot,PetscReal a,Mat A,Mat B,Userctx *user)
{
  PetscErrorCode ierr;
  SNES           snes;
  PetscScalar    atmp = (PetscScalar) a;
  PetscInt       i,row;

  PetscFunctionBegin;
  user->t = t;

  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = ResidualJacobian(snes,X,A,B,user);CHKERRQ(ierr);
  for (i=0;i < ngen;i++) {
    row = 9*i;
    ierr = MatSetValues(A,1,&row,1,&row,&atmp,ADD_VALUES);CHKERRQ(ierr);
    row  = 9*i+1;
    ierr = MatSetValues(A,1,&row,1,&row,&atmp,ADD_VALUES);CHKERRQ(ierr);
    row  = 9*i+2;
    ierr = MatSetValues(A,1,&row,1,&row,&atmp,ADD_VALUES);CHKERRQ(ierr);
    row  = 9*i+3;
    ierr = MatSetValues(A,1,&row,1,&row,&atmp,ADD_VALUES);CHKERRQ(ierr);
    row  = 9*i+6;
    ierr = MatSetValues(A,1,&row,1,&row,&atmp,ADD_VALUES);CHKERRQ(ierr);
    row  = 9*i+7;
    ierr = MatSetValues(A,1,&row,1,&row,&atmp,ADD_VALUES);CHKERRQ(ierr);
    row  = 9*i+8;
    ierr = MatSetValues(A,1,&row,1,&row,&atmp,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
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
  
  //Jacobian
  ierr = DMCreateMatrix(networkdm,&J);CHKERRQ(ierr);
 // Gets empty Jacobian for a DM
 // This properly preallocates the number of nonzeros in the sparse matrix so you do not need to do it yourself.
  ierr = MatSetOption(J,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);
  
  
  ierr = SetInitialGuess(networkdm, X, V0); CHKERRQ(ierr);
  //VecView(X,PETSC_VIEWER_STDOUT_WORLD);
  
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"Transient stability fault options","");CHKERRQ(ierr);
  {
    user.tfaulton  = 1.0;
    user.tfaultoff = 1.2;
    user.Rfault    = 0.0001;
    user.setisdiff = PETSC_FALSE;
    user.faultbus  = 8;
     neqs_gen   = 9*ngen; /* # eqs. for generator subsystem */
    neqs_net   = 2*nbus; /* # eqs. for network subsystem   */
    user.neqs_pgrid=9*3+2*9;
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
  ierr = TSSetApplicationContext(ts,&user);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetIFunction(ts,NULL, (TSIFunction) FormIFunction,&user);CHKERRQ(ierr);
  
  //Jacobian
   ierr = TSSetIJacobian(ts,J,J,(TSIJacobian)IJacobian,&user);CHKERRQ(ierr);
   //VecView(F,PETSC_VIEWER_STDOUT_WORLD);
  
  
  ierr = TSSetDuration(ts,1000,user.tfaulton);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,0.0,0.01);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  //ierr = TSSetPostStep(ts,SaveSolution);CHKERRQ(ierr);// do the save solution

  
  user.alg_flg = PETSC_FALSE;
  /* Prefault period */
   ierr = TSSolve(ts,X);CHKERRQ(ierr);
  // VecView(X,PETSC_VIEWER_STDOUT_WORLD);
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
  //VecView(X,PETSC_VIEWER_STDOUT_WORLD);
  
  
  
  
  /* Disturbance period */
  ierr = TSSetDuration(ts,1000,user.tfaultoff);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,user.tfaulton,.01);CHKERRQ(ierr);
  ierr = TSSetIFunction(ts,NULL, (TSIFunction) FormIFunction,&user);CHKERRQ(ierr);

  user.alg_flg = PETSC_TRUE;

  ierr = TSSolve(ts,X);CHKERRQ(ierr);
  //VecView(X,PETSC_VIEWER_STDOUT_WORLD);
  
  
  
  
  /* Remove the fault */

  ierr = SNESSetFunction(snes_alg,F_alg,AlgFunction,&user);CHKERRQ(ierr);
 //ierr = MatZeroEntries(J);CHKERRQ(ierr);
 //ierr = SNESSetJacobian(snes_alg,J,J,AlgJacobian,&user);CHKERRQ(ierr);
  ierr = SNESSetOptionsPrefix(snes_alg,"alg_");CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes_alg);CHKERRQ(ierr);
  


  user.alg_flg = PETSC_FALSE;
  /* Solve the algebraic equations */
  ierr = SNESSolve(snes_alg,NULL,X);CHKERRQ(ierr);
  //VecView(X,PETSC_VIEWER_STDOUT_WORLD);
  
  
  
  
  /* Post-disturbance period */
  ierr = TSSetDuration(ts,1000,user.tmax);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,user.tfaultoff,.01);CHKERRQ(ierr);
  ierr = TSSetIFunction(ts,NULL, (TSIFunction) FormIFunction,&user);CHKERRQ(ierr);

  user.alg_flg = PETSC_FALSE;

  ierr = TSSolve(ts,X);CHKERRQ(ierr);
  //VecView(X,PETSC_VIEWER_STDOUT_WORLD);
  
  
  
  

  
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
