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

#include "pf.h"
#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>

#define freq 60
#define w_s (2*PETSC_PI*freq)

PetscMPIInt rank;

/* Sizes and indices */
const PetscInt nbus    = 9; /* Number of network buses */
const PetscInt ngen    = 3; /* Number of generators */
const PetscInt nload   = 3; /* Number of loads */
const PetscInt gbus[3] = {0,1,2}; /* Buses at which generators are incident */
const PetscInt lbus[3] = {4,5,7}; /* Buses at which loads are incident */

/* Generator real and reactive powers (found via loadflow) */
const PetscScalar PG[3] = {0.716786142395021,1.630000000000000,0.850000000000000};
const PetscScalar QG[3] = {0.270702180178785,0.066120127797275,-0.108402221791588};
/* Generator constants */
const PetscScalar H[3]    = {23.64,6.4,3.01};   /* Inertia constant */
const PetscScalar Rs[3]   = {0.0,0.0,0.0}; /* Stator Resistance */
const PetscScalar Xd[3]   = {0.146,0.8958,1.3125};  /* d-axis reactance */
const PetscScalar Xdp[3]  = {0.0608,0.1198,0.1813}; /* d-axis transient reactance */
const PetscScalar Xq[3]   = {0.4360,0.8645,1.2578}; /* q-axis reactance Xq(1) set to 0.4360, value given in text 0.0969 */
const PetscScalar Xqp[3]  = {0.0969,0.1969,0.25}; /* q-axis transient reactance */
const PetscScalar Td0p[3] = {8.96,6.0,5.89}; /* d-axis open circuit time constant */
const PetscScalar Tq0p[3] = {0.31,0.535,0.6}; /* q-axis open circuit time constant */
PetscScalar M[3]; /* M = 2*H/w_s */
PetscScalar D[3]; /* D = 0.1*M */

PetscScalar TM[3]; /* Mechanical Torque */
/* Exciter system constants */
const PetscScalar KA[3] = {20.0,20.0,20.0};  /* Voltage regulartor gain constant */
const PetscScalar TA[3] = {0.2,0.2,0.2};     /* Voltage regulator time constant */
const PetscScalar KE[3] = {1.0,1.0,1.0};     /* Exciter gain constant */
const PetscScalar TE[3] = {0.314,0.314,0.314}; /* Exciter time constant */
const PetscScalar KF[3] = {0.063,0.063,0.063};  /* Feedback stabilizer gain constant */
const PetscScalar TF[3] = {0.35,0.35,0.35};    /* Feedback stabilizer time constant */
const PetscScalar k1[3] = {0.0039,0.0039,0.0039};
const PetscScalar k2[3] = {1.555,1.555,1.555};  /* k1 and k2 for calculating the saturation function SE = k1*exp(k2*Efd) */

PetscScalar Vref[3];
/* Load constants
  We use a composite load model that describes the load and reactive powers at each time instant as follows
  P(t) = \sum\limits_{i=0}^ld_nsegsp \ld_alphap_i*P_D0(\frac{V_m(t)}{V_m0})^\ld_betap_i
  Q(t) = \sum\limits_{i=0}^ld_nsegsq \ld_alphaq_i*Q_D0(\frac{V_m(t)}{V_m0})^\ld_betaq_i
  where
    ld_nsegsp,ld_nsegsq - Number of individual load models for real and reactive power loads
    ld_alphap,ld_alphap - Percentage contribution (weights) or loads
    P_D0                - Real power load
    Q_D0                - Reactive power load
    V_m(t)              - Voltage magnitude at time t
    V_m0                - Voltage magnitude at t = 0
    ld_betap, ld_betaq  - exponents describing the load model for real and reactive part

    Note: All loads have the same characteristic currently.
*/
const PetscScalar PD0[3] = {1.25,0.9,1.0};
const PetscScalar QD0[3] = {0.5,0.3,0.35};
const PetscInt    ld_nsegsp[3] = {3,3,3};
const PetscScalar ld_alphap[3] = {1.0,0.0,0.0};
const PetscScalar ld_betap[3]  = {2.0,1.0,0.0};
const PetscInt    ld_nsegsq[3] = {3,3,3};
const PetscScalar ld_alphaq[3] = {1.0,0.0,0.0};
const PetscScalar ld_betaq[3]  = {2.0,1.0,0.0};

/* #undef __FUNCT__
#define __FUNCT__ "GetListofEdges"
PetscErrorCode GetListofEdges(PetscInt nbranches, EDGEDATA branch,int edges[])
{
  PetscInt       i, fbus,tbus;

  PetscFunctionBegin;
  for (i=0; i < nbranches; i++) {
    fbus = branch[i].internal_i;
    tbus = branch[i].internal_j;
    edges[2*i]   = fbus;
    edges[2*i+1] = tbus;
  }
  PetscFunctionReturn(0);
} */

/* typedef struct{
  PetscScalar  Sbase;
}UserCtx; */




/* #undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction(SNES snes,Vec X, Vec F,void *appctx) // Q: Why is it used in such a way in SNESSetFunction (line 572) where the arguments cannot be supplied?  */

#undef __FUNCT__
#define __FUNCT__ "IFunction"
PetscErrorCode IFunction(TS ts,PetscReal t,Vec X,Vec Xdot,Vec F,void* ctx)
{
  PetscErrorCode ierr;
  DM             networkdm;
  UserCtx       *User=(UserCtx*)appctx;
  Vec           localX,localF;
  PetscInt      e;
  PetscInt      v,vStart,vEnd,vfrom,vto;
  const PetscScalar *xarr;
  PetscScalar   *farr;
  PetscInt      offsetfrom,offsetto,offset;
  DMNetworkComponentGenericDataType *arr; // arr is a pointer to the start of the local array of data components.  
  PetscFunctionBegin;
  ierr = SNESGetDM(snes,&networkdm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr); 
  ierr = DMGetLocalVector(networkdm,&localF);CHKERRQ(ierr);
  ierr = VecSet(F,0.0);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);  // From the global vector X, the local vector 'localX' are formed.
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr); // xarr is a pointer to localX, to access its contents using xarr[n]. 
  ierr = VecGetArray(localF,&farr);CHKERRQ(ierr);     // In VecGetArray, you can modify contents of the vector, in VecGetArrayRead, you cannot modify.

  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr); 
  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr); // Q: What is now inside arr? 

  for (v=vStart; v < vEnd; v++) { // Vertices contain all the info about each component (bus, generator, branch and load)
    PetscInt    i,j,offsetd,key;
    PetscScalar Vm;
    PetscScalar Sbase=User->Sbase;
    VERTEXDATA  bus=NULL;
    GEN         gen;
    LOAD        load;
    PetscBool   ghostvtex;
    PetscInt    numComps;

    ierr = DMNetworkIsGhostVertex(networkdm,v,&ghostvtex);CHKERRQ(ierr);  // Ghost vertices may be buses belonging to other processors whose info is needed by the local processor.
    ierr = DMNetworkGetNumComponents(networkdm,v,&numComps);CHKERRQ(ierr); // Buses, branches, gen and loads are the components here.
    ierr = DMNetworkGetVariableOffset(networkdm,v,&offset);CHKERRQ(ierr); // Get the offset for accessing the variable associated with the given vertex/edge from the local vector.
    for (j = 0; j < numComps; j++) {
      ierr = DMNetworkGetComponentTypeOffset(networkdm,v,j,&key,&offsetd);CHKERRQ(ierr);
      if (key == 1) {  // bus structure
        PetscInt       nconnedges;
	const PetscInt *connedges;

	bus = (VERTEXDATA)(arr+offsetd);  // Q: How does this directly go to the bus required?
	/* Handle reference bus constrained dofs */
	if (bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) {
	  farr[offset] = xarr[offset] - bus->va*PETSC_PI/180.0;  // bus->va, xarr[offset]: Phase angle, farr[offset]  = Real power equation
	  farr[offset+1] = xarr[offset+1] - bus->vm;  // bus->vm, xarr[offset+1]:Voltage magnitude, farr[offset]  = Reactive power equation
	  break;
	}

	if (!ghostvtex) {
	  Vm = xarr[offset+1];

	  /* Shunt injections, account for any shunt elements */
	  farr[offset] += Vm*Vm*bus->gl/Sbase; 
	  if(bus->ide != PV_BUS) farr[offset+1] += -Vm*Vm*bus->bl/Sbase; // Q: Why not for a PV bus?
	}

	ierr = DMNetworkGetSupportingEdges(networkdm,v,&nconnedges,&connedges);CHKERRQ(ierr);  // nconnedges: no. of branches connected to the bus, connedges: list of branches
	for (i=0; i < nconnedges; i++) { // loop for connected branches to the specific bus
	  EDGEDATA       branch;
	  PetscInt       keye;
          PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
          const PetscInt *cone;
          PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;

	  e = connedges[i];  // The current branch
	  ierr = DMNetworkGetComponentTypeOffset(networkdm,e,0,&keye,&offsetd);CHKERRQ(ierr);
	  branch = (EDGEDATA)(arr+offsetd);
	  if (!branch->status) continue;
	  Gff = branch->yff[0];  // ff, tt: diagonal elements, tf, ft: off-diagonal elements
	  Bff = branch->yff[1];
	  Gft = branch->yft[0];
	  Bft = branch->yft[1];
	  Gtf = branch->ytf[0];
	  Btf = branch->ytf[1];
	  Gtt = branch->ytt[0];
	  Btt = branch->ytt[1];

	  ierr = DMNetworkGetConnectedNodes(networkdm,e,&cone);CHKERRQ(ierr);
	  vfrom = cone[0];  // From bus voltage
	  vto   = cone[1];   // To bus voltage

	  ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
	  ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

	  thetaf = xarr[offsetfrom];
	  Vmf     = xarr[offsetfrom+1];
	  thetat = xarr[offsetto];
	  Vmt     = xarr[offsetto+1];
	  thetaft = thetaf - thetat; // phase angle difference
	  thetatf = thetat - thetaf;

	  if (vfrom == v) {
	    farr[offsetfrom]   += Gff*Vmf*Vmf + Vmf*Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
	    farr[offsetfrom+1] += -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
	  } else {
	    farr[offsetto]   += Gtt*Vmt*Vmt + Vmt*Vmf*(Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
	    farr[offsetto+1] += -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
	  }
	}
      } else if (key == 2) {  // generator structure
	if (!ghostvtex) {
	  gen = (GEN)(arr+offsetd);
	  if (!gen->status) continue;
	  farr[offset] += -gen->pg/Sbase;
	  farr[offset+1] += -gen->qg/Sbase;
	}
      } else if (key == 3) { // load structure
	if (!ghostvtex) {
	  load = (LOAD)(arr+offsetd);
	  farr[offset] += load->pl/Sbase;
	  farr[offset+1] += load->ql/Sbase;
	}
      }
    }
    if (bus && bus->ide == PV_BUS) { // For a PV bus, Vm is supplied, and real power is given
      farr[offset+1] = xarr[offset+1] - bus->vm;
    }
  }
  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(localF,&farr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localF);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian(SNES snes,Vec X, Mat J,Mat Jpre,void *appctx)
{
  PetscErrorCode ierr;
  DM            networkdm;
  UserCtx       *User=(UserCtx*)appctx;
  Vec           localX;
  PetscInt      e;
  PetscInt      v,vStart,vEnd,vfrom,vto;
  const PetscScalar *xarr;
  PetscInt      offsetfrom,offsetto,goffsetfrom,goffsetto;
  DMNetworkComponentGenericDataType *arr;
  PetscInt      row[2],col[8];
  PetscScalar   values[8];

  PetscFunctionBegin;
  ierr = MatZeroEntries(J);CHKERRQ(ierr);

  ierr = SNESGetDM(snes,&networkdm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

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
      ierr = DMNetworkGetComponentTypeOffset(networkdm,v,j,&key,&offsetd);CHKERRQ(ierr);
      if (key == 1) {
        PetscInt       nconnedges;
	const PetscInt *connedges;

	bus = (VERTEXDATA)(arr+offsetd);
	if (!ghostvtex) {
	  /* Handle reference bus constrained dofs */
	  if (bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) {
	    row[0] = goffset; row[1] = goffset+1;
	    col[0] = goffset; col[1] = goffset+1; col[2] = goffset; col[3] = goffset+1;
	    values[0] = 1.0; values[1] = 0.0; values[2] = 0.0; values[3] = 1.0;
	    ierr = MatSetValues(J,2,row,2,col,values,ADD_VALUES);CHKERRQ(ierr);
	    break;
	  }
	  
	  Vm = xarr[offset+1];
	  
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
	for (i=0; i < nconnedges; i++) {
	  EDGEDATA       branch;
	  VERTEXDATA     busf,bust;
	  PetscInt       offsetfd,offsettd,keyf,keyt;
          PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
          const PetscInt *cone;
          PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;

	  e = connedges[i];
	  ierr = DMNetworkGetComponentTypeOffset(networkdm,e,0,&key,&offsetd);CHKERRQ(ierr);
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
	  vfrom = cone[0];
	  vto   = cone[1];

	  ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
	  ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);
	  ierr = DMNetworkGetVariableGlobalOffset(networkdm,vfrom,&goffsetfrom);CHKERRQ(ierr);
	  ierr = DMNetworkGetVariableGlobalOffset(networkdm,vto,&goffsetto);CHKERRQ(ierr);

	  if (goffsetto < 0) goffsetto = -goffsetto - 1;

	  thetaf = xarr[offsetfrom];
	  Vmf     = xarr[offsetfrom+1];
	  thetat = xarr[offsetto];
	  Vmt     = xarr[offsetto+1];
	  thetaft = thetaf - thetat;
	  thetatf = thetat - thetaf;

	  ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&keyf,&offsetfd);CHKERRQ(ierr);
	  ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&keyt,&offsettd);CHKERRQ(ierr);
	  busf = (VERTEXDATA)(arr+offsetfd);
	  bust = (VERTEXDATA)(arr+offsettd);

	  if (vfrom == v) {
	    if (busf->ide != REF_BUS) {
	      /*    farr[offsetfrom]   += Gff*Vmf*Vmf + Vmf*Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));  */
	      row[0]  = goffsetfrom;
	      col[0]  = goffsetfrom; col[1] = goffsetfrom+1; col[2] = goffsetto; col[3] = goffsetto+1;
	      values[0] =  Vmf*Vmt*(Gft*-PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft)); /* df_dthetaf */    
	      values[1] =  2.0*Gff*Vmf + Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft)); /* df_dVmf */
	      values[2] =  Vmf*Vmt*(Gft*PetscSinScalar(thetaft) + Bft*-PetscCosScalar(thetaft)); /* df_dthetat */
	      values[3] =  Vmf*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft)); /* df_dVmt */
	      
	      ierr = MatSetValues(J,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
	    }
	    if (busf->ide != PV_BUS && busf->ide != REF_BUS) {
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
	    if (bust->ide != REF_BUS) {
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
	if (!ghostvtex && bus->ide == PV_BUS) {
	  row[0] = goffset+1; col[0] = goffset+1;
	  values[0]  = 1.0;
	  ierr = MatSetValues(J,1,row,1,col,values,ADD_VALUES);CHKERRQ(ierr);
	}
      }
    }
  }
  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetInitialValues"
PetscErrorCode SetInitialValues(DM networkdm,Vec X,void* appctx)
{
  PetscErrorCode ierr;
  VERTEXDATA     bus;
  GEN            gen;
  PetscInt       v, vStart, vEnd, offset;
  PetscBool      ghostvtex;
  Vec            localX;
  PetscScalar    *xarr;
  PetscInt       key,numComps,j,offsetd;
  DMNetworkComponentGenericDataType *arr;
  
  PetscFunctionBegin;
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
	bus = (VERTEXDATA)(arr+offsetd);
	xarr[offset] = bus->va*PETSC_PI/180.0;
	xarr[offset+1] = bus->vm;
      } else if(key == 2) {
	gen = (GEN)(arr+offsetd);
	if (!gen->status) continue;
	xarr[offset+1] = gen->vs; 
	break;
      }
    }
  }
  ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char ** argv)
{
  PetscErrorCode ierr;
  PetscInt       numEdges=9,numVertices=9;
  int            *edges = NULL;
  PetscInt       i;  
  DM             networkdm;
  PetscInt       componentkey[4];
  UserCtx        User;
  PetscLogStage  stage1;
  PetscMPIInt    size;
  PetscInt       eStart, eEnd, vStart, vEnd,j,neqs_gen,neqs_net,neqs_pgrid;//
 // PetscInt       genj,loadj;
  Vec            X,F,F_alg,Xdot,V0;//V0: initial real and imaginary voltage of all buses
  Mat            J,A;
 //SNES           snes;
  TS                ts;
  //SNES           snes_alg;
  PetscViewer    Xview,Ybusview,viewer;
  PetscInt       i,idx,*idx2,row_loc,col_loc;
  PetscScalar    *x,*mat,val,*amat;
  //Vec            vatol;
  Mat         Ybus; /* Network admittance matrix */
      
  neqs_gen   = 9*ngen; /* # eqs. for generator subsystem */
  neqs_net   = 2*nbus; /* # eqs. for network subsystem   */
  neqs_pgrid = neqs_gen + neqs_net;
  
  ierr = PetscInitialize(&argc,&argv,"petscoptions",help);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  


 /* Read initial voltage vector and Ybus */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"X.bin",FILE_MODE_READ,&Xview);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ybus.bin",FILE_MODE_READ,&Ybusview);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&V0);CHKERRQ(ierr);
  ierr = VecSetSizes(V0,PETSC_DECIDE,neqs_net);CHKERRQ(ierr);
  ierr = VecLoad(V0,Xview);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&Ybus);CHKERRQ(ierr);
  ierr = MatSetSizes(Ybus,PETSC_DECIDE,PETSC_DECIDE,neqs_net,neqs_net);CHKERRQ(ierr);
  ierr = MatSetType(Ybus,MATBAIJ);CHKERRQ(ierr);
  /*  ierr = MatSetBlockSize(Ybus,2);CHKERRQ(ierr); */
  ierr = MatLoad(Ybus,Ybusview);CHKERRQ(ierr);
  //MatView(Ybus, PETSC_VIEWER_STDOUT_SELF);
  VecView(V0, PETSC_VIEWER_STDOUT_SELF);
 
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  
   if (!rank) {
    /*    READ DATA */
    /* Only rank 0 reads the data */
    ierr = PetscOptionsGetString(NULL,NULL,"-pfdata",pfdata_file,PETSC_MAX_PATH_LEN-1,NULL);CHKERRQ(ierr);
    ierr = PetscNew(&pfdata);CHKERRQ(ierr);
    ierr = PFReadMatPowerData(pfdata,pfdata_file);CHKERRQ(ierr);  // The strtucture pfdata reads the input of the case file
    User.Sbase = pfdata->sbase;

    numEdges = pfdata->nbranch;
    numVertices = pfdata->nbus;

    ierr = PetscMalloc(2*numEdges*sizeof(int),&edges);CHKERRQ(ierr); // Q: Why allocate '2'*numEdges into the 'edges' structure?
    ierr = GetListofEdges(pfdata->nbranch,pfdata->branch,edges);CHKERRQ(ierr);
  }
  
  
}
// #undef __FUNCT__
// #define __FUNCT__ "main"
// int main(int argc,char ** argv)
// {
  // PetscErrorCode ierr;
  // char           pfdata_file[PETSC_MAX_PATH_LEN]="datafiles/case13659pegase.m";
  // PFDATA         *pfdata;
  // PetscInt       numEdges=0,numVertices=0;
  // int            *edges = NULL;
  // PetscInt       i;  
  // DM             networkdm;
  // PetscInt       componentkey[4];
  // UserCtx        User;
  // PetscLogStage  stage1,stage2;
  // PetscMPIInt    size;
  // PetscInt       eStart, eEnd, vStart, vEnd,j;
  // PetscInt       genj,loadj;
  // Vec            X,F;
  // Mat            J;
 //SNES           snes;
  // TS                ts;

  // PetscInitialize(&argc,&argv,"pfoptions",help);
  // ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  // /* Create an empty network object */
  // ierr = DMNetworkCreate(PETSC_COMM_WORLD,&networkdm);CHKERRQ(ierr);
  // /* Register the components in the network */
  // ierr = DMNetworkRegisterComponent(networkdm,"branchstruct",sizeof(struct _p_EDGEDATA),&componentkey[0]);CHKERRQ(ierr);
  // ierr = DMNetworkRegisterComponent(networkdm,"busstruct",sizeof(struct _p_VERTEXDATA),&componentkey[1]);CHKERRQ(ierr);
  // ierr = DMNetworkRegisterComponent(networkdm,"genstruct",sizeof(struct _p_GEN),&componentkey[2]);CHKERRQ(ierr);
  // ierr = DMNetworkRegisterComponent(networkdm,"loadstruct",sizeof(struct _p_LOAD),&componentkey[3]);CHKERRQ(ierr);

  // ierr = PetscLogStageRegister("Read Data",&stage1);CHKERRQ(ierr);
  // PetscLogStagePush(stage1);
  // /* READ THE DATA */
  // if (!rank) {
    // /*    READ DATA */
    // /* Only rank 0 reads the data */
    // ierr = PetscOptionsGetString(NULL,NULL,"-pfdata",pfdata_file,PETSC_MAX_PATH_LEN-1,NULL);CHKERRQ(ierr);
    // ierr = PetscNew(&pfdata);CHKERRQ(ierr);
    // ierr = PFReadMatPowerData(pfdata,pfdata_file);CHKERRQ(ierr);  // The strtucture pfdata reads the input of the case file
    // User.Sbase = pfdata->sbase;

    // numEdges = pfdata->nbranch;
    // numVertices = pfdata->nbus;

    // ierr = PetscMalloc(2*numEdges*sizeof(int),&edges);CHKERRQ(ierr); // Q: Why allocate '2'*numEdges into the 'edges' structure?
    // ierr = GetListofEdges(pfdata->nbranch,pfdata->branch,edges);CHKERRQ(ierr);
  // }
  // PetscLogStagePop();
  // ierr = MPI_Barrier(PETSC_COMM_WORLD);CHKERRQ(ierr);
  // ierr = PetscLogStageRegister("Create network",&stage2);CHKERRQ(ierr);
  // PetscLogStagePush(stage2);
  // /* Set number of nodes/edges */
  // ierr = DMNetworkSetSizes(networkdm,numVertices,numEdges,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  // /* Add edge connectivity */
  // ierr = DMNetworkSetEdgeList(networkdm,edges);CHKERRQ(ierr);
  // /* Set up the network layout */
  // ierr = DMNetworkLayoutSetUp(networkdm);CHKERRQ(ierr);

  // if (!rank) {
    // ierr = PetscFree(edges);CHKERRQ(ierr);
  // }
  // /* Add network components */
  
  // genj=0; loadj=0;
  // ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  // for (i = eStart; i < eEnd; i++) {
    // ierr = DMNetworkAddComponent(networkdm,i,componentkey[0],&pfdata->branch[i-eStart]);CHKERRQ(ierr);
  // }
  // ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  // for (i = vStart; i < vEnd; i++) {
    // ierr = DMNetworkAddComponent(networkdm,i,componentkey[1],&pfdata->bus[i-vStart]);CHKERRQ(ierr);
    // if (pfdata->bus[i-vStart].ngen) {
      // for (j = 0; j < pfdata->bus[i-vStart].ngen; j++) {
	// ierr = DMNetworkAddComponent(networkdm,i,componentkey[2],&pfdata->gen[genj++]);CHKERRQ(ierr);
      // }
    // }
    // if (pfdata->bus[i-vStart].nload) {
      // for (j=0; j < pfdata->bus[i-vStart].nload; j++) {
	// ierr = DMNetworkAddComponent(networkdm,i,componentkey[3],&pfdata->load[loadj++]);CHKERRQ(ierr);
      // }
    // }
    // /* Add number of variables */
    // ierr = DMNetworkAddNumVariables(networkdm,i,2);CHKERRQ(ierr);
  // }
  // /* Set up DM for use */
  // ierr = DMSetUp(networkdm);CHKERRQ(ierr);

  // if (!rank) {
    // ierr = PetscFree(pfdata->bus);CHKERRQ(ierr);
    // ierr = PetscFree(pfdata->gen);CHKERRQ(ierr);
    // ierr = PetscFree(pfdata->branch);CHKERRQ(ierr);
    // ierr = PetscFree(pfdata->load);CHKERRQ(ierr);
    // ierr = PetscFree(pfdata);CHKERRQ(ierr);
  // }


  // ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  // if (size > 1) {
    // DM distnetworkdm;
    // /* Network partitioning and distribution of data */
    // ierr = DMNetworkDistribute(networkdm,0,&distnetworkdm);CHKERRQ(ierr);
    // ierr = DMDestroy(&networkdm);CHKERRQ(ierr);
    // networkdm = distnetworkdm;
  // }

  // PetscLogStagePop();
  // ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  // ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  
// #if 0
  // PetscInt numComponents;
  // EDGEDATA edge;
  // PetscInt offset,key,kk;
  // DMNetworkComponentGenericDataType *arr;
  // VERTEXDATA     bus;
  // GEN            gen;
  // LOAD           load;
   
  // for (i = eStart; i < eEnd; i++) {
    // ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);
    // ierr = DMNetworkGetComponentTypeOffset(networkdm,i,0,&key,&offset);CHKERRQ(ierr);
    // edge = (EDGEDATA)(arr+offset);
    // ierr = DMNetworkGetNumComponents(networkdm,i,&numComponents);CHKERRQ(ierr);
    // ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d ncomps = %d Line %d ---- %d\n",rank,numComponents,edge->internal_i,edge->internal_j);CHKERRQ(ierr);
  // }    

  // for (i = vStart; i < vEnd; i++) {
    // ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);
    // ierr = DMNetworkGetNumComponents(networkdm,i,&numComponents);CHKERRQ(ierr);
    // for (kk=0; kk < numComponents; kk++) {
      // ierr = DMNetworkGetComponentTypeOffset(networkdm,i,kk,&key,&offset);CHKERRQ(ierr);
      // if (key == 1) {
	// bus = (VERTEXDATA)(arr+offset);
	// ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d ncomps = %d Bus %d\n",rank,numComponents,bus->internal_i);CHKERRQ(ierr);
      // } else if (key == 2) {
	// gen = (GEN)(arr+offset);
	// ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d Gen pg = %f qg = %f\n",rank,gen->pg,gen->qg);CHKERRQ(ierr);
      // } else if (key == 3) {
	// load = (LOAD)(arr+offset);
	// ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d Load pl = %f ql = %f\n",rank,load->pl,load->ql);CHKERRQ(ierr);
      // }
    // }
  // }  
// #endif  
  // /* Broadcast Sbase to all processors */
  // ierr = MPI_Bcast(&User.Sbase,1,MPIU_SCALAR,0,PETSC_COMM_WORLD);CHKERRQ(ierr);

  // ierr = DMCreateGlobalVector(networkdm,&X);CHKERRQ(ierr);
  // ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

  // ierr = DMCreateMatrix(networkdm,&J);CHKERRQ(ierr);
  // ierr = MatSetOption(J,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);

  // ierr = SetInitialValues(networkdm,X,&User);CHKERRQ(ierr);

  // /* HOOK UP SOLVER */
  // ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  // ierr = SNESSetDM(snes,networkdm);CHKERRQ(ierr);
  // ierr = SNESSetFunction(snes,F,FormFunction,&User);CHKERRQ(ierr);
  // ierr = SNESSetJacobian(snes,J,J,FormJacobian,&User);CHKERRQ(ierr);
  // ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  // ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);
  // VecView(X,PETSC_VIEWER_STDOUT_WORLD);
// VecView(F,PETSC_VIEWER_STDOUT_WORLD);
  // MatView(J,PETSC_VIEWER_STDOUT_WORLD);
  // ierr = VecDestroy(&X);CHKERRQ(ierr);
  // ierr = VecDestroy(&F);CHKERRQ(ierr);
  // ierr = MatDestroy(&J);CHKERRQ(ierr);

  // ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  // ierr = DMDestroy(&networkdm);CHKERRQ(ierr);

  // PetscFinalize();
  // return 0;
// }
