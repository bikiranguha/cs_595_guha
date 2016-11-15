// #undef __FUNCT__
// #define __FUNCT__ "SetInitialValues"
// PetscErrorCode SetInitialValues(DM networkdm,Vec X,void* appctx)
// {
  // PetscErrorCode ierr;
  // VERTEXDATA     bus;
  // GEN            gen;
  // PetscInt       v, vStart, vEnd, offset;
  // PetscBool      ghostvtex;
  // Vec            localX;
  // PetscScalar    *xarr;
  // PetscInt       key,numComps,j,offsetd;
  // DMNetworkComponentGenericDataType *arr;
  
  // PetscFunctionBegin;
  // ierr = DMNetworkGetVertexRange(networkdm,&vStart, &vEnd);CHKERRQ(ierr);

  // ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);

  // ierr = VecSet(X,0.0);CHKERRQ(ierr);
  // ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  // ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  // ierr = VecGetArray(localX,&xarr);CHKERRQ(ierr);
  // ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);
  // for (v = vStart; v < vEnd; v++) {
    // ierr = DMNetworkIsGhostVertex(networkdm,v,&ghostvtex);CHKERRQ(ierr);
    // if (ghostvtex) continue;
   
    // ierr = DMNetworkGetVariableOffset(networkdm,v,&offset);CHKERRQ(ierr);
    // ierr = DMNetworkGetNumComponents(networkdm,v,&numComps);CHKERRQ(ierr);
    // for (j=0; j < numComps; j++) {
      // ierr = DMNetworkGetComponentTypeOffset(networkdm,v,j,&key,&offsetd);CHKERRQ(ierr);
      // if (key == 1) {
	// bus = (VERTEXDATA)(arr+offsetd);
	// xarr[offset] = bus->va*PETSC_PI/180.0;
	// xarr[offset+1] = bus->vm;
      // } else if(key == 2) {
	// gen = (GEN)(arr+offsetd);
	// if (!gen->status) continue;
	// xarr[offset+1] = gen->vs; 
	// break;
      // }
    // }
  // }
  // ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);
  // ierr = DMLocalToGlobalBegin(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  // ierr = DMLocalToGlobalEnd(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  // ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);
  // PetscFunctionReturn(0);
// }

#undef __FUNCT__
#define __FUNCT__ "SetInitialGuess"
PetscErrorCode SetInitialGuess(DM networkdm, Vec X, Vec V0)
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
  PetscInt       i,idx=0;
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
	 bus = (bus)(arr+offsetd);
	 xarr[offset] = bus.vr;
	 xarr[offset+1] = bus.vi;//->
	 Vr = bus.vr;
	 Vi= bus.vi;
         } else if(key == 2) {
	gen = (GEN)(arr+offsetd);
    Vm  = PetscSqrtScalar(Vr*Vr + Vi*Vi); Vm2 = Vm*Vm;
    IGr = (Vr*gen.PG + Vi*gen.QG)/Vm2; // Real part of gen current
    IGi = (Vi*gen.PG - Vr*gen.QG)/Vm2; // Imaginary part of gen current

    
	/* Machine angle */
	delta = atan2(Vi+gen.Xq*IGr,Vr-gen.Xq*IGi);  // Link for atan2: https://en.wikipedia.org/wiki/Atan2 

    theta = PETSC_PI/2.0 - delta;

    Id = IGr*PetscCosScalar(theta) - IGi*PetscSinScalar(theta); /* d-axis stator current */
    Iq = IGr*PetscSinScalar(theta) + IGi*PetscCosScalar(theta); /* q-axis stator current */

    Vd = Vr*PetscCosScalar(theta) - Vi*PetscSinScalar(theta);
    Vq = Vr*PetscSinScalar(theta) + Vi*PetscCosScalar(theta);

    Edp = Vd + gen.Rs*Id - gen.Xqp*Iq; /* d-axis transient EMF */
    Eqp = Vq + gen.Rs*Iq + gen.Xdp*Id; /* q-axis transient EMF */

    gen.TM = gen.PG;
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
    Efd = Eqp + (gen.Xd - gen.Xdp)*Id;
    SE  = gen.k1*PetscExpScalar(gen.k2*Efd);
    VR  =  gen.KE*Efd + SE;
    RF  =  gen.KF*Efd/gen.TF;

    xarr[idx]   = Efd;
    xarr[idx+1] = RF;
    xarr[idx+2] = VR;

    gen.Vref = Vm + (VR/gen.KA);
	
	
      }
     }
   }
  ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
 }
  
  
  
  
  
  
  
  
  
  
  
  
  for (i=0; i < ngen; i++) {
    Vr  = xnet[2*gbus[i]]; /* Real part of generator terminal voltage */
    Vi  = xnet[2*gbus[i]+1]; /* Imaginary part of the generator terminal voltage */
    Vm  = PetscSqrtScalar(Vr*Vr + Vi*Vi); Vm2 = Vm*Vm;
    IGr = (Vr*PG[i] + Vi*QG[i])/Vm2; // Real part of gen current
    IGi = (Vi*PG[i] - Vr*QG[i])/Vm2; // Imaginary part of gen current

    
	/* Machine angle */
	delta = atan2(Vi+Xq[i]*IGr,Vr-Xq[i]*IGi);  // Link for atan2: https://en.wikipedia.org/wiki/Atan2 

    theta = PETSC_PI/2.0 - delta;

    Id = IGr*PetscCosScalar(theta) - IGi*PetscSinScalar(theta); /* d-axis stator current */
    Iq = IGr*PetscSinScalar(theta) + IGi*PetscCosScalar(theta); /* q-axis stator current */

    Vd = Vr*PetscCosScalar(theta) - Vi*PetscSinScalar(theta);
    Vq = Vr*PetscSinScalar(theta) + Vi*PetscCosScalar(theta);

    Edp = Vd + Rs[i]*Id - Xqp[i]*Iq; /* d-axis transient EMF */
    Eqp = Vq + Rs[i]*Iq + Xdp[i]*Id; /* q-axis transient EMF */

    TM[i] = PG[i];

    /* The generator variables are ordered as [Eqp,Edp,delta,w,Id,Iq] */
    xgen[idx]   = Eqp;
    xgen[idx+1] = Edp;
    xgen[idx+2] = delta;
    xgen[idx+3] = w_s;

    idx = idx + 4;

    xgen[idx]   = Id;
    xgen[idx+1] = Iq;

    idx = idx + 2;

    /* Exciter */
    Efd = Eqp + (Xd[i] - Xdp[i])*Id;
    SE  = k1[i]*PetscExpScalar(k2[i]*Efd);
    VR  =  KE[i]*Efd + SE;
    RF  =  KF[i]*Efd/TF[i];

    xgen[idx]   = Efd;
    xgen[idx+1] = RF;
    xgen[idx+2] = VR;

    Vref[i] = Vm + (VR/KA[i]);

    idx = idx + 3;
  }
  
  
  
  
  
  
  
  

//  ierr = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);

  /* Network subsystem initialization */
  ierr = VecCopy(user->V0,Xnet);CHKERRQ(ierr);

  /* Generator subsystem initialization */
  ierr = VecGetArray(Xgen,&xgen);CHKERRQ(ierr);
  ierr = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);

  

  ierr = VecRestoreArray(Xgen,&xgen);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xnet,&xnet);CHKERRQ(ierr);

  /* ierr = VecView(Xgen,0);CHKERRQ(ierr); */
  ierr = DMCompositeGather(user->dmpgrid,X,INSERT_VALUES,Xgen,Xnet);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



// #undef __FUNCT__
// #define __FUNCT__ "SetInitialGuess"
// PetscErrorCode SetInitialGuess(Vec X,Userctx *user)
// {
  // PetscErrorCode ierr;
  // Vec            Xgen,Xnet;
  // PetscScalar    *xgen,*xnet;
  // PetscInt       i,idx=0;
  // PetscScalar    Vr,Vi,IGr,IGi,Vm,Vm2;
  // PetscScalar    Eqp,Edp,delta;
  // PetscScalar    Efd,RF,VR; /* Exciter variables */
  // PetscScalar    Id,Iq;  /* Generator dq axis currents */
  // PetscScalar    theta,Vd,Vq,SE;

  // PetscFunctionBegin;
  // M[0] = 2*H[0]/w_s; M[1] = 2*H[1]/w_s; M[2] = 2*H[2]/w_s;
  // D[0] = 0.1*M[0]; D[1] = 0.1*M[1]; D[2] = 0.1*M[2];

  // ierr = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);

  // /* Network subsystem initialization */
  // ierr = VecCopy(user->V0,Xnet);CHKERRQ(ierr);

  // /* Generator subsystem initialization */
  // ierr = VecGetArray(Xgen,&xgen);CHKERRQ(ierr);
  // ierr = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);

  // for (i=0; i < ngen; i++) {
    // Vr  = xnet[2*gbus[i]]; /* Real part of generator terminal voltage */
    // Vi  = xnet[2*gbus[i]+1]; /* Imaginary part of the generator terminal voltage */
    // Vm  = PetscSqrtScalar(Vr*Vr + Vi*Vi); Vm2 = Vm*Vm;
    // IGr = (Vr*PG[i] + Vi*QG[i])/Vm2; // Real part of gen current
    // IGi = (Vi*PG[i] - Vr*QG[i])/Vm2; // Imaginary part of gen current

    
	// /* Machine angle */
	// delta = atan2(Vi+Xq[i]*IGr,Vr-Xq[i]*IGi);  // Link for atan2: https://en.wikipedia.org/wiki/Atan2 

    // theta = PETSC_PI/2.0 - delta;

    // Id = IGr*PetscCosScalar(theta) - IGi*PetscSinScalar(theta); /* d-axis stator current */
    // Iq = IGr*PetscSinScalar(theta) + IGi*PetscCosScalar(theta); /* q-axis stator current */

    // Vd = Vr*PetscCosScalar(theta) - Vi*PetscSinScalar(theta);
    // Vq = Vr*PetscSinScalar(theta) + Vi*PetscCosScalar(theta);

    // Edp = Vd + Rs[i]*Id - Xqp[i]*Iq; /* d-axis transient EMF */
    // Eqp = Vq + Rs[i]*Iq + Xdp[i]*Id; /* q-axis transient EMF */

    // TM[i] = PG[i];

    // /* The generator variables are ordered as [Eqp,Edp,delta,w,Id,Iq] */
    // xgen[idx]   = Eqp;
    // xgen[idx+1] = Edp;
    // xgen[idx+2] = delta;
    // xgen[idx+3] = w_s;

    // idx = idx + 4;

    // xgen[idx]   = Id;
    // xgen[idx+1] = Iq;

    // idx = idx + 2;

    // /* Exciter */
    // Efd = Eqp + (Xd[i] - Xdp[i])*Id;
    // SE  = k1[i]*PetscExpScalar(k2[i]*Efd);
    // VR  =  KE[i]*Efd + SE;
    // RF  =  KF[i]*Efd/TF[i];

    // xgen[idx]   = Efd;
    // xgen[idx+1] = RF;
    // xgen[idx+2] = VR;

    // Vref[i] = Vm + (VR/KA[i]);

    // idx = idx + 3;
  // }

  // ierr = VecRestoreArray(Xgen,&xgen);CHKERRQ(ierr);
  // ierr = VecRestoreArray(Xnet,&xnet);CHKERRQ(ierr);

  // /* ierr = VecView(Xgen,0);CHKERRQ(ierr); */
  // ierr = DMCompositeGather(user->dmpgrid,X,INSERT_VALUES,Xgen,Xnet);CHKERRQ(ierr);
  // ierr = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  // PetscFunctionReturn(0);
// }