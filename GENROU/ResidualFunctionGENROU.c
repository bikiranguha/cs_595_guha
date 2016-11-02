/* --------------------------------------------------
   Routine for GENROU Model Function Evaluation
   -------------------------------------------------- */

#include "../../include/TS3ph.h"
#include "../../include/blocks.h"
#undef __FUNCT__
#define __FUNCT__ "ResidualFunctionGENROU"
PetscErrorCode ResidualFunctionGENROU(TS3phDynCtx* Gen, TS3phNetCtx* Net, PetscInt i)
{

  PetscErrorCode ierr;
  
  PetscScalar    *xgen_arr, *fgen_arr, *v_arr, *ugen_arr, *xdot_arr;
  PetscScalar    *vgen, theta_a, theta_b, theta_c, MBoSB;
  
  PetscScalar	T_d0p,T_d0dp,T_q0p,T_q0dp,H,D,X_ddp,X_qdp;
  PetscScalar	Gaind1,Gaind2,Gaind3,Gaind4,Gaind5;
  PetscScalar	Gainq1,Gainq2,Gainq3,Gainq4,Gainq5,Gainq6;
  PetscScalar   Esat,Scoef;
  
  PetscInt      Sat_flag;
  
  PetscInt     Eqp_idx,Edp_idx,psi1d_idx,psi2q_idx,n_idx,delta_idx;
  PetscInt     EFD_idx,PM_idx,IFD_idx,pelec_idx,Id_idx,Iq_idx; // xgen_arr[GENROU->PM_idx[idx_gen]] in pu on MBase
  
  PetscScalar	Eqp,Edp,psi1d,psi2q,n,delta;
  PetscScalar	EFD,PM,IFD,pelec; // xgen_arr[GENROU->PM_idx[idx_gen]] in pu on MBase
  PetscScalar	TM,DIFD;
  PetscScalar	Id,Iq,Vd,Vq;
  PetscScalar   Psiddp,Psiqdp,Psidp;

  PetscInt     Eqp_dot,Edp_dot,psi1d_dot,psi2q_dot,n_dot,delta_dot;
  PetscInt     EFD_dot,PM_dot,IFD_dot,pelec_dot;
  
  Vec            Xgen = Gen->xdyn;
  Vec            Fgen = Gen->f_dyn;
  Vec            V = Net->V;
  Vec            Xdot = Gen->xdotdyn;
  Vec            Ugen = Gen->udyn;
  
  GENDATA        *Gdata = Gen->GenData;
  MAC_INFO_data  *MAC_INFO = Gdata->MAC_INFO;
  GENROU_MODEL   *GENROU = Gdata->GEN_MODEL->GENROU;
  GMAP_data      *GMAP = Gdata->GMAP;
  
  PetscInt       idx_gen, bus_index, seq_number;

  PetscFunctionBegin;
  
  /*External vectors*/
  ierr = VecGetArray(Xgen,&xgen_arr);CHKERRQ(ierr);
  ierr = VecGetArray(Fgen,&fgen_arr);CHKERRQ(ierr);
  ierr = VecGetArray(Ugen,&ugen_arr);CHKERRQ(ierr);
  ierr = VecGetArray(Xdot,&xdot_arr);CHKERRQ(ierr);
  ierr = VecGetArray(V,&v_arr);CHKERRQ(ierr);
  
  idx_gen = MAC_INFO->gen_idx[i];

  /*State Variables*/
  Eqp_idx = GENROU->Eqp_idx[idx_gen];
  Edp_idx = GENROU->Edp_idx[idx_gen];
  psi1d_idx = GENROU->psi1d_idx[idx_gen];
  psi2q_idx = GENROU->psi2q_idx[idx_gen];
  n_idx = GENROU->n_idx[idx_gen];
  delta_idx = GENROU->delta_idx[idx_gen];
  
  EFD_idx=GENROU->EFD_idx[idx_gen];
  PM_idx=GENROU->PM_idx[idx_gen];
  
  pelec_idx = GENROU->pelec_idx[idx_gen];
  Id_idx = GENROU->Id_idx[idx_gen];
  Iq_idx = GENROU->Iq_idx[idx_gen];
  IFD_idx = GENROU->IFD_idx[idx_gen];
  
  Eqp = xgen_arr[Eqp_idx];
  Edp = xgen_arr[Edp_idx];
  psi1d = xgen_arr[psi1d_idx];
  psi2q = xgen_arr[psi2q_idx];
  n = xgen_arr[n_idx];
  delta = xgen_arr[delta_idx];
  
  Id = xgen_arr[Id_idx];
  Iq = xgen_arr[Iq_idx];
  
  Eqp_dot = xdot_arr[Eqp_idx];
  Edp_dot = xdot_arr[Edp_idx];
  psi1d_dot = xdot_arr[psi1d_idx];
  psi2q_dot = xdot_arr[psi2q_idx];
  n_dot = xdot_arr[n_idx];
  delta_dot = xdot_arr[delta_idx];

  /*State Parameters*/
  T_d0p = GENROU->T_d0p[idx_gen];
  T_d0dp = GENROU->T_d0dp[idx_gen];
  T_q0p = GENROU->T_q0p[idx_gen];
  T_q0dp = GENROU->T_q0dp[idx_gen];
  H = GENROU->H[idx_gen];
  D = GENROU->D[idx_gen];
  X_ddp = GENROU->X_ddp[idx_gen];
  X_qdp = GENROU->X_ddp[idx_gen]; 
  MBoSB = GENROU->MBoSB[idx_gen]; 
  Gaind1=GENROU->Gaind1[idx_gen];
  Gaind2=GENROU->Gaind2[idx_gen];
  Gaind3=GENROU->Gaind3[idx_gen];
  Gaind4=GENROU->Gaind4[idx_gen];
  Gaind5=GENROU->Gaind5[idx_gen];
  Gainq1=GENROU->Gainq1[idx_gen];
  Gainq2=GENROU->Gainq2[idx_gen];
  Gainq3=GENROU->Gainq3[idx_gen];
  Gainq4=GENROU->Gainq4[idx_gen];
  Gainq5=GENROU->Gainq5[idx_gen];
  Gainq6=GENROU->Gainq6[idx_gen];
  
  Sat_flag=GENROU->Sat_flag[idx_gen];
  Esat=GENROU->Esat[idx_gen];
  Scoef=GENROU->Scoef[idx_gen];

  /*GENROU EQUATIONS*/
  bus_index  = MAC_INFO->gen_bus[i];
  seq_number = bus_index;
  
  vgen = v_arr + twonphase*(seq_number);

  theta_a = delta;
  theta_b = delta - twoPI_3;
  theta_c = delta + twoPI_3;
  PetscScalar     sin_theta[3] = {0.0}, cos_theta[3] = {0.0};
  sin_theta[0] = sin(theta_a); sin_theta[1] = sin(theta_b); sin_theta[2] = sin(theta_c);
  cos_theta[0] = cos(theta_a); cos_theta[1] = cos(theta_b); cos_theta[2] = cos(theta_c);

  Vd = (sin_theta[0]*vgen[0] + sin_theta[1]*vgen[1] + sin_theta[2]*vgen[2]
        - cos_theta[0]*vgen[3] - cos_theta[1]*vgen[4] - cos_theta[2]*vgen[5])
    / 3.0;
  Vq = (cos_theta[0]*vgen[0] + cos_theta[1]*vgen[1] + cos_theta[2]*vgen[2]
        + sin_theta[0]*vgen[3] + sin_theta[1]*vgen[4] + sin_theta[2]*vgen[5])
    /3.0;
	
  Psiddp = Gaind3*Eqp + Gaind2*psi1d;
  Psiqdp = -Gainq3*Edp + Gainq2*psi2q;
  if (Sat_flag==1) {
	  Psidp=sqrt(Psiddp*Psiddp+Psiqdp*Psiqdp);
  }
 
  if(MAC_INFO->tgov_model[i] == 0){
    PM = ugen_arr[2*i + 1]*MBoSB; // convert from MBASE to SBASE
  } else {
    PM = xgen_arr[PM_idx]*MBoSB; // convert from MBASE to SBASE
  }
  
  TM = (PM - D*n)/(1+n);

  if(MAC_INFO->exc_model[i] == 0){
    EFD = ugen_arr[2*i];
  } else {
    EFD = xgen_arr[EFD_idx];
  }

  /*Function mismatch construction*/
  
  /* psi1d */
  PetscScalar psi1din;
  psi1din=Eqp - Id*Gaind1 - psi1d;
  /* psi2q */
  PetscScalar psi2qin;
  psi2qin=-Edp - Iq*Gainq1 - psi2q;
  /* DIFD */
  if (Sat_flag==1) {
	  DIFD=Scoef*pow(Psidp-Esat,2);
  }
  /* IFD */
  PetscScalar IFDin;
  IFDin=Eqp+Gaind5*(Id+Gaind4*psi1din);
  if (Sat_flag==1) {
	  IFDin+=(Psiddp/Psidp)*DIFD;
  }
  if(IFD_idx != -1){
    IFD = xgen_arr[IFD_idx];
    fgen_arr[IFD_idx] = IFDin-IFD;
  }
  /* Eqp */
  PetscScalar Eqpin;
  Eqpin=EFD-IFDin;
  /* Edp */
  PetscScalar Edpin;
  Edpin=(Iq+psi2qin*Gainq4)*Gainq5-Edp;
  if (Sat_flag==1) {
	  Edpin+=(Psiqdp/Psidp)*Gainq6*DIFD;
  }
  /* pelec */
  PetscScalar pelecin;
  pelecin=Psiddp*Iq - Psiqdp*Id;
  pelec = xgen_arr[pelec_idx];
  fgen_arr[pelec_idx] = pelecin/MBoSB - pelec;
  /* Id & Iq */
  fgen_arr[Id_idx]=(Psiddp - Vq)/X_ddp-Id;
  fgen_arr[Iq_idx]=(Psiqdp + Vd)/X_qdp-Iq;
  /* n (SPEED DERIVATION) */
  PetscScalar nin;
  nin=(TM - pelec*MBoSB)/(2*H);
  /* delta */  
  PetscScalar deltain;
  deltain = 2*PI*60*n;

  if (Gen->net_mode!=2) {
	   fgen_arr[psi1d_idx]=psi1din/T_d0dp;
	   fgen_arr[psi2q_idx]=psi2qin/T_q0dp;
	   fgen_arr[Eqp_idx]=Eqpin/T_d0p;
	   fgen_arr[Edp_idx]=Edpin/T_q0p;
	   fgen_arr[n_idx]=nin;
	   fgen_arr[delta_idx]=deltain;
	   if (Gen->net_mode==0) {
		   fgen_arr[psi1d_idx]-=xdot_arr[psi1d_idx];
		   fgen_arr[psi2q_idx]-=xdot_arr[psi2q_idx];
		   fgen_arr[Eqp_idx]-=xdot_arr[Eqp_idx];
		   fgen_arr[Edp_idx]-=xdot_arr[Edp_idx];
		   fgen_arr[n_idx]-=xdot_arr[n_idx];
		   fgen_arr[delta_idx]-=xdot_arr[delta_idx];
	   }
  } else {
	  fgen_arr[psi1d_idx]=0.0;
	  fgen_arr[psi2q_idx]=0.0;
	  fgen_arr[Eqp_idx]=0.0;
	  fgen_arr[Edp_idx]=0.0;
	  fgen_arr[n_idx]=0.0;
	  fgen_arr[delta_idx]=0.0;
  }
  
  //VecView(Xgen,PETSC_VIEWER_DEFAULT);
  //VecView(Ugen,PETSC_VIEWER_DEFAULT);
  //VecView(Xdot,PETSC_VIEWER_DEFAULT);
  //VecView(Fgen,PETSC_VIEWER_DEFAULT);

  ierr = VecRestoreArray(Xgen,&xgen_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(Fgen,&fgen_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(Ugen,&ugen_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xdot,&xdot_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(V,&v_arr);CHKERRQ(ierr);

PetscReal norm, maxval, norm2, maxval2;
PetscInt N,maxind,N2,maxind2;
if(Gen->net_mode == 0 && i==1262){
/*  VecNorm(Fgen,NORM_2,&norm);
 * VecGetSize(Fgen,&N);
 * VecMax(Fgen,&maxind,&maxval);
 * printf("Fgen_size=%d \tNORM=%G \tmaxind=%d \tmaxval=%G\n",N,norm,maxind,maxval);
 */

  /*VecNorm(Xdot,NORM_2,&norm2);
  VecGetSize(Xdot,&N2);
  VecMax(Xdot,&maxind2,&maxval2);
  printf("Xdot_size=%d \tNORM=%G \t\tmaxind=%d \t\tmaxval=%G\n",N2,norm2,maxind2,maxval2);*/
}

  PetscFunctionReturn(0);
  
}
