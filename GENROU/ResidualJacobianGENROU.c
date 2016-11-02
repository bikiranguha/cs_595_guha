#include "../../include/TS3ph.h"
#include "../../include/jacobian.h"

/*
 * Module: GENROU Residual Jacobian 
 * --------------------
 * 
 * Description: builds the residual jacobian equations for the GENROU model.
 *
 */

#undef __FUNCT__
#define __FUNCT__ "ResidualJacobianGENROU"
PetscErrorCode ResidualJacobianGENROU(TS3phDynCtx* Gen, TS3phNetCtx* Net, Mat J, PetscInt i,PetscInt rstart)
{

  PetscErrorCode ierr;
  
  PetscScalar    *xgen_arr, *ugen_arr, *v_arr, *f_eval, *xdot_arr, *fgen_arr;
  PetscScalar    *vgen, theta_a, theta_b, theta_c, MBoSB;
  
  PetscScalar	T_d0p,T_d0dp,T_q0p,T_q0dp,H,D,X_ddp,X_qdp;
  PetscScalar	Gaind1,Gaind2,Gaind3,Gaind4,Gaind5;
  PetscScalar	Gainq1,Gainq2,Gainq3,Gainq4,Gainq5,Gainq6;
  PetscScalar   Esat,Scoef;
  
  PetscInt      Sat_flag;
  
  PetscInt     Eqp_idx,Edp_idx,psi1d_idx,psi2q_idx,n_idx,delta_idx;
  PetscInt     EFD_idx,PM_idx,IFD_idx,pelec_idx,Id_idx,Iq_idx;
  PetscInt     EFD_init_idx,PM_init_idx;
  
  PetscScalar	Eqp,Edp,psi1d,psi2q,n,delta;
  PetscScalar	EFD,PM,IFD,pelec; // xgen_arr[GENROU->PM_idx[idx_gen]] in pu on MBase
  PetscScalar	TM,DIFD;
  PetscScalar	Id,Iq,Vd,Vq;
  PetscScalar    Psiddp, Psiqdp,Psidp;
  
  PetscInt     Eqp_dot,Edp_dot,psi1d_dot,psi2q_dot,n_dot,delta_dot;
  PetscInt     EFD_dot,PM_dot,IFD_dot,pelec_dot;

  Vec            Xgen = Gen->xdyn;
  Vec            Fgen = Gen->f_dyn;
  Vec            V = Net->V;
  Vec            Xdot = Gen->xdotdyn;
  Vec            Ugen = Gen->udyn;

  GENDATA        *Gdata = Gen->GenData;
  GENROU_MODEL   *GENROU = Gdata->GEN_MODEL->GENROU;
  GMAP_data      *GMAP = Gdata->GMAP;
  MAC_INFO_data  *MAC_INFO = Gdata->MAC_INFO;
  
  PetscInt       idx_gen, bus_index, seq_number, j;
  
  PetscReal      a = Gen->a; 
  
  PetscFunctionBegin;
  
  ierr = VecGetArray(Xgen,&xgen_arr);CHKERRQ(ierr);
  ierr = VecGetArray(Xdot,&xdot_arr);CHKERRQ(ierr);
  ierr = VecGetArray(Ugen,&ugen_arr);CHKERRQ(ierr);
  ierr = VecGetArray(Fgen,&fgen_arr);CHKERRQ(ierr);
  ierr = VecGetArray(V,&v_arr);CHKERRQ(ierr);
  
  idx_gen = MAC_INFO->gen_idx[i];
  
  /*State Variables*/
  Eqp_idx = GENROU->Eqp_idx[idx_gen];
  Edp_idx = GENROU->Edp_idx[idx_gen];
  psi1d_idx = GENROU->psi1d_idx[idx_gen];
  psi2q_idx = GENROU->psi2q_idx[idx_gen];
  n_idx = GENROU->n_idx[idx_gen];
  delta_idx = GENROU->delta_idx[idx_gen];
  
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
  
  pelec = xgen_arr[pelec_idx];
  Id = xgen_arr[Id_idx];
  Iq = xgen_arr[Iq_idx];

  Eqp_idx += rstart;
  Edp_idx += rstart;
  psi1d_idx += rstart;
  psi2q_idx += rstart;
  n_idx += rstart;
  delta_idx += rstart;
  
  pelec_idx += rstart;
  Id_idx += rstart;
  Iq_idx += rstart;  
  IFD_idx += rstart;
       
  /* State Parameters*/
  T_d0p = GENROU->T_d0p[idx_gen];
  T_d0dp = GENROU->T_d0dp[idx_gen];
  T_q0p = GENROU->T_q0p[idx_gen];
  T_q0dp = GENROU->T_q0dp[idx_gen];
  H = GENROU->H[idx_gen];
  D = GENROU->D[idx_gen];
  X_ddp = GENROU->X_ddp[idx_gen];
  X_qdp = X_ddp;
  MBoSB = GENROU->MBoSB[idx_gen];
  Gaind1 = GENROU->Gaind1[idx_gen];
  Gaind2 = GENROU->Gaind2[idx_gen];
  Gaind3 = GENROU->Gaind3[idx_gen];
  Gaind4 = GENROU->Gaind4[idx_gen];
  Gaind5 = GENROU->Gaind5[idx_gen];
  Gainq1 = GENROU->Gainq1[idx_gen];
  Gainq2 = GENROU->Gainq2[idx_gen];
  Gainq3 = GENROU->Gainq3[idx_gen];
  Gainq4 = GENROU->Gainq4[idx_gen];
  Gainq5 = GENROU->Gainq5[idx_gen];
  Gainq6 = GENROU->Gainq6[idx_gen];
  Sat_flag = GENROU->Sat_flag[idx_gen];
  Esat = GENROU->Esat[idx_gen];
  Scoef = GENROU->Scoef[idx_gen];
  
  /*Sat_flag=0;
  if (S1>0 && S2>0 && S2*1.2>S1) {
		Z = sqrt(S1/(S2*1.2));
		Esat = (Z*1.2-1)/(Z-1);
		Scoef = pow((sqrt(S2*1.2)-sqrt(S1))*5, 2);
		//GENROU->S1[idx_gen]=Esat;
		//GENROU->S2[idx_gen]=Scoef;
		Sat_flag=1;
  } else {
	  if (S1>0 && S2<0){
		  printf("NO!!!");
	  }
  }
  
  Gaind1=X_dp-Xl;
  Gaind2=(X_dp-X_ddp)/Gaind1;
  Gaind3=(X_ddp-Xl)/Gaind1;
  Gaind4=Gaind2/Gaind1;
  Gaind5=X_d-X_dp;
  
  Gainq1=X_qp-Xl;
  Gainq2=(X_qp-X_qdp)/Gainq1;
  Gainq3=(X_qdp-Xl)/Gainq1;
  Gainq4=Gainq2/Gainq1;
  Gainq5=X_q-X_qp;
  Gainq6=(X_q-Xl)/(X_d-Xl);*/

  /* Jacobian indexing and filling */
  PetscScalar     val[30];
  PetscInt        col[30];
  PetscInt        row[1];
  
  PetscScalar     sin_theta[3] = {0.0}, cos_theta[3] = {0.0};
  
  PetscInt          nmode = Gen->net_mode;
  
  /*External vectors*/

  /*Auxiliary operations*/
  bus_index  = MAC_INFO->gen_bus[i];
  seq_number = bus_index;
  
  PetscInt idx_node = rstart + Gen->dyn_size + twonphase*seq_number;
  
  vgen = v_arr + twonphase*(seq_number);
  
  signl signalOne;signl signalTwo;signl signalThree;
  
  theta_a = delta;
  theta_b = delta - twoPI_3; 
  theta_c = delta + twoPI_3;
  sin_theta[0] = sin(theta_a); sin_theta[1] = sin(theta_b); sin_theta[2] = sin(theta_c);
  cos_theta[0] = cos(theta_a); cos_theta[1] = cos(theta_b); cos_theta[2] = cos(theta_c);

  Vd = (sin_theta[0]*vgen[0] + sin_theta[1]*vgen[1] + sin_theta[2]*vgen[2] - cos_theta[0]*vgen[3] - cos_theta[1]*vgen[4] - cos_theta[2]*vgen[5])/3.0;
  Vq = (cos_theta[0]*vgen[0] + cos_theta[1]*vgen[1] + cos_theta[2]*vgen[2] + sin_theta[0]*vgen[3] + sin_theta[1]*vgen[4] + sin_theta[2]*vgen[5])/3.0;
  
  Psiddp = Gaind3*Eqp + Gaind2*psi1d;
  Psiqdp = -Gainq3*Edp + Gainq2*psi2q;
  if (Sat_flag==1) {
	  Psidp=sqrt(Psiddp*Psiddp+Psiqdp*Psiqdp);
  }
  
  /* Psiddp */
  signl signalPsiddp;
  signalPsiddp.idx[0]=Eqp_idx;
  signalPsiddp.gain[0]=Gaind3;
  signalPsiddp.idx[1]=psi1d_idx;
  signalPsiddp.gain[1]=Gaind2;  
  signalPsiddp.n=2;
  
  /* Psiqdp */
  signl signalPsiqdp;
  signalPsiqdp.idx[0]=Edp_idx;
  signalPsiqdp.gain[0]=-Gainq3;
  signalPsiqdp.idx[1]=psi2q_idx;
  signalPsiqdp.gain[1]=Gainq2;  
  signalPsiqdp.n=2;
  
  /* Psidp */
  signl signalPsidp;
  if (Sat_flag==1) {
	signalPsidp.idx[0]=Eqp_idx;
	signalPsidp.gain[0]=Gaind3*Psiddp/Psidp;
	signalPsidp.idx[1]=psi1d_idx;
	signalPsidp.gain[1]=Gaind2*Psiddp/Psidp;  
	signalPsidp.idx[2]=Edp_idx;
	signalPsidp.gain[2]=-Gainq3*Psiqdp/Psidp;
	signalPsidp.idx[3]=psi2q_idx;
	signalPsidp.gain[3]=Gainq2*Psiqdp/Psidp;
	signalPsidp.n=4;
  }

  /* Id & Iq */
  signl signalId;
  signalId.idx[0]=Eqp_idx;
  signalId.gain[0]=Gaind3/X_ddp;
  signalId.idx[1]=psi1d_idx;
  signalId.gain[1]=Gaind2/X_ddp; 
  signalId.idx[2]=delta_idx;
  signalId.gain[2]=Vd/X_ddp; 
  signalId.n=3;
  
  signl signalIq;
  signalIq.idx[0]=Edp_idx;
  signalIq.gain[0]=-Gainq3/X_qdp;
  signalIq.idx[1]=psi2q_idx;
  signalIq.gain[1]=Gainq2/X_qdp;
  signalIq.idx[2]=delta_idx;
  signalIq.gain[2]=Vq/X_qdp;
  signalIq.n=3;
  
  if (Gen->net_mode!=1) {
	  
		signalIq.idx[3]=idx_node;
		signalIq.gain[3]=sin_theta[0]/(3.0*X_qdp);		
		signalIq.idx[4]=idx_node+1;
		signalIq.gain[4]=sin_theta[1]/(3.0*X_qdp);		
		signalIq.idx[5]=idx_node+2;
		signalIq.gain[5]=sin_theta[2]/(3.0*X_qdp);		
		signalIq.idx[6]=idx_node+3;
		signalIq.gain[6]=-cos_theta[0]/(3.0*X_qdp);		
		signalIq.idx[7]=idx_node+4;
		signalIq.gain[7]=-cos_theta[1]/(3.0*X_qdp);		
		signalIq.idx[8]=idx_node+5;
		signalIq.gain[8]=-cos_theta[2]/(3.0*X_qdp);		
		signalIq.n=9;
		
		signalId.idx[3]=idx_node;
		signalId.gain[3]=-cos_theta[0]/(3.0*X_ddp);
		signalId.idx[4]=idx_node+1;
		signalId.gain[4]=-cos_theta[1]/(3.0*X_ddp);
		signalId.idx[5]=idx_node+2;
		signalId.gain[5]=-cos_theta[2]/(3.0*X_ddp);		
		signalId.idx[6]=idx_node+3;
		signalId.gain[6]=-sin_theta[0]/(3.0*X_ddp);		
		signalId.idx[7]=idx_node+4;
		signalId.gain[7]=-sin_theta[1]/(3.0*X_ddp);		
		signalId.idx[8]=idx_node+5;
		signalId.gain[8]=-sin_theta[2]/(3.0*X_ddp);			
		signalId.n=9;
		
	}
	
  row[0] = Id_idx;  
  col[0] = Id_idx;
  val[0] = -1;
  ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
  for (j = 0; j < signalId.n; j++){
	  col[j] = signalId.idx[j];
      val[j] = signalId.gain[j]; // Coefficients of input U
  }  
  ierr = MatSetValues(J,1,row,signalId.n,col,val,ADD_VALUES);CHKERRQ(ierr);  
  
  row[0] = Iq_idx;  
  col[0] = Iq_idx;
  val[0] = -1;
  ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
  for (j = 0; j < signalIq.n; j++){
	  col[j] = signalIq.idx[j];
      val[j] = signalIq.gain[j]; // Coefficients of input U
  }  
  ierr = MatSetValues(J,1,row,signalIq.n,col,val,ADD_VALUES);CHKERRQ(ierr);

  signalId.idx[0]=Id_idx;
  signalId.gain[0]=1.0;
  signalId.n=1;
  
  signalIq.idx[0]=Iq_idx;
  signalIq.gain[0]=1.0;
  signalIq.n=1; 
  
  /********/

  if(MAC_INFO->exc_model[i] == 0){
	  if (Gen->net_mode==1) {
		    EFD = ugen_arr[2*i];
			EFD_init_idx=rstart +  Gen->dyn_size + 2*i;
	  }
  } else {
	EFD_idx = GENROU->EFD_idx[idx_gen];
	EFD = xgen_arr[EFD_idx];
    EFD_idx += rstart;
  }

  if(MAC_INFO->tgov_model[i] == 0){
	  if (Gen->net_mode==1) {
		    PM = ugen_arr[2*i + 1]*MBoSB; // convert from MBASE to SBASE
			PM_init_idx=rstart +  Gen->dyn_size + 2*i+1;
	  }
  } else {
	PM_idx = GENROU->PM_idx[idx_gen];
	PM = xgen_arr[PM_idx]*MBoSB; // convert from MBASE to SBASE
	PM_idx += rstart;
  }
  
  /* psi1d */
  signl signalpsi1din;
  signalpsi1din.idx[0]=psi1d_idx;
  signalpsi1din.gain[0]=-1.0;
  signalpsi1din.idx[1]=Eqp_idx;
  signalpsi1din.gain[1]=1.0;
  signalpsi1din.idx[2]=Id_idx;
  signalpsi1din.gain[2]=-Gaind1;
  signalpsi1din.n=3;
    
  row[0] = psi1d_idx;  
  if (Gen->net_mode!=2) {
	  
	for (j = 0; j < signalpsi1din.n; j++){
	  col[j] = signalpsi1din.idx[j];
      val[j] = signalpsi1din.gain[j]/T_d0dp; // Coefficients of input U
	} 
	ierr = MatSetValues(J,1,row, signalpsi1din.n,col,val,ADD_VALUES);CHKERRQ(ierr);	
	
	if (Gen->net_mode==0) {
		col[0] = psi1d_idx;  
		val[0] = -a;
		ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);	
	}
	
  } else {
	col[0] = psi1d_idx;  
	val[0] = 1;
	ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
  }

  /* psi2q */
  signl signalpsi2qin;
  signalpsi2qin.idx[0]=psi2q_idx;
  signalpsi2qin.gain[0]=-1.0;
  signalpsi2qin.idx[1]=Edp_idx;
  signalpsi2qin.gain[1]=-1.0;
  signalpsi2qin.idx[2]=Iq_idx;
  signalpsi2qin.gain[2]=-Gainq1;
  signalpsi2qin.n=3;  
  
  row[0] = psi2q_idx;  
  if (Gen->net_mode!=2) {
	  
	for (j = 0; j < signalpsi2qin.n; j++){
	  col[j] = signalpsi2qin.idx[j];
      val[j] = signalpsi2qin.gain[j]/T_q0dp; // Coefficients of input U
	} 
	ierr = MatSetValues(J,1,row, signalpsi2qin.n,col,val,ADD_VALUES);CHKERRQ(ierr);	
	
	if (Gen->net_mode==0) {
		col[0] = psi2q_idx;  
		val[0] = -a;
		ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);	
	}
	
  } else {
	col[0] = psi2q_idx;  
	val[0] = 1;
	ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
  }
   
  /* DIFD */
  PetscScalar Gain1,Gain2;
  signl signalDIFD;
  if (Sat_flag==1) {
	    Gain1=Scoef*2*(Psidp-Esat);
	    for (j = 0; j < signalPsidp.n; j++){
			signalDIFD.idx[j] = signalPsidp.idx[j];
			signalDIFD.gain[j] = signalPsidp.gain[j]*Gain1; // Coefficients of input U
		}
		signalDIFD.n=signalPsidp.n;
		DIFD=Scoef*pow(Psidp-Esat,2);
  }
  
  /* DIFDd */
  signl signalDIFDd;
  if (Sat_flag==1){
	  	for (j = 0; j < signalDIFD.n; j++){
			signalDIFDd.idx[j] = signalDIFD.idx[j];
			signalDIFDd.gain[j] = signalDIFD.gain[j]*Psiddp; // Coefficients of input U
		}
		signalDIFDd.n=signalDIFD.n;
		signalDIFDd.gain[0]+=DIFD*signalPsiddp.gain[0];
		signalDIFDd.gain[1]+=DIFD*signalPsiddp.gain[1];
		Gain1=-1/(Psidp*Psidp);
	  	for (j = 0; j < signalPsidp.n; j++){
			signalOne.idx[j] = signalPsidp.idx[j];
			signalOne.gain[j] = signalPsidp.gain[j]*Gain1; // Coefficients of input U
		}
		Gain1=1/Psidp;Gain2=DIFD*Psiddp;
	  	for (j = 0; j < signalDIFDd.n; j++){
			signalDIFDd.gain[j] *= Gain1; // Coefficients of input U
			signalDIFDd.gain[j] += signalOne.gain[j]*Gain2;
		}
  }
  
  /* DIFDq */
  signl signalDIFDq;
  if (Sat_flag==1){
	  	for (j = 0; j < signalDIFD.n; j++){
			signalDIFDq.idx[j] = signalDIFD.idx[j];
			signalDIFDq.gain[j] = signalDIFD.gain[j]*Psiqdp; // Coefficients of input U
		}
		signalDIFDq.n=signalDIFD.n;
		signalDIFDq.gain[2]+=DIFD*signalPsiqdp.gain[0];
		signalDIFDq.gain[3]+=DIFD*signalPsiqdp.gain[1];
		Gain1=Gainq6/Psidp;Gain2=DIFD*Psiqdp*Gainq6;
	  	for (j = 0; j < signalDIFDq.n; j++){
			signalDIFDq.gain[j] *= Gain1; // Coefficients of input U
			signalDIFDq.gain[j] += signalOne.gain[j]*Gain2;
		}
  }
   
  /* IFD */
  signl signalIFD;
  signalIFD.idx[0]=Eqp_idx;
  signalIFD.gain[0]=1.0+Gaind4*Gaind5;
  signalIFD.idx[1]=psi1d_idx;
  signalIFD.gain[1]=-Gaind4*Gaind5;
  signalIFD.idx[2]=Id_idx;
  signalIFD.gain[2]=Gaind5*(1.0-Gaind1*Gaind4);
  signalIFD.n=3;
       
  if (Sat_flag==1) {
	  /*signalIFD.gain[0]+=signalDIFDd.gain[0];
	  signalIFD.gain[1]+=signalDIFDd.gain[1];
	  signalIFD.idx[3]=signalDIFDd.idx[2];
	  signalIFD.gain[3]=signalDIFDd.gain[2];
	  signalIFD.idx[4]=signalDIFDd.idx[3];
	  signalIFD.gain[4]=signalDIFDd.gain[3];
	  signalIFD.n=5;*/
	  joinSignal(&signalIFD,&signalDIFDd);
  }
  
  if(GENROU->IFD_idx[idx_gen] != -1){
	    row[0] = IFD_idx;  
        col[0] = IFD_idx;
        val[0] = -1;
        for (j = 0; j < signalIFD.n; j++){
			col[1 + j] = signalIFD.idx[j];
			val[1 + j] = signalIFD.gain[j]; // Coefficients of input U
		}  
		ierr = MatSetValues(J,1,row,1 + signalIFD.n,col,val,ADD_VALUES);CHKERRQ(ierr);
	  
	  signalIFD.idx[0]=IFD_idx;
	  signalIFD.gain[0]=1.0;
	  signalIFD.n=1;
  }

  /* Eqp */
  row[0] = Eqp_idx;  
  if (Gen->net_mode!=2) {
	  
	for (j = 0; j < signalIFD.n; j++){
	  col[j] = signalIFD.idx[j];
      val[j] = -signalIFD.gain[j]/T_d0p; // Coefficients of input U
	} 
	ierr = MatSetValues(J,1,row, signalIFD.n,col,val,ADD_VALUES);CHKERRQ(ierr);	
	if (MAC_INFO->exc_model[i] != 0) {
		col[0] = EFD_idx;
		val[0] = 1/T_d0p;
	    ierr = MatSetValues(J,1,row, 1,col,val,ADD_VALUES);CHKERRQ(ierr);	
	} else {
	  if (Gen->net_mode==1) {
		    col[0] = EFD_init_idx;
			val[0] = 1/T_d0p;
			ierr = MatSetValues(J,1,row, 1,col,val,ADD_VALUES);CHKERRQ(ierr);
	  }
	}
	
	if (Gen->net_mode==0) {
		col[0] = Eqp_idx;  
		val[0] = -a;
		ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);	
	}
	
  } else {
	col[0] = Eqp_idx;  
	val[0] = 1;
	ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
  }
  
  /* Edp */
  signl signalEdpin;
  signalEdpin.idx[0]=Edp_idx;
  signalEdpin.gain[0]=-Gainq4*Gainq5-1.0;
  signalEdpin.idx[1]=psi2q_idx;
  signalEdpin.gain[1]=-Gainq4*Gainq5;
  signalEdpin.idx[2]=Iq_idx;
  signalEdpin.gain[2]=Gainq5*(1.0-Gainq1*Gainq4);
  signalEdpin.n=3;
       
  if (Sat_flag==1) {
	  /*signalEdpin.gain[0]+=signalDIFDq.gain[2];
	  signalEdpin.gain[1]+=signalDIFDq.gain[3];
	  signalEdpin.idx[3]=signalDIFDq.idx[0];
	  signalEdpin.gain[3]=signalDIFDq.gain[0];
	  signalEdpin.idx[4]=signalDIFDq.idx[1];
	  signalEdpin.gain[4]=signalDIFDq.gain[1];
	  signalEdpin.n=5;*/
	  joinSignal(&signalEdpin,&signalDIFDq); 
  }
  
  row[0] = Edp_idx;  
  if (Gen->net_mode!=2) {
	  
	for (j = 0; j < signalEdpin.n; j++){
	  col[j] = signalEdpin.idx[j];
      val[j] = signalEdpin.gain[j]/T_q0p; // Coefficients of input U
	} 
	ierr = MatSetValues(J,1,row, signalEdpin.n,col,val,ADD_VALUES);CHKERRQ(ierr);	
	
	if (Gen->net_mode==0) {
		col[0] = Edp_idx;  
		val[0] = -a;
		ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);	
	}
	
  } else {
	col[0] = Edp_idx;  
	val[0] = 1;
	ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
  } 

  /* pelec */
  row[0] = pelec_idx;  
  col[0] = pelec_idx;
  val[0] = -1;
  col[1] = Iq_idx;
  val[1] = Psiddp/MBoSB;
  col[2] = Id_idx;
  val[2] = -Psiqdp/MBoSB;  
  ierr = MatSetValues(J,1,row,3,col,val,ADD_VALUES);CHKERRQ(ierr);
  for (j = 0; j < signalPsiddp.n; j++){
	  col[j] = signalPsiddp.idx[j];
      val[j] = signalPsiddp.gain[j]*Iq/MBoSB; // Coefficients of input U
  }  
  ierr = MatSetValues(J,1,row,signalPsiddp.n,col,val,ADD_VALUES);CHKERRQ(ierr);  
  /*for (j = 0; j < signalIq.n; j++){
	  col[j] = signalIq.idx[j];
      val[j] = signalIq.gain[j]*Psiddp/MBoSB; // Coefficients of input U
  }  
  ierr = MatSetValues(J,1,row,signalIq.n,col,val,ADD_VALUES);CHKERRQ(ierr);*/
  for (j = 0; j < signalPsiqdp.n; j++){
	  col[j] = signalPsiqdp.idx[j];
      val[j] = -signalPsiqdp.gain[j]*Id/MBoSB; // Coefficients of input U
  }  
  ierr = MatSetValues(J,1,row,signalPsiqdp.n,col,val,ADD_VALUES);CHKERRQ(ierr);  
  /*for (j = 0; j < signalId.n; j++){
	  col[j] = signalId.idx[j];
      val[j] = -signalId.gain[j]*Psiqdp/MBoSB; // Coefficients of input U
  }  
  ierr = MatSetValues(J,1,row,signalId.n,col,val,ADD_VALUES);CHKERRQ(ierr);*/

  /* n (SPEED DERIVATION) */
  H=2*H;
  row[0] = n_idx;   
  if (Gen->net_mode!=2) {
	  
	col[0] = n_idx;  
	val[0] = -(PM+D)/(H*pow(1+n,2));
	col[1] = pelec_idx;
	val[1] = -MBoSB/H;
	ierr = MatSetValues(J,1,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
	
	if(MAC_INFO->tgov_model[i] == 0){
	  if (Gen->net_mode==1) {
		  	col[0] = PM_init_idx;
			val[0] = MBoSB/(H*(n+1));
			ierr = MatSetValues(J,1,row, 1,col,val,ADD_VALUES);CHKERRQ(ierr);
	  }
    } else {
		col[0] = PM_idx;
		val[0] = MBoSB/(H*(n+1));
	    ierr = MatSetValues(J,1,row, 1,col,val,ADD_VALUES);CHKERRQ(ierr);	
    }
  	
	if (Gen->net_mode==0) {
		col[0] = n_idx;  
		val[0] = -a;
		ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);	
	}
	
  } else {
	col[0] = n_idx;  
	val[0] = 1;
	ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
  }
  
  /* delta */
  row[0] = delta_idx;  
  if (Gen->net_mode!=2) {
	  
	col[0] = n_idx;  
	val[0] = 120*PI;
	ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	if (Gen->net_mode==0) {
		col[0] = delta_idx;  
		val[0] = -a;
		ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);	
	}
	
  } else {
	col[0] = delta_idx;  
	val[0] = 1;
	ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
  }

  /*network variables*/   
  if (Gen->net_mode != 1) {
	 
    row[0]  = idx_node;
	/*for (j = 0; j < signalId.n; j++){
	  col[j] = signalId.idx[j];
      val[j] = signalId.gain[j]*(-cos_theta[0]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalId.n,col,val,ADD_VALUES);CHKERRQ(ierr);
	for (j = 0; j < signalIq.n; j++){
	  col[j] = signalIq.idx[j];
      val[j] = signalIq.gain[j]*(sin_theta[0]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalIq.n,col,val,ADD_VALUES);CHKERRQ(ierr);
	col[0]  = delta_idx;
	val[0]  = Id*sin_theta[0]+Iq*cos_theta[0];
	ierr = MatSetValues(J,1,row, 1,col,val,ADD_VALUES);CHKERRQ(ierr);*/
	col[0]  = Id_idx;
	val[0]  = -cos_theta[0];
	col[1]  = Iq_idx;
	val[1]  = sin_theta[0];
	col[2]  = delta_idx;
	val[2]  = Id*sin_theta[0]+Iq*cos_theta[0];	
	ierr = MatSetValues(J,1,row, 3,col,val,ADD_VALUES);CHKERRQ(ierr);
	
    row[0]  = idx_node + 1;
	/*for (j = 0; j < signalId.n; j++){
	  col[j] = signalId.idx[j];
      val[j] = signalId.gain[j]*(-cos_theta[1]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalId.n,col,val,ADD_VALUES);CHKERRQ(ierr);  
	for (j = 0; j < signalIq.n; j++){
	  col[j] = signalIq.idx[j];
      val[j] = signalIq.gain[j]*(sin_theta[1]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalIq.n,col,val,ADD_VALUES);CHKERRQ(ierr);
	col[0]  = delta_idx;
	val[0]  = Id*sin_theta[1]+Iq*cos_theta[1];
	ierr = MatSetValues(J,1,row, 1,col,val,ADD_VALUES);CHKERRQ(ierr);*/
	col[0]  = Id_idx;
	val[0]  = -cos_theta[1];
	col[1]  = Iq_idx;
	val[1]  = sin_theta[1];
	col[2]  = delta_idx;
	val[2]  = Id*sin_theta[1]+Iq*cos_theta[1];
	ierr = MatSetValues(J,1,row, 3,col,val,ADD_VALUES);CHKERRQ(ierr);
	
	row[0]  = idx_node + 2;
	/*for (j = 0; j < signalId.n; j++){
	  col[j] = signalId.idx[j];
      val[j] = signalId.gain[j]*(-cos_theta[2]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalId.n,col,val,ADD_VALUES);CHKERRQ(ierr);  
	for (j = 0; j < signalIq.n; j++){
	  col[j] = signalIq.idx[j];
      val[j] = signalIq.gain[j]*(sin_theta[2]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalIq.n,col,val,ADD_VALUES);CHKERRQ(ierr);
	col[0]  = delta_idx;
	val[0]  = Id*sin_theta[2]+Iq*cos_theta[2];
	ierr = MatSetValues(J,1,row, 1,col,val,ADD_VALUES);CHKERRQ(ierr);*/
	col[0]  = Id_idx;
	val[0]  = -cos_theta[2];
	col[1]  = Iq_idx;
	val[1]  = sin_theta[2];
	col[2]  = delta_idx;
	val[2]  = Id*sin_theta[2]+Iq*cos_theta[2];
	ierr = MatSetValues(J,1,row, 3,col,val,ADD_VALUES);CHKERRQ(ierr);
	
    row[0]  = idx_node + 3;
	/*for (j = 0; j < signalId.n; j++){
	  col[j] = signalId.idx[j];
      val[j] = signalId.gain[j]*(sin_theta[0]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalId.n,col,val,ADD_VALUES);CHKERRQ(ierr);  
	for (j = 0; j < signalIq.n; j++){
	  col[j] = signalIq.idx[j];
      val[j] = signalIq.gain[j]*(cos_theta[0]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalIq.n,col,val,ADD_VALUES);CHKERRQ(ierr);
	col[0]  = delta_idx;
	val[0]  = Id*cos_theta[0]-Iq*sin_theta[0];
	ierr = MatSetValues(J,1,row, 1,col,val,ADD_VALUES);CHKERRQ(ierr);*/
	col[0]  = Id_idx;
	val[0]  = sin_theta[0];
	col[1]  = Iq_idx;
	val[1]  = cos_theta[0];
	col[2]  = delta_idx;
	val[2]  = Id*cos_theta[0]-Iq*sin_theta[0];
	ierr = MatSetValues(J,1,row, 3,col,val,ADD_VALUES);CHKERRQ(ierr);
	
	row[0]  = idx_node + 4;
	/*for (j = 0; j < signalId.n; j++){
	  col[j] = signalId.idx[j];
      val[j] = signalId.gain[j]*(sin_theta[1]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalId.n,col,val,ADD_VALUES);CHKERRQ(ierr);  
	for (j = 0; j < signalIq.n; j++){
	  col[j] = signalIq.idx[j];
      val[j] = signalIq.gain[j]*(cos_theta[1]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalIq.n,col,val,ADD_VALUES);CHKERRQ(ierr);
	col[0]  = delta_idx;
	val[0]  = Id*cos_theta[1]-Iq*sin_theta[1];
	ierr = MatSetValues(J,1,row, 1,col,val,ADD_VALUES);CHKERRQ(ierr);*/
	col[0]  = Id_idx;
	val[0]  = sin_theta[1];
	col[1]  = Iq_idx;
	val[1]  = cos_theta[1];
	col[2]  = delta_idx;
	val[2]  = Id*cos_theta[1]-Iq*sin_theta[1];
	ierr = MatSetValues(J,1,row, 3,col,val,ADD_VALUES);CHKERRQ(ierr);	
	
	row[0]  = idx_node + 5;
	/*for (j = 0; j < signalId.n; j++){
	  col[j] = signalId.idx[j];
      val[j] = signalId.gain[j]*(sin_theta[2]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalId.n,col,val,ADD_VALUES);CHKERRQ(ierr);  
	for (j = 0; j < signalIq.n; j++){
	  col[j] = signalIq.idx[j];
      val[j] = signalIq.gain[j]*(cos_theta[2]); // Coefficients of input U
	}  
	ierr = MatSetValues(J,1,row, signalIq.n,col,val,ADD_VALUES);CHKERRQ(ierr);
	col[0]  = delta_idx;
	val[0]  = Id*cos_theta[2]-Iq*sin_theta[2];
	ierr = MatSetValues(J,1,row, 1,col,val,ADD_VALUES);CHKERRQ(ierr);*/
	col[0]  = Id_idx;
	val[0]  = sin_theta[2];
	col[1]  = Iq_idx;
	val[1]  = cos_theta[2];
	col[2]  = delta_idx;
	val[2]  = Id*cos_theta[2]-Iq*sin_theta[2];
	ierr = MatSetValues(J,1,row, 3,col,val,ADD_VALUES);CHKERRQ(ierr);
	
  }
    
  /* Restore arrays */
  ierr = VecRestoreArray(Xgen,&xgen_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xdot,&xdot_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(Fgen,&fgen_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(Ugen,&ugen_arr);CHKERRQ(ierr);
  ierr = VecRestoreArray(V,&v_arr);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
  
}

/* LOG

TODO: optimize equations using python scripts


*/
