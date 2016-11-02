/* --------------------------------------------------
   Routine for GENROU Model Jacobian Preallocation
   -------------------------------------------------- */

#include "../../include/TS3ph.h"
#undef __FUNCT__
#define __FUNCT__ "PreallocateJacobianGENROU"
PetscErrorCode PreallocateJacobianGENROU(TS3phDynCtx* Gen, PetscInt* d_nnz, PetscInt i)
{
  PetscErrorCode ierr;
  
  GENDATA        *Gdata = Gen->GenData;
  MAC_INFO_data  *MAC_INFO = Gdata->MAC_INFO;
  GENROU_MODEL   *GENROU = Gdata->GEN_MODEL->GENROU;
  GMAP_data      *GMAP = Gdata->GMAP;
  
  PetscInt       row[1];
  
  PetscInt       idx_gen,idx_node, bus_index, seq_number;;      
  
  PetscInt     Eqp_idx,Edp_idx,psi1d_idx,psi2q_idx,n_idx,delta_idx;
  PetscInt     IFD_idx,pelec_idx,Id_idx,Iq_idx; 
  
  PetscFunctionBegin;

  idx_gen = MAC_INFO->gen_idx[i];
  
  /*State Variables indices*/
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

  bus_index  = MAC_INFO->gen_bus[i];
  seq_number = bus_index;
  
  idx_node    =Gen->dyn_size + twonphase*seq_number;  
  
  /* psi1d */ 
  d_nnz[psi1d_idx] += 3; 
  
  /* psi2q */
  d_nnz[psi2q_idx] += 3;
  
  /* IFD */  
  if(GENROU->IFD_idx[idx_gen] != -1){
	  d_nnz[IFD_idx] += 6;  
  }  
  
  /* Eqp */ 
  if (GENROU->IFD_idx[idx_gen] != -1) {
      d_nnz[Eqp_idx] += 1;
  } else {
      d_nnz[Eqp_idx] += 5;
  }
	  
  if (MAC_INFO->exc_model[i] != 0) {
	d_nnz[Eqp_idx] += 1;
  } else {
	  if (Gen->net_mode==1) {
		d_nnz[Eqp_idx] += 1;  
	  }
  }
   
  /* Edp */    
  d_nnz[Edp_idx] += 5;
  
  /* pelec */ 
  d_nnz[pelec_idx] += 7;  
  
  /* Id */
  d_nnz[Id_idx] += 4;  
  if(Gen->net_mode!=1){
    d_nnz[Id_idx] += 6;  
  } 

  /* Iq */  
  d_nnz[Iq_idx] += 4;  
  if(Gen->net_mode!=1){
    d_nnz[Iq_idx] += 6;  
  }  
  
  /* n (SPEED DERIVATION) */ 
  d_nnz[n_idx] += 2;
  if(MAC_INFO->tgov_model[i] != 0){
		d_nnz[n_idx] += 1;
  } else {
	  if (Gen->net_mode==1) {	  
		d_nnz[n_idx] += 1;  
	  }
  }
  
  /* delta */
  d_nnz[delta_idx] += 2;
  
  /*network variables*/
  if(Gen->net_mode == 1) {
    d_nnz[Gen->dyn_size + i*2] += 4;
    d_nnz[Gen->dyn_size + i*2 + 1] += 4;
  } else { 
	/*d_nnz[idx_node + 3] += 9;
	d_nnz[idx_node + 4] += 9;  
	d_nnz[idx_node + 5] += 9;    
	d_nnz[idx_node] += 9;
	d_nnz[idx_node + 1] += 9;
	d_nnz[idx_node + 2] += 9;*/
	
	d_nnz[idx_node + 3] += 3;
	d_nnz[idx_node + 4] += 3;  
	d_nnz[idx_node + 5] += 3;    
	d_nnz[idx_node] += 3;
	d_nnz[idx_node + 1] += 3;
	d_nnz[idx_node + 2] += 3;
	
	/*d_nnz[idx_node + 3] += 11;
	d_nnz[idx_node + 4] += 11;  
	d_nnz[idx_node + 5] += 11;    
	d_nnz[idx_node] += 11;
	d_nnz[idx_node + 1] += 11;
	d_nnz[idx_node + 2] += 11;*/
  }

  PetscFunctionReturn(0);
  
}
