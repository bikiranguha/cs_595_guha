/* --------------------------------------------------                        
   Routine for GENROU Model Jacobian Structure Set Up                           
   -------------------------------------------------- */

#include "../../include/TS3ph.h"
#undef __FUNCT__
#define __FUNCT__ "SetUpJacobianStructureGENROU"
PetscErrorCode SetUpJacobianStructureGENROU(TS3phDynCtx* Gen, Mat L, PetscInt i, PetscInt rstart)
{
  PetscErrorCode ierr;
  
  GENDATA        *Gdata = Gen->GenData;
  MAC_INFO_data  *MAC_INFO = Gdata->MAC_INFO;
  GENROU_MODEL   *GENROU = Gdata->GEN_MODEL->GENROU;
  GMAP_data      *GMAP = Gdata->GMAP;
  
  PetscInt       row[1], col[30];
  PetscScalar    values[30] = {0};
  
  PetscInt       idx_gen,idx_node, bus_index, seq_number;

  PetscInt     Eqp_idx,Edp_idx,psi1d_idx,psi2q_idx,n_idx,delta_idx;
  PetscInt     EFD_idx,PM_idx,IFD_idx,pelec_idx,Id_idx,Iq_idx; 
  
  PetscFunctionBegin;

  idx_gen = MAC_INFO->gen_idx[i];

  /*State Variables indices*/
  Eqp_idx = GENROU->Eqp_idx[idx_gen]+ rstart;
  Edp_idx = GENROU->Edp_idx[idx_gen]+ rstart;
  psi1d_idx = GENROU->psi1d_idx[idx_gen]+ rstart;
  psi2q_idx = GENROU->psi2q_idx[idx_gen]+ rstart;
  n_idx = GENROU->n_idx[idx_gen]+ rstart;
  delta_idx = GENROU->delta_idx[idx_gen]+ rstart;
  
  pelec_idx = GENROU->pelec_idx[idx_gen]+ rstart;
  Id_idx = GENROU->Id_idx[idx_gen]+ rstart;
  Iq_idx = GENROU->Iq_idx[idx_gen]+ rstart;
  IFD_idx = GENROU->IFD_idx[idx_gen]+ rstart;
 
  bus_index  = MAC_INFO->gen_bus[i];
  seq_number = bus_index;
  
  idx_node    = rstart + Gen->dyn_size + twonphase*seq_number;

  /* psi1d */
  row[0] = psi1d_idx;
  
  col[0] = psi1d_idx;
  col[1] = Eqp_idx;
  col[2] = Id_idx;
  ierr = MatSetValues(L,1,row,3,col,values,ADD_VALUES);CHKERRQ(ierr);  
  
  /* psi2q */
  row[0] = psi2q_idx;
  
  col[0] = psi2q_idx;
  col[1] = Edp_idx;
  col[2] = Iq_idx;
  ierr = MatSetValues(L,1,row,3,col,values,ADD_VALUES);CHKERRQ(ierr);   
   
  /* IFD */  
  if (GENROU->IFD_idx[idx_gen] != -1) {
	  
	    row[0] = IFD_idx; 
		
        col[0] = IFD_idx;
		col[1] = psi1d_idx;
		col[2] = Eqp_idx;
		col[3] = psi2q_idx;
		col[4] = Edp_idx;
		col[5] = Id_idx;		
		ierr = MatSetValues(L,1,row,6,col,values,ADD_VALUES);CHKERRQ(ierr);  
			
  }
  
  /* Eqp */
  row[0] = Eqp_idx;
  
  if (GENROU->IFD_idx[idx_gen] != -1) {
	  col[0] = IFD_idx;
	  ierr = MatSetValues(L,1,row,1,col,values,ADD_VALUES);CHKERRQ(ierr);
  } else {
	  col[0] = Eqp_idx;
	  col[1] = psi1d_idx;
	  col[2] = psi2q_idx;
	  col[3] = Edp_idx;
	  col[4] = Id_idx;
	  ierr = MatSetValues(L,1,row,5,col,values,ADD_VALUES);CHKERRQ(ierr);
  }
  
  if (MAC_INFO->exc_model[i] != 0) {
	EFD_idx = GENROU->EFD_idx[idx_gen];
    EFD_idx += rstart;
	col[0] = EFD_idx;
	ierr = MatSetValues(L,1,row,1,col,values,ADD_VALUES);CHKERRQ(ierr);
  } else {
	  if (Gen->net_mode==1) {
		col[0] = rstart + Gen->dyn_size + 2*i;
		ierr = MatSetValues(L,1,row,1,col,values,ADD_VALUES);CHKERRQ(ierr);
	  }
  }
   
  /* Edp */  
  row[0] = Edp_idx;
  
  col[0] = Edp_idx;
  col[1] = psi1d_idx;
  col[2] = Iq_idx;
  col[3] = psi2q_idx;
  col[4] = Eqp_idx;  
  ierr = MatSetValues(L,1,row,5,col,values,ADD_VALUES);CHKERRQ(ierr);  
 
  /* pelec */  
  row[0] = pelec_idx;
  
  col[0] = Edp_idx;
  col[1] = psi1d_idx;
  col[2] = psi2q_idx;
  col[3] = Eqp_idx;
  col[4] = pelec_idx;
  col[5] = Id_idx;
  col[6] = Iq_idx;
  ierr = MatSetValues(L,1,row,7,col,values,ADD_VALUES);CHKERRQ(ierr);
  
  /* Id */
  row[0] = Id_idx;
  
  col[0] = Id_idx;
  col[1] = psi1d_idx;
  col[2] = Eqp_idx;
  col[3] = delta_idx;
  ierr = MatSetValues(L,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
  
  if(Gen->net_mode!=1){
	col[0] = idx_node;
	col[1] = idx_node+1;
	col[2] = idx_node+2;
	col[3] = idx_node+3;
	col[4] = idx_node+4;
	col[5] = idx_node+5;
	ierr = MatSetValues(L,1,row,6,col,values,ADD_VALUES);CHKERRQ(ierr);  
  } 

  /* Iq */  
  row[0] = Iq_idx;
  
  col[0] = Iq_idx;
  col[1] = psi2q_idx;
  col[2] = Edp_idx;
  col[3] = delta_idx;
  ierr = MatSetValues(L,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
  
  if(Gen->net_mode!=1){
	col[0] = idx_node;
	col[1] = idx_node+1;
	col[2] = idx_node+2;
	col[3] = idx_node+3;
	col[4] = idx_node+4;
	col[5] = idx_node+5;
	ierr = MatSetValues(L,1,row,6,col,values,ADD_VALUES);CHKERRQ(ierr);  
  }

  /* n (SPEED DERIVATION) */
  row[0] = n_idx;
  
  col[0] = n_idx;
  col[1] = pelec_idx;
  ierr = MatSetValues(L,1,row,2,col,values,ADD_VALUES);CHKERRQ(ierr);
  
  if(MAC_INFO->tgov_model[i] != 0){
	PM_idx = GENROU->PM_idx[idx_gen];
	PM_idx += rstart;
	col[0] = PM_idx;
	ierr = MatSetValues(L,1,row,1,col,values,ADD_VALUES);CHKERRQ(ierr);
  } else {
	if (Gen->net_mode==1) {
		col[0] = rstart + Gen->dyn_size + 2*i+1;
		ierr = MatSetValues(L,1,row,1,col,values,ADD_VALUES);CHKERRQ(ierr);
	}
  }
  
  /* delta */
  row[0] = delta_idx;
  
  col[0] = delta_idx;
  col[1] = n_idx;
  ierr = MatSetValues(L,1,row,2,col,values,ADD_VALUES);CHKERRQ(ierr);  
  
  /*network variables*/ 
  if(Gen->net_mode == 1) {
	  
    /*row[0] = rstart + Gen->dyn_size + i*2;
	col[0] = rstart + Gen->dyn_size + i*2;
    col[1] = Id_idx;      
    col[2] = Iq_idx; 	
    col[3] = delta_idx;	
    ierr = MatSetValues(L,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);

    row[0] = rstart + Gen->dyn_size + i*2 + 1;
	col[0] = rstart + Gen->dyn_size + i*2 + 1;
    col[1] = Id_idx;      
    col[2] = Iq_idx;	
    col[3] = delta_idx;
    ierr = MatSetValues(L,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);*/
	
    row[0] = rstart + Gen->dyn_size + i*2;
	col[0] = rstart + Gen->dyn_size + i*2;
    col[1] = Eqp_idx;      
    col[2] = psi1d_idx; 	
    col[3] = delta_idx;	
    ierr = MatSetValues(L,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);

    row[0] = rstart + Gen->dyn_size + i*2 + 1;
	col[0] = rstart + Gen->dyn_size + i*2 + 1;
    col[1] = Edp_idx;      
    col[2] = psi2q_idx;	
    col[3] = delta_idx;
    ierr = MatSetValues(L,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
	
  } else { 
       
	col[0]  = Id_idx;
	col[1]  = Iq_idx;     
	col[2]  = delta_idx;
	/*col[3]  = idx_node;     
	col[4]  = idx_node + 1;   
	col[5]  = idx_node + 2;
	col[6]  = idx_node + 3;   
	col[7]  = idx_node + 4;   
	col[8] = idx_node + 5;*/
  
	row[0]  = idx_node + 3;
	ierr = MatSetValues(L,1,row,3,col,values,ADD_VALUES);CHKERRQ(ierr);  
  
	row[0]  = idx_node + 4;
	ierr = MatSetValues(L,1,row,3,col,values,ADD_VALUES);CHKERRQ(ierr);  

	row[0]  = idx_node + 5;
	ierr = MatSetValues(L,1,row,3,col,values,ADD_VALUES);CHKERRQ(ierr); 
  
	row[0]  = idx_node;
	ierr = MatSetValues(L,1,row,3,col,values,ADD_VALUES);CHKERRQ(ierr); 

	row[0]  = idx_node+1;
	ierr = MatSetValues(L,1,row,3,col,values,ADD_VALUES);CHKERRQ(ierr);  

	row[0]  = idx_node+2;
	ierr = MatSetValues(L,1,row,3,col,values,ADD_VALUES);CHKERRQ(ierr);
	
	/*col[0]  = Eqp_idx;      
	col[1]  = Edp_idx;     
	col[2]  = psi1d_idx;
	col[3]  = psi2q_idx;     
	col[4]  = delta_idx;
	col[5]  = idx_node;     
	col[6]  = idx_node + 1;   
	col[7]  = idx_node + 2;
	col[8]  = idx_node + 3;   
	col[9]  = idx_node + 4;   
	col[10] = idx_node + 5;
  
	row[0]  = idx_node + 3;
	ierr = MatSetValues(L,1,row,11,col,values,ADD_VALUES);CHKERRQ(ierr);  
  
	row[0]  = idx_node + 4;
	ierr = MatSetValues(L,1,row,11,col,values,ADD_VALUES);CHKERRQ(ierr);  

	row[0]  = idx_node + 5;
	ierr = MatSetValues(L,1,row,11,col,values,ADD_VALUES);CHKERRQ(ierr); 
  
	row[0]  = idx_node;
	ierr = MatSetValues(L,1,row,11,col,values,ADD_VALUES);CHKERRQ(ierr); 

	row[0]  = idx_node+1;
	ierr = MatSetValues(L,1,row,11,col,values,ADD_VALUES);CHKERRQ(ierr);  

	row[0]  = idx_node+2;
	ierr = MatSetValues(L,1,row,11,col,values,ADD_VALUES);CHKERRQ(ierr);*/
		
  }
  
  PetscFunctionReturn(0);
  
}
