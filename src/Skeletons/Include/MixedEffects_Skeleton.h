#ifndef __MIXEDEFFECTS_SKELETON_H__
#define __MIXEDEFFECTS_SKELETON_H__

#include "../../FdaPDE.h"
#include "../../Regression/Include/FPIRLS.h"
//#include "../../Regression/Include/FPIRLS_Factory.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"

template<typename InputHandler,UInt ORDER, UInt mydim, UInt ndim>
SEXP MixedEffects_skeleton(InputHandler & MixedEffectsData, OptimizationData & optimizationData, SEXP Rmesh/*, SEXP Rmu0, std::string family, SEXP RscaleParam*/)
{
  MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, MixedEffectsData.getSearch());

	// Factory:
	FPIRLS_MixedEffects<InputHandler, ORDER, mydim, ndim> fpirls(mesh, MixedEffectsData, optimizationData);

	fpirls.apply();
	
	const MatrixXv& solution = fpirls.getSolution();
  	const MatrixXr& dof = fpirls.getDOF();
  	const std::vector<std::vector<Real>>& J_value = fpirls.get_J();
  	const MatrixXv& fn_hat = fpirls.getFunctionEst();
  	const std::vector<std::vector<VectorXr>>& Sigma_b = fpirls.getSigma_b();
  	const std::vector<std::vector<std::vector<VectorXr>>>& b_hat = fpirls.get_b_hat();
  	const std::vector<std::vector<Real>>& GCV = fpirls.getGCV();

  	const UInt bestLambda = optimizationData.get_best_lambda_S();
  	const UInt lambdaS_len = fpirls.get_size_S();
  	const UInt lambdaT_len = 1;

	const std::vector<std::vector<UInt>>& n_iterations = fpirls.getIterations();

	MatrixXv beta;
  	if(MixedEffectsData.getCovariates()->rows()==0)
   	{
		beta.resize(1,1);
		beta(0,0).resize(1);
		beta(0,0)(0) = 10e20;
	}
	else
		 beta = fpirls.getBetaEst();

	const MatrixXr & barycenters = fpirls.getBarycenters();
	const VectorXi & elementIds = fpirls.getElementIds();

  	// COMPOSIZIONE SEXP result FOR RETURN

	//Copy result in R memory
	SEXP result = R_NilValue;
 	result = PROTECT(Rf_allocVector(VECSXP, 5+3+6+2+1));
  	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution(0).size(), solution.size()));
  	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, dof.size()));
  	SET_VECTOR_ELT(result, 2, Rf_allocMatrix(REALSXP, lambdaS_len, lambdaT_len));
  	SET_VECTOR_ELT(result, 3, Rf_allocVector(INTSXP, 1));
  	SET_VECTOR_ELT(result, 4, Rf_allocMatrix(REALSXP, beta(0).size(), beta.size()));

	//return solution
	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < solution.size(); j++)
	{
		for(UInt i = 0; i < solution(0).size(); i++)
			rans[i + solution(0).size()*j] = solution(j)(i);
	}

	//return DoF
	//std::cout << "rans1" << std::endl;
	Real *rans1 = REAL(VECTOR_ELT(result, 1));
	for(UInt i = 0; i < dof.size(); i++)
	{
		rans1[i] = dof(i);
	}

	//return GCV values
	//std::cout << "rans2" << std::endl;
  	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < lambdaS_len; i++)
		for(UInt j = 0; j < lambdaT_len; j++)
			rans2[i + lambdaS_len * j] = GCV[i][j];

	// Copy best lambda
	//std::cout << "rans3" << std::endl;
	UInt *rans3 = INTEGER(VECTOR_ELT(result, 3));
	rans3[0] = bestLambda;

	//return beta hat
	//std::cout << "rans4" << std::endl;
	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt j = 0; j < beta.size(); j++)
	{
		for(UInt i = 0; i < beta(0).size(); i++)
			rans4[i + beta(0).size()*j] = beta(j)(i);
	}

	//std::cout << "rans-entro-if" << std::endl;
	if(MixedEffectsData.getSearch()==2){

		//SEND TREE INFORMATION TO R
	//std::cout << "rans5" << std::endl;
		SET_VECTOR_ELT(result, 5, Rf_allocVector(INTSXP, 1)); //tree_header information
		int *rans5 = INTEGER(VECTOR_ELT(result, 5));
		rans5[0] = mesh.getTree().gettreeheader().gettreelev();

	//std::cout << "rans6" << std::endl;
		SET_VECTOR_ELT(result, 6, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
		Real *rans6 = REAL(VECTOR_ELT(result, 6));
		for(UInt i = 0; i < ndim*2; i++)
			rans6[i] = mesh.getTree().gettreeheader().domainorig(i);

	//std::cout << "rans7" << std::endl;
		SET_VECTOR_ELT(result, 7, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
		Real *rans7 = REAL(VECTOR_ELT(result, 7));
		for(UInt i = 0; i < ndim*2; i++)
			rans7[i] = mesh.getTree().gettreeheader().domainscal(i);


		UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
	///std::cout << "rans8" << std::endl;
		SET_VECTOR_ELT(result, 8, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
		int *rans8 = INTEGER(VECTOR_ELT(result, 8));
		for(UInt i = 0; i < num_tree_nodes; i++)
				rans8[i] = mesh.getTree().gettreenode(i).getid();

		for(UInt i = 0; i < num_tree_nodes; i++)
				rans8[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

		for(UInt i = 0; i < num_tree_nodes; i++)
				rans8[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

	//std::cout << "rans9" << std::endl;
		SET_VECTOR_ELT(result, 9, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
		Real *rans9 = REAL(VECTOR_ELT(result, 9));
		for(UInt j = 0; j < ndim*2; j++)
		{
			for(UInt i = 0; i < num_tree_nodes; i++)
				rans9[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
		}
	}
	
	//SEND BARYCENTER INFORMATION TO R
	//std::cout << "rans10" << std::endl;
	SET_VECTOR_ELT(result, 10, Rf_allocVector(INTSXP, elementIds.rows())); //element id of the locations point (vector)
	int *rans10 = INTEGER(VECTOR_ELT(result, 10));
	for(UInt i = 0; i < elementIds.rows(); i++)
		rans10[i] = elementIds(i);

	//std::cout << "rans11" << std::endl;
	SET_VECTOR_ELT(result, 11, Rf_allocMatrix(REALSXP, barycenters.rows(), barycenters.cols())); //barycenter information (matrix)
	Real *rans11 = REAL(VECTOR_ELT(result, 11));
	for(UInt j = 0; j < barycenters.cols(); j++)
	{
		for(UInt i = 0; i < barycenters.rows(); i++)
			rans11[i + barycenters.rows()*j] = barycenters(i,j);
	}

	// GAM PARAMETER ESTIMATIONS
	SET_VECTOR_ELT(result, 12, Rf_allocMatrix(REALSXP, fn_hat(0,0).size(), fn_hat.rows()*fn_hat.cols()));
	SET_VECTOR_ELT(result, 13, Rf_allocMatrix(REALSXP, lambdaS_len, lambdaT_len));
	SET_VECTOR_ELT(result, 14, Rf_allocMatrix(REALSXP, Sigma_b[0][0].size(), lambdaS_len*lambdaT_len));
	SET_VECTOR_ELT(result, 15, Rf_allocMatrix(REALSXP, b_hat[0][0][0].size()*MixedEffectsData.getGroupNumber(), lambdaS_len*lambdaT_len));
	

	//return fn hat
	//std::cout << "rans12" << std::endl;
	Real *rans12 = REAL(VECTOR_ELT(result, 12));
	for(UInt i = 0; i < fn_hat.rows(); i++){
		for(UInt j = 0; j < fn_hat.cols(); j++)
			for(UInt k = 0; k < fn_hat(0,0).size(); k++)
				rans12[k + fn_hat(0,0).size()*i + fn_hat(0,0).size()*fn_hat.rows()*j] = fn_hat.coeff(i,j)(k);
	}	
	for (UInt j = 0; j < fn_hat.size(); j++) {
        	for (UInt i = 0; i < fn_hat(0, 0).size(); i++)
            		rans12[i + fn_hat(0 , 0).size() * j] = fn_hat(j)(i);
    	}

	//return J_value
	//std::cout << "rans13" << std::endl;
  	Real *rans13 = REAL(VECTOR_ELT(result, 13));
  	for(UInt i = 0; i < lambdaS_len; i++){
  		for(UInt j = 0; j < lambdaT_len; j++)
			rans13[i + lambdaS_len*j] = J_value[i][j];
	}

	//return Sigma_b (estimated variance of random effects)
	//std::cout << "rans14" << std::endl;
	Real *rans14 = REAL(VECTOR_ELT(result, 14));
	for(UInt i = 0; i < lambdaS_len; i++){
		for(UInt j = 0; j < lambdaT_len; j++){
			for(UInt k = 0; k < Sigma_b[0][0].size(); k++){
				rans14[k + Sigma_b[0][0].size()*i + Sigma_b[0][0].size()*lambdaS_len*j] = Sigma_b[i][j].coeff(k);
			}
		}
	}

	// return b_i_hat (prediction of covariates for each group)
	//std::cout << "rans15" << std::endl;
	Real *rans15 = REAL(VECTOR_ELT(result, 15));
	for(UInt i = 0; i < lambdaS_len; i++){
		for(UInt j = 0; j < lambdaT_len; j++){
			for(UInt k = 0; k < b_hat[0][0][0].size(); k++){
				for(UInt h = 0; h < MixedEffectsData.getGroupNumber(); h++)
				rans15[k + b_hat[0][0][0].size()*h + b_hat[0][0][0].size()*MixedEffectsData.getGroupNumber()*i + b_hat[0][0][0].size()*MixedEffectsData.getGroupNumber()*lambdaS_len*j] = b_hat[i][j][h].coeff(k);
			}
		}
	}

	//return DoF
	//std::cout << "rans1" << std::endl;
	SET_VECTOR_ELT(result, 16, Rf_allocVector(REALSXP, n_iterations.size()));
	Real *rans16 = REAL(VECTOR_ELT(result, 16));
	for(UInt i = 0; i < n_iterations.size(); i++)
	{
		rans16[i] = n_iterations[i][0];
	}

	UNPROTECT(1);

	return(result);

}

#endif
