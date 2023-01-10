#ifndef __FPIRLS_IMP_H__
#define __FPIRLS_IMP_H__

#include "FPIRLS.h"


/*********** FPIRLS_base Methods ************/

// Constructor
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim> & mesh, InputHandler & inputData, OptimizationData & optimizationData):
  mesh_(mesh), inputData_(inputData), optimizationData_(optimizationData), regression_(inputData, optimizationData, mesh.num_nodes()), lenS_(optimizationData.get_size_S()), lenT_(optimizationData.get_size_T())
{
  //Pre-allocate memory for all quatities
  current_J_values.resize(lenS_, std::vector<std::array<Real, 2>>(lenT_));
  past_J_values.resize(lenS_, std::vector<std::array<Real, 2>>(lenT_));
  n_iterations.resize(lenS_, std::vector<UInt>(lenT_));
  _J_minima.resize(lenS_, std::vector<Real>(lenT_));
  _GCV.resize(lenS_, std::vector<Real>(lenT_, -1));
  
  //initialization of current_J_values and past_J_values.
  for(UInt i=0; i<optimizationData_.get_size_S() ; i++){
   for(UInt j=0; j<optimizationData_.get_size_T() ; j++){
    current_J_values[i][j] = std::array<Real,2>{1,1};
    past_J_values[i][j] = std::array<Real,2>{1,1};
   }
  }
};

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim> & mesh, const std::vector<Real>& mesh_time, InputHandler & inputData, OptimizationData & optimizationData):
  mesh_(mesh), mesh_time_(mesh_time), inputData_(inputData), optimizationData_(optimizationData), regression_(mesh_time, inputData, optimizationData, mesh.num_nodes()), lenS_(optimizationData.get_size_S()), lenT_(optimizationData.get_size_T())
{
  //Pre-allocate memory for all quatities
  current_J_values.resize(lenS_, std::vector<std::array<Real, 2>>(lenT_));
  past_J_values.resize(lenS_, std::vector<std::array<Real, 2>>(lenT_));
  n_iterations.resize(lenS_, std::vector<UInt>(lenT_));
  _J_minima.resize(lenS_, std::vector<Real>(lenT_));
  _GCV.resize(lenS_, std::vector<Real>(lenT_, -1));
  
  //initialization of mu, current_J_values and past_J_values.
  for(UInt i=0; i<optimizationData_.get_size_S() ; i++){
   for(UInt j=0; j<optimizationData_.get_size_T() ; j++){
    current_J_values[i][j] = std::array<Real,2>{1,1};
    past_J_values[i][j] = std::array<Real,2>{1,1};
   }
  }
};

// FPIRLS_base methods

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::apply( const ForcingTerm& u){
  // f-PRILS implementation

  // Initialize the outputs. The temporal dimension is not implemented, for this reason the 2nd dimension is set to 1.
  if( this->inputData_.getCovariates()->rows() > 0 )
  	_beta_hat.resize(lenS_, lenT_);
  _fn_hat.resize(lenS_, lenT_);
  _dof.resize(lenS_, lenT_);
  _solution.resize(lenS_,lenT_);

  if(isSpaceVarying)
  {
    FiniteElement<ORDER, mydim, ndim> fe;
  	Assembler::forcingTerm(mesh_, fe, u, forcingTerm);
  }

  for(UInt i=0 ; i < lenS_ ; i++){//for-cycle for each spatial penalization (lambdaS).
   for(UInt j=0 ; j < lenT_ ; j++){
    current_J_values[i][j][0] = past_J_values[i][j][0] + 2*inputData_.get_treshold();
    current_J_values[i][j][1] = past_J_values[i][j][1] + 2*inputData_.get_treshold();

    this->optimizationData_.setCurrentLambda(i, j); // set right lambda for the current iteration.

    //Rprintf("Start FPIRLS for the lambda number %d \n", i+1);

    // start the iterative method for the lambda index i
    while(stopping_criterion(i, j)){

      // STEP (1)
      // Compute the weights and other auxiliary quantities
      prepare_weighted_regression(i, j);

      // STEP (2)
      solve_weighted_regression(i, j);

      // STEP (3)
      update_parameters(i, j);

      // update J
      past_J_values[i][j] = current_J_values[i][j];
      current_J_values[i][j] = compute_J(i, j);

      if( !regression_.isMatrixNoCov_factorized() ) {

         Rprintf("WARNING: System matrix cannot be factorized for optimization parameters in position %d (Space) and  %d (Time). Try increasing optimization parameter.\n", i+1, j+1) ;
          break;
      }
      n_iterations[i][j]++;

    } //end while

    //Rprintf("\t n. iterations: %d\n \n", n_iterations[i]);

    _J_minima[i][j] = current_J_values[i][j][0]+current_J_values[i][j][1]; // compute the minimum value of the J fuctional

    if(this->optimizationData_.get_loss_function()=="GCV"){ // compute GCV if it is required

        if( !regression_.isMatrixNoCov_factorized() ){

            _GCV[i][j] = std::numeric_limits<double>::quiet_NaN();

        }else{

      	   compute_GCV(i,j);

        }

    }
  
  } // end time for
 }// end space for

  // Additional estimates
  additional_estimates();
}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::solve_weighted_regression(const UInt& lambdaS_index, const UInt& lambdaT_index){
  // performs step (2) of PIRLS. It requires pseudo data after step(1) and mimic regression skeleton behaviour

  // Here we have to solve a weighted regression problem.
  regression_.recomputeWTW(); // at each iteration of FPIRLS W is updated, so WTW has to be recomputed as well.
  regression_.preapply(this->mesh_);
  regression_.apply();

  // if the system matrix is correctly factorized
  if( regression_.isMatrixNoCov_factorized() ) {
  	const SpMat * Psi = regression_.getpsi_(); // get Psi matrix. It is used for the computation of fn_hat.

  	// get the solutions from the regression object.
  	_solution(lambdaS_index, lambdaT_index) = regression_.getSolution()(0,0);
	_dof(lambdaS_index, lambdaT_index) = regression_.getDOF()(0,0);

  	if(inputData_.getCovariates()->rows()>0){
    	_beta_hat(lambdaS_index, lambdaT_index) = regression_.getBeta()(0,0);
  	}
	_fn_hat(lambdaS_index, lambdaT_index) = (*Psi) *_solution(lambdaS_index,lambdaT_index).topRows(Psi->cols());
  }

}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
bool FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::stopping_criterion(const UInt& lambdaS_index, const UInt& lambdaT_index){
  // return true if the f-PIRLS has to perform another iteration, false if it has to be stopped

  bool do_stop_by_iteration = false;  // Do I need to stop becouse n_it > n_max?
  bool do_stop_by_treshold = false; // Do I need to stop becouse |J(k) - J(k+1)| < treshold?

  if(n_iterations[lambdaS_index][lambdaT_index] > inputData_.get_maxiter()){
    do_stop_by_iteration = true;
  }

  if(n_iterations[lambdaS_index][lambdaT_index] > 1){
    if(abs(past_J_values[lambdaS_index][lambdaT_index][0]+past_J_values[lambdaS_index][lambdaT_index][1] - current_J_values[lambdaS_index][lambdaT_index][0] - current_J_values[lambdaS_index][lambdaT_index][1]) < inputData_.get_treshold()){
        do_stop_by_treshold = true;
   }
  }

  return !(do_stop_by_iteration || do_stop_by_treshold );
}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
std::array<Real,2> FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_J(const UInt& lambdaS_index, const UInt& lambdaT_index){
	// compute the functional J: it is divided in parametric and non parametric part

	// Computation of parametric part: it is done in a pure virtual function since it is model-dependent
	Real parametric_value = compute_J_parametric(lambdaS_index, lambdaT_index);

	// Computation of nonparametric part: it is done here since it is PDE dependent
	Real non_parametric_value = 0;

	VectorXr Lf;

	Lf.resize(this->_solution(lambdaS_index,lambdaT_index).size()/2);
	for(UInt i=0; i< Lf.size(); i++){
	Lf(i) = this->_solution(lambdaS_index, lambdaT_index)(Lf.size() + i);
	}

	if(this->isSpaceVarying)
	{

	  if( this->inputData_.isSpaceTime() ){
	 
	   UInt M_ =  this->regression_.getM_();
	   UInt N_ = this->regression_.getN_();
	  
	   VectorXr forcingTerm_correction_;
	   forcingTerm_correction_.resize(M_*N_); // New size
	  
	   for(UInt i=0; i<N_; i++) // Update forcing term (i.e. repeat identically M_ times)
	   {
	   	for(UInt j=0; j<M_; j++)
		{
		 forcingTerm_correction_(i+j*N_) = this->forcingTerm(i);
		}
	   }
	  
	  	Lf = Lf - forcingTerm_correction_;
	  }
	  else{
	  	Lf = Lf - this->forcingTerm;
	  }
	}

	SpMat Int;
	std::array<Real, 2> lambdaST = {
	(*this->optimizationData_.get_LambdaS_vector())[lambdaS_index],
		(*this->optimizationData_.get_LambdaT_vector())[lambdaT_index]};
	if (this->inputData_.isSpaceTime()) {
	UInt correction_size = this->inputData_.getFlagParabolic() ? -1 : +2;
		VectorXr intcoef(this->mesh_time_.size() + correction_size);
		intcoef.setConstant(this->mesh_time_[1] - this->mesh_time_[0]);
		intcoef(0) *= 0.5;
		SpMat IN(this->mesh_.num_nodes(), this->mesh_.num_nodes());
		IN.setIdentity();
		SpMat tmp = intcoef.asDiagonal().toDenseMatrix().sparseView();
		tmp = kroneckerProduct(tmp, IN);
		Int.resize(tmp.rows(), tmp.cols());
		Int = lambdaST[0] * (*this->regression_.getR0_()) * tmp;
	} else {
		Int.resize(this->mesh_.num_nodes(), this->mesh_.num_nodes());
		Int = lambdaST[0] * (*this->regression_.getR0_());
	}

	non_parametric_value = Lf.transpose() * Int * Lf;

	std::array<Real,2> returnObject{parametric_value, non_parametric_value};

	return returnObject;
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_GCV(const UInt& lambdaS_index, const UInt& lambdaT_index){

        return;

}


/*********** FPIRLS apply template specialization ************/

// Laplace or Elliptic case
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<InputHandler,ORDER, mydim, ndim>::apply(){

  FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::apply(ForcingTerm());

}

// SpaceVarying case
template <UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<GAMDataEllipticSpaceVarying,ORDER, mydim, ndim>::apply(){

  this->isSpaceVarying = true;
  FPIRLS_Base<GAMDataEllipticSpaceVarying,ORDER, mydim, ndim>::apply(this->inputData_.getU());

}

/*********** FPIRLS_GAM Methods ************/

// Constructor
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::FPIRLS_GAM(const MeshHandler<ORDER,mydim,ndim> & mesh, InputHandler & inputData, OptimizationData & optimizationData,  VectorXr mu0, bool scale_parameter_flag, Real scale_param):
  	  FPIRLS<InputHandler,ORDER, mydim, ndim>(mesh, inputData, optimizationData), scale_parameter_flag_(scale_parameter_flag), _scale_param(scale_param)
{
  //Pre-allocate memory for all quatities
  WeightsMatrix_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
  mu_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
  pseudoObservations_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
  G_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
  
  //initialization of mu, current_J_values and past_J_values.
  for(UInt i=0; i<this->optimizationData_.get_size_S() ; i++){
   for(UInt j=0; j<this->optimizationData_.get_size_T() ; j++){
    mu_[i][j] = mu0;
   }
  }
};

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::FPIRLS_GAM(const MeshHandler<ORDER,mydim,ndim> & mesh, const std::vector<Real>& mesh_time, InputHandler & inputData, OptimizationData & optimizationData,  VectorXr mu0, bool scale_parameter_flag, Real scale_param):
	FPIRLS<InputHandler,ORDER, mydim, ndim>(mesh, mesh_time, inputData, optimizationData), scale_parameter_flag_(scale_parameter_flag), _scale_param(scale_param)
{
  //Pre-allocate memory for all quatities
  WeightsMatrix_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
  mu_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
  pseudoObservations_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
  G_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
  
  //initialization of mu, current_J_values and past_J_values.
  for(UInt i=0; i<this->optimizationData_.get_size_S() ; i++){
   for(UInt j=0; j<this->optimizationData_.get_size_T() ; j++){
    mu_[i][j] = mu0;
   }
  }
};


// FPIRLS_GAM methods

// STEP (1) methods

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::compute_G(const UInt& lambdaS_index, const UInt& lambdaT_index){
  // compute the G matrix as G_ii = diag( g'(mu_i))

  G_[lambdaS_index][lambdaT_index].resize(mu_[lambdaS_index][lambdaT_index].size());

  for(UInt i = 0; i<mu_[lambdaS_index][lambdaT_index].size(); i++){
    G_[lambdaS_index][lambdaT_index](i) = link_deriv(mu_[lambdaS_index][lambdaT_index](i));
  }

}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::compute_pseudoObs(const UInt& lambdaS_index, const UInt& lambdaT_index){
  // compute pseudodata observations

  VectorXr first_addendum; // G_ii( z_i - mu_i)
  VectorXr g_mu; // g( mu_i )

  const VectorXr * z = this->inputData_.getInitialObservations();

  first_addendum.resize(mu_[lambdaS_index][lambdaT_index].size());
  g_mu.resize(mu_[lambdaS_index][lambdaT_index].size());

  //compute the vecotr first_addendum and g_mu
  for(auto i=0; i < mu_[lambdaS_index][lambdaT_index].size(); i++){
    g_mu(i) = link(mu_[lambdaS_index][lambdaT_index](i));
    first_addendum(i) = G_[lambdaS_index][lambdaT_index](i)*((*z)(i)-mu_[lambdaS_index][lambdaT_index](i));
  }

  pseudoObservations_[lambdaS_index][lambdaT_index] = first_addendum + g_mu;

}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::compute_Weights(const UInt& lambdaS_index, const UInt& lambdaT_index){
  // computed W elementwise (it is a diagonal matrix)

  this->WeightsMatrix_[lambdaS_index][lambdaT_index].resize( mu_[lambdaS_index][lambdaT_index].size());

  for(auto i=0; i < mu_[lambdaS_index][lambdaT_index].size(); i++){
    this->WeightsMatrix_[lambdaS_index][lambdaT_index](i) = 1/(pow(G_[lambdaS_index][lambdaT_index](i),2)*(var_function( mu_[lambdaS_index][lambdaT_index](i))));
  }

}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::prepare_weighted_regression(const UInt& lambdaS_index, const UInt& lambdaT_index)
{
	// Compute weights and pseudo data
	compute_G(lambdaS_index, lambdaT_index);
	compute_pseudoObs(lambdaS_index, lambdaT_index);
	compute_Weights(lambdaS_index, lambdaT_index);
	
	// Store weights in RegressionData (it makes them visible to the solver for step (2) of f-PIRLS)
	this->inputData_.updatePseudodata(pseudoObservations_[lambdaS_index][lambdaT_index], this->WeightsMatrix_[lambdaS_index][lambdaT_index]);
	
}

// STEP (3) methods

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::compute_mu(const UInt& lambdaS_index, const UInt& lambdaT_index){
  //compute mu as mu_i = g-1( w_ii*beta + fn_hat)

  VectorXr W_beta = VectorXr::Zero(mu_[lambdaS_index][lambdaT_index].size()); // initialize the vector w_ii*beta

  if(this->inputData_.getCovariates()->rows()>0)
    W_beta = (*(this->inputData_.getCovariates()))*this->_beta_hat(lambdaS_index,lambdaT_index);

  for(UInt j=0; j < W_beta.size(); j++){
      mu_[lambdaS_index][lambdaT_index](j) = inv_link(W_beta[j] + this->_fn_hat(lambdaS_index,lambdaT_index)(j));
  }

}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::update_parameters(const UInt& lambdaS_index, const UInt& lambdaT_index)
{
	// Compute the vector mu
	compute_mu(lambdaS_index, lambdaT_index);
}

// Other methods

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
Real FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::compute_J_parametric(const UInt& lambdaS_index, const UInt& lambdaT_index){
  
  Real parametric_value = 0;
  Real tmp;

  const VectorXr * z = this->inputData_.getInitialObservations();

  for(UInt i=0; i < mu_.size(); i++){
    tmp =std::sqrt( var_function( mu_[lambdaS_index][lambdaT_index](i)) ) * ((*z)(i) - mu_[lambdaS_index][lambdaT_index](i)) ;
    parametric_value += tmp*tmp;
  }

  return parametric_value;
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::compute_GCV(const UInt& lambdaS_index, const UInt& lambdaT_index){

        if (this->optimizationData_.get_DOF_evaluation() != "not_required") //in this case surely we have already the ls
        { // is DOF_matrix to be computed?
        this->regression_.computeDegreesOfFreedom(0, 0, (*this->optimizationData_.get_LambdaS_vector())[lambdaS_index],
        					  							(*this->optimizationData_.get_LambdaT_vector())[lambdaT_index]);
        this->_dof(lambdaS_index, lambdaT_index) = this->regression_.getDOF()(0,0);
        }
        else this->_dof(lambdaS_index, lambdaT_index) = this->regression_.getDOF()(lambdaS_index, lambdaT_index);
 
        const VectorXr * y = this->inputData_.getInitialObservations();
        Real GCV_value = 0;

        for(UInt j=0; j < y->size();j++)
        	GCV_value += dev_function(mu_[lambdaS_index][lambdaT_index][j], (*y)[j]); //norm computation

        GCV_value *= y->size();

        GCV_value /= (y->size()-this->optimizationData_.get_tuning()*this->_dof(lambdaS_index,lambdaT_index))*(y->size()-this->optimizationData_.get_tuning()*this->_dof(lambdaS_index,lambdaT_index));

        this->_GCV[lambdaS_index][lambdaT_index] = GCV_value;

        //best lambda
        if(GCV_value < this->optimizationData_.get_best_value())
        {
        this->optimizationData_.set_best_lambda_S(lambdaS_index);
        this->optimizationData_.set_best_lambda_T(lambdaT_index);
        this->optimizationData_.set_best_value(GCV_value);
        }

}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::compute_variance_est(){
  Real phi;
  if(this->scale_parameter_flag_ && this->optimizationData_.get_loss_function()!="GCV"){// if scale param should be
    _variance_estimates.resize(this->lenS_, std::vector<Real>(this->lenT_, 0.0));
    const UInt n_obs = this->inputData_.getObservations()->size();

    //scale parameter computed as: mean((var.link(mu)*phi)/mu), and phi is computed as in Wood IGAM
    for(UInt i=0; i < this->lenS_;i++){
    	for(UInt j=0; j< this->lenT_; j++){ 
		phi = (this->scale_parameter_flag_ )? this->current_J_values[i][j][0]/(n_obs - this->_dof(i,j) ) : _scale_param;
		for(UInt k=0; k < this->mu_[i][j].size(); k++){
        	_variance_estimates[i][j] += phi* this->var_function(this->mu_[i][j](k))/this->mu_[i][j](k);
      }
      _variance_estimates[i][j] /= this->mu_[i][j].size();
    }
   }
  }else{
    _variance_estimates.resize(this->lenS_, std::vector<Real>(this->lenT_,-1));
  }
}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_GAM<InputHandler,ORDER, mydim, ndim>::additional_estimates()
{
	// Compute an estimate of the scale parameter when needed (only for proper distributions)
	compute_variance_est();
}

/*********** FPIRLS_MixedEffects Methods ************/

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::initialize_matrices()
{	
	// Compute and store Z_ and ZTZ_ for each group
	UInt starting_i = 0;
	for(UInt i=0; i<n_groups_; i++){
		// Store Z_ of group i
		Z_[i] = this->inputData_.getRandomEfectsCovariates()->block(starting_i, 0, group_sizes_[i], q_);
		starting_i += group_sizes_[i];
		
		// Compute and store Z_TZ_ of group i
		ZTZ_[i] = Z_[i].transpose() * Z_[i];
	}
	
	// Compute the initial guess of D_
	VectorXr D0(q_);
	for(UInt k=0; k<q_; k++){
		for(UInt i=0; i<n_groups_; i++){
			for(UInt j=0; j<group_sizes_[i]; j++){
				D0(k) += Z_[i](j,k) * Z_[i](j,k);
			}
		}
		D0(k) = std::sqrt( D0(k)/n_groups_ ) * 3 / 8 ;
	}	
	
	// Store it as the current D_ for each lambda to check
	for(UInt i=0; i<this->optimizationData_.get_size_S() ; i++){
		for(UInt j=0; j<this->optimizationData_.get_size_T() ; j++){
			D_[i][j] = D0;
		}
	}

	// Resize weights matrix
	WeightsMatrix_.resize(this->lenS_, std::vector<std::vector<MatrixXr>>(this->lenT_));
	for(UInt i=0; i<this->lenS_; i++){
		for(UInt j=0; j<this->lenT_; j++){
			//WeightsMatrix_[i][j].resize(n_groups_);
			for(UInt k=0; k<n_groups_; k++){
				WeightsMatrix_[i][j].push_back(MatrixXr(group_sizes_[k],group_sizes_[k]));
			}
		}
	}

	// Resize b_hat
	b_hat_.resize(this->lenS_, std::vector<std::vector<VectorXr>>(this->lenT_));
	for(UInt i=0; i<this->lenS_; i++){
		for(UInt j=0; j<this->lenT_; j++){
			b_hat_[i][j].resize(n_groups_, VectorXr(q_));
		}
	}
}

// Constructor
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::FPIRLS_MixedEffects(const MeshHandler<ORDER,mydim,ndim> & mesh, 
													InputHandler & inputData, OptimizationData & optimizationData):
		FPIRLS<InputHandler,ORDER, mydim, ndim>(mesh, inputData, optimizationData), 
		group_sizes_(*inputData.getGroupSizes()), n_groups_(inputData.getGroupNumber()), q_(inputData.get_q())
{
	// Pre-allocate memory for all quatities
	Z_.resize(n_groups_);
	ZTZ_.resize(n_groups_, MatrixXr(q_, q_));
	ZtildeTZtilde_.resize(n_groups_);
	D_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
	Sigma_b_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
	
	// Construct matrices ZTZ and initialize D_ for each lambda
	initialize_matrices();
	
};

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::FPIRLS_MixedEffects(const MeshHandler<ORDER,mydim,ndim> & mesh, const std::vector<Real>& mesh_time, 
													InputHandler & inputData, OptimizationData & optimizationData):
		FPIRLS<InputHandler,ORDER, mydim, ndim>(mesh, mesh_time, inputData, optimizationData),
		group_sizes_(*inputData.getGroupSizes()), n_groups_(inputData.getGroupNumber()), q_(inputData.get_q())
{
	// Pre-allocate memory for all quatities
	Z_.resize(n_groups_);
	ZTZ_.resize(n_groups_, MatrixXr(q_, q_));
	ZtildeTZtilde_.resize(n_groups_);
	D_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
	Sigma_b_.resize(this->lenS_, std::vector<VectorXr>(this->lenT_));
	
	// Construct matrices ZTZ and initialize D_ for each lambda
	initialize_matrices();
};


// FPIRLS_MixedEffects methods

// STEP (1) methods

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::compute_ZtildeTZtilde(const UInt& lambdaS_index, const UInt& lambdaT_index){

	for(auto i=0; i<n_groups_; i++){
		
		// Foe each group, consider Z^T Z
		MatrixXr ZtildeTZtilde_temp = ZTZ_[i];
		
		// Then add the diagonal elements stored in D_
		for(auto k=0; k<q_; k++){
			ZtildeTZtilde_temp(k,k) += D_[lambdaS_index][lambdaT_index](k) * D_[lambdaS_index][lambdaT_index](k);
		}	
		
		ZtildeTZtilde_[i].compute(ZtildeTZtilde_temp);
		
	}
}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::compute_Weights(const UInt& lambdaS_index, const UInt& lambdaT_index){
	// computed W blocktwise (it is a block diagonal matrix)
	
	//UInt starting_i = 0;
	for(auto i=0; i < n_groups_; i++){
		// Compute the current block with the Woodbury identity
		this->WeightsMatrix_[lambdaS_index][lambdaT_index][i] = - Z_[i] * ZtildeTZtilde_[i].solve( Z_[i].transpose() );
		
		// Add the identity matrix (1 on the diagonal)
		for(auto k=0; k < group_sizes_[i]; k++){
			this->WeightsMatrix_[lambdaS_index][lambdaT_index][i](k,k) =  1 + this->WeightsMatrix_[lambdaS_index][lambdaT_index][i](k, k);	
		}
	}

}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::prepare_weighted_regression(const UInt& lambdaS_index, const UInt& lambdaT_index)
{

	// Compute Ztilde^T Ztilde
	compute_ZtildeTZtilde(lambdaS_index, lambdaT_index);
	
	// Compute weights
	compute_Weights(lambdaS_index, lambdaT_index);
	
	// Store weights in RegressionData (it makes them visible to the solver for step (2) of f-PIRLS)
	this->inputData_.updateWeights(this->WeightsMatrix_[lambdaS_index][lambdaT_index]);
	
}



// STEP (3) methods

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::compute_bhat(const UInt& lambdaS_index, const UInt& lambdaT_index)
{

	VectorXr const * y = this->inputData_.getObservations();
	
	UInt starting_i = 0;
	for(auto i=0; i<n_groups_; i++){
		
		VectorXr y_i = y->segment(starting_i, group_sizes_[i]);
		VectorXr y_i_hat = this->_fn_hat(lambdaS_index,lambdaT_index).segment(starting_i, group_sizes_[i]);

		UInt p = this->inputData_.getCovariates()->cols();
		if( p > 0 ){
			MatrixXr X_i = this->inputData_.getCovariates()->block(starting_i, 0, group_sizes_[i], p);

			y_i_hat = y_i_hat + X_i*this->_beta_hat(lambdaS_index,lambdaT_index);
		}
		
		VectorXr res_i = y->segment(starting_i, group_sizes_[i]);
		res_i -= this->_fn_hat(lambdaS_index,lambdaT_index).segment(starting_i, group_sizes_[i]);
		if( p > 0 ){
			res_i -= this->inputData_.getCovariates()->block(starting_i, 0, group_sizes_[i], p)*this->_beta_hat(lambdaS_index,lambdaT_index);
		}

		b_hat_[lambdaS_index][lambdaT_index][i] = ZtildeTZtilde_[i].solve( Z_[i].transpose() * res_i );
		
		starting_i += group_sizes_[i];
	}
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::compute_sigma_sq_hat(const UInt& lambdaS_index, const UInt& lambdaT_index)
{

	this->regression_.computeDegreesOfFreedom(0, 0, (*this->optimizationData_.get_LambdaS_vector())[lambdaS_index],
							  (*this->optimizationData_.get_LambdaT_vector())[lambdaT_index]  );
	this->_dof(lambdaS_index, lambdaT_index) = this->regression_.getDOF()(0,0);

	sigma_sq_hat_ = 0;	

	const VectorXr * y = this->inputData_.getObservations();
	
	UInt starting_i = 0;
	for(auto i=0; i<n_groups_; i++){
		
		// log-likelihood of observations
		VectorXr res_i = y->segment(starting_i, group_sizes_[i]);
		
		//VectorXr f_n_hat_i = ;
		VectorXr f_n_hat_i = this->_fn_hat(lambdaS_index,lambdaT_index).segment(starting_i, group_sizes_[i]);
		res_i -= ( f_n_hat_i + Z_[i]*b_hat_[lambdaS_index][lambdaT_index][i] );

		UInt p = this->inputData_.getCovariates()->cols();
		if( p > 0 ){
			MatrixXr X_i = this->inputData_.getCovariates()->block(starting_i, 0, group_sizes_[i], p);
			res_i -= ( X_i*this->_beta_hat(lambdaS_index,lambdaT_index) );
		}
		sigma_sq_hat_ += res_i.dot(res_i);
		
		starting_i += group_sizes_[i];
	}
	
	sigma_sq_hat_ /= (y->size() - this->_dof(lambdaS_index, lambdaT_index));
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::build_LTL(const UInt& lambdaS_index, const UInt& lambdaT_index)
{

	MatrixXr LTL_temp = MatrixXr::Zero(q_, q_);
	
	for(auto i=0; i<n_groups_; i++){
	
		LTL_temp += b_hat_[lambdaS_index][lambdaT_index][i]*(b_hat_[lambdaS_index][lambdaT_index][i]).transpose() / sigma_sq_hat_;
		LTL_temp += ZtildeTZtilde_[i].solve(MatrixXr::Identity(q_, q_));
	}
	
	LTL_.compute(LTL_temp);
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::compute_A(const UInt& lambdaS_index, const UInt& lambdaT_index)
{

	MatrixXr A_temp = LTL_.matrixL();
	
	A_ = A_temp.triangularView<Eigen::Lower>().solve( MatrixXr::Identity(q_, q_) );
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::update_parameters(const UInt& lambdaS_index, const UInt& lambdaT_index)
{
	compute_bhat(lambdaS_index, lambdaT_index);
	
	compute_sigma_sq_hat(lambdaS_index, lambdaT_index);

	build_LTL(lambdaS_index, lambdaT_index);

	compute_A(lambdaS_index, lambdaT_index);

	for(auto i=0; i<q_; i++){
		D_[lambdaS_index][lambdaT_index](i) = A_(i,i) * std::sqrt(n_groups_);
	}
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
Real FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::compute_J_parametric(const UInt& lambdaS_index, const UInt& lambdaT_index) 
{
	Real parametric_value = 0;

	parametric_value += /*y->size() -*/ this->_dof(lambdaS_index, lambdaT_index);

	const VectorXr * y = this->inputData_.getObservations();

	parametric_value -= ( n_groups_*q_ - y->size() ) * std::log(sigma_sq_hat_);
	
	UInt starting_i = 0;
	for(auto i=0; i<n_groups_; i++){
		
		// log-likelihood of random effects	(completed outside the for cycle)
		VectorXr Db_i = D_[lambdaS_index][lambdaT_index].asDiagonal() * b_hat_[lambdaS_index][lambdaT_index][i];
		parametric_value -= sigma_sq_hat_ * ( Db_i ).dot( Db_i );
		
		starting_i += group_sizes_[i];
	}
	
	// Compute the determinant of D 	(NOTE: D is a diagonal matrix stored in a vector!)
	Real detD = 1;
	for(auto k=0; k<q_; k++){
		detD *= D_[lambdaS_index][lambdaT_index][k];
	}
	parametric_value += 2 * n_groups_ * std::log(detD);

	return parametric_value/2 ;
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::additional_estimates()
{
	for(UInt i=0; i < this->lenS_;i++){
    	for(UInt j=0; j< this->lenT_; j++){ 
			Sigma_b_[i][j] = D_[i][j];
			
			for(auto k=0; k<q_; k++){
				Sigma_b_[i][j](k) *= D_[i][j](k);
				Sigma_b_[i][j](k) = 1/Sigma_b_[i][j](k);
			}
		}
	}
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_MixedEffects<InputHandler,ORDER, mydim, ndim>::compute_GCV(const UInt& lambdaS_index, const UInt& lambdaT_index)
{
	
	const VectorXr * y = this->inputData_.getObservations();
	Real GCV_value = 0;
	
	UInt starting_i = 0;
	for(auto i=0; i<n_groups_; i++){
		
		VectorXr res_i = (*y).segment(starting_i, group_sizes_[i]);

		MatrixXr f_n_hat_i = this->_fn_hat(lambdaS_index,lambdaT_index).segment(starting_i, group_sizes_[i]);
		res_i -= ( f_n_hat_i + Z_[i]*b_hat_[lambdaS_index][lambdaT_index][i] );

		UInt p = this->inputData_.getCovariates()->cols();
		if( p > 0 ){
			MatrixXr X_i = this->inputData_.getCovariates()->block(starting_i, 0, group_sizes_[i], p);
			res_i = res_i - X_i*this->_beta_hat(lambdaS_index,lambdaT_index);
		}
		GCV_value += res_i.dot(res_i);
		
		starting_i += group_sizes_[i];
	}

	GCV_value *= y->size();

	GCV_value /= (y->size()-this->optimizationData_.get_tuning()*this->_dof(lambdaS_index,lambdaT_index))*(y->size()-this->optimizationData_.get_tuning()*this->_dof(lambdaS_index,lambdaT_index));

	this->_GCV[lambdaS_index][lambdaT_index] = GCV_value;

	//best lambda
	if(GCV_value < this->optimizationData_.get_best_value())
	{
		this->optimizationData_.set_best_lambda_S(lambdaS_index);
		this->optimizationData_.set_best_lambda_T(lambdaT_index);
		this->optimizationData_.set_best_value(GCV_value);
	}

}


#endif
