#ifndef __REGRESSION_DATA_IMP_H__
#define __REGRESSION_DATA_IMP_H__

// -- GAM CONSTRUCTORS and Utilities--
template<typename RegressionHandler>
void RegressionDataGAM<RegressionHandler>::initializeWeights()
{
	// Allocate required memory
	this->WeightsMatrix_.reserve(initialObservations_.size());
	this->WeightsMatrix_.resize(initialObservations_.size(), initialObservations_.size());

	// Initialize required elements
	for(auto i=0; i<initialObservations_.size(); i++){
		this->WeightsMatrix_.insert(i,i) = 1;
	}
}



template<typename RegressionHandler>
void RegressionDataGAM<RegressionHandler>::updatePseudodata(VectorXr& z_, VectorXr& P)
{
	this-> observations_ = z_; 
	for(auto i=0; i<P.size(); i++){
		this-> WeightsMatrix_.coeffRef(i,i) = P.coeff(i);
	}
}


// Laplace
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
	SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch,
	SEXP Rmax_num_iteration, SEXP Rthreshold):
	RegressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	initialObservations_ = this->observations_;
	initializeWeights();
	this->isGAM = true;
}

// PDE
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
	SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
	SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch, SEXP Rmax_num_iteration, SEXP Rthreshold):
	RegressionDataElliptic(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc,
		Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	initialObservations_ = this->observations_;
	initializeWeights();
	this->isGAM = true;
}

// PDE SpaceVarying
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
	SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
	SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch, SEXP Rmax_num_iteration, SEXP Rthreshold):
	RegressionDataEllipticSpaceVarying(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc, Ru,
		Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	initialObservations_ = this->observations_;
	initializeWeights();
	this->isGAM = true;
}

//Laplace time
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, 
	SEXP Rorder, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, 
	SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic,SEXP Rflag_iterative, 
	SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Ric, SEXP Rsearch, SEXP Rmax_num_iteration_pirls, SEXP Rthreshold_pirls):
	RegressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, Rcovariates, RBCIndices, RBCValues, 	
		RincidenceMatrix, RarealDataAvg, 
		Rflag_mass, Rflag_parabolic, Rflag_iterative, Rmax_num_iteration, Rthreshold, Ric, Rsearch) {
    max_num_iterations_ = INTEGER(Rmax_num_iteration_pirls)[0];
    threshold_ = REAL(Rthreshold_pirls)[0];
    initialObservations_ = this->observations_;
	initializeWeights();
    this->isGAM = true;
}

// PDE time
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, 
	SEXP Rorder, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix,
    	SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic,SEXP Rflag_iterative, SEXP Rmax_num_iteration, 
    	SEXP Rthreshold, SEXP Ric, SEXP Rsearch, SEXP Rmax_num_iteration_pirls, SEXP Rthreshold_pirls):
    		RegressionDataElliptic(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, RK, Rbeta, Rc, 
    				       Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, 
    				       Rflag_mass, Rflag_parabolic, Rflag_iterative, Rmax_num_iteration, Rthreshold, Ric, Rsearch) {
    max_num_iterations_ = INTEGER(Rmax_num_iteration_pirls)[0];
    threshold_ = REAL(Rthreshold_pirls)[0];
    initialObservations_ = this->observations_;
	initializeWeights();
    this->isGAM = true;
}

// PDE SpaceVarying time
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, 
	SEXP Rorder, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix,
	SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic,SEXP Rflag_iterative, SEXP Rmax_num_iteration, 
	SEXP Rthreshold, SEXP Ric, SEXP Rsearch, SEXP Rmax_num_iteration_pirls, SEXP Rthreshold_pirls): 
		RegressionDataEllipticSpaceVarying(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, RK, Rbeta, Rc, Ru, 
						   Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, 
						   Rflag_mass, Rflag_parabolic, Rflag_iterative, Rmax_num_iteration_pirls, Rthreshold_pirls, Ric,Rsearch) {
    max_num_iterations_ = INTEGER(Rmax_num_iteration_pirls)[0];
    threshold_ = REAL(Rthreshold_pirls)[0];
    initialObservations_ = this->observations_;
	initializeWeights();
    this->isGAM = true;
}




// -- MixedEffects CONSTRUCTORS and Utilities --
template<typename RegressionHandler>
void RegressionDataMixedEffects<RegressionHandler>::initializeWeights()
{
	// Allocate required memory
	UInt num_elems = 0;
	for(auto k=0; k<n_groups_; k++){	
		num_elems += group_sizes_[k]*group_sizes_[k];
	}
	this->WeightsMatrix_.reserve(num_elems);
	this->WeightsMatrix_.resize(this->n_, this->n_);

	// Initialize required elements
	UInt starting_ind = 0;
	for(auto k=0; k<n_groups_; k++){
		for(auto i=0; i<group_sizes_[k]; i++){
			for(auto j=0; j<group_sizes_[k]; j++){
				this->WeightsMatrix_.insert(starting_ind+i,starting_ind+j) = 1;
			}
		}
	starting_ind += group_sizes_[k];
	}
}



template<typename RegressionHandler>
void RegressionDataMixedEffects<RegressionHandler>::updateWeights(std::vector<MatrixXr> P)
{
	UInt starting_ind = 0;
	for(auto k=0; k<n_groups_; k++){
		for(auto i=0; i<group_sizes_[k]; i++){
			for(auto j=0; j<group_sizes_[k]; j++){
				this->WeightsMatrix_.coeffRef(starting_ind+i,starting_ind+j) = P[k].coeff(i,j);
			}
		}
	starting_ind += group_sizes_[k];
	}
}



template<typename RegressionHandler>
void RegressionDataMixedEffects<RegressionHandler>::setGroupSizes(SEXP Rgroup_sizes)
{
	group_sizes_.resize(n_groups_, 0);

	for(auto i=0;i<n_groups_;++i)
	{
		group_sizes_[i] = INTEGER(Rgroup_sizes)[i];	
	}
}


template<typename RegressionHandler>
void RegressionDataMixedEffects<RegressionHandler>::setRandomEffectsCovariates(SEXP Rrandom_effects_covariates)
{
	q_ = INTEGER(Rf_getAttrib(Rrandom_effects_covariates, R_DimSymbol))[1];
	UInt k=0;

	random_effects_covariates_.resize(this->n_, q_);

	const std::vector<UInt> observations_na = *this->getObservationsNA();

	for(auto i=0; i<this->n_; ++i)
	{
		for(auto j=0; j<q_ ; ++j)
		{
			if(observations_na.size()>k && i==observations_na[k])
			{
				random_effects_covariates_(i,j)=0;
				k++;
			}
			else
			{
				random_effects_covariates_(i,j)=REAL(Rrandom_effects_covariates)[i+ this->n_*j];
			}
		}
	}
}


// Laplace
template<typename RegressionHandler>
RegressionDataMixedEffects<RegressionHandler>::RegressionDataMixedEffects(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
	SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch,
	SEXP Rmax_num_iteration, SEXP Rthreshold, 
	SEXP Rrandom_effects_covariates, SEXP Rgroup_sizes, SEXP Rn_groups):
	RegressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	n_groups_ = INTEGER(Rn_groups)[0];
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	setGroupSizes(Rgroup_sizes);
	setRandomEffectsCovariates(Rrandom_effects_covariates);
	initializeWeights();
	this->isGAM = true;
}

// PDE
template<typename RegressionHandler>
RegressionDataMixedEffects<RegressionHandler>::RegressionDataMixedEffects(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
	SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
	SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch, 
	SEXP Rmax_num_iteration, SEXP Rthreshold, 
	SEXP Rrandom_effects_covariates, SEXP Rgroup_sizes, SEXP Rn_groups):
	RegressionDataElliptic(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc,
		Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	n_groups_ = INTEGER(Rn_groups)[0];
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	setGroupSizes(Rgroup_sizes);
	setRandomEffectsCovariates(Rrandom_effects_covariates);
	initializeWeights();
	this->isGAM = true;
}

// PDE SpaceVarying
template<typename RegressionHandler>
RegressionDataMixedEffects<RegressionHandler>::RegressionDataMixedEffects(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
	SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
	SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch, 
	SEXP Rmax_num_iteration, SEXP Rthreshold, 
	SEXP Rrandom_effects_covariates, SEXP Rgroup_sizes, SEXP Rn_groups):
	RegressionDataEllipticSpaceVarying(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc, Ru,
		Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	n_groups_ = INTEGER(Rn_groups)[0];
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	setGroupSizes(Rgroup_sizes);
	setRandomEffectsCovariates(Rrandom_effects_covariates);
	initializeWeights();
	this->isGAM = true;
}

//Laplace time
template <typename RegressionHandler>
RegressionDataMixedEffects<RegressionHandler>::RegressionDataMixedEffects(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, 
	SEXP Rorder, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, 
	SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic,SEXP Rflag_iterative, 
	SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Ric, SEXP Rsearch, 
	SEXP Rmax_num_iteration_pirls, SEXP Rthreshold_pirls, 
	SEXP Rrandom_effects_covariates, SEXP Rgroup_sizes, SEXP Rn_groups):
	RegressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, Rcovariates, RBCIndices, RBCValues, 	
		RincidenceMatrix, RarealDataAvg, 
		Rflag_mass, Rflag_parabolic, Rflag_iterative, Rmax_num_iteration, Rthreshold, Ric, Rsearch)
{
	n_groups_ = INTEGER(Rn_groups)[0];
	max_num_iterations_ = INTEGER(Rmax_num_iteration_pirls)[0];
	threshold_ = REAL(Rthreshold_pirls)[0];
	setGroupSizes(Rgroup_sizes);
	setRandomEffectsCovariates(Rrandom_effects_covariates);
	initializeWeights();
	this->isGAM = true;
}

// PDE time
template <typename RegressionHandler>
RegressionDataMixedEffects<RegressionHandler>::RegressionDataMixedEffects(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, 
	SEXP Rorder, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix,
    	SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic,SEXP Rflag_iterative, SEXP Rmax_num_iteration, 
    	SEXP Rthreshold, SEXP Ric, SEXP Rsearch, 
		SEXP Rmax_num_iteration_pirls, SEXP Rthreshold_pirls, 
		SEXP Rrandom_effects_covariates, SEXP Rgroup_sizes, SEXP Rn_groups):
    		RegressionDataElliptic(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, RK, Rbeta, Rc, 
    				       Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, 
    				       Rflag_mass, Rflag_parabolic, Rflag_iterative, Rmax_num_iteration, Rthreshold, Ric, Rsearch)
{
	n_groups_ = INTEGER(Rn_groups)[0];
	max_num_iterations_ = INTEGER(Rmax_num_iteration_pirls)[0];
	threshold_ = REAL(Rthreshold_pirls)[0];
	setGroupSizes(Rgroup_sizes);
	setRandomEffectsCovariates(Rrandom_effects_covariates);
	initializeWeights();
	this->isGAM = true;
}

// PDE SpaceVarying time
template <typename RegressionHandler>
RegressionDataMixedEffects<RegressionHandler>::RegressionDataMixedEffects(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, 
	SEXP Rorder, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix,
	SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic,SEXP Rflag_iterative, SEXP Rmax_num_iteration, 
	SEXP Rthreshold, SEXP Ric, SEXP Rsearch, 
	SEXP Rmax_num_iteration_pirls, SEXP Rthreshold_pirls, 
	SEXP Rrandom_effects_covariates, SEXP Rgroup_sizes, SEXP Rn_groups): 
		RegressionDataEllipticSpaceVarying(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, RK, Rbeta, Rc, Ru, 
						   Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, 
						   Rflag_mass, Rflag_parabolic, Rflag_iterative, Rmax_num_iteration_pirls, Rthreshold_pirls, Ric,Rsearch)
{
	n_groups_ = INTEGER(Rn_groups)[0];
	max_num_iterations_ = INTEGER(Rmax_num_iteration_pirls)[0];
	threshold_ = REAL(Rthreshold_pirls)[0];
	setGroupSizes(Rgroup_sizes);
	setRandomEffectsCovariates(Rrandom_effects_covariates);
	initializeWeights();
	this->isGAM = true;
	}

#endif
