CPP_smooth.MixedEffects.FEM.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                    covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL,
                                    areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, 
                                    max.steps.FPIRLS = 15, threshold.FPIRLS = 0.0002020,
                                    search, bary.locations, optim, lambdaS = NULL, lambdaT = NULL,
                                    DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                    DOF.matrix = NULL, GCV.inflation.factor = 1,
                                    lambda.optimization.tolerance = 0.05,
                                    rand.effects.covariates, group.ids, n.groups) 
{
  # Indexes in C++ starts from 0, in R from 1
  FEMbasis$mesh$triangles <- FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges <- FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] <-
    FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  group.ids = group.ids - 1

  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }

  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }

  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }

  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }

  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  FLAG_ITERATIVE<-as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE)<-"integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  rand.effects.covariates <- as.matrix(rand.effects.covariates)
  storage.mode(rand.effects.covariates) <-"double"
  group.ids <- as.matrix(group.ids)
  storage.mode(group.ids) <-"integer"
  storage.mode(n.groups) <- "integer"

  ## IC estimation for parabolic smoothing from the first column of observations
  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]

    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }

    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[notNAIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    # adjust Random Effects input for parabolic smoothing
    if (nrow(rand.effects.covariates) == 0) {
      rand.effects.covariatesIC <- rand.effects.covariates
      n.groupsIC <- n.groups
      group.idsIC <- group.ids
    } else {
      rand.effects.covariatesIC <- rand.effects.covariates[notNAIC, ]
      rand.effects.covariatesIC <- as.matrix(rand.effects.covariatesIC)
      storage.mode(rand.effects.covariatesIC) <- "double"
      group_idsIC = as.factor(group.ids[notNAIC,])
      group_idsIC_tab = table(group_idsIC)
      n.groupsIC = length(group_idsIC_tab)
      group.idsIC = as.integer(group_idsIC)
      group.idsIC = group.idsIC- 1
      group.idsIC = as.matrix(group.idsIC)
      storage.mode(group.idsIC) <-"integer"
      storage.mode(n.groupsIC) <- "integer"
    }
    
    

    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC

    ICsol <- .Call(
      "MixedEffects_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg,
      max.steps.FPIRLS, threshold.FPIRLS, search,
      optim, lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      rand.effects.covariatesIC, group.idsIC, n.groupsIC, 
      PACKAGE = "fdaPDE"
    )
    if (nrow(covariates) != 0) {
      betaIC = matrix(data=ICsol[[5]],nrow=ncol(covariatesIC),ncol=length(lambdaSIC))
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
    
    # fix Mixed Effects with parabolic smoothing
    if (nrow(rand.effects.covariates) != 0) {
      rand.effects.covariates <- rand.effects.covariates[(NobsIC + 1):nrow(rand.effects.covariates), ]
      rand.effects.covariates <- as.matrix(rand.effects.covariates)
    }
    group_ids = as.factor(group.ids[(NobsIC + 1):length(group.ids)])
    group_ids_tab = table(group_ids)
    n.groups = length(group_ids_tab)
    group.idsIC = as.integer(group_idsIC)
    group.idsIC = group.idsIC- 1
    group.idsIC = as.matrix(group.idsIC)
    storage.mode(group.idsIC) <-"integer"
    storage.mode(n.groupsIC) <- "integer"
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"

  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
    each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"

  ## Call C++ function
  bigsol <- .Call(
    "MixedEffects_Laplace_time", locations, bary.locations, time_locations,
    observations, FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim,
    covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
    FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, max.steps.FPIRLS, threshold.FPIRLS,
    search, optim, lambdaS, lambdaT,
    DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix,
    GCV.inflation.factor, lambda.optimization.tolerance,
    rand.effects.covariates, group.ids, n.groups, 
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}

CPP_smooth.MixedEffects.FEM.PDE.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                        covariates = NULL, PDE_parameters, ndim, mydim, BC = NULL,
                                        incidence_matrix = NULL, areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, 
                                        max.steps, threshold, IC, max.steps.FPIRLS = 15,
                                        threshold.FPIRLS = 0.0002020, search, bary.locations, optim, lambdaS = NULL,
                                        lambdaT = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                        DOF.matrix = NULL, GCV.inflation.factor = 1,
                                        lambda.optimization.tolerance = 0.05,
                                        rand.effects.covariates, group.ids, n.groups) 
{

  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation

  FEMbasis$mesh$triangles <- FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges <- FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] <-
    FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  group.ids = group.ids - 1

  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }


  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }

  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }

  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }

  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }

  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }

  if (is.null(time_locations)) {
    time_locations <- matrix(ncol = 0, nrow = 0)
  } else {
    time_locations <- as.matrix(time_locations)
  }

  if (is.null(time_mesh)) {
    time_mesh <- matrix(ncol = 0, nrow = 0)
  } else {
    time_mesh <- as.matrix(time_mesh)
  }


  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(PDE_parameters$K) <- "double"
  storage.mode(PDE_parameters$b) <- "double"
  storage.mode(PDE_parameters$c) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  FLAG_ITERATIVE<-as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE)<-"integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  rand.effects.covariates <- as.matrix(rand.effects.covariates)
  storage.mode(rand.effects.covariates) <-"double"
  group.ids <- as.matrix(group.ids)
  storage.mode(group.ids) <-"integer"
  storage.mode(n.groups) <- "integer"

  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]

    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }

    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[1:NobsIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    # adjust Random Effects input for parabolic smoothing
    if (nrow(rand.effects.covariates) == 0) {
      rand.effects.covariatesIC <- rand.effects.covariates
      n.groupsIC <- n.groups
      group.idsIC <- group.ids
    } else {
      rand.effects.covariatesIC <- rand.effects.covariates[1:NobsIC, ]
      rand.effects.covariatesIC <- as.matrix(rand.effects.covariatesIC)
      storage.mode(rand.effects.covariatesIC) <- "double"
      
      group_idsIC = as.factor(group.ids[1:NobsIC, ])
      group_idsIC_tab = table(group_idsIC)
      n.groupsIC = length(group_idsIC_tab)
      group.idsIC = as.integer(group_idsIC)
      group.idsIC = group.idsIC- 1
      group.idsIC = as.matrix(group.idsIC)
      storage.mode(group.idsIC) <-"integer"
      storage.mode(n.groupsIC) <- "integer"
    }

    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ICsol <- .Call(
      "MixedEffects_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg,
      max.steps.FPIRLS, threshold.FPIRLS, search,
      as.integer(c(0, 1, 1)), lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      rand.effects.covariatesIC, group.idsIC, n.groupsIC, 
      PACKAGE = "fdaPDE"
    )


    if (nrow(covariates) != 0) {
      betaIC <- ICsol[[15]]
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
    
    # fix Mixed Effects with parabolic smoothing
    if (nrow(rand.effects.covariates) != 0) {
      rand.effects.covariates <- rand.effects.covariates[(NobsIC + 1):nrow(rand.effects.covariates), ]
      rand.effects.covariates <- as.matrix(rand.effects.covariates)
    }
    group_ids = as.factor(group.ids[(NobsIC + 1):length(group.ids)])
    group_ids_tab = table(group_ids)
    n.groups = length(group_ids_tab)
    group.idsIC = as.integer(group_idsIC)
    group.idsIC = group.idsIC- 1
    group.idsIC = as.matrix(group.idsIC)
    storage.mode(group.idsIC) <-"integer"
    storage.mode(n.groupsIC) <- "integer"

    }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"

  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
    each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"

  ## Call C++ function
  bigsol <- .Call(
    "MixedEffects_PDE_time", locations, bary.locations, time_locations, observations,
    FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim, PDE_parameters$K,
    PDE_parameters$b, PDE_parameters$c, covariates, BC$BC_indices,
    BC$BC_values, incidence_matrix, areal.data.avg, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold,
    IC, max.steps.FPIRLS, threshold.FPIRLS,search,
    optim, lambdaS, lambdaT, DOF.stochastic.realizations, DOF.stochastic.seed,
    DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance,
    rand.effects.covariates, group.ids, n.groups, 
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}

CPP_smooth.MixedEffects.FEM.PDE.sv.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                           covariates = NULL, PDE_parameters, ndim, mydim, BC = NULL,
                                           incidence_matrix = NULL, areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, 
                                           max.steps, threshold, IC, max.steps.FPIRLS = 15,
                                           threshold.FPIRLS = 0.0002020, search, bary.locations, optim, lambdaS = NULL,
                                           lambdaT = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                           DOF.matrix = NULL, GCV.inflation.factor = 1,
                                           lambda.optimization.tolerance = 0.05,
                                           rand.effects.covariates, group.ids, n.groups) 
{

  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation
  FEMbasis$mesh$triangles <- FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges <- FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] <-
    FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  group.ids = group.ids - 1

  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }

  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }

  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }

  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }

  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }
  if(is.null(time_locations)){
    time_locations<-matrix(ncol=0, nrow=0)
  }else
  {
    time_locations <- as.matrix(time_locations)
  }

  if(is.null(time_mesh)){
    time_mesh<-matrix(ncol=0, nrow=0)
  }else
  {
    time_mesh <- as.matrix(time_mesh)
  }

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  rand.effects.covariates <- as.matrix(rand.effects.covariates)
  storage.mode(rand.effects.covariates) <-"double"
  group.ids <- as.matrix(group.ids)
  storage.mode(group.ids) <-"integer"
  storage.mode(n.groups) <- "integer"


  PDE_param_eval <- NULL
  points_eval <-
    matrix(CPP_get_evaluations_points(
      mesh = FEMbasis$mesh, order = FEMbasis$order
    ), ncol = 2)
  PDE_param_eval$K <- (PDE_parameters$K)(points_eval)
  PDE_param_eval$b <- (PDE_parameters$b)(points_eval)
  PDE_param_eval$c <- (PDE_parameters$c)(points_eval)
  PDE_param_eval$u <- (PDE_parameters$u)(points_eval)

  storage.mode(PDE_param_eval$K) <- "double"
  storage.mode(PDE_param_eval$b) <- "double"
  storage.mode(PDE_param_eval$c) <- "double"
  storage.mode(PDE_param_eval$u) <- "double"

  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]

    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }

    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[notNAIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    # adjust Random Effects input for parabolic smoothing
    if (nrow(rand.effects.covariates) == 0) {
      rand.effects.covariatesIC <- rand.effects.covariates
      n.groupsIC <- n.groups
      group.idsIC <- group.ids
    } else {
      rand.effects.covariatesIC <- rand.effects.covariates[notNAIC, ]
      rand.effects.covariatesIC <- as.matrix(rand.effects.covariatesIC)
      storage.mode(rand.effects.covariatesIC) <- "double"
      
      group_idsIC = as.factor(group.ids[notNAIC, ])
      group_idsIC_tab = table(group_idsIC)
      n.groupsIC = length(group_idsIC_tab)
      group.idsIC = as.integer(group_idsIC)
      group.idsIC = group.idsIC- 1
      group.idsIC = as.matrix(group.idsIC)
      storage.mode(group.idsIC) <-"integer"
      storage.mode(n.groupsIC) <- "integer"
    }

    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC

    ICsol <- .Call(
      "MixedEffects_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg,
      max.steps.FPIRLS, threshold.FPIRLS, search,
      as.integer(c(0, 1, 1)), lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      rand.effects.covariatesIC, group.idsIC, n.groupsIC,
      PACKAGE = "fdaPDE"
    )
    
  ## shifting the lambdas interval if the best lambda is the higher one and retry smoothing
    if (nrow(covariates) != 0) {
      betaIC = matrix(data=ICsol[[5]],nrow=ncol(covariatesIC),ncol=length(lambdaSIC))
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
    
    # fix Mixed Effects with parabolic smoothing
    if (nrow(rand.effects.covariates) != 0) {
      rand.effects.covariates <- rand.effects.covariates[(NobsIC + 1):nrow(rand.effects.covariates), ]
      rand.effects.covariates <- as.matrix(rand.effects.covariates)
    }
    group_ids = as.factor(group.ids[(NobsIC + 1):length(group.ids)])
    group_ids_tab = table(group_ids)
    n.groups = length(group_ids_tab)
    group.idsIC = as.integer(group_idsIC)
    group.idsIC = group.idsIC- 1
    group.idsIC = as.matrix(group.idsIC)
    storage.mode(group.idsIC) <-"integer"
    storage.mode(n.groupsIC) <- "integer"
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"

  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
    each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"


  ## Call C++ function
  bigsol <- .Call(
    "MixedEffects_PDE_space_varying_time", locations, bary.locations, time_locations,
    observations, FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim,
    PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u,
    covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
    FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, max.steps.FPIRLS, threshold.FPIRLS,
    search, optim, lambdaS, lambdaT,
    DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix,
    GCV.inflation.factor, lambda.optimization.tolerance,
    rand.effects.covariates, group.ids, n.groups, 
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}

CPP_smooth.manifold.MixedEffects.FEM.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                    covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL,
                                    areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC,
                                    max.steps.FPIRLS = 15, threshold.FPIRLS = 0.0002020,
                                    search, bary.locations, optim, lambdaS = NULL, lambdaT = NULL,
                                    DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                    DOF.matrix = NULL, GCV.inflation.factor = 1,
                                    lambda.optimization.tolerance = 0.05,
                                    rand.effects.covariates, group.ids, n.groups) 
{
  # Indexes in C++ starts from 0, in R from 1
  FEMbasis$mesh$triangles <- FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges <- FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] <-
    FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  group.ids = group.ids - 1
  
  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }
  
  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }
  
  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }
  
  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }
  
  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }
  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  FLAG_ITERATIVE<-as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE)<-"integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  rand.effects.covariates <- as.matrix(rand.effects.covariates)
  storage.mode(rand.effects.covariates) <-"double"
  group.ids <- as.matrix(group.ids)
  storage.mode(group.ids) <-"integer"
  storage.mode(n.groups) <- "integer"
  
  ## IC estimation for parabolic smoothing from the first column of observations
  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]
    
    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }
    
    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[notNAIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    # adjust Random Effects input for parabolic smoothing
    if (nrow(rand.effects.covariates) == 0) {
      rand.effects.covariatesIC <- rand.effects.covariates
      n.groupsIC <- n.groups
      group.idsIC <- group.ids
    } else {
      rand.effects.covariatesIC <- rand.effects.covariates[notNAIC, ]
      rand.effects.covariatesIC <- as.matrix(rand.effects.covariatesIC)
      storage.mode(rand.effects.covariatesIC) <- "double"
      
      group_idsIC = as.factor(group.ids[notNAIC, ])
      group_idsIC_tab = table(group_idsIC)
      n.groupsIC = length(group_idsIC_tab)
      group.idsIC = as.integer(group_idsIC)
      group.idsIC = group.idsIC- 1
      group.idsIC = as.matrix(group.idsIC)
      storage.mode(group.idsIC) <-"integer"
      storage.mode(n.groupsIC) <- "integer"
    }
    
    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC
    
    ICsol <- .Call(
      "MixedEffects_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg,
      max.steps.FPIRLS, threshold.FPIRLS, search,
      optim, lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      rand.effects.covariatesIC, group.idsIC, n.groupsIC, 
      PACKAGE = "fdaPDE"
    )
    
    if (nrow(covariates) != 0) {
      betaIC = matrix(data=ICsol[[5]],nrow=ncol(covariatesIC),ncol=length(lambdaSIC))
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
    
    # fix Mixed Effects with parabolic smoothing
    if (nrow(rand.effects.covariates) != 0) {
      rand.effects.covariates <- rand.effects.covariates[(NobsIC + 1):nrow(rand.effects.covariates), ]
      rand.effects.covariates <- as.matrix(rand.effects.covariates)
    }
    group_ids = as.factor(group.ids[(NobsIC + 1):length(group.ids)])
    group_ids_tab = table(group_ids)
    n.groups = length(group_ids_tab)
    group.idsIC = as.integer(group_idsIC)
    group.idsIC = group.idsIC- 1
    group.idsIC = as.matrix(group.idsIC)
    storage.mode(group.idsIC) <-"integer"
    storage.mode(n.groupsIC) <- "integer"
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  
  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
                       each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  
  
  ## Call C++ function
  bigsol <- .Call(
    "MixedEffects_Laplace_time", locations, bary.locations, time_locations,
    observations, FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim,
    covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
    FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, max.steps.FPIRLS, threshold.FPIRLS,
    search, optim, lambdaS, lambdaT,
    DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix,
    GCV.inflation.factor, lambda.optimization.tolerance,
    rand.effects.covariates, group.ids, n.groups, 
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}

CPP_smooth.volume.MixedEffects.FEM.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                             covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL,
                                             areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC,
                                             max.steps.FPIRLS = 15, threshold.FPIRLS = 0.0002020,
                                             search, bary.locations, optim, lambdaS = NULL, lambdaT = NULL,
                                             DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                             DOF.matrix = NULL, GCV.inflation.factor = 1,
                                             lambda.optimization.tolerance = 0.05,
                                             rand.effects.covariates, group.ids, n.groups) 
{
  # Indexes in C++ starts from 0, in R from 1
  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = 
    FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  group.ids = group.ids - 1
  
  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }
  
  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }
  
  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }
  
  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }
  
  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }
  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  FLAG_ITERATIVE<-as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE)<-"integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  rand.effects.covariates <- as.matrix(rand.effects.covariates)
  storage.mode(rand.effects.covariates) <-"double"
  group.ids <- as.matrix(group.ids)
  storage.mode(group.ids) <-"integer"
  storage.mode(n.groups) <- "integer"
  
  ## IC estimation for parabolic smoothing from the first column of observations
  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]
    
    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }
    
    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[notNAIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    # adjust Random Effects input for parabolic smoothing
    if (nrow(rand.effects.covariates) == 0) {
      rand.effects.covariatesIC <- rand.effects.covariates
      n.groupsIC <- n.groups
      group.idsIC <- group.ids
    } else {
      rand.effects.covariatesIC <- rand.effects.covariates[notNAIC, ]
      rand.effects.covariatesIC <- as.matrix(rand.effects.covariatesIC)
      storage.mode(rand.effects.covariatesIC) <- "double"
      
      group_idsIC = as.factor(group.ids[notNAIC, ])
      group_idsIC_tab = table(group_idsIC)
      n.groupsIC = length(group_idsIC_tab)
      group.idsIC = as.integer(group_idsIC)
      group.idsIC = group.idsIC- 1
      group.idsIC = as.matrix(group.idsIC)
      storage.mode(group.idsIC) <-"integer"
      storage.mode(n.groupsIC) <- "integer"
    }
    
    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC
    
    ICsol <- .Call(
      "MixedEffects_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg,
      max.steps.FPIRLS, threshold.FPIRLS, search,
      optim, lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      rand.effects.covariatesIC, group.idsIC, n.groupsIC, 
      PACKAGE = "fdaPDE"
    )
    
    if (nrow(covariates) != 0) {
      betaIC = matrix(data=ICsol[[5]],nrow=ncol(covariatesIC),ncol=length(lambdaSIC))
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
    
    # fix Mixed Effects with parabolic smoothing
    if (nrow(rand.effects.covariates) != 0) {
      rand.effects.covariates <- rand.effects.covariates[(NobsIC + 1):nrow(rand.effects.covariates), ]
      rand.effects.covariates <- as.matrix(rand.effects.covariates)
    }
    group_ids = as.factor(group.ids[(NobsIC + 1):length(group.ids)])
    group_ids_tab = table(group_ids)
    n.groups = length(group_ids_tab)
    group.idsIC = as.integer(group_idsIC)
    group.idsIC = group.idsIC- 1
    group.idsIC = as.matrix(group.idsIC)
    storage.mode(group.idsIC) <-"integer"
    storage.mode(n.groupsIC) <- "integer"
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  
  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
                       each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  
  
  ## Call C++ function
  bigsol <- .Call(
    "MixedEffects_Laplace_time", locations, bary.locations, time_locations,
    observations, FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim,
    covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
    FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, max.steps.FPIRLS, threshold.FPIRLS,
    search, optim, lambdaS, lambdaT,
    DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix,
    GCV.inflation.factor, lambda.optimization.tolerance,
    rand.effects.covariates, group.ids, n.groups, 
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}

CPP_smooth.volume.MixedEffects.FEM.PDE.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                        covariates = NULL, PDE_parameters, ndim, mydim, BC = NULL,
                                        incidence_matrix = NULL, areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, 
                                        max.steps, threshold, IC, max.steps.FPIRLS = 15,
                                        threshold.FPIRLS = 0.0002020, search, bary.locations, optim, lambdaS = NULL,
                                        lambdaT = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                        DOF.matrix = NULL, GCV.inflation.factor = 1,
                                        lambda.optimization.tolerance = 0.05,
                                        rand.effects.covariates, group.ids, n.groups) 
{

  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation
  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  group.ids = group.ids - 1
  
  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }


  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }

  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }

  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }

  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }

  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }

  if (is.null(time_locations)) {
    time_locations <- matrix(ncol = 0, nrow = 0)
  } else {
    time_locations <- as.matrix(time_locations)
  }

  if (is.null(time_mesh)) {
    time_mesh <- matrix(ncol = 0, nrow = 0)
  } else {
    time_mesh <- as.matrix(time_mesh)
  }


  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(PDE_parameters$K) <- "double"
  storage.mode(PDE_parameters$b) <- "double"
  storage.mode(PDE_parameters$c) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  FLAG_ITERATIVE<-as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE)<-"integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  rand.effects.covariates <- as.matrix(rand.effects.covariates)
  storage.mode(rand.effects.covariates) <-"double"
  group.ids <- as.matrix(group.ids)
  storage.mode(group.ids) <-"integer"
  storage.mode(n.groups) <- "integer"

  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]

    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }

    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[1:NobsIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    # adjust Random Effects input for parabolic smoothing
    if (nrow(rand.effects.covariates) == 0) {
      rand.effects.covariatesIC <- rand.effects.covariates
      n.groupsIC <- n.groups
      group.idsIC <- group.ids
    } else {
      rand.effects.covariatesIC <- rand.effects.covariates[1:NobsIC, ]
      rand.effects.covariatesIC <- as.matrix(rand.effects.covariatesIC)
      storage.mode(rand.effects.covariatesIC) <- "double"
      
      group_idsIC = as.factor(group.ids[1:NobsIC, ])
      group_idsIC_tab = table(group_idsIC)
      n.groupsIC = length(group_idsIC_tab)
      group.idsIC = as.integer(group_idsIC)
      group.idsIC = group.idsIC- 1
      group.idsIC = as.matrix(group.idsIC)
      storage.mode(group.idsIC) <-"integer"
      storage.mode(n.groupsIC) <- "integer"
    }

    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ICsol <- .Call(
      "MixedEffects_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg,
      max.steps.FPIRLS, threshold.FPIRLS, search,
      as.integer(c(0, 1, 1)), lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      rand.effects.covariatesIC, group.idsIC, n.groupsIC, 
      PACKAGE = "fdaPDE"
    )


    if (nrow(covariates) != 0) {
      betaIC <- ICsol[[15]]
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
    
    # fix Mixed Effects with parabolic smoothing
    if (nrow(rand.effects.covariates) != 0) {
      rand.effects.covariates <- rand.effects.covariates[(NobsIC + 1):nrow(rand.effects.covariates), ]
      rand.effects.covariates <- as.matrix(rand.effects.covariates)
    }
    group_ids = as.factor(group.ids[(NobsIC + 1):length(group.ids)])
    group_ids_tab = table(group_ids)
    n.groups = length(group_ids_tab)
    group.idsIC = as.integer(group_idsIC)
    group.idsIC = group.idsIC- 1
    group.idsIC = as.matrix(group.idsIC)
    storage.mode(group.idsIC) <-"integer"
    storage.mode(n.groupsIC) <- "integer"

    }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"

  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
    each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"

  ## Call C++ function
  bigsol <- .Call(
    "MixedEffects_PDE_time", locations, bary.locations, time_locations, observations,
    FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim, PDE_parameters$K,
    PDE_parameters$b, PDE_parameters$c, covariates, BC$BC_indices,
    BC$BC_values, incidence_matrix, areal.data.avg, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold,
    IC, max.steps.FPIRLS, threshold.FPIRLS, search,
    optim, lambdaS, lambdaT, DOF.stochastic.realizations, DOF.stochastic.seed,
    DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance,
    rand.effects.covariates, group.ids, n.groups, 
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}

CPP_smooth.volume.MixedEffects.FEM.PDE.sv.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                           covariates = NULL, PDE_parameters, ndim, mydim, BC = NULL,
                                           incidence_matrix = NULL, areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, 
                                           max.steps, threshold, IC, max.steps.FPIRLS = 15,
                                           threshold.FPIRLS = 0.0002020, search, bary.locations, optim, lambdaS = NULL,
                                           lambdaT = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                           DOF.matrix = NULL, GCV.inflation.factor = 1,
                                           lambda.optimization.tolerance = 0.05,
                                           rand.effects.covariates, group.ids, n.groups) 
{

  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation
  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  group.ids = group.ids - 1

  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }

  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }

  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }

  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }

  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }
  if(is.null(time_locations)){
    time_locations<-matrix(ncol=0, nrow=0)
  }else
  {
    time_locations <- as.matrix(time_locations)
  }

  if(is.null(time_mesh)){
    time_mesh<-matrix(ncol=0, nrow=0)
  }else
  {
    time_mesh <- as.matrix(time_mesh)
  }

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  rand.effects.covariates <- as.matrix(rand.effects.covariates)
  storage.mode(rand.effects.covariates) <-"double"
  group.ids <- as.matrix(group.ids)
  storage.mode(group.ids) <-"integer"
  storage.mode(n.groups) <- "integer"


  PDE_param_eval <- NULL
  points_eval <-
    matrix(CPP_get_evaluations_points(
      mesh = FEMbasis$mesh, order = FEMbasis$order
    ), ncol = 2)
  PDE_param_eval$K <- (PDE_parameters$K)(points_eval)
  PDE_param_eval$b <- (PDE_parameters$b)(points_eval)
  PDE_param_eval$c <- (PDE_parameters$c)(points_eval)
  PDE_param_eval$u <- (PDE_parameters$u)(points_eval)

  storage.mode(PDE_param_eval$K) <- "double"
  storage.mode(PDE_param_eval$b) <- "double"
  storage.mode(PDE_param_eval$c) <- "double"
  storage.mode(PDE_param_eval$u) <- "double"

  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]

    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }

    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[notNAIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    # adjust Random Effects input for parabolic smoothing
    if (nrow(rand.effects.covariates) == 0) {
      rand.effects.covariatesIC <- rand.effects.covariates
      n.groupsIC <- n.groups
      group.idsIC <- group.ids
    } else {
      rand.effects.covariatesIC <- rand.effects.covariates[notNAIC, ]
      rand.effects.covariatesIC <- as.matrix(rand.effects.covariatesIC)
      storage.mode(rand.effects.covariatesIC) <- "double"
      
      group_idsIC = as.factor(group.ids[notNAIC, ])
      group_idsIC_tab = table(group_idsIC)
      n.groupsIC = length(group_idsIC_tab)
      group.idsIC = as.integer(group_idsIC)
      group.idsIC = group.idsIC- 1
      group.idsIC = as.matrix(group.idsIC)
      storage.mode(group.idsIC) <-"integer"
      storage.mode(n.groupsIC) <- "integer"
    }

    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC

    ICsol <- .Call(
      "MixedEffects_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg,
      max.steps.FPIRLS, threshold.FPIRLS, search,
      as.integer(c(0, 1, 1)), lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      rand.effects.covariatesIC, group.idsIC, n.groupsIC, 
      PACKAGE = "fdaPDE"
    )
    
  ## shifting the lambdas interval if the best lambda is the higher one and retry smoothing
    if (nrow(covariates) != 0) {
      betaIC = matrix(data=ICsol[[5]],nrow=ncol(covariatesIC),ncol=length(lambdaSIC))
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
    
    # fix Mixed Effects with parabolic smoothing
    if (nrow(rand.effects.covariates) != 0) {
      rand.effects.covariates <- rand.effects.covariates[(NobsIC + 1):nrow(rand.effects.covariates), ]
      rand.effects.covariates <- as.matrix(rand.effects.covariates)
    }
    group_ids = as.factor(group.ids[(NobsIC + 1):length(group.ids)])
    group_ids_tab = table(group_ids)
    n.groups = length(group_ids_tab)
    group.idsIC = as.integer(group_idsIC)
    group.idsIC = group.idsIC- 1
    group.idsIC = as.matrix(group.idsIC)
    storage.mode(group.idsIC) <-"integer"
    storage.mode(n.groupsIC) <- "integer"
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"

  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
    each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"


  ## Call C++ function
  bigsol <- .Call(
    "MixedEffects_PDE_space_varying_time", locations, bary.locations, time_locations,
    observations, FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim,
    PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u,
    covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
    FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, max.steps.FPIRLS, threshold.FPIRLS,
    search, optim, lambdaS, lambdaT,
    DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix,
    GCV.inflation.factor, lambda.optimization.tolerance,
    rand.effects.covariates, group.ids, n.groups, 
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}

CPP_smooth.graph.MixedEffects.FEM.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                             covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL,
                                             areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC,
                                             max.steps.FPIRLS = 15, threshold.FPIRLS = 0.0002020,
                                             search, bary.locations, optim, lambdaS = NULL, lambdaT = NULL,
                                             DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                             DOF.matrix = NULL, GCV.inflation.factor = 1,
                                             lambda.optimization.tolerance = 0.05,
                                             rand.effects.covariates, group.ids, n.groups) 
{
  # Indexes in C++ starts from 0, in R from 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  num_sides = 2*dim(FEMbasis$mesh$edges)[1] 
  for(i in 1:num_sides){
    if( dim(FEMbasis$mesh$neighbors[[i]] )[1] > 0)
      FEMbasis$mesh$neighbors[[i]] = FEMbasis$mesh$neighbors[[i]] - 1
  }
  
  group.ids = group.ids - 1
  
  
  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }
  
  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }
  
  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }
  
  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }
  
  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }
  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  for(i in 1:num_sides)
    storage.mode(FEMbasis$mesh$neighbors[[i]]) <- "integer" 
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  FLAG_ITERATIVE<-as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE)<-"integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  rand.effects.covariates <- as.matrix(rand.effects.covariates)
  storage.mode(rand.effects.covariates) <-"double"
  group.ids <- as.matrix(group.ids)
  storage.mode(group.ids) <-"integer"
  storage.mode(n.groups) <- "integer"
  
  ## IC estimation for parabolic smoothing from the first column of observations
  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]
    
    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }
    
    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[notNAIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    # adjust Random Effects input for parabolic smoothing
    if (nrow(rand.effects.covariates) == 0) {
      rand.effects.covariatesIC <- rand.effects.covariates
      n.groupsIC <- n.groups
      group.idsIC <- group.ids
    } else {
      rand.effects.covariatesIC <- rand.effects.covariates[notNAIC, ]
      rand.effects.covariatesIC <- as.matrix(rand.effects.covariatesIC)
      storage.mode(rand.effects.covariatesIC) <- "double"
      
      group_idsIC = as.factor(group.ids[notNAIC, ])
      group_idsIC_tab = table(group_idsIC)
      n.groupsIC = length(group_idsIC_tab)
      group.idsIC = as.integer(group_idsIC)
      group.idsIC = group.idsIC- 1
      group.idsIC = as.matrix(group.idsIC)
      storage.mode(group.idsIC) <-"integer"
      storage.mode(n.groupsIC) <- "integer"
    }
    
    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC
    
    ICsol <- .Call(
      "MixedEffects_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg,
      max.steps.FPIRLS, threshold.FPIRLS, search,
      optim, lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      rand.effects.covariatesIC, group.idsIC, n.groupsIC, 
      PACKAGE = "fdaPDE"
    )
    
    if (nrow(covariates) != 0) {
      betaIC = matrix(data=ICsol[[5]],nrow=ncol(covariatesIC),ncol=length(lambdaSIC))
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
    
    # fix Mixed Effects with parabolic smoothing
    if (nrow(rand.effects.covariates) != 0) {
      rand.effects.covariates <- rand.effects.covariates[(NobsIC + 1):nrow(rand.effects.covariates), ]
      rand.effects.covariates <- as.matrix(rand.effects.covariates)
    }
    group_ids = as.factor(group.ids[(NobsIC + 1):length(group.ids)])
    group_ids_tab = table(group_ids)
    n.groups = length(group_ids_tab)
    group.idsIC = as.integer(group_idsIC)
    group.idsIC = group.idsIC- 1
    group.idsIC = as.matrix(group.idsIC)
    storage.mode(group.idsIC) <-"integer"
    storage.mode(n.groupsIC) <- "integer"
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  
  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
                       each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  
  
  ## Call C++ function
  bigsol <- .Call(
    "MixedEffects_Laplace_time", locations, bary.locations, time_locations,
    observations, FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim,
    covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
    FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, max.steps.FPIRLS, threshold.FPIRLS,
    search, optim, lambdaS, lambdaT,
    DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix,
    GCV.inflation.factor, lambda.optimization.tolerance,
    rand.effects.covariates, group.ids, n.groups, 
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}


