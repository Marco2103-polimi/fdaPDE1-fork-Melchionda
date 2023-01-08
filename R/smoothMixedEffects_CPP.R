CPP_smooth.MixedEffects.FEM<-function(locations, observations, FEMbasis, covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL, areal.data.avg = FALSE, rand.effects.covariates, group.sizes, n.groups, max.steps.FPIRLS = 15, threshold.FPIRLS = 0.0002020, search, bary.locations, optim, lambda = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05)
{
  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  max.steps.FPIRLS = max.steps.FPIRLS - 1
  
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(DOF.matrix))
  {
    DOF.matrix<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
  }
  
  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }
  
  if(is.null(lambda))
  {
    lambda<-vector(length=0)
  }
  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
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
  storage.mode(BC$BC_values) <-"double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"
  rand.effects.covariates <- as.matrix(rand.effects.covariates)
  storage.mode(rand.effects.covariates) <-"double"
  group.sizes <- as.matrix(group.sizes)
  storage.mode(group.sizes) <-"integer"
  storage.mode(n.groups) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"  
  storage.mode(lambda) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  
  ## Call C++ function
  bigsol <- .Call("MixedEffects_Laplace", locations, bary.locations, observations, FEMbasis$mesh, FEMbasis$order,
                  mydim, ndim, covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
                  max.steps.FPIRLS, threshold.FPIRLS, search,
                  optim, lambda, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance,
                  rand.effects.covariates, group.sizes, n.groups, PACKAGE = "fdaPDE")
  
  return(bigsol)
}