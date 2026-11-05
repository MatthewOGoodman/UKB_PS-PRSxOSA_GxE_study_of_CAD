glm_function = function(
  
  formula, 
  family, 
  subset.vec, 
  data, 
  Q.cvts = NULL,
  Q.list = NULL,
  q.list = NULL,
  CI.alpha = 0.05,
  save_model_data = FALSE
  
){
  
  start = Sys.time()
  
  glm_fit = 
    try({
      glm(formula = formula, 
          family  = family, 
          subset  = subset.vec,
          data    = data,
          model   = save_model_data # whether to return the model frame
          
      )
    })
  
  n_glm =  length(glm_fit$linear.predictors)
  
  if(class(glm_fit)[1] != "try-error" ){
    
    glm_summary = summary(glm_fit)
    
    glm_fit$data = NULL  # original UKB_pheno_PGRS data frame too large to save
    # sapply(glm_fit, object.size)/object.size(glm_fit$data)

    if(!save_model_data){
      glm_fit$residuals = 
        glm_fit$fitted.values = 
        glm_fit$linear.predictors = 
        glm_fit$weights      = 
        glm_fit$prior.weights = 
        glm_fit$y = 
        glm_fit$qr$qr =
        glm_summary$deviance.resid = NULL }
    
  } else {
    glm_summary = NULL
  }
  
  # names(glm_fit)
  # names(glm_summary)
  
  if(!is.null(Q.cvts)  & !is.null(glm_summary)){
    
    glm_coef = glm_summary$coefficients
    glm_cov = glm_summary$cov.scaled
    
    # rownames(glm_summary$coef)
    # rownames(glm_summary$cov.scaled)
    # cbind(glm_summary$coefficients[,1]/sqrt(diag(glm_summary$cov.unscaled)), glm_summary$coefficients[,3], sqrt(diag(glm_summary$cov.unscaled)), glm_summary$coefficients[,2])
    cvt.coef = glm_coef[Q.cvts,"Estimate", drop = FALSE]
    cvt.cov  = glm_cov[Q.cvts,Q.cvts]
    
    Q.thetas = sapply(Q.list, function(x) t(x) %*% cvt.coef)
    Q.theta.diffs = Q.thetas - unlist(q.list)
    Q.thetas.se = sqrt(sapply(Q.list, function(x) t(x) %*% cvt.cov %*% x))
    Q.z.vals = Q.thetas / Q.thetas.se # glm_coef[Q.cvts,]
    Q.p.vals = 2*pnorm(-abs(Q.z.vals))
    Q.thetas.CI.lo = Q.thetas - qnorm(1-CI.alpha/2)*Q.thetas.se
    Q.thetas.CI.hi = Q.thetas + qnorm(1-CI.alpha/2)*Q.thetas.se
    
    Q = list(
      Q.thetas       = Q.thetas, 
      Q.thetas.se    = Q.thetas.se,
      Q.thetas.CI.lo = Q.thetas.CI.lo,
      Q.thetas.CI.hi = Q.thetas.CI.hi,
      Q.theta.diffs  = Q.theta.diffs,
      Q.diff.z.vals  = Q.z.vals, 
      Q.diff.p.vals  = Q.p.vals
    )
    
  } else {
    Q = list(
      Q.thetas       = NULL, 
      Q.thetas.se    = NULL,
      Q.thetas.CI.lo = NULL,
      Q.thetas.CI.hi = NULL,
      Q.theta.diffs  = NULL,
      Q.diff.z.vals  = NULL, 
      Q.diff.p.vals  = NULL
    )
  }
  
  # class(glm_fit)
  # names(glm_fit)
  # class(glm_summary)
  # names(summary(glm_fit))
  names.both = intersect( names(glm_fit), names(glm_summary))
  # lapply(names.both, function(x) all.equal(glm_fit[[x]], summary(glm_fit)[[x]]))
  # cbind(glm_fit$coefficients, glm_summary$coefficients) #glm_fit$coefficients is redundant, the others are equal
  lapply(names.both, function(x) glm_fit[[x]] <<- NULL) # set redundant variables to NULL
  
  if(class(glm_fit)[1] != "try-error" ){
    glm_res = c(glm_fit, glm_summary,  Q, list(Q.mat = as.matrix(as.data.frame(Q))) )
    print(paste("number of glm subjects", n_glm ))
  } else {
    glm_res = glm_fit
    print(paste("try-error"))
  }

  print(paste("glm time =", Sys.time() - start))
  
  return(glm_res)
  
}

####
####
####

gam_function = function(
  
  formula, 
  family, 
  subset.vec, 
  data, 
  Q.cvts = NULL,
  Q.list = NULL,
  q.list = NULL,
  CI.alpha = 0.05,
  save_model_data = FALSE, 
  gam_fit_fn = "gam"
  
){
  
  start = Sys.time()
  
  if(gam_fit_fn == "gam"){
    gam_fn = mgcv::gam #mgcv::gam, default gam fit model
  } else if(gam_fit_fn == "bam"){
    gam_fn = mgcv::bam #mgcv::bam, supposed to be big data efficient
  }
  
  
  gam_fit = 
    try({
      gam_fn(
        formula = formula, 
        family  = family, 
        subset  = subset.vec,
        data    = data
      )
    })
  
  n_gam = length(gam_fit$linear.predictors)

  if(class(gam_fit)[1] != "try-error" ){
    
    gam_summary = summary(gam_fit)
    if(!save_model_data){
      gam_fit$model =   # model.frame too large to save
        gam_fit$smooth =  # data for smooths basis
        gam_fit$residuals = 
        gam_fit$fitted.values = 
        gam_fit$linear.predictors = 
        gam_fit$weights      = 
        gam_fit$working.weights = 
        gam_fit$prior.weights = 
        gam_fit$y = 
        gam_fit$hat =
        gam_fit$offset = NULL
    }
    
  } else {
    gam_summary = NULL
  }
  

  # names(glm_fit)
  # names(glm_summary)
  
  if(!is.null(Q.cvts) & !is.null(gam_summary)){
    
    #spline coefs
    gam_coef.s = gam_summary$s.table
    #parametric coefs
    gam_coef = gam_summary$p.table
    gam_cov = gam_summary$cov.scaled
    
    # rownames(gam_summary$coef)
    # rownames(gam_summary$cov.scaled)
    # cbind(gam_summary$coefficients[,1]/sqrt(diag(gam_summary$cov.unscaled)), gam_summary$coefficients[,3], sqrt(diag(gam_summary$cov.unscaled)), gam_summary$coefficients[,2])
    cvt.coef = gam_coef[Q.cvts,"Estimate", drop = FALSE]
    cvt.cov  = gam_cov[Q.cvts,Q.cvts]
    
    Q.thetas = sapply(Q.list, function(x) t(x) %*% cvt.coef)
    Q.theta.diffs = Q.thetas - unlist(q.list)
    Q.thetas.se = sqrt(sapply(Q.list, function(x) t(x) %*% cvt.cov %*% x))
    Q.z.vals = Q.theta.diffs / Q.thetas.se 
    Q.p.vals = 2*pnorm(-abs(Q.z.vals))
    # gam_coef[Q.cvts,]
    Q.thetas.CI.lo = Q.thetas - qnorm(1-CI.alpha/2)*Q.thetas.se
    Q.thetas.CI.hi = Q.thetas + qnorm(1-CI.alpha/2)*Q.thetas.se
    
    Q = list(
      Q.thetas       = Q.thetas, 
      Q.thetas.se    = Q.thetas.se,
      Q.thetas.CI.lo = Q.thetas.CI.lo,
      Q.thetas.CI.hi = Q.thetas.CI.hi,
      Q.theta.diffs  = Q.theta.diffs,
      Q.diff.z.vals  = Q.z.vals, 
      Q.diff.p.vals  = Q.p.vals
    )
    
  } else {
    Q = list(
      Q.thetas       = NULL, 
      Q.thetas.se    = NULL,
      Q.thetas.CI.lo = NULL,
      Q.thetas.CI.hi = NULL,
      Q.theta.diffs  = NULL,
      Q.diff.z.vals  = NULL, 
      Q.diff.p.vals  = NULL
    )
  }
  
  names.both = intersect( names(gam_fit), names(gam_summary))
  # lapply(names.both, function(x) all.equal(gam_fit[[x]], summary(gam_fit)[[x]]))
  # cbind(gam_fit$edf, gam_summary$edf) #gam_fit$coefficients is redundant, the others are equal
  gam_summary$edf = NULL # set redundant variables to NULL
  
  if(class(gam_fit)[1] != "try-error" ){
    gam_res = c(gam_fit, gam_summary, Q, list(Q.mat = as.matrix(as.data.frame(Q))) )
    print(paste("number of gam subjects", n_gam ))
  } else {
    gam_res = gam_fit
    print(paste("try-error"))
  }
  
  print(paste("gam time =", Sys.time() - start))
  
  return(gam_res)
  
}

####
####
####

