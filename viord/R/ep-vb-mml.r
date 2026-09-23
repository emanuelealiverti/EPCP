# Routines for maximum aproximate marginal likelihood
# 1] alternating optimisation

#' @noRd
check_spd_matrix = function(x, name, p){
	if(!is.matrix(x) || !is.numeric(x) || any(dim(x) != c(p, p))){
		stop(sprintf("%s must be a numeric %d x %d matrix.", name, p, p))
	}
	if(any(!is.finite(x))){
		stop(sprintf("%s must contain only finite values.", name))
	}
	if(!isSymmetric(x, tol = sqrt(.Machine$double.eps))){
		stop(sprintf("%s must be symmetric.", name))
	}
	is_spd = tryCatch({
		chol(x)
		TRUE
	}, error = function(e) FALSE)
	if(!is_spd){
		stop(sprintf("%s must be positive definite.", name))
	}
	invisible(TRUE)
}

#' @noRd
check_random_effects = function(Z = NULL, Z_group = NULL, n){
	if(is.null(Z)){
		if(!is.null(Z_group)){
			stop("Z_group was provided, but Z is NULL.")
		}
		return(list(Z = matrix(0, n, 1), Z_group = 0L,
			    has_random = FALSE, labels = NULL))
	}

	if(!is.matrix(Z) || !is.numeric(Z)){
		stop("Z must be a numeric matrix.")
	}
	if(NROW(Z) != n){
		stop("Z must have the same number of rows as X.")
	}

	q = NCOL(Z)
	if(q <= 1){
		return(list(Z = matrix(0, n, 1), Z_group = 0L,
			    has_random = FALSE, labels = NULL))
	}

	if(is.null(Z_group)){
		stop("Z_group must be provided when Z has columns.")
	}
	if(length(Z_group) != q){
		stop("Z_group must have length equal to ncol(Z).")
	}
	if(any(is.na(Z_group))){
		stop("Z_group cannot contain NA values.")
	}

	labels = unique(as.character(Z_group))
	Z_group = as.integer(match(as.character(Z_group), labels) - 1L)
	list(Z = Z, Z_group = Z_group, has_random = TRUE, labels = labels)
}

#' @noRd
check_prior = function(prior, p, algorithm, has_random = FALSE){
	if(!is.list(prior)){
		stop("prior must be a list.")
	}

	check_mu0 = function(){
		if(is.null(prior$mu0)){
			stop(sprintf("For %s, prior must contain mu0.", algorithm))
		}
		if(!is.numeric(prior$mu0) || length(prior$mu0) != p){
			stop(sprintf("For %s, prior$mu0 must be a numeric vector of length %d.", algorithm, p))
		}
		if(any(!is.finite(prior$mu0))){
			stop(sprintf("For %s, prior$mu0 must contain only finite values.", algorithm))
		}
	}

	check_mu0()

	if(algorithm == "VB_prior"){
		if(is.null(prior$a0) || is.null(prior$b0)){
			stop(sprintf("For %s, prior must contain mu0, a0, and b0.", algorithm))
		}
		if(!is.numeric(prior$a0) || length(prior$a0) != 1 || !is.finite(prior$a0) || prior$a0 <= 0){
			stop(sprintf("For %s, prior$a0 must be a strictly positive finite number.", algorithm))
		}
		if(!is.numeric(prior$b0) || length(prior$b0) != 1 || !is.finite(prior$b0) || prior$b0 <= 0){
			stop(sprintf("For %s, prior$b0 must be a strictly positive finite number.", algorithm))
		}
		if(algorithm == "VB_prior" && has_random){
			if(is.null(prior$au0) || is.null(prior$bu0)){
				stop("For VB_prior with random effects, prior must contain au0 and bu0.")
			}
			if(!is.numeric(prior$au0) || length(prior$au0) != 1 || !is.finite(prior$au0) || prior$au0 <= 0){
				stop("For VB_prior, prior$au0 must be a strictly positive finite number.")
			}
			if(!is.numeric(prior$bu0) || length(prior$bu0) != 1 || !is.finite(prior$bu0) || prior$bu0 <= 0){
				stop("For VB_prior, prior$bu0 must be a strictly positive finite number.")
			}
		}
	} else if(algorithm == "PMF_mixed"){
		if(!has_random){
			stop("PMF_mixed requires a random-effects design Z (with at least two columns) and Z_group.")
		}
		if(is.null(prior$Q0)){
			stop("For PMF_mixed, prior must contain mu0, Q0, and s_sigma.")
		}
		check_spd_matrix(prior$Q0, "prior$Q0 for PMF_mixed", p)
		if(!is.numeric(prior$s_sigma) || length(prior$s_sigma) != 1 || !is.finite(prior$s_sigma) || prior$s_sigma <= 0){
			stop("For PMF_mixed, prior$s_sigma must be a strictly positive finite number.")
		}
	} else {
		if(is.null(prior$S0) || is.null(prior$Q0)){
			stop(sprintf("For %s, prior must contain mu0, S0, and Q0.", algorithm))
		}
		check_spd_matrix(prior$S0, sprintf("prior$S0 for %s", algorithm), p)
		check_spd_matrix(prior$Q0, sprintf("prior$Q0 for %s", algorithm), p)
	}

	invisible(TRUE)
}

#' @noRd
optim_ep_ml = function(Y,X,prior, maxit = 100, conv_tr = 1e-6){

	check_prior(prior, NCOL(X), "EP")

	
	# optimizer for alpha
	optim_alpha = function(response, lin_pred,...) {
		suppressWarnings(ordinal::clm.fit(y=response,offset=lin_pred, link = 'probit')$alpha)
	}
	# initial estimate
	conv = F
	ll_old = Inf
	it = 0
	lp = numeric(NROW(X))
	while(!conv & it < maxit) {
		alpha = optim_alpha(Y, lp)
		tmp = ep_ordinal(Y = as.numeric(Y),
				 X = X,
				 alpha = c(-Inf, alpha, Inf), 
				 mu0 = prior$mu0,
				 S0 = prior$S0,
				 Q0 = prior$Q0,
				 maxit = maxit, 
				 full_out = T)
		ll = tmp$logZ
		lp = drop(X %*% tmp$m)
		# relative tolerance: the log marginal likelihood scales with n
		conv = abs(ll - ll_old) < conv_tr * (1 + abs(ll))
		ll_old = ll

		it = it + 1
	}
	out = list('est' = tmp, 'alpha' = c(-Inf, alpha, Inf))
	return(out)
}



#' @noRd
optim_vb_ml = function(Y,X,prior, 
		       maxit = 100, conv_tr = 1e-6,
		       method = 'grad', vb_factor = "MF", full_path = F,
		       Z = NULL, Z_group = NULL, warm_start = TRUE,
		       conv_crit = c("elbo", "coef")){

	conv_crit = match.arg(conv_crit)

	vb_factor = match.arg(vb_factor, c("MF", "PMF", "VB_prior", "PMF_mixed"))
	random = check_random_effects(Z = Z, Z_group = Z_group, n = NROW(X))
	if(random$has_random && !(vb_factor %in% c("VB_prior", "PMF_mixed"))){
		stop("Z and Z_group are currently supported only with algorithm = 'VB_prior' or 'PMF_mixed'.")
	}
	check_prior(prior, NCOL(X), vb_factor, has_random = random$has_random)

	# optimizer for alpha (Newton Rapson)
	optim_alpha = function(response, lin_pred, ...) {
		       suppressWarnings(ordinal::clm.fit(y=response,offset=lin_pred, link = 'probit',...)$alpha)
	}
	conv   =  F
	conv_warm = F
	ll_old =  Inf
	it     =  0
	lp     =  rep(0, NROW(X))
	# warm start: state carried from one threshold iteration to the next
	init   =  NULL
	state_old = NULL # used by the "coef" convergence criterion
	warm_state = function(tmp) {
		switch(vb_factor,
		       MF        = tmp["m"],
		       PMF       = tmp["meanZ"],
		       VB_prior  = tmp[c("m_joint", "sigma_b2_b", "sigma_u2_b")],
		       PMF_mixed = tmp[c("meanZ", "sigma_u2_b", "a_u_b")])
	}

	while(!conv & it < maxit) {
		alpha = optim_alpha(response = Y, lin_pred = lp)
		if(vb_factor == "VB_prior"){
			au0 = if(random$has_random) prior$au0 else NA_real_
			bu0 = if(random$has_random) prior$bu0 else NA_real_
			tmp = vb_ordinal_prior(Y = as.numeric(Y),
					X = X,
					alpha = c(-Inf, alpha, Inf),
					mu0 = prior$mu0,
					a0 = prior$a0,
					b0 = prior$b0,
					Z = random$Z,
					Z_group = random$Z_group,
					au0 = au0,
					bu0 = bu0,
					maxit = maxit,
					conv_crit = conv_crit,
					full_out = T,
					init = init)
		} else if(vb_factor == "PMF_mixed"){
			tmp = pmf_ordinal_mixed(Y = as.numeric(Y),
					X = X,
					alpha = c(-Inf, alpha, Inf),
					mu0 = prior$mu0,
					Q0 = prior$Q0,
					Z = random$Z,
					Z_group = random$Z_group,
					s_sigma = prior$s_sigma,
					maxit = maxit,
					conv_crit = conv_crit,
					full_out = T,
					init = init)
		} else {
			tmp = switch(vb_factor,
				     MF = vb_ordinal(Y = as.numeric(Y),
						      X = X,
						      alpha = c(-Inf, alpha, Inf),
						      mu0 = prior$mu0,
						      S0 = prior$S0,
						      Q0 = prior$Q0,
						      maxit = maxit,
						      conv_crit = conv_crit,
						      full_out = T,
						      init = init),
				     PMF = pmf_ordinal(Y = as.numeric(Y),
							X = X,
							alpha = c(-Inf, alpha, Inf),
							mu0 = prior$mu0,
							S0 = prior$S0,
							Q0 = prior$Q0,
							maxit = maxit,
							conv_crit = conv_crit,
							full_out = T,
							init = init))
		}
		ll = tmp$elbo
		lp = drop(X %*% tmp$m)
		if(vb_factor %in% c("VB_prior", "PMF_mixed") && random$has_random){
			lp = lp + drop(random$Z %*% tmp$m_u)
		}
		if(warm_start) init = warm_state(tmp)
		if(conv_crit == "coef"){
			# monitor thresholds and coefficients rather than the ELBO
			state = c(alpha, tmp$m, if(is.null(tmp$m_u)) NULL else tmp$m_u)
			conv = !is.null(state_old) &&
				max(abs(state - state_old) / (1 + abs(state))) < conv_tr
			state_old = state
		} else {
			# relative tolerance: the ELBO scales with n
			conv = abs(ll - ll_old) < conv_tr * (1 + abs(ll))
		}
		it = it + 1
		ll_old = ll

	}
	out = list('est' = tmp, 'alpha' = c(-Inf, alpha, Inf))
	return(out)
}

# Crude optimization with optim (unefficient)

#' @noRd
optim_vb_brute = function(Y,X,prior){

	K = length(levels(Y))
	# optimize via optim creating a dedicated function
	fn = function(pars) {
		aa = t2a(pars)
		tmp = vb_ordinal(Y = as.numeric(Y), 
				 X = X, 
				 alpha = aa, 
				 mu0 = prior$mu0,
				 S0 = prior$S0,
				 Q0 = prior$Q0,
				 maxit = 100, 
				 full_out = F)
		-tmp$elbo
	}
	margL = try(optim(rep(0,K-2), fn))
	aa = t2a(margL$par)
	tmp = vb_ordinal(Y = as.numeric(Y), 
			 X = X, 
			 alpha = aa, 
			 mu0 = prior$mu0,
			 S0 = prior$S0,
			 Q0 = prior$Q0,
			 maxit = 500, 
			 full_out = F)
	return(list('est' = tmp, 'alpha' = aa))
}



#' @noRd
optim_ep_brute = function(Y,X,prior){

	K = length(levels(Y))
	# optimize via optim creating a dedicated function
	fn = function(pars) {
		aa = t2a(pars)
		tmp = ep_ordinal(Y = as.numeric(Y), 
				 X = X, 
				 alpha = aa, 
				 mu0 = prior$mu0,
				 S0 = prior$S0,
				 Q0 = prior$Q0,
				 maxit = 100, 
				 full_out = F)
		-tmp$logZ
	}
	margL = try(optim(rep(0,K-2), fn))
	aa = t2a(margL$par)
	tmp = ep_ordinal(Y = as.numeric(Y), 
			 X = X, 
			 alpha = aa, 
			 mu0 = prior$mu0,
			 S0 = prior$S0,
			 Q0 = prior$Q0,
			 maxit = 500, 
			 full_out = F)
	return(list('est' = tmp, 'alpha' = aa))
}
