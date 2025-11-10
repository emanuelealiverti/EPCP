# Routines for maximum aproximate marginal likelihood
# 1] alternating optimisation

#' @noRd
optim_ep_ml = function(Y,X,prior, maxit = 100, conv_tr = 1e-6){

	
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
		conv = abs(ll - ll_old) < conv_tr
		ll_old = ll

		it = it + 1
	}
	out = list('est' = tmp, 'alpha' = c(-Inf, alpha, Inf))
	return(out)
}



#' @noRd
optim_vb_ml = function(Y,X,prior, 
		       maxit = 100, conv_tr = 1e-6,
		       method = 'grad', vb_factor = "MF",full_path = F){

	# mean-field factorization
	vb_method = function(...){
		switch(vb_factor,
		MF = vb_ordinal(...),
		PMF = pmf_ordinal(...),
		)
	}

	# optimizer for alpha (Newton Rapson)
	optim_alpha = function(response, lin_pred, ...) {
		       suppressWarnings(ordinal::clm.fit(y=response,offset=lin_pred, link = 'probit',...)$alpha)
	}
	conv   =  F
	conv_warm = F
	ll_old =  Inf
	it     =  0
	lp     =  rep(0, NROW(X))
	# quantities used by PMF
	sigmaZ =  rep(1,NROW(X))
	muZ    =  rep(0, NROW(X))
	# initialize V and H

	while(!conv & it < maxit) {
		alpha = optim_alpha(response = Y, lp = lp, sigmaZ = sigmaZ, muZ = muZ)
		tmp = vb_method(Y = as.numeric(Y),
				X = X,
				 alpha = c(-Inf, alpha, Inf), 
				 mu0 = prior$mu0,
				 S0 = prior$S0,
				 Q0 = prior$Q0,
				 maxit = maxit, 
				 full_out = T)
		ll = tmp$elbo
		lp = drop(X %*% tmp$m)
		conv = abs(ll - ll_old) < conv_tr
		it = it + 1
		ll_old = ll
#		cat(it)
		if(vb_factor %in% c("PMF", "MIX")){
			muZ = lp
			sigmaZ = tmp$sigmaZ
		}


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




