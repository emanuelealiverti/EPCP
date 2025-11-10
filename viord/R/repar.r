# transform free parameters (tau \in R^{K-2}) into constrained tresholds (alpha, ordered positve tresholds. -Inf and Inf set for identifiability. Inf as the last class is y_i = K iif z_i > alpha_k)
t2a = function(b) c(-Inf, cumsum(exp(b)), Inf)
# invese of t2a (remove -Inf and Inf)
a2t = function(a) log(diff(a[-c(1, length(a))]))

