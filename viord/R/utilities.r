print_hash_row = function(l){
	cat(paste(rep('#',l), collapse=''),'\n')
}

print_msg = function(s, l = cli::console_width()-1){
	print_hash_row(l)
	cat(s, "\n")
	print_hash_row(l)

}

print_time = function() {
	time_str = paste("# Current time:", format(Sys.time(), "%H:%M:%S #"))
	print_hash_row(nchar(time_str))
	cat(time_str,'\n')
	print_hash_row(nchar(time_str))
}

# relevel ordinal factor removing missing
fct_drop_ordinal = function(x) {
	if(!is.ordered(x)) stop("input should be an ordered factor")
	K = length(unique(x))
	xd = factor(rank(x, ties.method = 'min'), ordered = T, labels = 1:K)
	return(xd)
}
