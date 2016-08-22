#' Write parameter file
#'
#' Write simulation parameters to a file
#'
#' This function writes paramater assignments from a named list of parameters
#' to a given file. Each element in the list is written as a line that
#' can be read into R (using \code{\link{source}}) such the name of the list
#' element becomes the object name. The function only works for vectors, lists
#' with one level, and functions. It then appends two lines defining 
#' niche parameter lists for hosts and symbionts.
#' 
#' @param parm_list (required) named list of parameters to record
#' @param file_name (required) path to file where parameters should be written
#' @param file_dir directory where file should be saved, if not included in 
#' 	\code{file_name}. Defaults to current directory.
#' @return nothing
#'
#' @export
write_parms = function(parm_list, file_name, file_dir='./'){
	# Function for writing a new line containing a parameter assignment
	make_newline = function(varname, value){
		if(length(value)==1){
			if(is.numeric(value)) new_line = paste(varname, value, sep='=')
			if(is.character(value)) new_line = paste0(varname, "='", value, "'")
		} else {
			if(is.numeric(value)) new_line = paste0(varname, '=c(', paste(value, collapse=','), ')')
			if(is.character(value)) new_line = paste0(varname, '=c(', paste0("'", value, "'", collapse=','), ')')
		}
		new_line
	}

	# Open a new file and write each parameter contained in the list
	this_file = file(file.path(file_dir, file_name), open='w')
	for(i in 1:length(parm_list)){
		varname = names(parm_list)[i]
		value = parm_list[[i]]
	
		if(class(value)=='function'){
			new_line = c(paste0(varname,'='), deparse(value))	
		} else {	

		if(class(value)=='list'){
			line_list = paste(sapply(1:length(value), function(j){
				subvalue = value[[j]]
				make_newline(names(value)[j], subvalue)
			}), collapse=',')
			new_line = paste0(varname, '=list(', line_list, ')')
		} else {
			new_line = make_newline(varname, value)
		}}
				
		writeLines(new_line, this_file)
	}

	# Create lists of niche parameters
	writeLines('nicheparms_a = list(mu = c(mu_a1, mu_a2), rho = rho_a, sigma = c(sigma_a1, sigma_a2), alpha = c(alpha_a1, alpha_a2), r=r_a)', this_file)
	writeLines('nicheparms_b = list(mu = c(mu_b1, mu_b2), rho = rho_b, sigma = c(sigma_b1, sigma_b2), alpha = c(alpha_b1, alpha_b2), r=r_b)', this_file)

	close(this_file)
}
