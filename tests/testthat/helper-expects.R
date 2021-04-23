# helper-expects


expect_warnings <- function(object, times=1){
	if (times<=1){
		expect_warning(object)
	} else {
		expect_warning(expect_warnings(object,times-1))
	}
}