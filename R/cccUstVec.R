#' @title Wrapper around cccUst to use vector input
cccUstVec <- function(x,y){
  fit_tbl <- data.frame(ry = c(x,y), 
                        rmet = rep(c("a", "b"), each = length(x)))
  if(sd(x) == 0 || sd(y) == 0) return(NA)
  cccUst_obj = suppressWarnings(cccrm::cccUst(dataset = fit_tbl,
                                              ry = "ry",
                                              rmet = "rmet",
                                              cl = 0.95))
  cccUst_obj[[1]]
}
