fnc_expectgee_sigparam <- function(res, unexpect, xncol){
  if (unexpect==0.5){
    out = sum(res^2)/(length(res) - xncol)
  } else {
    phi <- fnc_expectgee_asymweight(res, unexpect)
    out = sum( (phi*res)^2 )/(length(res) - xncol)
  }
  return(out)
}