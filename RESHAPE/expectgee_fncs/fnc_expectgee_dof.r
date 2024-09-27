fnc_expectgee_dof <- function(tpslengvec,  xncol, corstr){

  if (corstr=="ar1") out = sum(tpslengvec-1) - xncol

  if (corstr=="exchangeable"){
    out = 0
    for (i in 1:length(tpslengvec)){
      mi = tpslengvec[i]
      out = out + mi*(mi-1)
    }
    out = out - xncol # out/2 - xncol
  }

  # if (corstr=="ind") out = 1
  return(out)
}