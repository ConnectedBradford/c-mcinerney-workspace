fnc_expectgee_ivcov <- function(scalparam, rhoparam, iobsize, corstr, scaled){


  if (corstr=="ar1") {
    ivcovmat = matrix(NA, ncol=iobsize, nrow=iobsize)
    for ( t in 1:iobsize){
      for ( s in 1:iobsize ){
        ivcovmat[t, s] = rhoparam^(abs(t-s))
        ivcovmat[s, t] = ivcovmat[t, s]
      }
    }
  }

  if (corstr=="exchangeable"){
    ivcovmat = matrix(rhoparam, ncol=iobsize, nrow=iobsize)
    diag(ivcovmat) = 1
  }

  # if (corstr=="ind") ivcovmat = diag(iobsize)

  if (scaled=='TRUE')   out = scalparam*ivcovmat
  if (scaled=='FALSE')  out = ivcovmat

  #
  return(out)

}