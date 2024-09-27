fnc_expectgee_listCov <- function(x, eres, unexpect, scaleparam, rhoparam, ltps, sampsize, corstr, scaled){

  if (corstr=="ind"){
    if (unexpect == 0.5) {
      ephi = rep(1, length(eres))
    } else {
      ephi = fnc_expectgee_asymweight(eres, unexpect)
    }
    Memat =  0
    Brmat = 0
    MematMod =  0
    for(i in 1:sampsize ){
      itps = ltps[[i]]
      ires = eres[itps]
      iphi = ephi[itps]
      xi = x[itps, ]
      sigmahat = tcrossprod(ires*iphi)
      Memat = Memat + t(xi) %*% sigmahat %*% xi
      Brmat = Brmat + t(xi) %*% diag(iphi) %*% xi
      MematMod = MematMod + scaleparam * t(xi) %*% xi
    }

    invBrmat = MASS::ginv(Brmat)
    sandCovhat = invBrmat %*% Memat %*% invBrmat # sandCovhat
    ModCovhat = invBrmat %*% MematMod %*% invBrmat

  } else {
    if (unexpect == 0.5) {
      ephi = rep(1, length(eres))
    } else {
      ephi = fnc_expectgee_asymweight(eres, unexpect)
    }
    Memat =  0
    Brmat = 0
    MematMod =  0
    for(i in 1:sampsize ){
      itps = ltps[[i]]
      ires = eres[itps]
      iphi = ephi[itps]
      xi = x[itps, ]
      sigmahat = tcrossprod(ires*iphi)
      ivcovmat = fnc_expectgee_ivcov(scaleparam, rhoparam, length(itps), corstr, scaled)
      ivcovmat = MASS::ginv(ivcovmat)
      corstxi  =  ivcovmat %*%  xi
      Memat = Memat + t(corstxi) %*% sigmahat %*% corstxi
      Brmat = Brmat + t(corstxi) %*% diag(iphi) %*% xi
      MematMod = MematMod + scaleparam * t(xi) %*% xi
    }
    invBrmat = MASS::ginv(Brmat)
    sandCovhat = invBrmat %*% Memat %*% invBrmat # sandCovhat
    ModCovhat = invBrmat %*% MematMod %*% invBrmat
  }

  return(list(sandCovhat=sandCovhat, ModCovhat=ModCovhat))

}