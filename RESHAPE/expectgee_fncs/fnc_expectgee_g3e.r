fnc_expectgee_g3e <- function(id, y, x, expectile, intercept='FALSE', corstr, scaled){

  if (intercept == 'TRUE') x = cbind(1, x)

  ltps = split(1:length(id), id)
  tpslengvec = unlist(lapply(ltps, length))
  sampsize = length(ltps)
  xncol = ncol(x)
  expectleng = length(expectile)

  listbeta = list()
  leres = list()
  lsandCovhat = list()
  lmodCovhat = list()
  lrho = list()
  lsig = list()
  lQlike = list()
  lQIC = list()

  if (corstr=="ind"){

    for (k in 1:expectleng){
      # elmfit = lm(y~x-1)
      elmfit <- geepack::geeglm(y~x[, -1], id=id, corstr=corstr)
      eres = elmfit$residuals
      betaOld = elmfit$coef
      deltaParam = 1
      iter = 1
      while( deltaParam > 1e-05) {
        firstDeriv = 0
        scoreFun = 0
        for(i in 1:sampsize ){
          itps = ltps[[i]]
          ires = eres[itps]
          xi = x[itps, ]
          iphi = fnc_expectgee_asymweight(ires, expectile[k])
          txiPhi = t(xi) %*% diag(iphi)
          firstDeriv = firstDeriv + txiPhi %*% xi
          scoreFun = scoreFun + txiPhi %*% ires
        }
        betaNew = betaOld + MASS::ginv(firstDeriv) %*% scoreFun
        eres = y - x %*% betaNew
        betaDiff = abs(betaNew-betaOld)
        deltaParam = max(abs((betaNew - betaOld)/(betaOld + .Machine$double.eps)))  #  max(betaDiff)
        betaOld = betaNew
        iter = iter + 1
      }

      listbeta[[k]] = c(betaNew)
      leres[[k]] = c(eres)
      sigparamhat = fnc_expectgee_sigparam(eres, expectile[k], xncol)
      rhoparam = 0
      lsig[[k]] = sigparamhat
      VarCovlist = listCov(x, eres, expectile[k], sigparamhat, rhoparam, ltps, sampsize, corstr, scaled)
      lsandCovhat[[k]] = VarCovlist$sandCovhat
      lmodCovhat[[k]] = VarCovlist$ModCovhat
    }

    out = list(id=id, y=y, x=x, expectile=expectile, listbeta=listbeta, lmodCovhat=lmodCovhat,
               lsandCovhat=lsandCovhat, lsig=lsig, lrho=lrho, leres=leres)


  } else {

    vdof = fnc_expectgee_dof(tpslengvec,  xncol, corstr)

    for (k in 1:expectleng){

      # elmfit = lm(y~x-1)
      elmfit <- geepack::geeglm(y~x[ ,-1], id=id, corstr=corstr)
      eres = elmfit$residuals
      sigparamhat = 1 # sum(eres^2)/(length(eres) - xncol)
      rhoparamOld = 0.5 # (length(eres) - xncol)/vdof
      betaOld = elmfit$coef
      deltaParam = 1
      iter = 1

      while( deltaParam > 1e-05) {
        firstDeriv = 0
        scoreFun = 0

        for(i in 1:sampsize ){
          itps = ltps[[i]]
          ires = eres[itps]
          xi = x[itps, ]
          iphi = fnc_expectgee_asymweight(ires, expectile[k])
          ivcovmat = fnc_expectgee_ivcov(sigparamhat, rhoparamOld, length(itps), corstr, scaled)
          ivcovmat = MASS::ginv(ivcovmat)
          corEprx = t(xi) %*% (ivcovmat %*% diag(iphi))
          firstDeriv = firstDeriv + corEprx %*% xi
          scoreFun = scoreFun + corEprx %*% ires
        }

        betaNew = betaOld + MASS::ginv(firstDeriv) %*% scoreFun
        eres = y - x %*% betaNew
        sigparamhat = fnc_expectgee_sigparam(eres, expectile[k], xncol)
        rhoparamNew = sum(unlist(lapply(1:sampsize, function(i)irhoparam(eres[ltps[[i]]], expectile[k], corstr))))/(vdof * sigparamhat)
        betaDiff = abs((betaNew - betaOld)/(betaOld + .Machine$double.eps))
        rhoDiff = abs((rhoparamNew-rhoparamOld)/(rhoparamOld + .Machine$double.eps))
        deltaParam = max(c(betaDiff, rhoDiff))  #  max(betaDiff) #max(c(betaDiff, rhoDiff))
        betaOld = betaNew
        rhoparamOld = rhoparamNew
        iter = iter + 1

      }
      lrho[[k]] = rhoparamOld
      lsig[[k]] = sigparamhat
      listbeta[[k]] = c(betaNew)
      leres[[k]] = c(eres)
      VarCovlist = fnc_expectgee_listCov(x, eres, expectile[k], sigparamhat, rhoparamNew, ltps, sampsize, corstr, scaled)
      lsandCovhat[[k]] = VarCovlist$sandCovhat
      lmodCovhat[[k]] = VarCovlist$ModCovhat
    }

    out = list(id=id, y=y, x=x, expectile=expectile, listbeta=listbeta, lmodCovhat=lmodCovhat,
               lsandCovhat=lsandCovhat, lsig=lsig, lrho=lrho, leres=leres)

  }


  out

}