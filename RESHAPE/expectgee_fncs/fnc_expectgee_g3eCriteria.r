fnc_expectgee_g3eCriteria = function(object){

  outInd = fnc_expectgee_g3e(object$id, object$y, object$x, object$expectile, intercept='FALSE', corstr="ind", scaled='TRUE')
  lSigmaIndhatHH = outInd$lmodCovhat

  expectleng = length(object$expectile)
  leres = object$leres
  lSigmaRhat = object$lsandCovhat

  lQlike = lapply(1:expectleng, function(k) {out = - 0.5 * sum(leres[[k]]^2); out} )

  lCIC = lapply(1:expectleng, function(k) sum(diag(MASS::ginv(lSigmaIndhatHH[[k]]) %*% lSigmaRhat[[k]])) )

  lQIC = lapply(1:expectleng, function(k) { out = -2 * lQlike[[k]] + 2 * lCIC[[k]] ; out })

  out = rbind(unlist(lQIC), unlist(lCIC))

  colnames(out) = paste0('tau', 1:length(object$expectile))
  rownames(out) = c('asymQIC', 'asymCIC')

  return(out)
}