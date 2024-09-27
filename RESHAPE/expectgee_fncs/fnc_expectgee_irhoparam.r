fnc_expectgee_irhoparam <- function(ires, unexpect, corstr){

  iobsize <- length(ires)

  if (unexpect == 0.5){
    if (corstr=="ar1") {
      irhoval = 0
      for ( t in 1:(iobsize-1) ){
        irhoval = irhoval + ires[t]*ires[t+1]
      }
      irhoval = irhoval #/((iobsize-1))
    }

    if (corstr=="exc") {
      irhoval = 0
      for ( t in 1:iobsize ){
        for ( s in 1:iobsize){
          irhoval = irhoval + ires[t]*ires[s]
        }
      }
      ires2sum = sum(ires^2)
      irhoval = (irhoval - ires2sum) #/(iobsize*(iobsize-1))
    }

  } else {
    iphi <- fnc_expectgee_asymweight(ires, unexpect)
    if (corstr=="ar1") {
      irhoval = 0
      for ( t in 1:(iobsize-1) ){
        irhoval = irhoval + iphi[t]*ires[t]*iphi[t+1]*ires[t+1]
      }
      irhoval = irhoval #/(iobsize-1)
    }

    if (corstr=="exc") {
      irhoval = 0
      for ( t in 1:iobsize ){
        for ( s in 1:iobsize){
          irhoval = irhoval + iphi[t]*ires[t]*iphi[s]*ires[s]
        }
      }
      ires2sum = sum((iphi*ires)^2)
      irhoval = (irhoval - ires2sum) #/(iobsize*(iobsize-1))
    }

  }

  return(irhoval)

}