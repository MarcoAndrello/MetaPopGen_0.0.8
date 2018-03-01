settler.survival <-
  function(S,kappa0) {
    return( (1 / (1 + (1/kappa0) * S ) ))
  }