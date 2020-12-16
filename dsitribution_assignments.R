


compute_assignments <- function(m, p, K){
  N = p * K
  if (m>N){
    print("Error: m must be < than N")
    return(NA)
  }
  if(m==0){
    return(list(rep(0,K)))
  }
  if(K ==1){
    return(m)
  }
  else{
    res = c()
    it = 1
    for (i in 0:min(p,m)){
      new_entries = compute_assignments(m-i,p, K-1)
      for (e in new_entries){
        res[[it]] <- c(i, e)
        it <- it + 1
        
      }
    }
    return(res)
  }
}




