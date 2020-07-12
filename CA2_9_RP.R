RP <- function(X1, X2, L, R){
  if (X1 < R){
    RP = 0
  }
  if (R < X1 < L){
    RP = X1 - R + X2*(X1 - R)/X1
  }
  if (X1 > L){
    RP = L - R + X2*(L-R)/L
  }
}
