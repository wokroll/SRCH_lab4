n_GARM <- 8
w_GRAN <- 1200
N <- 1024

generate_signal <- function(n_GARM = 8, w_GRAN = 1200, N){
  
  #амплітуда і зсув
  A <- runif(n_GARM,0,1)
  PHI <- runif(n_GARM,0,6.28)
  #список зі змінними омега
  w <- seq(w_GRAN/n_GARM,w_GRAN,w_GRAN/n_GARM)
  my_x <- seq(1,n_GARM)
  my_y = as.data.frame(matrix(0, nrow = N, ncol = my_x))
  
  for (i in my_x) for (j in seq(1,N)){
    my_y[j,i] <- A[i]*sin(w[i]*j+PHI[i])
  }
  return(rowSums(my_y))
}  


get_FFT <- function(x){
  freal11 <- numeric(N/2)
  freal12 <- numeric(N/2)
  freal1 <- numeric(N)
  fimage11 <- numeric(N/2)
  fimage12 <- numeric(N/2)
  fimage1 <- numeric(N)
  f1 <- numeric(N)
  
  for(p in 1:(N/2)){
    for(m in 1:(N/2)){
      freal11[p] <- freal11[p] + x[2 * m ] * cos(4 * pi / N * p * m)
      fimage11[p] <- fimage11[p] + x[2 * m ] * sin(4 * pi / N* p * m)
      freal12[p] <- freal12[p] + x[2 * m - 1] *  cos(4 * pi / N * p * m)
      fimage12[p] <- fimage12[p] + x[2 * m - 1] * sin(4 * pi / N * p * m)
    }
    freal1[p] <-  freal12[p] + freal11[p] * cos(2 * pi / N * p) - fimage11[p] * sin(2 * pi / N * p)
    fimage1[p] <-  fimage12[p] + fimage11[p] * cos(2 * pi / N * p) + freal11[p] * sin(2 * pi / N * p)
    freal1[p + N / 2] <-  freal12[p] - (freal11[p] * cos(2 * pi / N * p) - fimage11[p] * sin(2 * pi / N * p))
    fimage1[p + N / 2] <-  fimage12[p] - (fimage11[p] * cos(2 * pi / N * p) + freal11[p] * sin(2 * pi / N * p))
    f1[p] <-  sqrt(freal1[p] ** 2 + fimage1[p] ** 2)
    f1[p + N / 2] <-  sqrt(freal1[p + N / 2] ** 2 + fimage1[p + N / 2] ** 2)
  }
  return(f1)
}

sig <- generate_signal(N = N)

my_fft <- get_FFT(sig)

my_plot <- plot(seq(1,N), my_fft, type = "l")