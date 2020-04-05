n_GARM <- 8
w_GRAN <- 1200
N <- 2048

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

get_table_dft <- function(sig){
  p = N
  k = N
  w_table_real <- matrix(0, p, k)
  w_table_image <- matrix(0, p ,k)
  w_table <- matrix(0, p ,k)
  
  for (i in 1:p) for(j in 1:k){
    w_table_real[i,j] <- cos(2*pi/N*i*j)
    w_table_image[i,j] <- sin(2*pi/N*i*j)
  }
  
  w_table <- list(w_table_real, w_table_image)
  
  Fp2_real <- numeric(N)
  Fp2_image <- numeric(N)

  for (i in 1:nrow(w_table[[1]])) for (j in  1:ncol(w_table[[1]])){
    Fp2_real[i] <- Fp2_real[i] + sig[j] * w_table[[1]][i,j]
    Fp2_image[i] <- Fp2_image[i] + sig[j] * w_table[[2]][i,j]
  }
  data <- data.frame(Fp2_real, Fp2_image)

  return(apply(data, 1, function(x) sqrt(x[[1]]**2+x[[2]]**2)))
}

sig <- generate_signal(N = N)

start <- Sys.time()
my_fft <- get_FFT(sig)
end <- Sys.time()
time_fft <- as.numeric(end - start)

start <- Sys.time()
table_dft <- get_table_dft(sig)
end <- Sys.time()
time_table_dft <- as.numeric(end - start)

result <- sprintf("Result for FFT is %.2f, result for Table DFT is %.2f, FFT is %.2fx faster", time_fft,time_table_dft,time_table_dft/time_fft)
