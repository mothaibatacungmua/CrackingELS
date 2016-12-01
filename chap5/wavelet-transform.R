father <- function(x){
  result = rep(0, length(x))
  result[which((x >=0) & (x <= 1))] = 1
  
  return(result)
}

haar_father <- function(x, j, k){
  return(2^(j/2)*(father(2^j*x - k)))
}

mother <- function(x){
  return(father(2*x) - father(2*x-1))  
}

haar_mother <- function(x, j, k){
  return(2^(j/2)*(mother(2^j*x - k)))
}
# plot some haar basis
par(mfrow=c(3,3))

x = seq(-2,2, length.out = 8192)
b = haar_father(x, 0, 0)
plot(x, b, type = 'l', main = 'Haar(0,0)')

b = haar_mother(x, 0,0)
plot(x, b, type = 'l', main = 'Haar(1,1)')

b = haar_mother(x, 1,1)
plot(x, b, type = 'l', main = 'Haar(1,1)')

b = haar_mother(x, 1, 2)
plot(x, b, type = 'l', main = 'Haar(1,2)')

b = haar_mother(x, 2, 1)
plot(x, b, type = 'l', main = 'Haar(2,1)')

b = haar_mother(x, 2, 3)
plot(x, b, type = 'l', main = 'Haar(2,3)')

b = haar_mother(x, 3, 1)
plot(x, b, type = 'l', main = 'Haar(3,1)')

b = haar_mother(x, 3, 2)
plot(x, b, type = 'l', main = 'Haar(3,2)')

b = haar_mother(x, 3, 3)
plot(x, b, type = 'l', main = 'Haar(3,3)')

b = haar_mother(x, 8, 1)
plot(x, b, type = 'l', main = 'Haar(8,1)')


# code for wavelet transforming
fcof = c(1/sqrt(2),1/sqrt(2))
mcof = c(1/sqrt(2), -1/sqrt(2))

N = 8
J = log(N, 2)

# using dynamic programming to do the transformation
basis_coeff <- function(y){
  raw_coeffs = rep(0, N*2)  
  raw_coeffs[(N+1):(N*2)] = y[1:N]
  ret_coeffs = rep(0, N)
  
  for(j in (J-1):0){
    for(k in 0:(2^j-1)){
        raw_coeffs[2^j + 1 + k] = fcof[1]*(raw_coeffs[2^(j+1)+1 + 2*k]) + 
                                  fcof[2]*(raw_coeffs[2^(j+1)+1 + 2*k + 1])
        
        ret_coeffs[2^j + 1 + k] = mcof[1]*(raw_coeffs[2^(j+1)+1 + 2*k]) + 
                                  mcof[2]*(raw_coeffs[2^(j+1)+1 + 2*k + 1])
    }
  }
  
  ret_coeffs[1] = raw_coeffs[2]
  
  return(ret_coeffs*2^(-J/2))
}

transform_func <- function(x, coeffs){
  ret = coeffs[1] * father(x)
  
  for(j in 0:(J-1)){
    for(k in 0:(2^j-1)){
      ret = ret + coeffs[2^j + 1 + k] * haar_mother(x, j, k)
    }
  }
  
  return(ret)
}

# testing case
N = 8
J = log(N,2)
x = seq(0,1,length.out = 8)
y = c(1, 0, -3, 2, 1, 0, 1, 2)

#coeffs should be c(1/2, -1/2, 1/(2*sqrt(2)), -1/(2*sqrt(2)), 1/4, -5/4, 1/4, -1/4)
coeffs = basis_coeff(y)

# using wavelet to reconstruct
N = 256
J = log(N, 2)
x = seq(0,1, length.out = N)
y = sin(50*x) + rnorm(N, 0, sd = 0.5)
plot(x,y, type = 'l')

coeffs = basis_coeff(y)
y_hat = transform_func(x, coeffs)
points(x, y_hat, type='l', col='red')
