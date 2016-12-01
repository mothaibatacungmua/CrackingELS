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

x = seq(0,2, length.out = 64)
b = haar_father(x, 0, 0)
plot(x, b, type = 'l', main = 'Haar(0,0)')

b = haar_mother(x, 1,1)
plot(x, b, type = 'l', main = 'Haar(2,1)')

b = haar_mother(x, 1, 2)
plot(x, b, type = 'l', main = 'Haar(2,2)')

b = haar_mother(x, 2, 1)
plot(x, b, type = 'l', main = 'Haar(2,2)')

b = haar_mother(x, 2, 3)
plot(x, b, type = 'l', main = 'Haar(2,3)')

b = haar_mother(x, 3, 1)
plot(x, b, type = 'l', main = 'Haar(3,1)')

b = haar_mother(x, 3, 2)
plot(x, b, type = 'l', main = 'Haar(3,2)')

b = haar_mother(x, 3, 3)
plot(x, b, type = 'l', main = 'Haar(3,3)')

b = haar_mother(x, 4, 1)
plot(x, b, type = 'l', main = 'Haar(4,1)')


# code for wavelet transforming
father_coefF_expand = 1/(sqrt(2))
mother_coeff_expand <- function(n){
  t = n %% 2
  ret = (-1)^t*father_coefF_expand
  return(ret)
}

N = 256
J = log(N, 2)

# using dynamic programming to do the transformation
basis_coeff <- function(y){
  raw_coeffs = rep(0, N*2)  
  raw_coeffs[(N+1):(N*2)] = y[1:N]
  ret_coeffs = rep(0, N)
  
  for(j in (J-1):0){
    raw_coeffs[2^j] = father_coefF_expand*raw_coeffs[2^(j+1)]
    ret_coeffs[2^j] = mother_coeff_expand(2^j) * raw_coeffs[2^(j+1)]
    for(k in 1:(2^j-1)){
        raw_coeffs[2^j + k] = father_coefF_expand*sum(raw_coeffs[2^(j+1):(2^(j+1) + 2*k - 1)])  
        h = mother_coeff_expand(seq(0, 2*k-1))
        ret_coeffs[2^j + k] = sum(h * raw_coeffs[2^(j+1):(2^(j+1) + 2*k - 1)])
    }
  }
  
  return(ret_coeffs)
}

transform_func <- function(x, coeffs){
  ret = coeffs[1] * father(x)
  
  for(j in 1:(J-1)){
    for(k in 0:(2^j-1)){
      ret = ret + coeffs[2^j + k] * haar_mother(x, j, k)
    }
  }
  
  return(ret/(sqrt(N)))
}

# testing
x = seq(0,10, length.out = N)
y = sin(x) + rnorm(N, 0, sd = 0.5)
plot(x,y)
