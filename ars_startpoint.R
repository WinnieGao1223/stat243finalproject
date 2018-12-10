#calculate the derivative of a function
Deriv_sub <- function(x, fun=h, a, b){
  if (x == a) {return ((fun(x + 1e-9)-fun(x))/1e-9)}
  if (x == b) {return ((fun(x) - fun(x - 1e-9))/1e-9)}
  if (a <= x && x <= b) {return((fun(x + 1e-9)-fun(x - 1e-9))/2e-9)}
}

init_mat <- function(a,b,fun,deriv){
  #this mat stores x,f(x) and df(x)
  mat <- matrix(nrow = 2,ncol = 3)
  mat[1,1] <- a
  mat[1,2] <- fun(a)
  mat[1,3] <- deriv(a,fun,a,b)
  mat[2,1] <- b
  mat[2,2] <- fun(b)
  mat[2,3] <- deriv(b,fun,a,b)
  return(mat)}
calc_z <-function(mat){
  return((mat[2, 2] - mat[1, 2] -mat[2, 1]*mat[2, 3] +
            mat[1, 1]*mat[1, 3])/(mat[1, 3] - mat[2, 3]))
}

init_z <- function(a,b,mat){
  z <- c()
  z[1] <- a
  z[2] <- calc_z(mat)
  z[3] <- b
  return(z)
}

update_z <- function(z,x_star,mat){
  index <- which(mat[,1] == x_star)
  z <- c(z,NA)
  tmp_mat1 <- mat[(index-1):index,]
  tmp_mat2 <- mat[index:(index+1),]
  z1 <- calc_z(tmp_mat1)
  z2 <- calc_z(tmp_mat2)
  z[(index+2):length(z)] <- z[(index+1):(length(z)-1)]
  z[index] <- z1
  z[index+1] <- z2
  return(z)
}
update_mat <- function(mat,x_star,fun,deriv,a,b){
  index <- sum(mat[,1] < x_star) #find the location to insert x*
  tmp <- matrix(NA,1,3)
  mat <- rbind(mat,tmp)
  mat[(index+2):nrow(mat),] <- mat[(index+1):(nrow(mat)-1),]
  mat[(index+1),] <- c(x_star,fun(x_star),deriv(x_star,fun,a,b))
  return(mat)
}

sample_x <- function(z, mat, n=1){
  uk1 <- exp(mat[,2]+(z[-1]-mat[,1])*mat[,3])
  uk2 <- exp(mat[,2]+(z[-length(z)]-mat[,1])*mat[,3])
  p <- rep(NA,nrow(mat))
  p <- (uk1-uk2)/mat[,3]
  
  # normalize p
  p_norm <- p/sum(p)
  w <- runif(n)
  # determine the range we want to sample from
  i <- sum(cumsum(p_norm) < w) + 1
  
  # sample x using inverse cdf
  samp_x <- (log(p[i]*mat[i,3]*runif(n)+exp(mat[i,2]+(z[i]-mat[i,1])*mat[i,3]))-mat[i,2])/mat[i,3]+mat[i,1]
  return(samp_x)
}

upper_bound <- function(z, mat, samp_x){
  #compare samp_x with z to find which segment of line to calculate
  index <- which(z > samp_x)[1] - 1
  u <- mat[index,2] + (samp_x - mat[index,1])*mat[index,3]
  return(u)
}

lower_bound <- function(mat, samp_x){
  index <- rev(which(mat[,1] < samp_x))[1]
  l <- ((mat[(index+1),1]-samp_x)*mat[index,2]+(samp_x-mat[index,1])*mat[(index+1),2])/(mat[(index+1),1]-mat[index,1])
  return(l)
}

check_concave <- function(mat){
  #check if h'(x) decrease monotonically
  return(prod(mat[,3][-1] <= mat[,3][-nrow(mat)])==1)
}

ars <- function(samp_n, g, a=-Inf, b=Inf){
  
  if(class(g) != "function"){
    stop('g should be a function')
  }
  
  ## Take log of input function ##
  h <- function(x){
    return(log(g(x)))
  }
  

  #initialize return value
  final_x <- rep(NA, samp_n)
  counter <- 0
  
  ##Check starting point
  ct <- 0
  step_width <- 0.5
  start_val <- 0
  
  if(a != -Inf && b != Inf){
    mat <- init_mat(a, b, h, Deriv_sub)
    z <- init_z(a, b, mat)
  }
  
  if(a == -Inf && b != Inf){
    if(start_val > b){
      start_val <- b - step_width
    }
    a_star <- start_val
    test_de <- Deriv_sub(a_star,fun=h,a,b)
    while( -Inf < test_de && test_de <= 0 && ct <= 50){
      a_star <- a_star - step_width
      test_de <- Deriv_sub(a_star,fun=h,a,b)
      ct <- ct + 1
    }
    mat <- init_mat(a_star, b, h, Deriv_sub)
    z <- init_z(a_star, b, mat) ## wrong? a_star or a?
  }
  
  if(a != -Inf && b == Inf){
    if(start_val < a){
      start_val <- a + step_width
    }
    b_star <- start_val
    test_de <- Deriv_sub(b_star, fun = h, a, b)
    while( 0 <= test_de && test_de < Inf && ct <= 50){
      b_star <- b_star + step_width
      test_de <- Deriv_sub(b_star, fun = h, a, b)
      ct <- ct + 1
    }
    mat <- init_mat(a, b_star, h, Deriv_sub)
    z <- init_z(a, b_star, mat)
  }
  
  if(a == -Inf && b == Inf){
    a_star <- start_val - step_width
    b_star <- start_val + step_width
    test_de_a <- Deriv_sub(a_star, fun = h, a, b)
    test_de_b <- Deriv_sub(b_star, fun = h, a, b)
    while(-Inf < test_de_a && test_de_a <= 0 && ct <= 50){
      a_star <- a_star - step_width
      test_de_a <- Deriv_sub(a_star, fun = h, a, b)
      ct <- ct + 1
    }
    while(0 <= test_de_b && test_de_b < Inf && ct <= 100){
      b_star <- b_star + step_width
      test_de_b <- Deriv_sub(b_star, fun = h, a, b)
      ct <- ct + 1
    }
    mat <- init_mat(a_star, b_star, h, Deriv_sub)
    z <- init_z(a_star, b_star, mat)
  }
  
  while(counter < samp_n){
    if (!check_concave(mat)){
      stop('Input must be a log-concave function')
    }
    samp_x_n = 1
    samp_x <- sample_x(z, mat, samp_x_n)
    w <- runif(samp_x_n)
    
    u <- upper_bound(z, mat, samp_x)
    l <- lower_bound(mat, samp_x)
    if(w <= exp(l-u)){
      counter = counter + 1
      final_x[counter] = samp_x
    }
    else{
      if (w <= exp(h(samp_x)-u)){
        counter = counter + 1
        final_x[counter] = samp_x
      }
      mat <- update_mat(mat, samp_x, h, Deriv_sub, a, b)
      z <- update_z(z, samp_x, mat)
    }
  }
  return(final_x)
}

dnor <- function(x){
  return((1/(sqrt(2*pi)))*exp(-(x^2)/2))
}

a <- ars(10000, dnor, a = -10, b = 10)
plot(density(a))
qqnorm(a)

b <- ars(10000, dnor, a = -Inf, b = 10)
plot(density(b))
c <- ars(10000, dnor, a = -Inf, b = Inf)
plot(density(c))
d <- ars(10000, dnor, a = -10, b = Inf)
plot(density(d))
