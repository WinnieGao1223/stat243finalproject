
#calculate the derivative of a function numerically
#x is the point to evaluate the derivate at
#fun is the log-concave function h
#a is the lower bound, b is the upper bound
#return the derivative approximately h'(x) at x
Deriv_sub <- function(x, fun = h, a, b){
  if (x == a) {return ((fun(x + 1e-9)-fun(x))/1e-9)}
  if (x == b) {return ((fun(x) - fun(x - 1e-9))/1e-9)}
  if (a <= x && x <= b) {return((fun(x + 1e-9)-fun(x - 1e-9))/2e-9)}
}

#initialize the info matrix giving lower bound and upper bound
#a is the lower bound, b is the upper bound
#fun is the log-concave function h
#deriv is the derivative function
#return a matrix including initial x, h(x), h'(x)
init_mat <- function(a,b,fun,deriv){
  mat <- matrix(nrow = 2,ncol = 3)
  mat[1,1] <- a
  mat[1,2] <- fun(a)
  if(!is.finite(mat[1,2])){
    stop('a is not defined on function h', call. = FALSE)
  }
  mat[1,3] <- deriv(a)
  mat[2,1] <- b
  mat[2,2] <- fun(b)
  if(!is.finite(mat[2,2])){
    stop('b is not defined on function h', call. = FALSE)
  }
  mat[2,3] <- deriv(b)
  return(mat)
}

#calculate z value given two adjacent values of x
calc_z <-function(mat){
  return((mat[2, 2] - mat[1, 2] -mat[2, 1]*mat[2, 3] +
            mat[1, 1]*mat[1, 3])/(mat[1, 3] - mat[2, 3]))
}

#mat is the initial infor matrix
#return initialize vector z given the initial infor matrix
#z_0 is lower bound, z_k is the upper bound
init_z <- function(a,b,mat){
  z <- c()
  z[1] <- mat[1,1]
  z[3] <- mat[2,1]
  #if the slopes at a and b are equal, set z as the mean of a and b
  if (abs(mat[2,3]-mat[1,3]) < 1e-8){
    z[2] <- (mat[2,1] + mat[1,1])/2
  } else {
    z[2] <- calc_z(mat)
  }
  return(z)
}

#mat is the ordered infor matrix
#x_star is the sample needed to be updated
#fun is the log-concave function h
#deriv is the derivative function
#retun updated infor matrix with ascending order of x
update_mat <- function(mat,x_star,fun,deriv){
  index <- sum(mat[,1] < x_star) #find the location to insert x*
  tmp <- matrix(NA,1,3)
  mat <- rbind(mat,tmp)
  mat[(index+2):nrow(mat),] <- mat[(index+1):(nrow(mat)-1),]
  mat[(index+1),] <- c(x_star,fun(x_star),deriv(x_star))
  return(mat)
}


update_z <- function(z,x_star,mat){
  if (nrow(mat) <= 2){
    stop("something wrong with info matrix", call. = FALSE)
  }
  index <- which(mat[,1] == x_star)
  if ((index > 1) && (index < nrow(mat))) {
    
  }
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
