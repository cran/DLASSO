.onAttach <- function( lib, pkg ) {
  packageStartupMessage(
    paste0(
      '\n ========================================================================',
      '\n If you have any question about this package and corresponded paper use ',
      '\n      hamedhaseli@gmail.com  or visit www.hamedhaseli.webs.com           ',
      '\n ========================================================================'
    ),
    domain = NULL,  appendLF = TRUE )
}
#############
#############
fmax = function(x) {
  if (x >= 0) {
    r = max(x)
  }  else{
    r = min(x)
  }
  return(r)
}
lseq <- function(from = 1,
                 to = 5,
                 length.out = 6 ,
                 adj = 1) {
  r = exp(seq(log(from) / adj, log(to), length.out = length.out))
  return(r)
}
compare <- function(x1, x0) {
  po = which(x1 > 0)
  ne = which(x1 < 0)
  x1[po] = pmin(x1[po], x0[po])
  x1[ne] = pmax(x1[ne], x0[ne])
  return(x1)

}
dec = function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else{
    return(0)
  }
}

kfold <- function(Nobs, K = 5) {
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs * seq(1, K - 1) / K)
  k <- matrix(c(0, rep(k, each = 2), Nobs), ncol = 2, byrow = TRUE)
  k[, 1] <- k[, 1] + 1
  l <- lapply(seq.int(K), function(x, k, d)
    list(train = d[!(seq(d) %in% seq(k[x, 1], k[x, 2]))],
         test = d[seq(k[x, 1], k[x, 2])]), k = k, d = id)
  return(l)
}
erfa3 <- function(x = 1, s = 1) {
  x = x / s
  y = tanh(39 * x / sqrt(4 * pi) - 111 / 2 * atan(35 * x / 111 / sqrt(pi)))
  return(y)
}


Sigma1  =  function(b, s, c, digit) {
  lb   = length(b)
  bsc  = (b / s) ^ c
  r    =  (2 * pnorm(bsc, 0, 1 / sqrt(2)) - 1 + 2 * c * bsc * dnorm(bsc, 0, 1 /
                                                                      sqrt(2))) / b
  if (length(r) > 1) {
    r = diag(as.vector(r), nrow = lb, ncol = lb)
  }
  return(r)
}

Sigma2  =  function(b, s = 1, c = 1, digit) {
  lb    = length(b)
  bsc   = (b / s) ^ c
  r     = -2 * c / b * bsc * (2 * c * bsc ^ 2 - c - 1) * (dnorm(bsc, 0, 1 /
                                                                  sqrt(2)))
  r[is.nan(r)] = 0
  if (length(r) > 1) {
    r = diag(as.vector(r), nrow = lb, ncol = lb)
  }
  return(r)
}

Sigma3  =  function(b, s, c, digit) {
  lb   = length(b)
  bsc  = (b / s) ^ c
  r    =  2/pi*(c*bsc/(1+bsc^2)+atan(bsc))/b
  if (length(r) > 1) {
    r = diag(as.vector(r), nrow = lb, ncol = lb)
  }
  return(r)
}


Sigma4  =  function(b, s = 1, c = 1, digit) {
  lb    = length(b)
  bsc   = (b / s) ^ c
  r     = 2*c*bsc*(-(c-1)*(bsc^2)+c+1)/(pi*b*(bsc^2+1)^2)
  r[is.nan(r)] = 0
  if (length(r) > 1) {
    r = diag(as.vector(r), nrow = lb, ncol = lb)
  }
  return(r)
}


tr <- function (m)
{
  total_sum <- 0
  if (is.matrix(m))
  {
    row_count <- nrow(m)
    col_count <- ncol(m)
    if (row_count == col_count)
    {
      total_sum <- sum(diag(m))
      total_sum
    }
    else
    {
      message ('Matrix is not square')
    }
  }
  else
  {
    message('Object is not a matrix')

  }
}
