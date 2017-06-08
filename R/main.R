dlasso  <-  function(x,
                     y,
                     s  =  1           ,
                     intercept = FALSE ,
                     c  = 1            ,
                     adp = TRUE        ,
                     lambda = NULL     ,
                     split  = 50       ,
                     maxIter = 500     ,
                     adj = 1.1         ,
                     lowlambda = 10^-3 ,
                     digit = 5         ,
                     cauchy = FALSE    ,
                     force = 'auto'    ,
                     trace = FALSE)
{
  if (intercept)
    x   = cbind(1, x)

  n     =  length(y)
  ncx   =  ncol(x)
  if(force == 'auto'){
    force = sqrt(n)>2*log(ncx)
  }
  #print(sqrt(n)/2/log(ncx))

  if (cauchy) {
    fn1  = Sigma3
    fn2  = Sigma4
  } else{
    fn1  = Sigma1
    fn2  = Sigma2
  }
  if (all(lambda <= 0) | all(is.null(lambda))) {
    lmax   =  max(abs(cov(y, x))) * adj * n
    lambda =  lseq(lowlambda , lmax  , length.out  =  split)
  } else{
    split  = length(lambda)
  }
  s  = ifelse(adp, 1, s)
  ls =  length(s)
  ##### initializing ...
  b  =  sqrt(n) * MASS::ginv(t(x) %*% x) %*% t(x) %*% y
  tmp2 = tmp = b
  outputMat  =  matrix(0, ncol  =  ncx + 6 , nrow =  split * ls)
  counter = 1
  if (trace)
    cat('\r\n\r\n Ind|Max \t lambda \t s \t AICc \t\t GIC \t BIC \t GCV \t Converge \t Error')
  time0 = Sys.time()
  for (ss in s) {
    active  =  as.vector(abs(b) > 10 ^ -digit)
    tmp = b#/2
    for (rep in seq(split)) {
      l =  lambda[rep]
      iter = err = 1
      while (any(active)  && (err >  10 ^ -(digit + 0))) {
        newx =  x[, active, drop = FALSE]
        if (adp) {
          inV  =  (t(newx) %*% newx)  + n * l * fn1(
            b  =  b[active] ,
            s  =  (abs(tmp2[active]) /
                     sqrt(n)) ^ c,
            c = c,
            digit = digit
          )
        } else{
          inV  =  (t(newx) %*% newx)  + n * l * fn1(
            b  =  b[active] ,
            s  =  ss ,
            c = c,
            digit = digit
          )
        }
        xInv =  chol2inv(chol(inV, tol = 10 ^ -25))
        b[active]  =  (xInv %*% t(newx) %*% y)
        ###
        if (counter > 1 && force) {
          b[active] = compare(b[active], tmp)
        }
        ###
        active  =  as.vector(abs(b) > 10 ^ -digit)
        #b[!active] = 0
        err  =  max(abs(tmp - b))
        tmp  =  b
        iter =  iter + 1
        if (iter >=  maxIter) {
          cat('\n not converged in ' ,
              maxIter,
              ' iterations. err=',
              err,
              '\n')
          break
        }
      }
      e   = sum((y - x[, active, drop = FALSE] %*% b[active]) ^ 2)
      if (sum(active) > 0) {
        nxx = x[, active, drop = FALSE]
        if (adp) {
          inVg = (t(nxx) %*% nxx) + n * l * fn2(
            b  =  b[active] ,
            s  =  abs(tmp2[active] /
                        sqrt(n)) ^ c ,
            c = c  ,
            digit = digit
          )
        } else{
          inVg = (t(nxx) %*% nxx) + n * l * fn2(
            b  =  b[active] ,
            s  =  ss ,
            c = c,
            digit = digit
          )
        }
        df  =  sum(diag(nxx %*% solve(inVg, tol = 10 ^ -25)  %*% t(nxx)))
      } else{
        df  = 0
      }
      tmp2 = b
      tau  = 1
      ll   = n * log(2 * pi * tau) + e / tau
      aicc = ll + 2 *  n * df * (df + 1) / (n - df - 1) + 2 * df
      bic  = ll + log(n) * df
      gcv  = ll / ((1 - df / n) ^ 2)
      gic  = ll  + 2 * df

      if (trace) {
        cat(
          '\r\n',
          format(counter, digits = 4),
          '|',
          format(split * ls, digits = 3, nsmall = 0),
          '\t ',
          format(l  , digits  = 4, nsmall = 2),
          '\t',
          format(ss , digits  = 4, nsmall = 2),
          '\t',
          format(aicc, digits  = 4, nsmall = 2),
          '\t',
          format(gic, digits  = 4, nsmall = 2),
          '\t',
          format(bic, digits  = 4, nsmall = 2),
          '\t',
          format(gcv, digits  = 4, nsmall = 2),
          '\t',
          format(iter, digits = 4, nsmall = 0),
          ' \t ',
          format(err, digits  = 4, nsmall = 2)
        )
      } else{
        cat('\r ', counter, '|', split * ls, rep(' ', 20))
      }
      outputMat[counter,]  =  c(l,  ss, aicc, gic, bic, gcv,
                                round(b, digit - 1))
      colnames(outputMat)  =  c(
        'lambda',
        's',
        'AICc',
        'GIC',
        'BIC',
        'GCV',
        colnames(x, do.NULL = FALSE, prefix = 'X.')
      )
      counter  = counter + 1
    }
  }
  cat('\n Executed in ', Sys.time() - time0, '(s) \n')
  output = outputMat
  class(output) = 'dlasso'
  return(output)
}


plot.dlasso <- function (x,
                         label = FALSE,
                         cex.lab = 1,
                         all = TRUE,
                         ...) {
  #x = x[-1,]
  object = x
  label = abs(label)
  if (all) {
    # AICc
    plot(
      object[, 1, drop = FALSE],
      object[, 3, drop = FALSE],
      type = 'l',
      ylab = 'AICc',
      xlab = 'Lambda',
      col = 'gray'
    )
    abline(
      v = object[mean(which(object[, 3, drop = FALSE] == min(object[, 3, drop =
                                                                      FALSE]))), 1] ,
      col = 2,
      lty = 2,
      lwd = 2
    )
    # GIC
    plot(
      object[, 1, drop = FALSE],
      object[, 4, drop = FALSE],
      type = 'l',
      ylab = 'GIC',
      xlab = 'Lambda',
      col = 'gray'
    )
    abline(
      v = object[mean(which(object[, 4, drop = FALSE] == min(object[, 4, drop =
                                                                      FALSE]))), 1] ,
      col = 3,
      lty = 3,
      lwd = 2
    )
    # BIC
    plot(
      object[, 1, drop = FALSE],
      object[, 5, drop = FALSE],
      type = 'l',
      ylab = 'BIC',
      xlab = 'Lambda',
      col = 'gray'
    )
    abline(
      v = object[mean(which(object[, 5, drop = FALSE] == min(object[, 5, drop =
                                                                      FALSE]))), 1] ,
      col = 4,
      lty = 4,
      lwd = 2
    )
    # GCV
    plot(
      object[, 1, drop = FALSE],
      object[, 6, drop = FALSE],
      type = 'l',
      ylab = 'GCV',
      xlab = 'Lambda',
      col = 'gray'
    )
    abline(
      v = object[mean(which(object[, 6, drop = FALSE] == min(object[, 6, drop =
                                                                      FALSE]))), 1] ,
      col = 5,
      lty = 5,
      lwd = 2
    )
  }
  coef1 = (object[, -c(1:6)])
  s1 <- apply(abs(coef1), 1, sum)
  if (max(s1) != 0)
    ms1 = max(s1)
  else
    ms1 = 1
  s1 = s1 / ms1
  matplot(
    s1,
    coef1,
    type = 'l',
    ylab = 'Coefficients',
    xlab = '|beta|/max|beta|',
    xlim = c(0, 1 + label) ,
    col = c(1:(ncol(object) - 4))[-2],
    ...
  )
  abline(h = 0, lty = 1, col = 'gray')
  if (label)
    text(
      x = 1 ,
      #+ label/2,
      y = object[1, -c(1:6)],
      labels = colnames(object)[-c(1:6)],
      cex = cex.lab,
      xpd = TRUE,
      pos = 4
    )
  if (1 == 1) {
    abline(
      v = s1[which(object[, 4, drop = FALSE] == min(object[, 4, drop = FALSE]))][1],
      col = 2,
      lty = ifelse(all, 2, 1),
      lwd = 1
    )
    if (all) {
      abline(
        v = s1[which(object[, 3, drop = FALSE] == min(object[, 3, drop = FALSE]))][1],
        col = 3,
        lty = 3,
        lwd = 1
      )
      abline(
        v = s1[which(object[, 5, drop = FALSE] == min(object[, 5, drop = FALSE]))][1],
        col = 4,
        lty = 4,
        lwd = 1
      )
      abline(
        v = s1[which(object[, 6, drop = FALSE] == min(object[, 6, drop = FALSE]))][1],
        col = 5,
        lty = 5,
        lwd = 1
      )
      legend(
        'topleft',
        legend = c('min.GIC', 'min.AICc', 'min.BIC', 'min.GCV'),
        lty = 2:5,
        col = 2:5,
        lwd = 1
      )
    }
  }
  invisible()
}


coef.dlasso <- function(object , ...) {
  #object = object[-1,]
  B_AICc  = object[mean(which(object[, 3, drop = FALSE] == min(object[, 3, drop =
                                                                        FALSE]))) + 0  , -c(1:6)]
  B_GIC  = object[mean(which(object[, 4, drop = FALSE] == min(object[, 4, drop =
                                                                       FALSE]))) + 0  , -c(1:6)]
  B_BIC  = object[mean(which(object[, 5, drop = FALSE] == min(object[, 5, drop =
                                                                       FALSE]))) + 0  , -c(1:6)]
  B_GCV  = object[mean(which(object[, 6, drop = FALSE] == min(object[, 6, drop =
                                                                       FALSE]))) + 0  , -c(1:6)]
  #print(apply(object[,3:6],2,function(x){min(x)[1]}))
  result = cbind(
    coef.AICc = B_AICc,
    coef.GIC = B_GIC,
    coef.BIC = B_BIC,
    coef.GCV = B_GCV
  )
  return(result)
}
