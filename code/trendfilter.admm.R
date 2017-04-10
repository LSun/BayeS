tfprox.admm = function (z, w, D, penalty, pen.par = NULL, lambda, beta.init = NULL, rho.init = NULL, h = NULL, max.iter = 100, objtol = 1e-6) {
  if (is.null(h)) {
    h = 1 / max(svd(D %*% diag(1 / w) %*% t(D))$d)
  }
  if (is.null(pen.par)) {pen.par = 1}
  prox_c = lambda / h
  W = diag(w)
  WhDtD = W + h * t(D) %*% D
  WhDtDinv = solve(WhDtD)
  WhDtDinvWz = WhDtDinv %*% W %*% z
  hWhDtDinvDt = h * WhDtDinv %*% t(D)
# hWinvDt = h * diag(1 / w) %*% t(D)
  conv = FALSE
  iter = 0
  obj = c()
# beta <- z
# obj = c(obj, tfobj(beta = z, z, w, D, penalty, pen.par, lambda))
  if (is.null(beta.init)) {
    beta <- WhDtDinvWz
  } else {
    beta <- beta.init
  }
  if (is.null(rho.init)) {
    rho <- rep(0, nrow(D))
  } else {
    rho <- rho.init
  }
  obj <- c(obj, tfobj(beta, z, w, D, penalty, pen.par, lambda))
  rho <- rep(0, nrow(D))
  while ((!conv & (iter < max.iter)) | (iter <= 10)) {
    iter <- iter + 1
#   gamma <- D %*% beta
#   gamma_new <- prox(z = gamma + rho, w = prox_c, penalty, pen.par)
#   rho <- rho + (gamma_new - gamma)
#   beta <- z - hWinvDt %*% rho
    gamma <- prox(z = rho + D %*% beta, w = prox_c, penalty, pen.par)
    beta <- WhDtDinvWz - hWhDtDinvDt %*% (rho - gamma)
    rho <- rho + (D %*% beta - gamma)
    obj <- c(obj, tfobj(beta, z, w, D, penalty, pen.par, lambda))
    conv <- abs((obj[iter + 1] - obj[iter]) / obj[iter]) < objtol
  }
  return(list(beta = beta, objval = obj, iter = iter, convergence = conv, z = z, w = w, D = D, penalty = penalty, pen.par = pen.par, lambda = lambda))
}

# compute the trend filter objective function
tfobj = function (beta, z, w, D, penalty, pen.par, lambda) {
  sum(w * (z - beta)^2) / 2 + lambda * phi(D, beta, penalty, pen.par)
}


# compute the penalty function $phi(Dbeta)$
phi = function (D, beta, penalty, pen.par) {
  if (penalty == "l1") {
    return(sum(abs(D %*% beta)))
  } else {
    if (penalty == "double-pareto") {
      a = pen.par[1]
      return(sum(log(1 + abs(D %*% beta) / a)))
    } else {
      if (penalty == "lq") {
        q = pen.par[1]
        return(sum(abs(D %*% beta)^q))
      } else{
        if (penalty == "l0") {
          return(sum((D %*% beta) != 0))
        } else {
          stop ("invalid penalty")
        }
      }
    }
  }
}


# proximal operator for common penalty functions
# options include "l1," "double-pareto," "lq", "l0"
prox = function(z, w, penalty, pen.par) {
  if (penalty == "l1") {
    prox_z = sign(z) * pmax(abs(z) - w, 0)
    return(prox_z)
  } else {
    if (penalty == "double-pareto") {

      a = pen.par[1]
      prox_z = rep(0, length(z))
      I1 = (abs(z) > (2 * sqrt(w) - a))
      I2 = (abs(z) >= w / a)
      zI1I2 = z[I1 & I2]
      zI1I2.abs = abs(zI1I2)
      prox_zI1I2 = (zI1I2.abs - a + sqrt((zI1I2.abs + a)^2 - 4 * w)) / 2
      prox_z[I1 & I2] = prox_zI1I2 * sign(zI1I2)

      if (sqrt(w) > a) {
        zI1nI2 = z[I1 & (!I2)]
        zI1nI2.abs = abs(zI1nI2)
        prox_zI1nI2 = (zI1nI2.abs - a + sqrt((zI1nI2.abs + a)^2 - 4 * w)) / 2
        fprox_zI1nI2 = (prox_zI1nI2 - zI1nI2.abs)^2 / 2 + w * log(1 + prox_zI1nI2 / a)
        f0 = zI1nI2.abs^2 / 2
        prox_zI1nI2 = (fprox_zI1nI2 < f0) * prox_zI1nI2
        prox_z[I1 & !I2] = prox_zI1nI2 * sign(zI1nI2)
      }

      return(prox_z)

    } else {
      if (penalty == "lq") {
        q = pen.par[1]
        bwq = (2 * w * (1 - q))^(1 / (2 - q))
        hwq = bwq + w * q * bwq^(q - 1)

        # threshold rule for lq when 0 < q < 1
        prox_z = rep(0, length(z))
        I = (abs(z) > hwq)

        gamma = abs(z[I])
        get_gamma = SQUAREM::squarem(par = gamma, fixptfn = lqfixpt, z = gamma, w = w, q = q)
        gamma = get_gamma$par

        prox_z[I] = gamma * sign(z[I])

        return(prox_z)
      } else {
        if (penalty == "l0") {
          prox_z = (abs(z) > sqrt(2 * w)) * z
          return(prox_z)
        } else {
          stop ("invalid penalty")
        }
      }
    }
  }
}
