tfprox.fixpt = function (z, w = NULL, D, dual.init = NULL, lambda, penalty, pen.par = NULL, objtol = 1e-10, eta = 0.2, v_tol = 1e-10, max_iter = 1000) {
  if (is.null(pen.par)) {pen.par = 1}

  n = length(z)
  m = nrow(D)
  if (is.null(w)) {
    w = rep(1, n)
  }
  WiDt = 1 / w * t(D)
  DWiDt = D %*% WiDt

  snDWiDt = svd(DWiDt)$d[1]
  h = 1 / snDWiDt

  hDWiDt = h * DWiDt
  # IcDWiDt = diag(1, m) - cDWiDt
  hWiDt = h * WiDt

  Dz = D %*% z
  prox_c = lambda / h

  if (is.null(dual.init)) {
    v = rep(0, m)
  } else {
    v = dual.init / h
  }

  beta.mat = c()

  conv = FALSE
  iter = 1
  obj = c()

  v.init = v
  beta.init = z - hWiDt %*% v.init
  obj.init = tfobj(beta.init, z, w, D, penalty, pen.par, lambda)
  objc = obj.init

  while (!conv & (iter <= max.iter)) {
    tmp = v - hDWiDt %*% v + Dz
    v = eta * v + (tmp - prox(z = tmp, w = prox_c, penalty = penalty, pen.par = pen.par)) * (1 - eta)
    beta = z - hWiDt %*% v
    objnew = tfobj(beta, z, w, D, penalty, pen.par, lambda)
    conv = abs(objnew - objc) / objc < objtol
    obj[iter] = objnew
    objc = objnew
    beta.mat = cbind(beta.mat, beta)
    iter <- iter + 1
    #   gamma <- D %*% beta
    #   gamma_new <- prox(z = gamma + rho, w = prox_c, penalty, pen.par)
    #   rho <- rho + (gamma_new - gamma)
    #   beta <- z - hWinvDt %*% rho
    # gamma <- prox(z = rho + D %*% beta, w = prox_c, penalty, pen.par)
    # beta <- WhDtDinvWz - hWhDtDinvDt %*% (rho - gamma)
    # rho <- rho + (D %*% beta - gamma)
    # obj <- c(obj, tfobj(beta, z, w, D, penalty, pen.par, lambda))
    # conv <- abs((obj[iter + 1] - obj[iter]) / obj[iter]) < objtol
  }


#
#   vfxpt = function (v) {
#     tmp = v - cDWiDt %*% v + Dz
#     vnew = eta * v + (tmp - prox(z = tmp, w = w, penalty = penalty, pen.par = pen.par)) * (1 - eta)
#     return(vnew)
#   }
#
#   v_iter <- SQUAREM::squarem(par = rep(0, m), fixptfn = vfxpt, control = list(maxiter = max_iter))
#   v = v_iter$par
#   converged = v_iter$convergence
#   iter_count = v_iter$fpevals

  # if (is.null(v.init)) {v = rep(0, m)}
  # v = v.init
  # betahat = z
  # iter_count = 0
  # converged = FALSE
  # obj_val = lambda * sum(abs(D %*% betahat))
  #
  # while(!converged && iter_count < max_iter) {
  # tmp = v - cDWiDt %*% v + Dz
  # vnew = eta * v + (tmp - prox(z = tmp, w = w, penalty = penalty, pen.par = pen.par)) * (1 - eta)
  # betahat_new = z - cWiDt %*% vnew
  # obj_val = c(obj_val, sum((z - betahat_new)^2) / 2 + lambda * sum(abs(D %*% betahat_new)))
  # betahat_diff = mean((betahat_new - betahat)^2) / mean(betahat_new^2)
  # converged = betahat_diff <= v_tol
  # betahat = betahat_new
  # # converged = max(abs(v - vnew)) <= v_tol
  # # converged = sum((v - vnew)^2) / (sum(vnew^2) + v_tol) <= v_tol
  # v = vnew
  # iter_count = iter_count + 1
  # }

  # beta_hat = z - cWiDt %*% v

  return(list(beta = beta, v = v, beta.mat = beta.mat, obj.vec = obj, iter = iter - 1, converged = conv, beta.init = beta.init, v.init = v.init, obj.init = obj.init))
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


prox = function(z, w, penalty, pen.par) {
  if (penalty == "l1") {
    prox_z = sign(z) * pmax(abs(z) - w, 0)
    return(prox_z)
  } else {
    if (penalty == "dp") {

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
