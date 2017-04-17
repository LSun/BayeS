reglm.prox = function (y, X, lambda, penalty, par = NULL, beta_init = NULL, max_iter = 1e4, vtol = 1e-10, eta = 0.1) {
  n = nrow(X)
  p = ncol(X)
  XtX = t(X) %*% X
  XtXi = solve(XtX)
  c = 1 / (svd(X)$d[p])^2
  cXtXi = c * XtXi
  beta0 = .lm.fit(x = X, y = y)$coef
  if (is.null(beta_init)) {v = rep(0, p)} else {v = XtX %*% (beta0 - beta_init) / c}
  w = lambda / c

  iter_count = 0
  converged = FALSE
  beta = beta0
  obj_vec = obj_val = obj_fun(beta, y, X, lambda, penalty, par)

  while(!converged && iter_count < max_iter) {
  tmp = v - cXtXi %*% v + beta0
  v_new = eta * v + (tmp - prox(z = tmp, w, penalty = penalty, pen.par = par)) * (1 - eta)
  beta_new = beta0 - cXtXi %*% v_new
  obj_val_new = obj_fun(beta_new, y, X, lambda, penalty, par)
  converged = sqrt(sum((obj_val - obj_val_new)^2) / sum(obj_val_new^2)) < vtol
  obj_vec = c(obj_vec, obj_val_new)
  obj_val = obj_val_new
  beta = beta_new
  v = v_new
  iter_count = iter_count + 1
  }

  return(list(beta = beta, obj_vec = obj_vec, obj_val = obj_val, iter_count = iter_count, converged = converged))
}

rss_LS = function (beta, y, X) {
  return(sum((y - X %*% beta)^2))
}

obj_fun = function (beta, y, X, lambda, penalty, par) {
  lmloss = rss_LS(beta, y, X)/2
  reg = lambda * phi(beta, penalty, par)
  obj_val = lmloss + reg
  return(obj_val = obj_val)
}

phi = function (beta, penalty, par) {
  if (penalty == "l1") {
    return(sum(abs(beta)))
  } else {
    if (penalty == "dp") {
      return(sum(log(1 + abs(beta) / par)))
    } else {
      if (penalty == "lq") {
        return(sum(abs(beta)^q))
      } else {
        if (penalty == "l0") {
          return(sum(beta != 0))
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
