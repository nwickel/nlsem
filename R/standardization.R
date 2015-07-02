# standardization.R
#
# created Jun/30/2015, NU
# last mod Jul/01/2015, NU

#--------------- main functions ---------------

phi_mix <- function(object, direct=TRUE) {

  n <- object$info$num.xi
  nl <- factorial(n + 2 - 1) / (factorial(2)*factorial(n - 1))

  if (direct) {
  
    phi.l <- list()
    for (class in paste0("class", seq_len(object$info$num.classes))) {

      phi <- matrix(nrow=n + nl, ncol=n + nl)

      phi[seq_len(n), seq_len(n)] <- second_moments_group(object, class)
      phi[(n+1):nrow(phi), seq_len(n)] <- cov_xyz(object,
        direct=direct)[[class]]
      phi[seq_len(n), (n+1):nrow(phi)] <- t(phi[(n+1):nrow(phi), seq_len(n)])
      tmp <- matrix(nrow=nl, ncol=nl)
      tmp[lower.tri(tmp, diag=T)] <- cov_xy(object, direct=direct)[[class]]
      phi[(n+1):nrow(phi), (n+1):nrow(phi)] <- fill_symmetric(tmp)

      phi.l[[class]] <- phi
    }
    phi.l

  } else {

    phi <- matrix(nrow=n + nl, ncol=n + nl)

    phi[seq_len(n), seq_len(n)] <- second_moments(object)
    phi[(n+1):nrow(phi), seq_len(n)] <- cov_xyz(object, direct=direct)
    phi[seq_len(n), (n+1):nrow(phi)] <- t(phi[(n+1):nrow(phi), seq_len(n)])
    tmp <- matrix(nrow=nl, ncol=nl)
    tmp[lower.tri(tmp, diag=T)] <- cov_xy(object, direct=direct)
    phi[(n+1):nrow(phi), (n+1):nrow(phi)] <- fill_symmetric(tmp)

    phi
  }
}

standardize <- function(object) {

  pars <- coef(object)

  # direct approach
  if (is.list(pars)) {

    gamma.dot <- list()
    for (class in paste0("class", seq_len(object$info$num.classes))) {
      l <- pars[[class]][grep("Gamma", names(pars[[class]]))]
      nl <- pars[[class]][grep("Omega", names(pars[[class]]))]
      gamma <- c(l, nl)

      phi <- phi_mix(object, direct=TRUE)[[class]]

      if (length(gamma) == nrow(phi)) {

        phi00 <- var_eta(parameters=gamma, phi=phi, psi=pars[[class]][grep("Psi",
          names(pars[[class]]))])

        Omega <- matrix(nrow=object$info$num.xi, ncol=object$info$num.xi)
        Omega[lower.tri(Omega, diag=TRUE)] <- nl
        Omega <- fill_symmetric(Omega)

        Omega2 <- Omega
        diag(Omega2) <- 2*diag(Omega)

        mu.i <- mu_group(object, class)

        # linear effects
        phi.l <- phi[seq_len(object$info$num.xi), seq_len(object$info$num.xi)]

        gamma.l <- (l + Omega2 %*% mu.i) * sqrt(diag(phi.l) / phi00)

        # nonlinear effects
        gamma.nl <- matrix(nrow=object$info$num.xi, ncol=object$info$num.xi)
        for (i in seq_len(object$info$num.xi)) {
          for (j in seq_len(object$info$num.xi)) {

            gamma.nl[i,j] <- Omega[i,j] * sqrt(phi.l[i,i]*phi[j,j] / phi00)
          }
        }

        gamma.all <- c(gamma.l, gamma.nl[upper.tri(gamma.nl, diag=TRUE)])
        names(gamma.all) <- names(gamma)
        gamma.dot[[class]] <- gamma.all

      } else if (object$model.class == "singleClass") {

      # TODO need to extend phi_mix for this one!

      } else {



      }

    }


  # indirect approach
  } else {

    # pars <- pars[[1]]
    # TOTHINK Depending on how I implement it this needs to be added.

    l <- pars[grep("Gamma", names(pars))]
    nl <- pars[grep("Omega", names(pars))]
    gamma <- c(l, nl)

    phi <- phi_mix(object, direct=FALSE)

    if (length(gamma) == nrow(phi)) {

      phi00 <- var_eta(parameters=gamma, phi=phi, psi=pars[grep("Psi",
        names(pars))])

      Omega <- matrix(nrow=object$info$num.xi, ncol=object$info$num.xi)
      Omega[lower.tri(Omega, diag=TRUE)] <- nl
      Omega <- fill_symmetric(Omega)

      Omega2 <- Omega
      diag(Omega2) <- 2*diag(Omega)

      mu.i <- pars[grep("tau", names(pars))]

      # linear effects
      phi.l <- phi[seq_len(object$info$num.xi), seq_len(object$info$num.xi)]

      gamma.l <- (l + Omega2 %*% mu.i) * sqrt(diag(phi.l) / phi00)

      # nonlinear effects
      gamma.nl <- matrix(nrow=object$info$num.xi, ncol=object$info$num.xi)
      for (i in seq_len(object$info$num.xi)) {
        for (j in seq_len(object$info$num.xi)) {

          gamma.nl[i,j] <- Omega[i,j] * sqrt(phi.l[i,i]*phi[j,j] / phi00)
        }
      }

      gamma.dot <- c(gamma.l, gamma.nl[upper.tri(gamma.nl, diag=TRUE)])
      names(gamma.dot) <- names(gamma)

    } else {


    }
  }

  gamma.dot
}

#--------------- helper functions ---------------

# first moments
mu_group <- function(object, class) {

  mu <- coef(object)[[class]][grep("tau", names(coef(object)[[class]]))]
  mu

}

mu <- function(object) {

  mu.class <- matrix(nrow=object$info$num.xi, ncol=object$info$num.classes)
  for (i in seq_len(object$info$num.classes)) {
    mu.class[,i] <- object$info$w[i] * mu_group(object, paste0("class",i))
  }
  sm <- apply(mu.class, 1, sum)
  sm
}

## central moments
nu_group <- function(object, class) {

  nu <- matrix(nrow=object$info$num.xi, ncol=object$info$num.xi)
  nu[lower.tri(nu, diag=TRUE)] <- coef(object)[[class]][grep("Phi", names(coef(object)[[class]]))]
  nu <- fill_symmetric(nu)
  nu

}

## noncentral class-specific moments

# Eq. 30
second_moments_group <- function(object, class="class1") {

  mu.i <- mu_group(object, class)
  nu.ij <- nu_group(object, class)
  
  mu.ij <- matrix(nrow=object$info$num.xi, ncol=object$info$num.xi)
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      mu.ij[i,j] <- mu.i[i]*mu.i[j] + nu.ij[i, j]
    }
  }
  mu.ij

}

# Eq. 31
third_moments_group <- function(object, class="class1") {

  mu.i <- mu_group(object, class)
  mu.ij <- second_moments_group(object, class)

  mu.ijk <- array(dim=rep(object$info$num.xi, 3))
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        mu.ijk[i,j,k] <- mu.ij[i, j]*mu.i[k] + mu.ij[i, k]*mu.i[j] +
          mu.ij[j, k]*mu.i[i] - 2*mu.i[i]*mu.i[j]*mu.i[k]
      }
    }
  }
  mu.ijk

}

# Eq. 32
fourth_moments_group <- function(object, class="class1") {

  nu.ij <- nu_group(object, class)

  # Eq. 34
  nu.ijkl <- array(dim=rep(object$info$num.xi, 4))
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        for (l in seq_len(object$info$num.xi)) {
          nu.ijkl[i,j,k,l] <- nu.ij[i,j]*nu.ij[k,l] + nu.ij[i,k]*nu.ij[j,l] +
            nu.ij[i,l]*nu.ij[j,k]
        }
      }
    }
  }

  mu.i <- mu_group(object, class)
  mu.ij <- second_moments_group(object, class)
  mu.ijk <- third_moments_group(object, class)

  mu.ijkl <- array(dim=rep(object$info$num.xi, 4))
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        for (l in seq_len(object$info$num.xi)) {
          mu.ijkl[i,j,k,l] <- nu.ijkl[i,j,k,l] + mu.ijk[i,j,k]*mu.i[l] +
            mu.ijk[i,j,l]*mu.i[k] + mu.ijk[i,k,l]*mu.i[j] +
            mu.ijk[j,k,l]*mu.i[i] - mu.ij[i,j]*mu.i[k]*mu.i[l] -
            mu.ij[i,k]*mu.i[j]*mu.i[l] - mu.ij[i,l]*mu.i[j]*mu.i[k] -
            mu.ij[j,k]*mu.i[i]*mu.i[l] - mu.ij[j,l]*mu.i[i]*mu.i[k] -
            mu.ij[k,l]*mu.i[i]*mu.i[j] + 3*mu.i[i]*mu.i[j]*mu.i[k]*mu.i[l]
        }
      }
    }
  }
  mu.ijkl

}

## noncentral moments of mixture variables

# Eq. 35
second_moments <- function(object) {

  mu.class <- array(dim=c(object$info$num.xi, object$info$num.xi,
                    object$info$num.classes))
  for (i in seq_len(object$info$num.classes)) {
    mu.class[,,i] <- object$info$w[i] * second_moments_group(object, paste0("class",i))
  }
  sm <- apply(mu.class, 1:2, sum)
  sm
}

third_moments <- function(object) {

  mu.class <- array(dim=c(object$info$num.xi, object$info$num.xi,
                    object$info$num.xi, object$info$num.classes))
  for (i in seq_len(object$info$num.classes)) {
    mu.class[,,,i] <- object$info$w[i] * third_moments_group(object, paste0("class",i))
  }
  sm <- apply(mu.class, 1:3, sum)
  sm
}

fourth_moments <- function(object) {

  mu.class <- array(dim=c(object$info$num.xi, object$info$num.xi,
                    object$info$num.xi, object$info$num.xi,
                    object$info$num.classes))
  for (i in seq_len(object$info$num.classes)) {
    mu.class[,,,,i] <- object$info$w[i] * fourth_moments_group(object, paste0("class",i))
  }
  sm <- apply(mu.class, 1:4, sum)
  sm
}

## central moments of mixture variables

# Eq. 36
second_moments_central <- function(object) {

  mu.i <- mu(object)
  mu.ij <- second_moments(object)

  nu.ij <- matrix(nrow=object$info$num.xi, ncol=object$info$num.xi)
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      nu.ij[i, j] <- mu.ij[i, j] - mu.i[i]*mu.i[j]
    }
  }
  nu.ij
}


# Eq. 37
third_moments_central <- function(object) {

  mu.i <- mu(object)
  mu.ij <- second_moments(object)
  mu.ijk <- third_moments(object)

  nu.ijk <- array(dim=c(object$info$num.xi, object$info$num.xi,
                    object$info$num.xi))
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        nu.ijk[i,j,k] <- mu.ijk[i,j,k] - mu.ij[i,j]*mu.i[k] -
          mu.ij[i,k]*mu.i[j] - mu.ij[j,k]*mu.i[i] + 2*mu.i[i]*mu.i[j]*mu.i[k]
      }
    }
  }
  nu.ijk

}

# Eq. 38
fourth_moments_central <- function(object) {

  mu.i <- mu(object)
  mu.ij <- second_moments(object)
  mu.ijk <- third_moments(object)
  mu.ijkl <- fourth_moments(object)

  nu.ijkl <- array(dim=c(object$info$num.xi, object$info$num.xi,
                    object$info$num.xi, object$info$num.xi))
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        for (l in seq_len(object$info$num.xi)) {
          nu.ijkl[i,j,k,l] <- mu.ijkl[i,j,k,l] - mu.ijk[i,j,k]*mu.i[l] -
            mu.ijk[i,j,l]*mu.i[k] - mu.ijk[i,k,l]*mu.i[j] -
            mu.ijk[j,k,l]*mu.i[i] + mu.ij[i,j]*mu.i[k]*mu.i[l] +
            mu.ij[i,k]*mu.i[j]*mu.i[l] + mu.ij[i,l]*mu.i[j]*mu.i[k] +
            mu.ij[j,k]*mu.i[i]*mu.i[l] + mu.ij[j,l]*mu.i[i]*mu.i[k] +
            mu.ij[k,l]*mu.i[i]*mu.i[j] - 3*mu.i[i]*mu.i[j]*mu.i[k]*mu.i[l]
        }
      }
    }
  }
  nu.ijkl
}

# Eq. 28: covariance of two product terms
cov_xy <- function(object, direct=TRUE) {

  if (direct) {

    cov.xy.l <- list()

    for (class in paste0("class", seq_len(object$info$num.classes))) {
      mu.i <- mu_group(object, class)
      nu.ij <- second_moments_group(object, class)
      nu.ijk <- third_moments_group(object, class)
      nu.ijkl <- fourth_moments_group(object, class)

      cov.xy <- NULL
      for (i in seq_len(object$info$num.xi)) {
        for (j in seq_len(object$info$num.xi)) {
          for (k in seq_len(object$info$num.xi)) {
            for (l in seq_len(object$info$num.xi)) {

              if ((i == j & j == k & k == l) || (i == j & k == l & i < k) || (j
              == k & k == l & i < j) || (i == l & j == k & i < j) || (i == j &
              j == k & i < l)) {

                cov.xy <- c(cov.xy, mu.i[i]*mu.i[k]*nu.ij[j,l] +
                  mu.i[i]*mu.i[l]*nu.ij[j,k] + mu.i[j]*mu.i[k]*nu.ij[i,l] +
                  mu.i[j]*mu.i[l]*nu.ij[i,k] - nu.ij[i,j]*nu.ij[k,l] +
                  mu.i[i]*nu.ijk[j,k,l] + mu.i[j]*nu.ijk[i,k,l] +
                  mu.i[k]*nu.ijk[i,j,l] + mu.i[l]*nu.ijk[i,j,k] +
                  nu.ijkl[i,j,k,l])
              }
            }
          }
        }
      }
      cov.xy.l[[class]] <- cov.xy
    }
    cov.xy.l

  } else {

    mu.i <- mu(object)
    nu.ij <- second_moments_central(object)
    nu.ijk <- third_moments_central(object)
    nu.ijkl <- fourth_moments_central(object)

    cov.xy <- NULL
    for (i in seq_len(object$info$num.xi)) {
      for (j in seq_len(object$info$num.xi)) {
        for (k in seq_len(object$info$num.xi)) {
          for (l in seq_len(object$info$num.xi)) {

            if ((i == j & j == k & k == l) || (i == j & k == l & i < k) || (j
            == k & k == l & i < j) || (i == l & j == k & i < j) || (i == j &
            j == k & i < l)) {

              cov.xy <- c(cov.xy, mu.i[i]*mu.i[k]*nu.ij[j,l] +
                mu.i[i]*mu.i[l]*nu.ij[j,k] + mu.i[j]*mu.i[k]*nu.ij[i,l] +
                mu.i[j]*mu.i[l]*nu.ij[i,k] - nu.ij[i,j]*nu.ij[k,l] +
                mu.i[i]*nu.ijk[j,k,l] + mu.i[j]*nu.ijk[i,k,l] +
                mu.i[k]*nu.ijk[i,j,l] + mu.i[l]*nu.ijk[i,j,k] +
                nu.ijkl[i,j,k,l])
            }
          }
        }
      }
    }
    cov.xy
  }
}

# Eq. 29: covariance of product term with third variable
cov_xyz <- function(object, direct=TRUE) {

  if (direct) {

    cov.xyz.l <- list()
      for (class in paste0("class", seq_len(object$info$num.classes))) {

      mu.i <- mu_group(object, class)
      nu.ij <- second_moments_group(object, class)
      nu.ijk <- third_moments_group(object, class)

      cov.ijk <- array(dim=c(object$info$num.xi, object$info$num.xi,
        object$info$num.xi))
      for (i in seq_len(object$info$num.xi)) {
        for (j in seq_len(object$info$num.xi)) {
          for (k in seq_len(object$info$num.xi)) {
            cov.ijk[i,j,k] <- mu.i[i]*nu.ij[j,k] + mu.i[j]*nu.ij[i,k] + nu.ijk[i,j,k]
          }
        }
      }
      cov.ijk
      cov.xyz.l[[class]] <- c(apply(cov.ijk, 3, function(x) x[lower.tri(x, diag=TRUE)]))
    }
    cov.xyz.l

  } else {

    mu.i <- mu(object)
    nu.ij <- second_moments_central(object)
    nu.ijk <- third_moments_central(object)

    cov.ijk <- array(dim=c(object$info$num.xi, object$info$num.xi,
      object$info$num.xi))
    for (i in seq_len(object$info$num.xi)) {
      for (j in seq_len(object$info$num.xi)) {
        for (k in seq_len(object$info$num.xi)) {
          cov.ijk[i,j,k] <- mu.i[i]*nu.ij[j,k] + mu.i[j]*nu.ij[i,k] + nu.ijk[i,j,k]
        }
      }
    }
    cov.ijk
    out <- c(apply(cov.ijk, 3, function(x) x[lower.tri(x, diag=TRUE)]))
    out
  }
}

# Eq. 39: variance of eta

var_eta <- function(parameters, phi, psi) {

  phi_00 <- t(parameters) %*% phi %*% parameters + psi
  phi_00

}

