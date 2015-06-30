# standardization.R
#
# created Jun/30/2015, NU
# last mod Jun/30/2015, NU

#--------------- main function ---------------

phi_mix <- function(object, direct=TRUE) {

  n <- object$info$num.xi
  nl <- factorial(n + 2 - 1) /
    (factorial(2)*factorial(n - 1))

  phi <- matrix(nrow=n + nl, ncol=n + nl)

  phi[seq_len(n), seq_len(n)] <- second_moments(object)
  phi[(n+1):nrow(phi), seq_len(n)] <- cov_xyz_indirect(object)
  phi[seq_len(n), (n+1):nrow(phi)] <- t(phi[(n+1):nrow(phi), seq_len(n)])
  tmp <- matrix(nrow=nl, ncol=nl)
  tmp[lower.tri(tmp, diag=T)] <- cov_xy_indirect(object)
  phi[(n+1):nrow(phi), (n+1):nrow(phi)] <- fill_symmetric(tmp)
  
  phi

}

#--------------- helper functions ---------------

# first moments
mu_group <- function(object, class) {

  mu <- coef(object)[[class]][grep("tau", names(coef(object)[[class]]))]
  mu

}

mu <- function(object) {

  mu_class <- matrix(nrow=object$info$num.xi, ncol=object$info$num.classes)
  for (i in seq_len(object$info$num.classes)) {
    mu_class[,i] <- object$info$w[i] * mu_group(object, paste0("class",i))
  }
  sm <- apply(mu_class, 1, sum)
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

  mu_i <- mu_group(object, class)
  nu_ij <- nu_group(object, class)
  
  mu_ij <- NULL
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      mu_ij <- c(mu_ij, mu_i[i]*mu_i[j] + nu_ij[i, j])
    }
  }
  mu_ij <- matrix(mu_ij, nrow=object$info$num.xi, ncol=object$info$num.xi)
  # TODO: Check if this works for num.xi > 2
  mu_ij

}

# Eq. 31
third_moments_group <- function(object, class="class1") {

  mu_i <- mu_group(object, class)
  mu_ij <- second_moments_group(object, class)

  mu_ijk <- NULL
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        mu_ijk <- c(mu_ijk, mu_ij[i, j]*mu_i[k] + mu_ij[i, k]*mu_i[j] +
          mu_ij[j, k]*mu_i[i] - 2*mu_i[i]*mu_i[j]*mu_i[k])
      }
    }
  }
  mu_ijk <- array(mu_ijk, rep(object$info$num.xi, 3))
  mu_ijk

}

# Eq. 32
fourth_moments_group <- function(object, class="class1") {

  nu_ij <- nu_group(object, class)

  nu_ijkl <- NULL
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        for (l in seq_len(object$info$num.xi)) {
          nu_ijkl <- c(nu_ijkl, nu_ij[i,j]*nu_ij[k,l] + nu_ij[i,k]*nu_ij[j,l] +
            nu_ij[i,l]*nu_ij[j,k])
        }
      }
    }
  }
  nu_ijkl <- array(nu_ijkl, rep(object$info$num.xi, 4))

  mu_i <- mu_group(object, class)
  mu_ij <- second_moments_group(object, class)
  mu_ijk <- third_moments_group(object, class)

  mu_ijkl <- NULL
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        for (l in seq_len(object$info$num.xi)) {
          mu_ijkl <- c(mu_ijkl, nu_ijkl + mu_ijk[i,j,k]*mu_i[l] +
            mu_ijk[i,j,l]*mu_i[k] + mu_ijk[i,k,l]*mu_i[j] +
            mu_ijk[j,k,l]*mu_i[i] - mu_ij[i,j]*mu_i[k]*mu_i[l] -
            mu_ij[i,k]*mu_i[j]*mu_i[l] - mu_ij[i,l]*mu_i[j]*mu_i[k] -
            mu_ij[j,k]*mu_i[i]*mu_i[l] - mu_ij[j,l]*mu_i[i]*mu_i[k] -
            mu_ij[k,l]*mu_i[i]*mu_i[j] + 3*mu_i[i]*mu_i[j]*mu_i[k]*mu_i[l])
        }
      }
    }
  }
  mu_ijkl <- array(mu_ijkl, rep(object$info$num.xi, 4))
  # TOTHINK: How important is the order in these arrays?
  mu_ijkl

}

## noncentral moments of mixture variables

# Eq. 35
second_moments <- function(object) {

  mu_class <- array(dim=c(object$info$num.xi, object$info$num.xi,
                    object$info$num.classes))
  for (i in seq_len(object$info$num.classes)) {
    mu_class[,,i] <- object$info$w[i] * second_moments_group(object, paste0("class",i))
  }
  sm <- apply(mu_class, 1:2, sum)
  sm
}

third_moments <- function(object) {

  mu_class <- array(dim=c(object$info$num.xi, object$info$num.xi,
                    object$info$num.xi, object$info$num.classes))
  for (i in seq_len(object$info$num.classes)) {
    mu_class[,,,i] <- object$info$w[i] * third_moments_group(object, paste0("class",i))
  }
  sm <- apply(mu_class, 1:3, sum)
  sm
}

fourth_moments <- function(object) {

  mu_class <- array(dim=c(object$info$num.xi, object$info$num.xi,
                    object$info$num.xi, object$info$num.xi,
                    object$info$num.classes))
  for (i in seq_len(object$info$num.classes)) {
    mu_class[,,,,i] <- object$info$w[i] * fourth_moments_group(object, paste0("class",i))
  }
  sm <- apply(mu_class, 1:4, sum)
  sm
}

## central moments of mixture variables

# Eq. 36
second_moments_central <- function(object) {

  mu_i <- mu(object)
  mu_ij <- second_moments(object)

  nu_ij <- matrix(nrow=object$info$num.xi, ncol=object$info$num.xi)
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      nu_ij[i, j] <- mu_ij[i, j] - mu_i[i]*mu_i[j]
    }
  }
  nu_ij
}


# Eq. 37
third_moments_central <- function(object) {

  mu_i <- mu(object)
  mu_ij <- second_moments(object)
  mu_ijk <- third_moments(object)

  nu_ijk <- array(dim=c(object$info$num.xi, object$info$num.xi,
                    object$info$num.xi))
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        nu_ijk[i,j,k] <- mu_ijk[i,j,k] - mu_ij[i,j]*mu_i[k] -
          mu_ij[i,k]*mu_i[j] - mu_ij[j,k]*mu_i[i] + 2*mu_i[i]*mu_i[j]*mu_i[k]
      }
    }
  }
  nu_ijk

}

# Eq. 38
fourth_moments_central <- function(object) {

  mu_i <- mu(object)
  mu_ij <- second_moments(object)
  mu_ijk <- third_moments(object)
  mu_ijkl <- fourth_moments(object)

  nu_ijkl <- array(dim=c(object$info$num.xi, object$info$num.xi,
                    object$info$num.xi, object$info$num.xi))
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        for (l in seq_len(object$info$num.xi)) {
          nu_ijkl[i,j,k,l] <- mu_ijkl[i,j,k,l] - mu_ijk[i,j,k]*mu_i[l] -
            mu_ijk[i,j,l]*mu_i[k] - mu_ijk[i,k,l]*mu_i[j] -
            mu_ijk[j,k,l]*mu_i[i] + mu_ij[i,j]*mu_i[k]*mu_i[l] +
            mu_ij[i,k]*mu_i[j]*mu_i[l] + mu_ij[i,l]*mu_i[j]*mu_i[k] +
            mu_ij[j,k]*mu_i[i]*mu_i[l] + mu_ij[j,l]*mu_i[i]*mu_i[k] +
            mu_ij[k,l]*mu_i[i]*mu_i[j] - 3*mu_i[i]*mu_i[j]*mu_i[k]*mu_i[l]
        }
      }
    }
  }
  nu_ijkl
}

## indirect approach

# Eq. 27: variance of product terms
var_xy_indirect <- function(object) {

  mu_i <- mu(object)
  nu_ij <- second_moments_central(object)
  nu_ijk <- third_moments_central(object)
  nu_ijkl <- fourth_moments_central(object)

  var_xy <- matrix(nrow=object$info$num.xi, ncol=object$info$num.xi)
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      var_xy[i, j] <- mu_i[i]^2*nu_ij[j,j] + mu_i[j]^2*nu_ij[i,i] +
        2*mu_i[i]*mu_i[j]*nu_ij[i,j] - nu_ij[i,j]^2 +
        2*mu_i[i]*nu_ijk[i,j,j] + 2*mu_i[j]*nu_ijk[i,i,j] +
        nu_ijkl[i,i,j,j]
    }
  }
  out <- var_xy[lower.tri(var_xy, diag=TRUE)]
  out
}
## --> function not needed since it is part of cov_xy_indirect, of course!

# Eq. 28: covariance of two product terms
# cov_xy_indirect <- function(object) {
# 
#   mu_i <- mu(object)
#   nu_ij <- second_moments_central(object)
#   nu_ijk <- third_moments_central(object)
#   nu_ijkl <- fourth_moments_central(object)
# 
#   cov_xy <- array(dim=c(object$info$num.xi, object$info$num.xi,
#     object$info$num.xi, object$info$num.xi))
#   for (i in seq_len(object$info$num.xi)) {
#     for (j in seq_len(object$info$num.xi)) {
#       for (k in seq_len(object$info$num.xi)) {
#         for (l in seq_len(object$info$num.xi)) {
#           cov_xy[i,j,k,l] <- mu_i[i]*mu_i[k]*nu_ij[j,l] +
#             mu_i[i]*mu_i[l]*nu_ij[j,k] + mu_i[j]*mu_i[k]*nu_ij[i,l] +
#             mu_i[j]*mu_i[l]*nu_ij[i,k] - nu_ij[i,j]*nu_ij[k,l] +
#             mu_i[i]*nu_ijk[j,k,l] + mu_i[j]*nu_ijk[i,k,l] +
#             mu_i[k]*nu_ijk[i,j,l] + mu_i[l]*nu_ijk[i,j,k] +
#             nu_ijkl[i,j,k,l]
#         }
#       }
#     }
#   }
#   m <- matrix(c(cov_xy), nrow=object$info$num.xi^2,
#     ncol=object$info$num.xi^2)
#   m_ii <- matrix(diag(m), nrow=object$info$num.xi,
#     ncol=object$info$num.xi)
#   m_ii <- m_ii[lower.tri(m_ii, diag=TRUE)]
#   m_i <- m[upper.tri(m)]
#   cov_xy
# }

cov_xy_indirect <- function(object) {

  mu_i <- mu(object)
  nu_ij <- second_moments_central(object)
  nu_ijk <- third_moments_central(object)
  nu_ijkl <- fourth_moments_central(object)

  cov_xy <- NULL
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        for (l in seq_len(object$info$num.xi)) {

          if ((i == j & j == k & k == l) || (i == j & k == l & i < k) || (j
          == k & k == l & i < j) || (i == l & j == k & i < j) || (i == j &
          j == k & i < l)) {

            cov_xy <- c(cov_xy, mu_i[i]*mu_i[k]*nu_ij[j,l] +
              mu_i[i]*mu_i[l]*nu_ij[j,k] + mu_i[j]*mu_i[k]*nu_ij[i,l] +
              mu_i[j]*mu_i[l]*nu_ij[i,k] - nu_ij[i,j]*nu_ij[k,l] +
              mu_i[i]*nu_ijk[j,k,l] + mu_i[j]*nu_ijk[i,k,l] +
              mu_i[k]*nu_ijk[i,j,l] + mu_i[l]*nu_ijk[i,j,k] +
              nu_ijkl[i,j,k,l])
          }
        }
      }
    }
  }
  cov_xy
}

# Eq. 29: covariance of product term with third variable
cov_xyz_indirect <- function(object) {

  mu_i <- mu(object)
  nu_ij <- second_moments_central(object)
  nu_ijk <- third_moments_central(object)

  cov_ijk <- array(dim=c(object$info$num.xi, object$info$num.xi,
    object$info$num.xi))
  for (i in seq_len(object$info$num.xi)) {
    for (j in seq_len(object$info$num.xi)) {
      for (k in seq_len(object$info$num.xi)) {
        cov_ijk[i,j,k] <- mu_i[i]*nu_ij[j,k] + mu_i[j]*nu_ij[i,k] + nu_ijk[i,j,k]
      }
    }
  }
  cov_ijk
  out <- c(apply(cov_ijk, 3, function(x) x[lower.tri(x, diag=TRUE)]))
  out
}


