# sbr = function (A, z, lambda, control = list(VERB = 0, qinit = c(), explore = "all", K = M))
#
# Minimization of the cost function J(x) = ||z-Ax||_2^2 + lambda ||x||_0
# by the Single Best Replacement (SBR) algorithm
#
# This R program is modified from the Matlab program provided by Charles Soussen at
# http://w3.cran.univ-lorraine.fr/perso/charles.soussen/software.html
#
# Contact information: sunl@uchicago.edu
#
# The Matlab program is a supplementary material to the paper:
#
# From Bernoulli-Gaussian deconvolution to sparse signal restoration,
# by Charles Soussen, Jerome Idier, David Brie, and Junbo Duan.
# IEEE Trans. Signal Process., vol. 59, no. 10, pp. 4572-4584, Oct. 2011.
#
# Contact information: Charles.Soussen@univ-lorraine.fr
#----------------------------------------------------------------------
#
# Version which STORES THE GRAM MATRIX A'A and updates the Cholesky
# factorization of (A'*A)(active_columns,active_columns). The
# insertion/removal tests are done in parallel
#
#
# INPUTS
#
# - z: data signal / regression coefficients, vector of size n
# - A: observation matrix / design matrix of size n x M (if A is sparse, then A'A is sparse)
# - lambda: sparsity parameter (scalar)
# - Optional / control argument VERB = 1 for debug (many prints on screen), 0 by default
# - Optional / control argument qinit to choose a non-zero initialization,
#   vector of indices of non-zero components in xhat for initialization
#   (qinit = c() by default)
# - Optional / control argument explore: strategy of exploration of the culumns a_i
#
#   explore = 'all' (default): both insertions and removals are tested
#   explore = 'remove': version which favors the removals: when a removal
#                       decreases J(x), the "best removal" is selected
#                       and no insertion test is done.
#   The 'remove' option may improve the approximation found with the
#   default option in some cases but the improvement does not occur
#   systematically.
#
# - Optional / control argument K: SBR is stopped when an iterate has K non-zero
#   components.
#   * When lambda > 0, this option shall not be activated since SBR
#     terminates without any stopping condition.
#   * When lambda = 0, SBR coincides with Orthogonal Least Squares (OLS).
#     Here, the K stopping condition is useful to stop OLS.
#
#
# OUTPUTS:
#
#   - x (SBR output after full convergence): sparse vector of size M,
#     with M = ncol(A)
#   - x.sparse: a sparse representation of x, with only the position
#     and value of non-zero components displayed
#   - Jit: vector of size (M+1) or (K+1), with M = ncol(A)
#     Save of the value of J(x) for each best iterate:
#     Jit[i] = the least J-value obtained with an iterate with i-1
#     entries (i >= 1).
#     NOTE: in the version with qinit != c(), Jit is not a full vector (there
#     may be "holes" for early cardinalities since some values Jit[i] are
#     lacking).
#   - FLAG: exit flag.
#     Equal to 0 if the execution terminates correctly,
#     Equal to -1 if an anomaly occurs (numerical instability), i.e.
#     if J(x) < 0 or if there are numerical issues when computing the
#     Cholesky factor of (t(A) %*% A)[active_columns, active_columns]
#
#---------------------------------------------------------------------

require(Matrix)

sbr = function (A, z, lambda, control = list()) {
  lambda = 2 * lambda

	M = ncol(A) # size of vector x

	# Optional argument
	control.default = list(
	VERB = 0,
	qinit = c(),
	explore = "all",
	K = M
	)
	control = modifyList(control.default, control)
	VERB = control$VERB
	qinit = control$qinit
	explore = control$explore
	K = control$K
	K = min(K, min(dim(A)))

	# WE ASSUME THAT THE GRAM MATRIX A'A CAN BE COMPUTED AND STORED ONCE FOR ALL
	AA = t(A) %*% A
	Az = t(A) %*% z # < h_i,z>
	h2 = diag(AA) # ||h_i||^2
	z2 = sum(z^2) # ||z||^2
	Jit = rep(0, length = K + 1)

	# zero-valued initial solution
	q = rep(0, length = M) # sparse vector of size Mx1, initialized to 0
	vecm1 = c() # support = active indices
	lvecm1 = 0       # size of support
	lvecm0 = M       # complement
#	R = c()           # Cholesky factor of the Gram matrix:
                  # R is the upper triangular matrix such that
                  #    AA(active_columns,active_columns) = R'*R
                  # With the notations of [Soussen et al, 2011], R is the
                  # transpose of L_Q (see Appendix C of [Soussen et al, 2011])
	crit = z2        # Jinit = J(0) = ||z||^2
	FLAG = 1         # Termination of SBR if FLAG = 0 (normal) or -1 (anomaly)
	iter = 0
	Jit[1] = crit
  sol_cour = 0    # R'\Az(vecm1);

# USER INITIALIZATION
if (!is.null(qinit)) {
    q[qinit] = 1
    vecm1 = qinit
    lvecm1 = length(vecm1)
    lvecm0 = M - lvecm1
    if (lvecm1) {
        R = chol(AA[vecm1, vecm1])   # factorisation AA = R'R
        sol_cour = solve(t(R), Az[vecm1])
        crit = z2 - t(sol_cour) %*% sol_cour + lambda * lvecm1
        Jit[lvecm1 + 1] = crit
    }
}

if (VERB) {cat("0:", "\t", "\t", crit, "\n")}

# EXPLORATION STRATEGY: ALL COLUMNS OR NOT

# all the columns are explored, starting by
                             # insertion tests

if (explore == 'all') {
	while ((FLAG == 1) & (lvecm1 < K)) {
		# The current value of the cost function J(x) is already stored in crit
        dcritopt = 0;       # J_new - J_old, dcritopt is always the lowest obj value
                            # in current loop
#        mopt = [];          # position of the new column to insert
#        indopt = [];        # location of the index to remove in vecm1
        mouv = 0;           # 1 if insertion, -1 if removal
        vecm0 = which(q == 0); # updated at the beginning of the loop contrary
                            # to vecm1

        #--------------------------------------------------------------
        #-------- I. INSERTION/REMOVAL TESTS     ----------------------
        #--------------------------------------------------------------
        # INSERTION TRIALS h_m (m-th column) for all m
# non empty support
	if (lvecm1) {
        	tmp_vect = solve(t(R), AA[vecm1, vecm0, drop = FALSE])  # size lvecm1 x lvecm0
            F22_vect = h2[vecm0] - colSums(tmp_vect^2);
            dcritnew_vect = lambda - (t(tmp_vect) %*% sol_cour - Az[vecm0])^2 / F22_vect;
        }        else {
        	dcritnew_vect = lambda - Az^2 / h2;
        }   # empty support

	val_min = min(dcritnew_vect);
	ind = which.min(dcritnew_vect);

        if (val_min < dcritopt) {
        	            mouv = 1;
            mopt = vecm0[ind];   # the best insertion
            dcritopt = val_min;  # deltaJ(insert m) = the least deltaJ
        }
        # END INSERTION TESTS

        # In some cases, it is not worth testing removals:
        if ((lambda > 0) & (dcritopt > -lambda) & (lvecm1 > 2)) {
            # REMOVAL TESTS
            # Inversion of R matrix
            R_inv = solve(R);
            # Computation of current amplitudes
            x_loc = R_inv %*% sol_cour;
            dcritnew_vect = x_loc^2 / rowSums(R_inv^2) - lambda;
            	val_min = min(dcritnew_vect);
	ind = which.min(dcritnew_vect);

            if(val_min < dcritopt) {
                mouv = -1;
                indopt = ind;       # index in VECM1 (not in q) of best removal
                dcritopt = val_min;
            }
            # END REMOVAL TESTS
        }

        #--------------------------------------------------------------
        #-------- II. DO INSERTION, REMOVAL OR NOTHING    -------------
        #--------------------------------------------------------------
        if (dcritopt < 0) {       # Modification of the support
            iter = iter + 1;
            crit = crit + dcritopt;

            if ((mouv == 1) & (lvecm1 > 0)) {  # instability test
                tmp = solve(t(R), AA[vecm1, mopt, drop = FALSE]);
                FLAG =  1 - 2 * (h2[mopt] < sum(tmp^2));
            }

            if ((crit < -1e-16) | (FLAG == -1)) {
                cat("ABORT - NUMERICAL INSTABILITY", "\t", "FLAG =", FLAG)
#                disp(sprintf('*** ABORT - NUMERICAL INSTABILITY, FLAG = #d ***\n',FLAG));
                FLAG = -1;
                # Amplitude computation
                x = rep(0, M);
                x[vecm1] = solve(R, sol_cour);
                stop;
            }

            if (mouv == 1) {  # INSERTION
                # Update of R
                if (lvecm1 > 0) {
                    R = cbind(R, tmp);
                    R = rbind(R, c(rep(0, lvecm1), sqrt(h2[mopt] - sum(tmp^2))));
                } else {
                    R = sqrt(h2[mopt]);
                }
                vecm1 = c(vecm1, mopt); # insertion of mopt at the last position
                lvecm0 = lvecm0 - 1;
                lvecm1 = lvecm1 + 1;
                #  note: vecm0 is updated in the beginning of the loop
                if (VERB) {cat(iter, ":", "\t", "+", mopt, "\t", crit, "\n")}
            } else {          # REMOVAL
                # Update of R
              if (indopt < nrow(R)) {
                FF = R[(indopt + 1) : nrow(R), (indopt + 1) : ncol(R)];
                E = R[indopt, (indopt + 1) : ncol(R)];
                R = R[-indopt, ];
                R = R[, -indopt];
                R[indopt : nrow(R), indopt : ncol(R)] = chol(t(FF) %*% FF + E %*% t(E));
              } else {
                R = R[-indopt, ];
                R = R[, -indopt];
              }
                mopt = vecm1[indopt];
                vecm1 = vecm1[-indopt];
                # vecm0 is updated in the beginning of the loop
                lvecm0 = lvecm0 + 1;
                lvecm1 = lvecm1 - 1;
                if (VERB) {cat(iter, ":", "\t", "-", mopt, "\t", crit, "\n")}
            }
            q[mopt] = 1 - q[mopt];   # 1-> 0 or 0 -> 1
            sol_cour = solve(t(R), Az[vecm1]);
            Jit[lvecm1 + 1] = crit;
        } else {
            FLAG = 0;   # Stop of SBR (no more decrease)
        }
	}
}	else if (explore == 'remove') { # the removals are first explored.
                                    # If one of them yields a decrease of
                                    # the cost, insertions are not tested
		    while ((FLAG == 1) & (lvecm1 < K)) {
        # The current value of the cost function J(x) is already stored in crit
        dcritopt = 0;       # J_new - J_old
#        mopt = [];          # position of the new column to insert
#        indopt = [];        # location of the index to remove in vecm1
        mouv = 0;           # 1 if insertion, -1 if removal
        vecm0 = which(q == 0); # updated at the beginning of the loop contrary
                            # to vecm1

        #--------------------------------------------------------------
        #-------- I. INSERTION/REMOVAL TESTS     ----------------------
        #--------------------------------------------------------------
        if ((lambda > 0) & (lvecm1 > 2)) {
            # REMOVAL TESTS
            # Inversion of R matrix
            R_inv = solve(R);
            # Computation of current amplitudes
            x_loc = R_inv %*% sol_cour;
            dcritnew_vect = x_loc^2 / rowSums(R_inv^2) - lambda;

            val_min = min(dcritnew_vect);
			ind = which.min(dcritnew_vect);

            if (val_min < dcritopt) {
                mouv = -1;
                indopt = ind;  # index in VECM1 (not in q) of best removal
                dcritopt = val_min;
            }
            # END REMOVAL TESTS
        }

        # If mouv=-1 the best removal is selected, else insertions are explored
        if (mouv == 0) {
            # INSERTION TRIALS h_m (m-th column) for all m
            if (lvecm1) {           # support of at least 2 elements
                tmp_vect = solve(t(R), AA[vecm1, vecm0, drop = FALSE]);  # size lvecm1 x lvecm0
                F22_vect = h2[vecm0] - colSums(tmp_vect^2);
                dcritnew_vect = lambda - (t(tmp_vect) %*% sol_cour - Az[vecm0])^2 / F22_vect;
            } else {   # empty support
                dcritnew_vect = lambda - Az^2 / h2;
            }

                       val_min = min(dcritnew_vect);
			ind = which.min(dcritnew_vect);

            if(val_min < dcritopt) {
                mouv = 1;
                mopt = vecm0[ind];  # the best insertion
                dcritopt = val_min; # deltaJ(insert m) = the least deltaJ
            }
            # END INSERTION TESTS
        }

        #--------------------------------------------------------------
        #-------- II. DO INSERTION, REMOVAL OR NOTHING    -------------
        #--------------------------------------------------------------
        if (dcritopt < 0) {       # Modification of the support
            iter = iter + 1;
            crit = crit + dcritopt;
            if ((mouv == 1) & (lvecm1 > 0)) {  # instability test
                tmp = solve(t(R), AA[vecm1, mopt, drop = FALSE]);
                FLAG = 1 - 2 * (h2[mopt] < sum(tmp^2));
            }

            if ((crit < -1e-16) | (FLAG == -1)) {
                cat("ABORT - NUMERICAL INSTABILITY", "\t", "FLAG =", FLAG);
                FLAG = -1;
                # Amplitude computation
                x = rep(0, M);
                x[vecm1] = solve(R, sol_cour);
                stop;
            }

            if (mouv == 1) {  # INSERTION
                # Update of R
                if (lvecm1 > 0) {
                    R = cbind(R, tmp);
                    R = rbind(R, c(rep(0, lvecm1), sqrt(h2[mopt] - sum(tmp^2))));
                } else {
                    R = sqrt(h2[mopt]);
                }
                vecm1 = c(vecm1, mopt); # insertion of mopt at the last position
                lvecm0 = lvecm0 - 1;
                lvecm1 = lvecm1 + 1;
                #  note: vecm0 is updated in the beginning of the loop
                if (VERB) {cat(iter, ":", "\t", "+", mopt, "\t", crit, "\n")}
            } else {          # REMOVAL
                # Update of R
              if (indopt < nrow(R)) {
                FF = R[(indopt + 1) : nrow(R), (indopt + 1) : ncol(R)];
                E = R[indopt, (indopt + 1) : ncol(R)];
                R = R[-indopt, ];
                R = R[, -indopt];
                R[indopt : nrow(R), indopt : ncol(R)] = chol(t(FF) %*% FF + E %*% t(E));
              } else {
                R = R[-indopt, ];
                R = R[, -indopt];
              }
              mopt = vecm1[indopt];
                vecm1 = vecm1[-indopt];
                # vecm0 is updated in the beginning of the loop
                lvecm0 = lvecm0 + 1;
                lvecm1 = lvecm1 - 1;
                if (VERB) {cat(iter, ":", "\t", "-", mopt, "\t", crit, "\n")}
            }
            q[mopt] = 1 - q[mopt];   # 1-> 0 or 0 -> 1
            sol_cour = solve(t(R), Az[vecm1]);
            Jit[lvecm1 + 1] = crit;
        } else {
            FLAG = 0;   # Stop of SBR (no more decrease)
        }
    } # WHILE
}

# Final amplitude computation
#
x = rep(0, M);
if (!all(sol_cour == 0)) {
  x[vecm1] = solve(R, sol_cour);
}
x.sparse = cbind(which(x != 0), x[which(x != 0)])
colnames(x.sparse) = c("position", "value")


return(list(x = x, x.sparse = x.sparse, Jit = Jit, FLAG = FLAG))
}
