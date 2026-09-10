#ifndef _HMMTMB_LIKELIHOOD_
#define _HMMTMB_LIKELIHOOD_

#include <algorithm>

//' Forward algorithm, exact or banded
//'
//' The exact forward algorithm accumulates log-likelihood contributions that
//' are each conditional on the entire process history, through the scaled
//' forward variable. The Hessian with respect to anything entering those
//' contributions time point by time point -- the weights of a latent Gaussian
//' field, say -- is therefore dense, which the automatic Laplace approximation
//' cannot afford.
//'
//' The banded version (Fischer, 2026) splits the series into blocks of length
//' bw and starts each block's forward variable from a fixed vector rho at the
//' beginning of the *previous* block. Each block then sees only the 2 * bw
//' observations nearest to it, cross-derivatives beyond that lag vanish, and
//' the Hessian is banded. Every block is traversed twice, so this costs about
//' twice the exact algorithm, and the error decays geometrically in bw because
//' an ergodic Markov chain forgets its initial condition exponentially fast.
//' This mirrors LaMa::forward_g(), so that the two can be checked against
//' each other.
//'
//' The exact algorithm is the special case of one block spanning the whole
//' time series, which is how bw < 2 is handled below.
//'
//' @param prob n x n_states matrix of state-dependent probabilities
//' @param tpm_array Transition probability matrices; entry i governs the
//' transition out of time step i
//' @param delta0 n_ID x n_states matrix of initial distributions
//' @param ID Vector of time series IDs, in blocks of equal value
//' @param bw Bandwidth; < 2 means the exact algorithm
//' @return Log-likelihood
template<class Type>
Type forward_alg(const matrix<Type>& prob,
                 const vector<matrix<Type> >& tpm_array,
                 const matrix<Type>& delta0,
                 const vector<Type>& ID,
                 int bw) {
  int n = prob.rows();
  matrix<Type> phi(1, prob.cols()), rho(1, prob.cols());
  rho.setConstant(Type(1) / Type(prob.cols()));
  Type llk = 0;

  // One time series at a time
  int id = 0;
  for(int start = 0; start < n; id++) {
    int end = start;
    while(end + 1 < n && ID(end + 1) == ID(end)) end++;
    int k = (bw < 2) ? (end - start + 1) : bw;

    // First block, from the model's own initial distribution
    phi = (delta0.row(id).array() * prob.row(start).array()).matrix();
    llk += log(phi.sum());
    phi /= phi.sum();
    for(int i = start + 1; i <= std::min(start + k - 1, end); i++) {
      phi = ((phi * tpm_array(i - 1)).array() * prob.row(i).array()).matrix();
      llk += log(phi.sum());
      phi /= phi.sum();
    }

    // Later blocks, each warmed up from rho over the preceding block
    for(int b = start + k; b <= end; b += k) {
      phi = (rho.array() * prob.row(b - k).array()).matrix();
      phi /= phi.sum();
      for(int i = b - k + 1; i < b; i++) {
        phi = ((phi * tpm_array(i - 1)).array() * prob.row(i).array()).matrix();
        phi /= phi.sum();
      }
      for(int i = b; i <= std::min(b + k - 1, end); i++) {
        phi = ((phi * tpm_array(i - 1)).array() * prob.row(i).array()).matrix();
        llk += log(phi.sum());
        phi /= phi.sum();
      }
    }

    start = end + 1;
  }

  return llk;
}

//' Negative log-density of the penalised coefficients of all smooth terms
//'
//' A smooth carries one or more penalty matrices, combined through possibly
//' fewer smoothing parameters via mgcv's L convention, log(lambda) = L * theta.
//' Penalties of one smooth are consecutive and share a block of coefficients,
//' which is how they are grouped below, and the blocks of S follow in the same
//' order.
//'
//' Two cases, chosen by gmrf. A fixed penalty matrix scaled by one parameter
//' (gmrf = 0) has log|lambda S|+ = Sn * log(lambda) + log|S|+ in closed form,
//' with log|S|+ supplied as data. A precision that depends on the parameters
//' through more than an overall scale (gmrf = 1), as an SPDE field's does,
//' needs its log-determinant differentiated, which density::GMRF does with a
//' sparse Cholesky.
//'
//' @param S Block diagonal matrix holding every penalty, one block each
//' @param coeff_re Coefficients of all smooth terms
//' @param theta Log smoothing parameters
//' @param ncol_re 2 x n_pen; first and last coefficient each penalty applies to
//' @param L n_pen x n_theta matrix with log(lambda) = L * theta
//' @param log_det_S Log generalised determinant of each penalty matrix
//' @param gmrf For each penalty, 1 if its log-determinant is computed here
//' @return Negative log-density of the penalised coefficients
template<class Type>
Type smooth_penalty(const Eigen::SparseMatrix<Type>& S,
                    const vector<Type>& coeff_re,
                    const vector<Type>& theta,
                    const matrix<int>& ncol_re,
                    const matrix<Type>& L,
                    const vector<Type>& log_det_S,
                    const vector<int>& gmrf) {
  int n_pen = ncol_re.cols();
  vector<Type> log_lambda = L * theta;
  Type pen = 0;

  // One smooth at a time, i.e. one run of penalties sharing a coefficient block
  int off = 0;
  for(int i = 0; i < n_pen; ) {
    int last = i;
    while(last + 1 < n_pen && ncol_re(0, last + 1) == ncol_re(0, i)) last++;
    int Sn = ncol_re(1, i) - ncol_re(0, i) + 1;
    vector<Type> this_coeff_re = coeff_re.segment(ncol_re(0, i) - 1, Sn);

    if(gmrf(i) == 0) {
      for(int j = i; j <= last; j++, off += Sn) {
        Eigen::SparseMatrix<Type> this_S = S.block(off, off, Sn, Sn);
        pen += Type(0.5) * Sn * log(2*M_PI) -
          Type(0.5) * log_det_S(j) -
          Type(0.5) * Sn * log_lambda(j) +
          Type(0.5) * exp(log_lambda(j)) * density::GMRF(this_S).Quadform(this_coeff_re);
      }
    } else {
      // Precision of the field, summed over its penalty matrices
      Eigen::SparseMatrix<Type> Q = exp(log_lambda(i)) *
        Eigen::SparseMatrix<Type>(S.block(off, off, Sn, Sn));
      for(int j = i + 1, o = off + Sn; j <= last; j++, o += Sn) {
        Eigen::SparseMatrix<Type> this_S = S.block(o, o, Sn, Sn);
        Q = Q + exp(log_lambda(j)) * this_S;
      }
      off += (last - i + 1) * Sn;
      pen += density::GMRF(Q)(this_coeff_re);
    }

    i = last + 1;
  }

  return pen;
}

#endif
