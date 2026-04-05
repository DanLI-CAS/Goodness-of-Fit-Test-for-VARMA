#include <RcppEigen.h>
#include <Eigen/Eigenvalues> 
#include <cmath>
#include <vector> 
#include <algorithm>
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace Rcpp;
using namespace std; 


// =========================================================================
// 计算 VARMA(p,q) 伴随矩阵的最大特征值模长
// =========================================================================
double getMaxEigenvalue(const Map<const MatrixXd>& Mat, int d, int order) {
  if (order == 0) return 0.0;
  
  MatrixXd comp_Mat = MatrixXd::Zero(d * order, d * order);
  comp_Mat.block(0, 0, d, d * order) = Mat; 
  
  if (order > 1) {
    comp_Mat.block(d, 0, d * (order - 1), d * (order - 1)) = MatrixXd::Identity(d * (order - 1), d * (order - 1));
  }
  
  EigenSolver<MatrixXd> es(comp_Mat, false);
  return es.eigenvalues().cwiseAbs().maxCoeff();
}

struct OptResult {
  double llk;
  VectorXd grad;
};

// =========================================================================
//  VARMA(p,q) MLE 
// =========================================================================
OptResult zzr_llk_grad_varma_adjoint_cpp(const Map<MatrixXd> &zt, const Map<VectorXd> &beta, int p, int q) {
  const int n = zt.cols();
  const int d = zt.rows();
  const int d2 = d * d;
  const int pqmax = std::max(p, q);
  
  Map<const MatrixXd> PH(beta.data(), d, d * p);
  Map<const MatrixXd> TH(beta.data() + p * d2, d, d * q);
  

  double max_eval_PH = getMaxEigenvalue(PH, d, p);
  double max_eval_TH = getMaxEigenvalue(TH, d, q);
  const double threshold = 0.999;
  
  if (max_eval_PH >= threshold || max_eval_TH >= threshold) {
    OptResult res;
    double violation = std::max(max_eval_PH, max_eval_TH) - threshold;
    res.llk = 1e10 + 1e5 * violation; 
    res.grad = VectorXd::Zero((p + q) * d2);
    return res;
  }
  

  MatrixXd at = zt; 
  for (int t = pqmax; t < n; ++t) {
    VectorXd curr_a = zt.col(t);
    for (int i = 1; i <= p; ++i) curr_a.noalias() -= PH.block(0, (i - 1) * d, d, d) * zt.col(t - i);
    for (int j = 1; j <= q; ++j) curr_a.noalias() += TH.block(0, (j - 1) * d, d, d) * at.col(t - j);
    at.col(t) = curr_a;
  }
  
  MatrixXd Sigma = (at * at.transpose()) / double(n);
  LLT<MatrixXd> llt(Sigma);
  
  OptResult res;
  res.grad = VectorXd::Zero((p + q) * d2);
  
  if (llt.info() != Success) {
    res.llk = 1e10; 
    return res;
  }
  
  MatrixXd SigmaInv = llt.solve(MatrixXd::Identity(d, d));
  double logDetSigma = 2.0 * llt.matrixLLT().diagonal().array().log().sum();
  res.llk = 0.5 * (n * logDetSigma + n * d + n * d * log(2.0 * M_PI));
  

  MatrixXd grad_PH_mat = MatrixXd::Zero(d, d * p);
  MatrixXd grad_TH_mat = MatrixXd::Zero(d, d * q);
  MatrixXd Lambda = MatrixXd::Zero(d, n);
  
  for (int t = n - 1; t >= pqmax; --t) {
    VectorXd lambda_t = SigmaInv * at.col(t); 
    for (int j = 1; j <= q; ++j) {
      if (t + j < n) lambda_t.noalias() += TH.block(0, (j - 1) * d, d, d).transpose() * Lambda.col(t + j);
    }
    Lambda.col(t) = lambda_t;
    
    for (int i = 1; i <= p; ++i) grad_PH_mat.block(0, (i - 1) * d, d, d).noalias() -= lambda_t * zt.col(t - i).transpose();
    for (int j = 1; j <= q; ++j) grad_TH_mat.block(0, (j - 1) * d, d, d).noalias() += lambda_t * at.col(t - j).transpose();
  }
  
  if (p > 0) res.grad.head(p * d2) = Map<VectorXd>(grad_PH_mat.data(), p * d2);
  if (q > 0) res.grad.tail(q * d2) = Map<VectorXd>(grad_TH_mat.data(), q * d2);
  
  return res;
}

// =========================================================================
//  加权 Bootstrap  VARMA(p,q) 
// =========================================================================
OptResult zzr_llk_grad_varma_boot_adjoint_cpp(const Map<MatrixXd> &zt, const Map<VectorXd> &beta, const Map<VectorXd> &w, int p, int q) {
  const int n = zt.cols();
  const int d = zt.rows();
  const int d2 = d * d;
  const int pqmax = std::max(p, q);
  
  Map<const MatrixXd> PH(beta.data(), d, d * p);
  Map<const MatrixXd> TH(beta.data() + p * d2, d, d * q);
  
  double max_eval_PH = getMaxEigenvalue(PH, d, p);
  double max_eval_TH = getMaxEigenvalue(TH, d, q);
  const double threshold = 0.999;
  
  if (max_eval_PH >= threshold || max_eval_TH >= threshold) {
    OptResult res;
    double violation = std::max(max_eval_PH, max_eval_TH) - threshold;
    res.llk = 1e10 + 1e5 * violation; 
    res.grad = VectorXd::Zero((p + q) * d2);
    return res;
  }
  
  MatrixXd at = zt; 
  MatrixXd Sigma_unweighted = at.col(0) * at.col(0).transpose();
  MatrixXd Sigma_w = w[0] * Sigma_unweighted;
  double W_sum = w.sum();
  

  for (int t = pqmax; t < n; ++t) {
    VectorXd curr_a = zt.col(t);
    for (int i = 1; i <= p; ++i) curr_a.noalias() -= PH.block(0, (i - 1) * d, d, d) * zt.col(t - i);
    for (int j = 1; j <= q; ++j) curr_a.noalias() += TH.block(0, (j - 1) * d, d, d) * at.col(t - j);
    at.col(t) = curr_a;
    
    MatrixXd outer_prod = curr_a * curr_a.transpose();
    Sigma_unweighted += outer_prod;
    Sigma_w += w[t] * outer_prod;
  }
  
  Sigma_unweighted /= double(n);
  LLT<MatrixXd> llt(Sigma_unweighted);
  
  OptResult res;
  res.grad = VectorXd::Zero((p + q) * d2);
  
  if (llt.info() != Success) {
    res.llk = 1e10; 
    return res;
  }
  
  MatrixXd SigmaInv = llt.solve(MatrixXd::Identity(d, d));
  double logDetSigma = 2.0 * llt.matrixLLT().diagonal().array().log().sum();
  double tr_term = (SigmaInv * Sigma_w).trace();
  res.llk = 0.5 * W_sum * (logDetSigma + d * log(2.0 * M_PI)) + 0.5 * tr_term;
  
  MatrixXd M = W_sum * SigmaInv - SigmaInv * Sigma_w * SigmaInv;
  

  MatrixXd grad_PH_mat = MatrixXd::Zero(d, d * p);
  MatrixXd grad_TH_mat = MatrixXd::Zero(d, d * q);
  MatrixXd Lambda = MatrixXd::Zero(d, n);
  
  for (int t = n - 1; t >= pqmax; --t) {
    VectorXd a_bar_t = (M / double(n) + w[t] * SigmaInv) * at.col(t);
    VectorXd lambda_t = a_bar_t;
    
    for (int j = 1; j <= q; ++j) {
      if (t + j < n) lambda_t.noalias() += TH.block(0, (j - 1) * d, d, d).transpose() * Lambda.col(t + j);
    }
    Lambda.col(t) = lambda_t;
    
    for (int i = 1; i <= p; ++i) grad_PH_mat.block(0, (i - 1) * d, d, d).noalias() -= lambda_t * zt.col(t - i).transpose();
    for (int j = 1; j <= q; ++j) grad_TH_mat.block(0, (j - 1) * d, d, d).noalias() += lambda_t * at.col(t - j).transpose();
  }
  
  if (p > 0) res.grad.head(p * d2) = Map<VectorXd>(grad_PH_mat.data(), p * d2);
  if (q > 0) res.grad.tail(q * d2) = Map<VectorXd>(grad_TH_mat.data(), q * d2);
  
  return res;
}

// =========================================================================
// Rcpp 导出
// =========================================================================

// [[Rcpp::export]]
List zzr_varma_adjoint_list_rcpp(const Map<MatrixXd> &zt, const Map<VectorXd> &beta, int p, int q) {
  OptResult res = zzr_llk_grad_varma_adjoint_cpp(zt, beta, p, q);
  return List::create(Named("objective") = res.llk, Named("gradient") = res.grad);
}

// [[Rcpp::export]]
List zzr_varma_boot_adjoint_list_rcpp(const Map<MatrixXd> &zt, const Map<VectorXd> &beta, const Map<VectorXd> &w, int p, int q) {
  OptResult res = zzr_llk_grad_varma_boot_adjoint_cpp(zt, beta, w, p, q);
  return List::create(Named("objective") = res.llk, Named("gradient") = res.grad);
}





// [[Rcpp::export]]
MatrixXd zzrvarmaResiduals_cpp(const Map<MatrixXd> &zt,
                                  const Map<MatrixXd> &PH,
                                  const Map<MatrixXd> &TH,
                                  int p,
                                  int q) {
  
  const int d = zt.rows();
  const int nT = zt.cols();
  

  const int pqmax = std::max(p, q);
  

  MatrixXd at = zt;
  

  for (int t = pqmax; t < nT; t++) {
    

    VectorXd current_a = zt.col(t);
    

    for (int i = 1; i <= p; ++i) {

      current_a.noalias() -= PH.block(0, (i - 1) * d, d, d) * zt.col(t - i);
    }
    

    for (int j = 1; j <= q; ++j) {

      current_a.noalias() += TH.block(0, (j - 1) * d, d, d) * at.col(t - j);
    }
    

    at.col(t) = current_a;
  }
  
  return at;
}


// [[Rcpp::export]]
NumericVector zzr_boot_cov_cpp(const Map<MatrixXd>& X_boot, const Map<VectorXd>& w, int k_max) {
  int d = X_boot.rows();
  int n = X_boot.cols();
  int d2 = d * d;
  

  NumericVector out(d2 + k_max * d2);
  

  MatrixXd X_boot2 = X_boot.array().rowwise() * w.transpose().array();
  

  MatrixXd sigma = (X_boot * X_boot.transpose()) / double(n);

  Map<VectorXd>(&out[0], d2) = Map<VectorXd>(sigma.data(), d2);
  

  MatrixXd sm(d, d);
  for (int k = 1; k <= k_max; ++k) {
    

    sm.noalias() = (X_boot2.rightCols(n - k) * X_boot.leftCols(n - k).transpose()) / double(n);
    

    int start = d2 + (k - 1) * d2;
    Map<VectorXd>(&out[start], d2) = Map<VectorXd>(sm.data(), d2);
  }
  
  return out;
}