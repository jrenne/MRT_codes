
// ---- Enable C++11: ----
// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <cmath>
#include <vector>

//[[Rcpp::depends(RcppEigen)]] 

using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::export]]
Rcpp::List delta_function(const Eigen::MatrixXd & Y,
                        const Rcpp::List      & StateSpace,
                        const Eigen::VectorXd & indic_pos_z) {
  
  int time = Y.rows();
  Rcpp::List Omega_t              = Rcpp::List(time);
  Rcpp::List M_t              = Rcpp::List(time);
  Rcpp::List M_0t              = Rcpp::List(time);
  
  for(int t = 0; t < time; t++){
  
  Eigen::MatrixXd M_all     = StateSpace("M") ; // add for new delta
  
  Eigen::VectorXd M_row = M_all.row(t); // add for new delta
  M_0t = M_row;
  Eigen::MatrixXd M = M_row.asDiagonal(); // add for new delta
  M_t = M;
  Omega_t(t) = M * M.transpose();
  
  }
  
  return List::create(
    Named("Omega_t")        = Omega_t,
    Named("M_t")        = M_t,
    Named("M_0t")        = M_0t
  );
  
  
}


// [[Rcpp::export]]
MatrixXd Q_function(List Model, MatrixXd X) {
  
  int n = Model["n"]; 
  int q = Model["q"]; 
  MatrixXd Theta = Model["Theta"];
  MatrixXd nu = Model["nu"];
  MatrixXd mu_z = Model["mu.z"];
  MatrixXd Phi_z = Model["Phi.z"];
  MatrixXd Gamma_z0 = Model["Gamma.z0"];
  MatrixXd Gamma_z1 = Model["Gamma.z1"];
  MatrixXd Gamma_Y1 = Model["Gamma.Y1"];
  MatrixXd Gamma_Y0 = Model["Gamma.Y0"];
  
  // Identify z latent variables
  MatrixXd z = X.block(n, 0, q, X.cols());
  MatrixXd diag_vec_1 = (Gamma_z0 + Gamma_z1.transpose()*z).asDiagonal();
  MatrixXd diag_vec_2 = (Gamma_Y0 + Gamma_Y1.transpose()*(mu_z + Phi_z*z)).asDiagonal();
  
  // MatrixXd diag_mu_3 = mu.array().cube().matrix().asDiagonal().toDenseMatrix();
  
  // Compute covariance matrix of transition-eq innovations
  /// Sigma_11 = Theta*diag(Gamma.z0 + Gamma.z1'z_t-1)*Theta' + diag(Gamma.z0 + Gamma.z1'(mu_z+ z_t-1))
  MatrixXd A = Theta*diag_vec_1*Theta.transpose() + diag_vec_2;
  
  //Sigma_22 = diag(Gamma.z0 + Gamma.z1'z_t-1)
  MatrixXd B = diag_vec_1;
  
  //Sigma_12 = Theta*Sigma_22
  MatrixXd C = Theta*B;
  MatrixXd Q1 = MatrixXd::Zero(A.rows(), A.cols() + C.cols());
  Q1 << A, C;
  // Concatenate C_transpose and B horizontally
  MatrixXd Q2 = MatrixXd::Zero(C.transpose().rows(), C.transpose().cols() + B.cols());
  Q2 << C.transpose(), B;
  MatrixXd Q = MatrixXd::Zero(Q1.rows() + Q2.rows(), Q1.cols());
  Q << Q1, Q2;
  
  return Q;
  
}

// [[Rcpp::export]]
Rcpp::List KF_filter_cpp(const Eigen::MatrixXd & Y,
                         const Rcpp::List      & StateSpace,
                         const Eigen::VectorXd & indic_pos_z){
  
  // Weight matrices are (# [states] * # [states])
  // Kalman filter:
  // Measurement:
  // y_t   = mu_t+G*rho_t+M*eps_t
  // Transition:
  // rho_t = nu_t + H*rho_t-1 +N*xi_t
  
  
  Eigen::MatrixXd nu_t  = StateSpace("nu_t") ;
  Eigen::MatrixXd H     = StateSpace("H") ;
  //Eigen::MatrixXd N     = StateSpace("N") ;
  List sigma = StateSpace("N") ;
  Eigen::MatrixXd mu_t  = StateSpace("mu_t") ;
  Eigen::MatrixXd G     = StateSpace("G") ;
  Eigen::MatrixXd M_all     = StateSpace("M") ; // add for new delta
  //Eigen::MatrixXd M     = StateSpace("M") ; // comment for new delta
  Eigen::MatrixXd P0    = StateSpace("Sigma_0") ;
  Eigen::MatrixXd W0    = StateSpace("rho_0") ;

  // sample length:
  int time = Y.rows();
  
  // Number of latent variables:
  int n_W = G.cols();
  
  // Number of observed variables:
  int n_obs = G.rows();
  
  // Building objects to retrieve
  //-----------------------------
  double pi = M_PI;//3.141592653589793;
  
  // W_t_t is xi_tt
  Eigen::MatrixXd W_t_t           = Eigen::MatrixXd::Zero(n_W, time);
  // W_t_t_nomodif is orginal W_t_t
  Eigen::MatrixXd W_t_t_nomodif   = Eigen::MatrixXd::Zero(n_W, time);
  // W_t_t1 is xi_ttm1
  Eigen::MatrixXd W_t_t1          = Eigen::MatrixXd::Zero(n_W, time);
  // Obs_t_t1 are the observables predicted
  Eigen::MatrixXd Obs_t_t1        = Eigen::MatrixXd::Zero(n_obs, time);
  // Obs_t_t are the observables predicted
  Eigen::MatrixXd Obs_t_t        = Eigen::MatrixXd::Zero(n_obs, time);
  // residuals_t_t1 is Yt - Obs_pred
  Eigen::MatrixXd residuals_t_t1  = Eigen::MatrixXd::Zero(n_obs, time);
  
  Rcpp::List P_t_t                = Rcpp::List(time);
  Rcpp::List P_t_t1               = Rcpp::List(time);
  Rcpp::List Omega_t              = Rcpp::List(time);
  Rcpp::List Gain_t              = Rcpp::List(time);  
  Rcpp::List N_t              = Rcpp::List(time);  
  //Rcpp::List R_t                  = Rcpp::List(time);
  Rcpp::List M_t_t1               = Rcpp::List(time);
  
  Eigen::VectorXd loglik          = Eigen::VectorXd::Zero(time);
  Eigen::VectorXd penalty_neg_vec = Eigen::VectorXd::Zero(time);
  double penalty_neg_factor = 100000;
  
  // Initialize the Filter
  //----------------------
  Eigen::VectorXd W = W0;
  Eigen::MatrixXd P = P0;
  
  Eigen::VectorXd W_pred          = Eigen::VectorXd::Zero(n_W);
  Eigen::VectorXd Obs_pred        = Eigen::VectorXd::Zero(n_obs);
  Eigen::VectorXd Observables_day = Eigen::VectorXd::Zero(n_obs);
  
  Eigen::MatrixXd P_pred          = Eigen::MatrixXd::Zero(n_W, n_W);
  Eigen::MatrixXd M_pred          = Eigen::MatrixXd::Zero(n_obs, n_obs);
  
  Eigen::MatrixXd Condi_P_pred    = Eigen::MatrixXd::Zero(n_W, n_W);
  
  
  // LAUNCH THE RECURSIONS
  //======================
  for(int t = 0; t < time; t++){
    //for(int t = 0; t < 2; t++){
    
    // intialize penalty_neg_t to 0 at each step
    double penalty_neg_t = 0; 
    
    Eigen::VectorXd nut  = nu_t.row(t) ;
    Eigen::VectorXd mut  = mu_t.row(t) ;
    Eigen::VectorXd Yt   = Y.row(t) ;

    // PREDICTION STEP
    //----------------
    //1  Forecast latent variables for next iteration
    // xi.t|t-1 = mu + Fxi.t-1|t-1
    // Mean of latent factors
    W_pred = nut + H * W ;
    
    // Update Q and R for the next iteration
    MatrixXd N;
    if(t==0){
      N = Q_function(sigma, W0);
    } else{
      N = Q_function(sigma, W);
    }
    
    N_t(t) = N;
    
    // Adjustment we have to make to the filter pertains to the fact that factors
    // z_t are non negative. For this purpose, after each updating step of the 
    // algorithm, negative entries in the z_t estimate are replaced by 0.
    if (indic_pos_z.sum()>0) {
      MatrixXd xi_ttm1_aux = W_pred;
      MatrixXd xi_ttm1_aux_aux =  W_pred;
      
      for (int i = 0; i < xi_ttm1_aux.rows(); ++i) {
        if (indic_pos_z(i) == 1 && xi_ttm1_aux(i) < 0) {
          xi_ttm1_aux_aux(i) = std::max(xi_ttm1_aux(i), 0.0);
        }
      }
      W_pred = xi_ttm1_aux_aux;
    }
    
    
    //3 Forecast observed variables for next iteration
    //  P.t|t-1 = F P.t-1|t-1 F'  + Q
    P_pred = H * P * H.transpose() + N ;//N * N.transpose()
    
    // Predicting the observables
    //2 Forecast of y_t: y.t|t-1 = A'X.t + H'xi.t|t-1 
    Obs_pred  = mut + G * W_pred;
    //4 h.t = H'P.t|t-1 H + R 
    Eigen::VectorXd M_row = M_all.row(t); // add for new delta
    Eigen::MatrixXd M = M_row.asDiagonal(); // add for new delta
    Omega_t(t) = M * M.transpose();
    M_pred    = G * P_pred * G.transpose() + M * M.transpose();
    
    // UPDATING STEP
    //--------------
    //6 Forecast error: eta.t = y.t - y.t|t-1
    Eigen::VectorXd residual = Yt - Obs_pred ;
    
    
    // Determining which components are observed
    Array<bool,Dynamic,1> bool_na(Yt.array() == Yt.array());
    
    int n_obs_modif = bool_na.count();
    
    // Fill the new components
    Eigen::VectorXd resid_modif = Eigen::VectorXd::Zero(n_obs_modif);
    Eigen::MatrixXd M_modif     = Eigen::MatrixXd::Zero(n_obs_modif, n_obs_modif);
    Eigen::MatrixXd G_modif     = Eigen::MatrixXd::Zero(n_obs_modif, n_W);
    
    // Control for unobserved components
    //----------------------------------
    int counter_rows = 0;
    for (int i = 0; i < n_obs; i++){
      
      if(bool_na(i)==TRUE){
        
        resid_modif(counter_rows)  = residual(i);
        G_modif.row(counter_rows)  = G.row(i);
        
        int counter_cols = 0 ;
        
        for (int j = 0; j < n_obs; j++){
          
          if (bool_na(j) == TRUE){
            
            M_modif(counter_rows, counter_cols) = M_pred(i,j);
            counter_cols += 1;
          }
        }
        counter_rows += 1;
      }
      
    }// End of the NA loop
    
    
    // Perform the update step and the loglik
    //---------------------------------------
    double det_M_pred      = M_modif.determinant();
    
    // if (det_M_pred < 0 ||det_M_pred > 1e80) {
    //
    //   loglik(t) = -1e20;
    //   //std::cout << det_M_pred << std::endl;
    //   break;
    //
    // }
    
    // Add by AT
    Eigen::MatrixXd M_inv;
    // Eigen::MatrixXd M_inv = M_modif.inverse();
    
    if (det_M_pred < 0 ||det_M_pred > 1e80) {
      M_inv = MatrixXd::Zero(M_modif.rows(), M_modif.cols());
    } else{
     M_inv = M_modif.inverse();
    }
    
    //5 Gain equation: K.t = [P.t|t-1*H + S']h_t^-1
    Eigen::MatrixXd Gain  = P_pred * G_modif.transpose() * M_inv;
    Gain_t(t) = Gain;
    
    //7 Final prediction of the latent variables:
    //  xi.t|t =  Fxi.t|t-1 - K.t*eta.t
    Eigen::VectorXd W_up  = W_pred + Gain * resid_modif;
    
    
    // Adjustment we have to make to the filter pertains to the fact that factors
    // z_t are non negative. For this purpose, after each updating step of the 
    // algorithm, negative entries in the z_t estimate are replaced by 0.
    MatrixXd xi_tt_aux =  W_up;
    MatrixXd xi_tt_aux_aux =  W_up;
    
    if (indic_pos_z.sum()>0) {
      
      // Condtion for penalty_neg_t
      for (int i = 0; i < xi_tt_aux.rows(); ++i) {
        if (indic_pos_z(i) == 1 && xi_tt_aux(i) < 0) {
          penalty_neg_t += xi_tt_aux(i) * xi_tt_aux(i);
          xi_tt_aux_aux(i) = std::max(xi_tt_aux(i), 0.0);
        }
      }
      penalty_neg_t = penalty_neg_t * penalty_neg_factor;
      penalty_neg_vec(t) = penalty_neg_t; 
    }
    W_up = xi_tt_aux_aux;
    
    
    //8 Final prediction for the variance of latent variables
    // P.t|t = P.t|t-1 - K.t(H'P.t|t-1 + S)
    Eigen::MatrixXd P_up  = P_pred - Gain * G_modif * P_pred;
    
    double lik_val;
    if (std::isnan(det_M_pred) || std::isinf(det_M_pred) || det_M_pred <= 0) { // determinant not correct => penalize function
      lik_val         = -.5*(n_obs_modif * log(2*pi) +
        resid_modif.transpose() * M_inv * resid_modif)  - 70000000 - penalty_neg_t;  
    } else {
        lik_val         = -.5*(n_obs_modif * log(2*pi) + log(det_M_pred) +
                                resid_modif.transpose() * M_inv * resid_modif) - penalty_neg_t;
    }
    
    //double lik_val         = -.5*(n_obs_modif * log(2*pi) + log(det_M_pred) +
    //                              resid_modif.transpose() * M_inv * resid_modif);
    
    // Update the values before looping back
    //--------------------------------------
    W = W_up;
    P = P_up;
   
    // STORING THE VALUES
    //---------------------
    W_t_t.col(t)          = W_up;
    W_t_t1.col(t)         = W_pred;
    Obs_t_t1.col(t)       = Obs_pred;
    residuals_t_t1.col(t) = residual;
    
    P_t_t(t)  = P_up;
    P_t_t1(t) = P_pred;
    M_t_t1(t) = M_pred;
    
    loglik(t) = lik_val;
    
  }
  
  // Compute fitted values for observables:
  //  y.t|t =  A*X + H*xi.t|t-1
  Obs_t_t  = mu_t.transpose() + G * W_t_t;
  
  double logl = loglik.sum();
  
  // Sending results back
  //---------------------
  return List::create(
    Named("loglik")         = logl,
    Named("loglik.vector")  = loglik,
    Named("W.updated")      = W_t_t,
    Named("P.updated")      = P_t_t,
    Named("W.predicted")    = W_t_t1,
    Named("P.predicted")    = P_t_t1,
    Named("Obs.predicted")  = Obs_t_t1,
    Named("Obs.updated")  = Obs_t_t,
    Named("M.predicted")    = M_t_t1,
    Named("W.0")            = W0,
    Named("P.0")            = P0,
    Named("residual")       = residuals_t_t1,
    Named("Omega_t")        = Omega_t,
    Named("Gain_t")        = Gain_t,
    Named("N_t")        = N_t,
    //Named("R_t")            = R_t,
    Named("penalty_neg_vec")        = penalty_neg_vec
  );
  
  
}


// [[Rcpp::export]]
Eigen::MatrixXd ginv_cpp(const Eigen::MatrixXd& matrix, double tol=1e-30) {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
  const auto& singular_values = svd.singularValues();
  Eigen::MatrixXd singular_values_inv(matrix.cols(), matrix.rows());
  singular_values_inv.setZero();
  for (Eigen::Index i = 0; i < singular_values.size(); ++i) {
    if (singular_values(i) > tol)
      singular_values_inv(i, i) = 1.0 / singular_values(i);
  }
  return svd.matrixV() * singular_values_inv * svd.matrixU().transpose();
}

// [[Rcpp::export]]
Rcpp::List Kalman_filter_cpp(List all_parameters, MatrixXd Y,MatrixXd X, 
                             MatrixXd xi_00, MatrixXd P_00, VectorXd indic_pos_z){
  
  // Extract dimensions of interest.
  int r = xi_00.rows();
  int T = Y.rows();
  int n = Y.cols();
  
  //  Extract matrices from the list "all.parameters
  MatrixXd mu = all_parameters["mu"];
  MatrixXd F = all_parameters["F"];
  List sigma = all_parameters["sigma"];
  MatrixXd A = all_parameters["A"];
  MatrixXd H = all_parameters["H"];
  MatrixXd delta = all_parameters["delta"];
  
  // Initialize Q and R for t=0
  MatrixXd Q = Q_function(sigma, xi_00); //Qfunction defined below (transition eq.)
  List list_Q_t(T);
  list_Q_t = Q;
  //List list_Q_t = List::create(Named("Q") = Q);
  MatrixXd R = delta*delta.transpose();
  List list_R_t(T);
  list_R_t = R;
  //List list_R_t = List::create(Named("R") = R);
  
  // Create S
  MatrixXd S = MatrixXd::Zero(n, r);
  
  // Outputs:
  // xi.tt is r x 1
  // Create a matrix matrix_xi_tt with dimensions T x r filled with NA values
  MatrixXd matrix_xi_tt = MatrixXd::Constant(T, r, std::numeric_limits<double>::quiet_NaN());
  MatrixXd matrix_xi_tt_aux = MatrixXd::Constant(T, r, std::numeric_limits<double>::quiet_NaN());
  MatrixXd matrix_xi_ttm1 = MatrixXd::Constant(T+1, r, std::numeric_limits<double>::quiet_NaN());
  List list_xi_tt(T);
  //List list_xi_tt = List::create(MatrixXd::Constant(T, r, std::numeric_limits<double>::quiet_NaN()));
  
  // P.tt is r x r
  //  matrix.P.tt  is of dimension T x (r^2)
  MatrixXd matrix_P_tt = MatrixXd::Constant(T, r*r, std::numeric_limits<double>::quiet_NaN());
  MatrixXd matrix_P_ttm1 = MatrixXd::Constant(T+1, r*r, std::numeric_limits<double>::quiet_NaN());
  List list_P_tt(T);
  //List list_P_tt = List::create(MatrixXd::Constant(T, r*r, std::numeric_limits<double>::quiet_NaN()));
  
  // Intermediate variables of T elements (filled with NA)
  Rcpp::List list_xi_ttm1(T+1);
  List list_P_ttm1(T+1);
  //List list_xi_ttm1 = List::create(MatrixXd::Constant(T+1, r, std::numeric_limits<double>::quiet_NaN()));
  //List list_P_ttm1 = List::create(MatrixXd::Constant(T+1, r*r, std::numeric_limits<double>::quiet_NaN()));
  
  // Initialize xi.1|0 and P.1|0 for the loop (Equations 1 and 3)
  //(avoids having the initial condition in the final results)
  list_xi_ttm1[0] = mu + F * xi_00; // = xi.1|0
  matrix_xi_ttm1.row(0) = (mu + F * xi_00).transpose();
  list_P_ttm1[0] = F* P_00*F.transpose() + Q ; // = P.1|0
  matrix_P_ttm1.row(0) = (F* P_00*F.transpose() + Q).transpose();
  
  //  Initialize gains to zero (K.t)
  List list_K_t(T);
  
  //  Initialize log likelihood to zero
  double log_lhd = 0;
  VectorXd log_lhd_vec = Eigen::VectorXd::Zero(T);
  VectorXd penalty_neg_vec = Eigen::VectorXd::Zero(T);
  double penalty_neg_factor = 100000;
  
  // Loop that calculate the filter
  for (int t = 0; t <= T-1; t++) {

    // intialize penalty_neg_t to 0 at each step
    double penalty_neg_t = 0; 
    //2 Forecast of y_t: y.t|t-1 = A'X.t + H'xi.t|t-1 A.transpose() * X.row(0).transpose() + H.transpose() *
    MatrixXd y_t =  A.transpose() * X.row(t).transpose() + H.transpose() * Rcpp::as<Eigen::MatrixXd>(list_xi_ttm1[t]);
    
    //4 h.t = H'P.t|t-1H + R + 2H'S'
    MatrixXd h_t = H.transpose() *  Rcpp::as<Eigen::MatrixXd>(list_P_ttm1[t]) * H + R + 2*H.transpose()*S.transpose();
    
    //6 Forecast error: eta.t = y.t - y.t|t-1
    MatrixXd eta_t = Y.row(t).transpose() - y_t;
    
    // Updating step for missing variables 
    // Determining which components are observed
    Array<bool,Dynamic,1> indic_notNaN(Y.row(t).array() == Y.row(t).array());
    
    int n_obs_modif = indic_notNaN.count();
    MatrixXd K_t;
    MatrixXd K_t_all;
    MatrixXd h_t_star;
    MatrixXd h_t_inverse;
    MatrixXd H_star;
    MatrixXd R_star;
    MatrixXd S_star;
    MatrixXd eta_t_star;
    

    if (n_obs_modif == 0) {// There is no observation on this date.
      
      // Define n.star to 0 as no observed value
      int n_star = 0;

      //4 h.t = H'P.t|t-1H + R + 2H'S'
      // No need to redefine it in this case.
      
      //5 Gain equation: K.t = [P.t|t-1*H + S']h_t^-1
      // Use equation 5 to export K.t results.
      //h_t_inverse = ginv_cpp(h_t);
      
      
      double determinant = h_t.determinant();
      //Eigen::JacobiSVD<Eigen::MatrixXd> svd(h_t);
      //double conditionNumber = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
      
      // Add by AT
      Eigen::MatrixXd h_t_inverse;
      //Eigen::MatrixXd M_inv = M_modif.inverse();
      
      if (determinant <= 0 ||determinant > 1e80) {
        h_t_inverse = MatrixXd::Zero(h_t.rows(), h_t.cols());
        
      } else{
        h_t_inverse = h_t.inverse();
      }
      
      // if (determinant == 0) { // Matrix is singular.
      //   h_t_inverse = MatrixXd::Zero(h_t.rows(), h_t.cols()); //h_t.inverse();
      // } else if (conditionNumber > 1e20) { // Matrix is close to being singular.
      //   h_t_inverse = MatrixXd::Zero(h_t.rows(), h_t.cols()); //h_t.inverse();
      // } else { // Matrix is not singular.
      //   h_t_inverse = h_t.completeOrthogonalDecomposition().pseudoInverse(); //h_t.inverse();
      // }
      
      //h_t_inverse = h_t.completeOrthogonalDecomposition().pseudoInverse(); //h_t.inverse();
      K_t_all = (Rcpp::as<Eigen::MatrixXd>(list_P_ttm1[t]) * H + S.transpose()) * h_t_inverse;
     
      for (int i = 0; i < K_t_all.rows(); i++)
      {
        for (int j = 0; j < K_t_all.cols(); j++)
        {
          if (!indic_notNaN(j))
            K_t_all(i, j) = std::numeric_limits<double>::quiet_NaN();
        }
      }    
      
      //7 Final prediction of the latent variables:
      //  xi.t|t =  Fxi.t|t-1 - K.t*eta.t
      list_xi_tt[t] = list_xi_ttm1[t];
      matrix_xi_tt.row(t) = Rcpp::as<Eigen::MatrixXd>(list_xi_tt[t]).transpose();
      matrix_xi_tt_aux.row(t) = Rcpp::as<Eigen::MatrixXd>(list_xi_tt[t]).transpose();
      
      //8 Final prediction for the variance of latent variables
      // P.t|t = P.t|t-1 - K.t(H'P.t|t-1 + S)
      list_P_tt[t] = list_P_ttm1[t];
      MatrixXd P_tt_aux_mat = Rcpp::as<Eigen::MatrixXd>(list_P_tt[t]);
      matrix_P_tt.row(t) = Eigen::Map<Eigen::VectorXd>(P_tt_aux_mat.data(), P_tt_aux_mat.size()).transpose();
      
      // log-likelihood update:
      // Past log.lik + -0.5*n*log(2*pi) - 0.5*log(det(h.t)) - 0.5*eta.t' h.t^-1 eta.t
      // Note: to calculate the log likelihood we do it incrementally (for each t).
      // First part of the log likelihood : -0.5*n*log(2*pi)
      // Log.lik = -0.5*n*T*log(2*pi) - 0.5*sum(det(h.t)) - 0.5*sum(eta.t' h.t eta.t)
      // n.star = eta.t.star.size()
      double logl_aux = -0.5 * n_star * log(2 * M_PI);
      
      // Second part of the log likelihood (no update as h.t and eta.t.star are empty)
      log_lhd += logl_aux;
      log_lhd_vec(t) = logl_aux; 

    } else {// There is at least one observation on this date.
      
      int n_star = n_obs_modif;
      
      // Modify H
      H_star = MatrixXd::Zero(r, n_star); //r=H.rows()
      int star_idx = 0;
      for (int i = 0; i < H.cols(); ++i) {
        if (indic_notNaN(i)) {
          H_star.col(star_idx) = H.col(i);
          star_idx++;
        }
      }
      // Modify R
      R_star = MatrixXd::Zero(n_star, n_star);
      star_idx = 0;
      for (int i = 0; i < R.rows(); ++i) {
        if (indic_notNaN(i)) {
          int col_star_idx = 0;
          for (int j = 0; j < R.cols(); ++j) {
            if (indic_notNaN(j)) {
              R_star(star_idx, col_star_idx) = R(i, j);
              col_star_idx++;
            }
          }
          star_idx++;
        }
      }
      // Modify S
      S_star = MatrixXd::Zero(n_star, r);
      star_idx = 0;
      for (int i = 0; i < S.rows(); ++i) {
        if (indic_notNaN(i)) {
          S_star.row(star_idx) = S.row(i);
          star_idx++;
        }
      }
      
      // Modify eta_t
      eta_t_star = MatrixXd::Zero(n_star, 1);
      star_idx = 0;
      for (int i = 0; i < eta_t.size(); ++i) {
        if (indic_notNaN(i)) {
          eta_t_star(star_idx) = eta_t(i);
          star_idx++;
        }
      }   
      
      //4 h.t = H'P.t|t-1H + R + 2H'S'
      // h.t.star <- t(H.star) %*% list.P.ttm1[[t]] %*% H.star + R.star + 2*t(H.star)%*%t(S.star).
      // Modify R
      h_t_star = MatrixXd::Zero(n_star, n_star);
      star_idx = 0;
      for (int i = 0; i < h_t.rows(); ++i) {
        if (indic_notNaN(i)) {
          int col_star_idx = 0;
          for (int j = 0; j < h_t.cols(); ++j) {
            if (indic_notNaN(j)) {
              h_t_star(star_idx, col_star_idx) = h_t(i, j);
              col_star_idx++;
            }
          }
          star_idx++;
        }
      }
      
      //5 Gain equation: K.t = [P.t|t-1*H + S']h_t^-1
      // Use equation 5 to export K.t results.
      //h_t_inverse = ginv_cpp(h_t);
      
      double determinant = h_t_star.determinant();
      // Eigen::JacobiSVD<Eigen::MatrixXd> svd(h_t_star);
      // double conditionNumber = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);

      // Add by AT
      Eigen::MatrixXd h_t_inverse;
      //Eigen::MatrixXd M_inv = M_modif.inverse();
      
      if (determinant <= 0 ||determinant > 1e80) {
        h_t_inverse = MatrixXd::Zero(h_t_star.rows(), h_t_star.cols());
        
      } else{
        h_t_inverse = h_t_star.inverse();
      }
      
      // if (determinant == 0) { // Matrix is singular.
      //   h_t_inverse = MatrixXd::Zero(h_t_star.rows(), h_t_star.cols()); //h_t.inverse();
      // } else if (conditionNumber > 1e20) { // Matrix is close to being singular.
      //   h_t_inverse = MatrixXd::Zero(h_t_star.rows(), h_t_star.cols()); //h_t.inverse();
      // } else { // Matrix is not singular.
      //   h_t_inverse = h_t_star.completeOrthogonalDecomposition().pseudoInverse(); //h_t.inverse();
      // }
      
      //h_t_inverse = h_t_star.completeOrthogonalDecomposition().pseudoInverse();
      K_t = (Rcpp::as<Eigen::MatrixXd>(list_P_ttm1[t]) * H_star + S_star.transpose()) * h_t_inverse;
      K_t_all = (Rcpp::as<Eigen::MatrixXd>(list_P_ttm1[t]) * H + S.transpose()) * h_t_inverse;

      // Replace with NA when indic_notNaN is false
      for (int i = 0; i < K_t_all.rows(); i++)
      {
        for (int j = 0; j < K_t_all.cols(); j++)
        {
          if (!indic_notNaN(j))
            K_t_all(i, j) = std::numeric_limits<double>::quiet_NaN();
        }
      }    
          
      //7 Final prediction of the latent variables:
      //  xi.t|t =  Fxi.t|t-1 - K.t*eta.t
      list_xi_tt[t] = Rcpp::as<Eigen::MatrixXd>(list_xi_ttm1[t]) + K_t * eta_t_star;
        
      // Adjustment we have to make to the filter pertains to the fact that factors
      // z_t are non negative. For this purpose, after each updating step of the 
      // algorithm, negative entries in the z_t estimate are replaced by 0.
      MatrixXd xi_tt_aux =  Rcpp::as<Eigen::MatrixXd>(list_xi_tt[t]);
      MatrixXd xi_tt_aux_aux =  Rcpp::as<Eigen::MatrixXd>(list_xi_tt[t]);
      
      if (indic_pos_z.sum()>0) {
        
        // Condtion for penalty_neg_t
        for (int i = 0; i < xi_tt_aux.rows(); ++i) {
          if (indic_pos_z(i) == 1 && xi_tt_aux(i) < 0) {
            penalty_neg_t += xi_tt_aux(i) * xi_tt_aux(i);
            xi_tt_aux_aux(i) = std::max(xi_tt_aux(i), 0.0);
          }
        }
        penalty_neg_t = penalty_neg_t *penalty_neg_factor;
        penalty_neg_vec(t) = penalty_neg_t; 
      }
      list_xi_tt[t] = xi_tt_aux_aux;
      matrix_xi_tt.row(t) = Rcpp::as<Eigen::MatrixXd>(list_xi_tt[t]).transpose();
      matrix_xi_tt_aux.row(t) = xi_tt_aux.transpose();
      
      //8 Final prediction for the variance of latent variables
      // P.t|t = P.t|t-1 - K.t(H'P.t|t-1 + S)
      list_P_tt[t] = Rcpp::as<Eigen::MatrixXd>(list_P_ttm1[t]) - K_t * (H_star.transpose() * Rcpp::as<Eigen::MatrixXd>(list_P_ttm1[t]));
      MatrixXd P_tt_aux_mat = Rcpp::as<Eigen::MatrixXd>(list_P_tt[t]);
      matrix_P_tt.row(t) = Eigen::Map<Eigen::VectorXd>(P_tt_aux_mat.data(), P_tt_aux_mat.size()).transpose();
      // 
      // log-likelihood update:
      // Past log.lik + -0.5*n*log(2*pi) - 0.5*log(det(h.t)) - 0.5*eta.t' h.t^-1 eta.t
      // Note: to calculate the log likelihood we do it incrementally (for each t).
      // First part of the log likelihood : -0.5*n*log(2*pi)
      // Log.lik = -0.5*n*T*log(2*pi) - 0.5*sum(det(h.t)) - 0.5*sum(eta.t' h.t eta.t)
      // n.star = eta.t.star.size()
      double logl_aux = -0.5 * n_star * log(2 * M_PI);
      
      // Second part of the leg likelihood
      //log.lhd <- log.lhd + logl.aux  - (1/2*log(det(h.t))) -
      //  (1/2*t(eta.t.star)%*%ginv(h.t.star)%*%(eta.t.star))
      double det_h_t; // Variable to store the determinant
      
      if (h_t_star.size() == 1) {
        det_h_t = h_t_star.coeff(0);
      } else {
        det_h_t = h_t_star.determinant();
      }
      
       if (std::isnan(det_h_t) || std::isinf(det_h_t) || det_h_t <= 0) { //std::isnan(det_h_t) || std::isinf(det_h_t) || std::isinf(det_h_t) || det_h_t < 0
        log_lhd = log_lhd + logl_aux - 70000000 - penalty_neg_t - (0.5 * eta_t_star.transpose() * h_t_inverse * eta_t_star).coeff(0);
        log_lhd_vec(t) = logl_aux - 70000000 - penalty_neg_t - (0.5 * eta_t_star.transpose() * h_t_inverse * eta_t_star).coeff(0);
       } else {
        log_lhd = log_lhd + logl_aux - (0.5* log(det_h_t)) - penalty_neg_t - (0.5 * eta_t_star.transpose() * h_t_star.inverse() * eta_t_star).coeff(0);
        log_lhd_vec(t) = logl_aux - (0.5* log(det_h_t)) - penalty_neg_t - (0.5 * eta_t_star.transpose() * h_t_star.inverse() * eta_t_star).coeff(0);
      }
      
    }
  
  //1  Forecast latent variables for next iteration
  // xi.t|t-1 =  Fxi.t-1|t-1
  list_xi_ttm1[t+1] = mu + F * Rcpp::as<Eigen::MatrixXd>(list_xi_tt[t]);
    
  // Adjustment we have to make to the filter pertains to the fact that factors
  // z_t are non negative. For this purpose, after each updating step of the 
  // algorithm, negative entries in the z_t estimate are replaced by 0.
  if (indic_pos_z.sum()>0) {
    MatrixXd xi_ttm1_aux = Rcpp::as<Eigen::MatrixXd>(list_xi_ttm1[t+1]);
    MatrixXd xi_ttm1_aux_aux =  Rcpp::as<Eigen::MatrixXd>(list_xi_ttm1[t+1]);
    
    for (int i = 0; i < xi_ttm1_aux.rows(); ++i) {
      if (indic_pos_z(i) == 1 && xi_ttm1_aux(i) < 0) {
        xi_ttm1_aux_aux(i) = std::max(xi_ttm1_aux(i), 0.0);
      }
    }
    list_xi_ttm1[t+1] = xi_ttm1_aux_aux;
  }
  
  // Update Q and R for the next iteration
  MatrixXd Q_t = Q_function(sigma, as<Eigen::MatrixXd>(list_xi_tt[t])); //Qfunction defined below (transition eq.)
  list_Q_t[t+1] = Q_t;
  MatrixXd R_t = delta*delta.transpose();
  list_R_t[t+1] = R_t;
  
  //3 Forecast observed variables for next iteration
  //  P.t|t-1 = F P.t-1|t-1 F'  + Q
  list_P_ttm1[t+1] = F * Rcpp::as<Eigen::MatrixXd>(list_P_tt[t]) * F.transpose() + Q_t;

  // Prepare xi.ttm1 and P.ttm1 for extraction
  matrix_xi_ttm1.row(t+1) = as<Eigen::MatrixXd>(list_xi_ttm1[t+1]).transpose();
  MatrixXd P_ttm1_aux_mat = Rcpp::as<Eigen::MatrixXd>(list_P_ttm1[t+1]);
  matrix_P_tt.row(t+1) = Eigen::Map<Eigen::VectorXd>(P_ttm1_aux_mat.data(), P_ttm1_aux_mat.size()).transpose();

  
    
  } // End of the loop
    
   MatrixXd fitted_obs = X * A + matrix_xi_tt * H ; // Fitted observables
      
  return List::create(Named("xi.tt")=matrix_xi_tt,
                      Named("xi.tt.aux")=matrix_xi_tt_aux,
                      Named("P.tt")=matrix_P_tt,
                      Named("P.tt-1")=matrix_P_ttm1,
                      Named("P.tt-1.list")=list_P_ttm1,
                      Named("xi.tt-1")=matrix_xi_ttm1,
                      Named("xi.tt-1.list")=list_xi_ttm1,
                      Named("K.t")=list_K_t,
                      Named("Q.t")=list_Q_t,
                      Named("R.t")=list_R_t,
                      Named("log.lik")=log_lhd,
                      Named("log.lik.vec")=log_lhd_vec,
                      Named("penalty.neg.vec")=penalty_neg_vec,
                      Named("fitted.obs") = fitted_obs);
    
    

  
  // return List::create(Named("H_star") = H_star, Named("S_star") = S_star,
  //                     Named("R_star") = R_star, Named("h_t_star") = h_t_star,
  //                     Named("eta_t_star") = eta_t_star,Named("indic_notNaN") = indic_notNaN, Named("n_obs_modif") = n_obs_modif);
  
}

