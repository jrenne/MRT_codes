#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <array>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
//#include <unsupported/Eigen/CXX11/Tensor>
//#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Rcpp;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::Map;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
MatrixXd a_1(MatrixXd gamma, List Model) {
  // Extract elements of interest in "Model"
  MatrixXd mu_X = Model["mu.X"];
  
  // Compute a_1
  MatrixXd result = gamma.transpose()*mu_X;
  return result;
}

// [[Rcpp::export]]
MatrixXd b_1(MatrixXd gamma, List Model) {
  // Extract elements of interest in "Model"
  MatrixXd Phi_X = Model["Phi.X"];
  
  // Compute b_1
  MatrixXd result = Phi_X.transpose() * gamma;
  return result;
}

// [[Rcpp::export]]
MatrixXd alpha_1(MatrixXd gamma, List Model, MatrixXd S_n, MatrixXd S_q) {
  // Extract elements of interest in "Model"
  int n = Model["n"];
  int q = Model["q"];
  MatrixXd Theta = Model["Theta"];
  MatrixXd Gamma_z0 = Model["Gamma.z0"];
  MatrixXd Gamma_Y0 = Model["Gamma.Y0"];
  MatrixXd Gamma_Y1 = Model["Gamma.Y1"];
  MatrixXd mu_z = Model["mu.z"];
  
  // Create gamma_Y, gamma_z and Id_q
  MatrixXd gamma_Y = gamma.block(0, 0, n, 1);
  MatrixXd gamma_z = gamma.block(n, 0, q, 1);
  MatrixXd Id_q = MatrixXd::Identity(q, q);
  
  // Compute alpha_1
  MatrixXd result = kroneckerProduct(gamma_Y,gamma_Y).transpose()*(kroneckerProduct(Theta,Theta)*S_q*Gamma_z0 +
    S_n*Gamma_Y0 + S_n*Gamma_Y1.transpose()*mu_z) + kroneckerProduct(gamma_z,gamma_z).transpose()*S_q*Gamma_z0 +
    2*kroneckerProduct(gamma_z,gamma_Y).transpose()*kroneckerProduct(Id_q, Theta)*S_q*Gamma_z0;

  return result;
}

// [[Rcpp::export]]
MatrixXd beta_1(MatrixXd gamma, List Model, MatrixXd S_n, MatrixXd S_q) {
  // Extract elements of interest in "Model"
  int n = Model["n"];
  int q = Model["q"];
  MatrixXd Theta = Model["Theta"];
  MatrixXd Gamma_z1 = Model["Gamma.z1"];
  MatrixXd Gamma_Y0 = Model["Gamma.Y0"];
  MatrixXd Gamma_Y1 = Model["Gamma.Y1"];
  MatrixXd Phi_z = Model["Phi.z"];
  
  // Create gamma_Y, gamma_z and Id_q
  MatrixXd gamma_Y = gamma.block(0, 0, n, 1);
  MatrixXd gamma_z = gamma.block(n, 0, q, 1);
  MatrixXd Id_q = MatrixXd::Identity(q, q);
  
  // Compute beta_1
  MatrixXd result_inter = (kroneckerProduct(gamma_Y,gamma_Y).transpose()*(kroneckerProduct(Theta,Theta)*S_q*Gamma_z1.transpose() +
    S_n *Gamma_Y1.transpose()* Phi_z ) + kroneckerProduct(gamma_z,gamma_z).transpose()*S_q*Gamma_z1.transpose() +
    2*kroneckerProduct(gamma_z,gamma_Y).transpose()*kroneckerProduct(Id_q, Theta)*S_q*Gamma_z1.transpose()).transpose();
  
  // Correct the dimension to depend on X_t
  MatrixXd result((n + q) , 1); // create a new matrix with the desired size
  result << MatrixXd::Zero(n , 1), result_inter; // set the first n columns to zeros and the remaining columns to the matrix M1_z
  
  return result;
}


// [[Rcpp::export]]
MatrixXd alpha_dot_dot_1(MatrixXd gamma, List Model, MatrixXd S_q) {
  // Extract elements of interest in "Model"
  int n = Model["n"];
  int q = Model["q"];
  MatrixXd Theta = Model["Theta"];
  MatrixXd nu = Model["nu"];
  MatrixXd mu = Model["mu"];
  MatrixXd Gamma_z0 = Model["Gamma.z0"];
  MatrixXd Gamma_Y1 = Model["Gamma.Y1"];
  
  // Create gamma_Y and gamma_z
  MatrixXd gamma_Y = gamma.block(0, 0, n, 1);
  MatrixXd gamma_z = gamma.block(n, 0, q, 1);
  
  // Compute alpha_dot_dot_1
  MatrixXd result = 2*(gamma_Y.transpose()*Theta + gamma_z.transpose()).array().cube().matrix()*(nu.array()*mu.array().cube()).matrix() +
    3*kroneckerProduct((Gamma_Y1*(gamma_Y.array().square().matrix())).transpose(), gamma_Y.transpose()*Theta + gamma_z.transpose())*S_q*Gamma_z0;
    
  return result;
}

// [[Rcpp::export]]
MatrixXd beta_dot_dot_1(MatrixXd gamma, List Model, MatrixXd S_q) {
  // Extract elements of interest in "Model"
  int n = Model["n"];
  int q = Model["q"];
  MatrixXd Theta = Model["Theta"];
  MatrixXd nu = Model["nu"];
  MatrixXd mu = Model["mu"];
  MatrixXd phi = Model["phi"];
  MatrixXd Gamma_z1 = Model["Gamma.z1"];
  MatrixXd Gamma_Y1 = Model["Gamma.Y1"];
  
  // Create gamma_Y and gamma_z
  MatrixXd gamma_Y = gamma.block(0, 0, n, 1);
  MatrixXd gamma_z = gamma.block(n, 0, q, 1);

  // Create diag(mu^3)
  MatrixXd diag_mu_3 = mu.array().cube().matrix().asDiagonal().toDenseMatrix();
  
  // Compute alpha_dot_dot_1
  MatrixXd result_inter = (6*(gamma_Y.transpose()*Theta + gamma_z.transpose()).array().cube().matrix()*(diag_mu_3*phi.transpose()) +
    3*kroneckerProduct((Gamma_Y1*(gamma_Y.array().square().matrix())).transpose(), gamma_Y.transpose()*Theta + gamma_z.transpose())*S_q*Gamma_z1.transpose()).transpose();
  
  // Correct the dimension to depend on X_t
  MatrixXd result((n + q) , 1); // create a new matrix with the desired size
  result << MatrixXd::Zero(n , 1), result_inter; // set the first n columns to zeros and the remaining columns to the matrix M1_z
  
  return result;
}

// [[Rcpp::export]]
MatrixXd alpha_dot_dot_dot_1(MatrixXd gamma, List Model, MatrixXd S_q, MatrixXd S_q_tilde) {
  // Extract elements of interest in "Model"
  int n = Model["n"];
  int q = Model["q"];
  MatrixXd Theta = Model["Theta"];
  MatrixXd nu = Model["nu"];
  MatrixXd mu = Model["mu"];
  MatrixXd Gamma_z0 = Model["Gamma.z0"];
  MatrixXd Gamma_Y1 = Model["Gamma.Y1"];
  
  // Create gamma_Y and gamma_z
  MatrixXd gamma_Y = gamma.block(0, 0, n, 1);
  MatrixXd gamma_z = gamma.block(n, 0, q, 1);
  
  // Compute alpha_dot_dot_dot_1
  MatrixXd result = 6*(gamma_Y.transpose()*Theta + gamma_z.transpose()).array().pow(4).matrix()*(nu.array()*mu.array().pow(4)).matrix() +
    3*kroneckerProduct(Gamma_Y1*(gamma_Y.array().square().matrix()), Gamma_Y1*(gamma_Y.array().square().matrix())).transpose()*S_q*Gamma_z0 +
    12*kroneckerProduct(gamma_Y.transpose()*Theta + gamma_z.transpose(),kroneckerProduct((Gamma_Y1*(gamma_Y.array().square().matrix())).transpose(), gamma_Y.transpose()*Theta + gamma_z.transpose()))*S_q_tilde*(nu.array()*mu.array().cube()).matrix();

  return result;
}

// [[Rcpp::export]]
MatrixXd beta_dot_dot_dot_1(MatrixXd gamma, List Model, MatrixXd S_q, MatrixXd S_q_tilde) {
  // Extract elements of interest in "Model"
  int n = Model["n"];
  int q = Model["q"];
  MatrixXd Theta = Model["Theta"];
  MatrixXd nu = Model["nu"];
  MatrixXd mu = Model["mu"];
  MatrixXd phi = Model["phi"];
  MatrixXd Gamma_z1 = Model["Gamma.z1"];
  MatrixXd Gamma_Y1 = Model["Gamma.Y1"];
  
  // Create gamma_Y and gamma_z
  MatrixXd gamma_Y = gamma.block(0, 0, n, 1);
  MatrixXd gamma_z = gamma.block(n, 0, q, 1);
  
  // Fourth cumulants of z use mu^4; mixed fourth-cumulant terms use the
  // third cumulants of z and therefore mu^3.
  MatrixXd diag_mu_4 = mu.array().pow(4).matrix().asDiagonal().toDenseMatrix();
  MatrixXd diag_mu_3 = mu.array().cube().matrix().asDiagonal().toDenseMatrix();
  
  // Compute alpha_dot_dot_dot_1
  MatrixXd result_inter = (24*(gamma_Y.transpose()*Theta + gamma_z.transpose()).array().pow(4).matrix()*(diag_mu_4*phi.transpose()) +
    3*kroneckerProduct(Gamma_Y1*(gamma_Y.array().square().matrix()), Gamma_Y1*(gamma_Y.array().square().matrix())).transpose()*S_q*Gamma_z1.transpose() +
    36*kroneckerProduct(gamma_Y.transpose()*Theta + gamma_z.transpose(),kroneckerProduct((Gamma_Y1*(gamma_Y.array().square().matrix())).transpose(), gamma_Y.transpose()*Theta + gamma_z.transpose()))*S_q_tilde*(diag_mu_3*phi.transpose())).transpose();
  
  // Correct the dimension to depend on X_t
  MatrixXd result((n + q) , 1); // create a new matrix with the desired size
  result << MatrixXd::Zero(n , 1), result_inter; // set the first n columns to zeros and the remaining columns to the matrix M1_z
  
  return result;
}

// [[Rcpp::export]]
NumericVector arrayC(NumericVector input, IntegerVector dim) { 
  input.attr("dim") = dim;
  return input;}

// [[Rcpp::export]]
NumericVector matrixXdToNumericVector(const MatrixXd& mat) {
  // Create a Map object that shares the data with the MatrixXd object
  Map<const MatrixXd> mappedMat(const_cast<double*>(mat.data()), mat.rows(), mat.cols());
  
  // Convert the mapped MatrixXd object to a NumericVector
  return wrap(mappedMat);
}


// [[Rcpp::export]]
MatrixXd test_1(MatrixXd Gamma, List Model, int freq = 12) {
  
  int max_H = 20;
  int n = 19;
  int q = 5;
  MatrixXd Gamma_2nd_case = MatrixXd::Zero(max_H, n + q);
  Gamma_2nd_case.row(max_H - 1) = Map<RowVectorXd>(Gamma.data(), Gamma.size());
  
  
  // 1st case
  MatrixXd Gamma_1st_case = Gamma_2nd_case;
  /// Extract the number elements to be display in aux_index and find the min value to put.
  int nb_elements = max_H/freq+1;
  int min_seq = max_H - (nb_elements-1)*freq;
  
  /// Create aux_index.
  VectorXd aux_index = VectorXd::LinSpaced(nb_elements, max_H, min_seq); // first element is the number of element of the vector, second is the min, last is the max 
  /// Create a matrix of ones with dimensions length(aux.index) x 1
  MatrixXd ones = MatrixXd::Ones(nb_elements, 1);
  MatrixXd replacement = ones * Gamma.transpose(); 
  
  for (int i = 0; i < aux_index.size(); i++) {
    Gamma_1st_case.row(aux_index(i)-1) = replacement.row(i);
  }
  
  //Gamma_1st_case.block(aux_index(nb_elements-1)-1, 0, nb_elements-1, n + q) = ones * Gamma.transpose();
  //Gamma_1st_case.row(aux_index(0)-1) = ones * Gamma.transpose();
  //Gamma_1st_case.row(aux_index(1)-1) = ones * Gamma.transpose();
  
  // 3rd case
  MatrixXd Gamma_3rd_case = MatrixXd::Zero(max_H, n + q);

  VectorXd aux_index_3rd = VectorXd::LinSpaced(4, max_H,max_H-(freq-freq/4));
  /// Create a matrix of ones with dimensions length(aux.index) x 1
  MatrixXd replacement_3rd = MatrixXd::Ones(4, 1) * Gamma.transpose() * 1/4; 
  
  /// Replace at the right positions
  for (int i = 0; i < aux_index_3rd.size(); i++) {
    Gamma_3rd_case.row(aux_index_3rd(i)-1) = replacement_3rd.row(i);
  } 
  
  // 4th case
  MatrixXd Gamma_4th_case = MatrixXd::Zero(max_H, n + q);
  
  VectorXd aux_index_4th = VectorXd::LinSpaced(freq, max_H,max_H-(freq-1));
  /// Create a matrix of ones with dimensions length(aux.index) x 1
  MatrixXd replacement_4th = MatrixXd::Ones(freq, 1) * Gamma.transpose() * 1/freq; 
  
  /// Replace at the right positions
  for (int i = 0; i < aux_index_4th.size(); i++) {
    Gamma_4th_case.row(aux_index_4th(i)-1) = replacement_4th.row(i);
  } 

  MatrixXd  a_h_1st_case = a_1(Gamma_1st_case.row(max_H-1).transpose(), Model); 

//return a_h_1st_case;
return aux_index_3rd;
}



// [[Rcpp::export]]
Rcpp::List compute_loadings(MatrixXd Gamma, List Model, VectorXd H, MatrixXd S_n, MatrixXd S_q, MatrixXd S_q_tilde, int indic_5y_in_Xy = 0, int freq = 12) {
  
  int n = Model["n"];
  int q = Model["q"];
  
  // Extract the length of H
  int nb_horizons = H.size();
  // Get maximum value of H
  int max_H = H.maxCoeff();
  
  // Extract Lambda.0, Lambda.1, M.0, M.1.
  MatrixXd Lambda_0 = Model["Lambda.0"];
  MatrixXd Lambda_1 = Model["Lambda.1"];
  MatrixXd M_0 = Model["M.0"];
  MatrixXd M_1 = Model["M.1"];
  
  // Create matrix for the loadings
  // 1st case
  MatrixXd a_h_1st_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd b_h_1st_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_h_1st_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_h_1st_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_dot_dot_h_1st_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_dot_dot_h_1st_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_dot_dot_dot_h_1st_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_dot_dot_dot_h_1st_case_mat = MatrixXd::Zero(max_H, n+q);
  
  // 2nd case
  MatrixXd a_h_2nd_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd b_h_2nd_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_h_2nd_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_h_2nd_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_dot_dot_h_2nd_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_dot_dot_h_2nd_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_dot_dot_dot_h_2nd_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_dot_dot_dot_h_2nd_case_mat = MatrixXd::Zero(max_H, n+q);
  
  // 3rd case
  MatrixXd a_h_3rd_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd b_h_3rd_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_h_3rd_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_h_3rd_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_dot_dot_h_3rd_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_dot_dot_h_3rd_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_dot_dot_dot_h_3rd_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_dot_dot_dot_h_3rd_case_mat = MatrixXd::Zero(max_H, n+q);
  
  // 4th case
  MatrixXd a_h_4th_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd b_h_4th_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_h_4th_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_h_4th_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_dot_dot_h_4th_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_dot_dot_h_4th_case_mat = MatrixXd::Zero(max_H, n+q);
  MatrixXd alpha_dot_dot_dot_h_4th_case_mat = MatrixXd::Zero(max_H, 1);
  MatrixXd beta_dot_dot_dot_h_4th_case_mat = MatrixXd::Zero(max_H, n+q);
  
  // Create the Gamma
  // 2nd case
  MatrixXd Gamma_2nd_case = MatrixXd::Zero(max_H, n + q);
  Gamma_2nd_case.row(max_H - 1) = Map<RowVectorXd>(Gamma.data(), Gamma.size());
  
  // 1st case
  MatrixXd Gamma_1st_case = Gamma_2nd_case;
  /// Extract the number elements to be display in aux_index and find the min value to put.
  int nb_elements = max_H/freq; //max_H/freq+1
  int min_seq = max_H - (nb_elements-1)*freq;
  
  /// Create aux_index.
  VectorXd aux_index = VectorXd::LinSpaced(nb_elements, min_seq, max_H); // VectorXd::LinSpaced(nb_elements, max_H, min_seq) first element is the number of element of the vector, second is the min, last is the max 
  /// Create a matrix of ones with dimensions length(aux.index) x 1
  MatrixXd ones = MatrixXd::Ones(nb_elements, 1);
  MatrixXd replacement = ones * Gamma.transpose(); 
  
  /// Replace at the right positions
  for (int i = 0; i < aux_index.size(); i++) {
    Gamma_1st_case.row(aux_index(i)-1) = replacement.row(i);
  }  
  
  // 3rd case
  MatrixXd Gamma_3rd_case = MatrixXd::Zero(max_H, n + q);

  VectorXd aux_index_3rd = VectorXd::LinSpaced(4, max_H-(freq-freq/4),max_H);//VectorXd::LinSpaced(4, max_H,max_H-(freq-freq/4))
  /// Create a matrix of ones with dimensions length(aux.index) x 1
  MatrixXd replacement_3rd = MatrixXd::Ones(4, 1) * Gamma.transpose() * 1/4; 
  
  /// Replace at the right positions
  for (int i = 0; i < aux_index_3rd.size(); i++) {
    Gamma_3rd_case.row(aux_index_3rd(i)-1) = replacement_3rd.row(i);
  } 
  
  // 4th case
  MatrixXd Gamma_4th_case = MatrixXd::Zero(max_H, n + q);

  VectorXd aux_index_4th = VectorXd::LinSpaced(freq, max_H-(freq-1), max_H); //VectorXd::LinSpaced(freq, max_H,max_H-(freq-1))
  /// Create a matrix of ones with dimensions length(aux.index) x 1
  MatrixXd replacement_4th = MatrixXd::Ones(freq, 1) * Gamma.transpose() * 1/freq; 
  
  /// Replace at the right positions
  for (int i = 0; i < aux_index_4th.size(); i++) {
    Gamma_4th_case.row(aux_index_4th(i)-1) = replacement_4th.row(i);
  } 
  
  
  if (indic_5y_in_Xy == 1 && Gamma_1st_case.rows() > freq*5) {
    Gamma_1st_case.block(0, 0, Gamma_1st_case.rows() - freq*5, Gamma_1st_case.cols()).setZero();
  }
  
  int count = 0;
  
  // Define all matrices
  
  MatrixXd a_Xd = MatrixXd::Zero(4,nb_horizons);
  //MatrixXd b_Xd = MatrixXd::Zero((n+q)*4,nb_horizons);
  MatrixXd b_Xd = MatrixXd::Zero((n+q)*nb_horizons,4);
  MatrixXd alpha_Xd = MatrixXd::Zero(4,nb_horizons);
  MatrixXd beta_Xd = MatrixXd::Zero((n+q)*nb_horizons,4);;
  MatrixXd alpha_dot_dot_Xd = MatrixXd::Zero(4,nb_horizons);
  MatrixXd beta_dot_dot_Xd = MatrixXd::Zero((n+q)*nb_horizons,4);;
  MatrixXd alpha_dot_dot_dot_Xd = MatrixXd::Zero(4,nb_horizons);
  MatrixXd beta_dot_dot_dot_Xd = MatrixXd::Zero((n+q)*nb_horizons,4);;
  
  //  1st case
  MatrixXd a_h_1st_case;
  MatrixXd b_h_1st_case;
  MatrixXd alpha_h_1st_case;
  MatrixXd beta_h_1st_case;
  MatrixXd alpha_dot_dot_h_1st_case;
  MatrixXd beta_dot_dot_h_1st_case;
  MatrixXd alpha_dot_dot_dot_h_1st_case;
  MatrixXd beta_dot_dot_dot_h_1st_case;
  
  //  2nd case
  MatrixXd a_h_2nd_case;
  MatrixXd b_h_2nd_case;
  MatrixXd alpha_h_2nd_case;
  MatrixXd beta_h_2nd_case;
  MatrixXd alpha_dot_dot_h_2nd_case;
  MatrixXd beta_dot_dot_h_2nd_case;
  MatrixXd alpha_dot_dot_dot_h_2nd_case;
  MatrixXd beta_dot_dot_dot_h_2nd_case;
  
  //  3rd case
  MatrixXd a_h_3rd_case;
  MatrixXd b_h_3rd_case;
  MatrixXd alpha_h_3rd_case;
  MatrixXd beta_h_3rd_case;
  MatrixXd alpha_dot_dot_h_3rd_case;
  MatrixXd beta_dot_dot_h_3rd_case;
  MatrixXd alpha_dot_dot_dot_h_3rd_case;
  MatrixXd beta_dot_dot_dot_h_3rd_case;
  
  //  4th case
  MatrixXd a_h_4th_case;
  MatrixXd b_h_4th_case;
  MatrixXd alpha_h_4th_case;
  MatrixXd beta_h_4th_case;
  MatrixXd alpha_dot_dot_h_4th_case;
  MatrixXd beta_dot_dot_h_4th_case;
  MatrixXd alpha_dot_dot_dot_h_4th_case;
  MatrixXd beta_dot_dot_dot_h_4th_case;
  
  for (int h_ind = 1; h_ind <= max_H; h_ind++) {
    
    if (h_ind == 1) {

      // Create the a_1
      //  1st case
      a_h_1st_case = a_1(Gamma_1st_case.row(max_H-1).transpose(), Model);
      b_h_1st_case = b_1(Gamma_1st_case.row(max_H-1).transpose(), Model);
      alpha_h_1st_case = alpha_1(Gamma_1st_case.row(max_H-1).transpose(),Model,S_n,S_q);
      beta_h_1st_case  = beta_1(Gamma_1st_case.row(max_H-1).transpose(),Model,S_n,S_q);
      alpha_dot_dot_h_1st_case = alpha_dot_dot_1(Gamma_1st_case.row(max_H-1).transpose(),Model,S_q);
      beta_dot_dot_h_1st_case = beta_dot_dot_1(Gamma_1st_case.row(max_H-1).transpose(),Model,S_q);
      alpha_dot_dot_dot_h_1st_case = alpha_dot_dot_dot_1(Gamma_1st_case.row(max_H-1).transpose(),Model,S_q, S_q_tilde);
      beta_dot_dot_dot_h_1st_case = beta_dot_dot_dot_1(Gamma_1st_case.row(max_H-1).transpose(),Model,S_q, S_q_tilde);
        
      // 2nd case
      a_h_2nd_case = a_h_1st_case;
      b_h_2nd_case = b_h_1st_case;
      alpha_h_2nd_case = alpha_h_1st_case;
      beta_h_2nd_case  = beta_h_1st_case;
      alpha_dot_dot_h_2nd_case = alpha_dot_dot_h_1st_case;
      beta_dot_dot_h_2nd_case  = beta_dot_dot_h_1st_case;
      alpha_dot_dot_dot_h_2nd_case = alpha_dot_dot_dot_h_1st_case;
      beta_dot_dot_dot_h_2nd_case  = beta_dot_dot_dot_h_1st_case;
       
      //  3rd case
      a_h_3rd_case = a_1(Gamma_3rd_case.row(max_H-1).transpose(), Model); 
      b_h_3rd_case = b_1(Gamma_3rd_case.row(max_H-1).transpose(), Model);
      alpha_h_3rd_case = alpha_1(Gamma_3rd_case.row(max_H-1).transpose(),Model,S_n,S_q);
      beta_h_3rd_case  = beta_1(Gamma_3rd_case.row(max_H-1).transpose(),Model,S_n,S_q);
      alpha_dot_dot_h_3rd_case = alpha_dot_dot_1(Gamma_3rd_case.row(max_H-1).transpose(),Model,S_q);
      beta_dot_dot_h_3rd_case = beta_dot_dot_1(Gamma_3rd_case.row(max_H-1).transpose(),Model,S_q);
      alpha_dot_dot_dot_h_3rd_case = alpha_dot_dot_dot_1(Gamma_3rd_case.row(max_H-1).transpose(),Model,S_q, S_q_tilde);
      beta_dot_dot_dot_h_3rd_case = beta_dot_dot_dot_1(Gamma_3rd_case.row(max_H-1).transpose(),Model,S_q, S_q_tilde);
      
      //  4th case
      a_h_4th_case = a_1(Gamma_4th_case.row(max_H-1).transpose(), Model); 
      b_h_4th_case = b_1(Gamma_4th_case.row(max_H-1).transpose(), Model);
      alpha_h_4th_case = alpha_1(Gamma_4th_case.row(max_H-1).transpose(),Model,S_n,S_q);
      beta_h_4th_case  = beta_1(Gamma_4th_case.row(max_H-1).transpose(),Model,S_n,S_q);
      alpha_dot_dot_h_4th_case = alpha_dot_dot_1(Gamma_4th_case.row(max_H-1).transpose(),Model,S_q);
      beta_dot_dot_h_4th_case = beta_dot_dot_1(Gamma_4th_case.row(max_H-1).transpose(),Model,S_q);
      alpha_dot_dot_dot_h_4th_case = alpha_dot_dot_dot_1(Gamma_4th_case.row(max_H-1).transpose(),Model,S_q, S_q_tilde);
      beta_dot_dot_dot_h_4th_case = beta_dot_dot_dot_1(Gamma_4th_case.row(max_H-1).transpose(),Model,S_q, S_q_tilde);
      
    } else{
      // Create the a_h for h > 1
      // 1st case
      /// Alphas
      a_h_1st_case = a_h_1st_case + a_1(Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case, Model);
      alpha_h_1st_case = alpha_1(Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case,Model,S_n,S_q) +
        alpha_h_1st_case + a_1(beta_h_1st_case,Model);
      alpha_dot_dot_h_1st_case = alpha_dot_dot_h_1st_case + a_1(beta_dot_dot_h_1st_case,Model) +
        alpha_dot_dot_1(Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case,Model,S_q) +
        3*kroneckerProduct(beta_h_1st_case, Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case).transpose()*Lambda_0;
      alpha_dot_dot_dot_h_1st_case = alpha_dot_dot_dot_h_1st_case + a_1(beta_dot_dot_dot_h_1st_case,Model) +
        alpha_dot_dot_dot_1(Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case,Model,S_q, S_q_tilde) +
        4*kroneckerProduct(beta_dot_dot_h_1st_case, Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case).transpose()*Lambda_0 +
        3*alpha_1(beta_h_1st_case,Model,S_n, S_q) +
        6*kroneckerProduct((Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case), kroneckerProduct(beta_h_1st_case, (Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case))).transpose()*M_0;
      
      /// Betas 
      beta_dot_dot_dot_h_1st_case =  b_1(beta_dot_dot_dot_h_1st_case,Model) +
        beta_dot_dot_dot_1(Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case,Model,S_q,S_q_tilde) +
        4*Lambda_1.transpose()*kroneckerProduct(beta_dot_dot_h_1st_case, Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case) +
        3*beta_1(beta_h_1st_case,Model,S_n, S_q) +
        6*M_1.transpose()*kroneckerProduct(Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case, kroneckerProduct(beta_h_1st_case, Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case));
      beta_dot_dot_h_1st_case =  b_1(beta_dot_dot_h_1st_case,Model) +
        beta_dot_dot_1(Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case,Model,S_q) +
        3*Lambda_1.transpose()*kroneckerProduct(beta_h_1st_case, Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case);
      beta_h_1st_case = beta_1(Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case,Model,S_n,S_q) +
        b_1(beta_h_1st_case,Model);
      b_h_1st_case = b_1(Gamma_1st_case.row(max_H-h_ind).transpose() + b_h_1st_case,Model);
      
      // 2nd case
      /// Alphas
      a_h_2nd_case = a_h_2nd_case + a_1(Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case, Model);
      alpha_h_2nd_case = alpha_1(Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case,Model,S_n,S_q) +
        alpha_h_2nd_case + a_1(beta_h_2nd_case,Model);
      alpha_dot_dot_h_2nd_case = alpha_dot_dot_h_2nd_case + a_1(beta_dot_dot_h_2nd_case,Model) +
        alpha_dot_dot_1(Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case,Model,S_q) +
        3*kroneckerProduct(beta_h_2nd_case, Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case).transpose()*Lambda_0;
      alpha_dot_dot_dot_h_2nd_case = alpha_dot_dot_dot_h_2nd_case + a_1(beta_dot_dot_dot_h_2nd_case,Model) +
        alpha_dot_dot_dot_1(Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case,Model,S_q, S_q_tilde) +
        4*kroneckerProduct(beta_dot_dot_h_2nd_case, Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case).transpose()*Lambda_0 +
        3*alpha_1(beta_h_2nd_case,Model,S_n, S_q) +
        6*kroneckerProduct((Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case), kroneckerProduct(beta_h_2nd_case, (Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case))).transpose()*M_0;
      
      /// Betas 
      beta_dot_dot_dot_h_2nd_case =  b_1(beta_dot_dot_dot_h_2nd_case,Model) +
        beta_dot_dot_dot_1(Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case,Model,S_q,S_q_tilde) +
        4*Lambda_1.transpose()*kroneckerProduct(beta_dot_dot_h_2nd_case, Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case) +
        3*beta_1(beta_h_2nd_case,Model,S_n, S_q) +
        6*M_1.transpose()*kroneckerProduct(Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case, kroneckerProduct(beta_h_2nd_case, Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case));
      beta_dot_dot_h_2nd_case =  b_1(beta_dot_dot_h_2nd_case,Model) +
        beta_dot_dot_1(Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case,Model,S_q) +
        3*Lambda_1.transpose()*kroneckerProduct(beta_h_2nd_case, Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case);
      beta_h_2nd_case = beta_1(Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case,Model,S_n,S_q) +
        b_1(beta_h_2nd_case,Model);
      b_h_2nd_case = b_1(Gamma_2nd_case.row(max_H-h_ind).transpose() + b_h_2nd_case,Model);
  
      // 3rd case
      /// Alphas
      a_h_3rd_case = a_h_3rd_case + a_1(Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case, Model);
      alpha_h_3rd_case = alpha_1(Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case,Model,S_n,S_q) +
        alpha_h_3rd_case + a_1(beta_h_3rd_case,Model);
      alpha_dot_dot_h_3rd_case = alpha_dot_dot_h_3rd_case + a_1(beta_dot_dot_h_3rd_case,Model) +
        alpha_dot_dot_1(Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case,Model,S_q) +
        3*kroneckerProduct(beta_h_3rd_case, Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case).transpose()*Lambda_0;
      alpha_dot_dot_dot_h_3rd_case = alpha_dot_dot_dot_h_3rd_case + a_1(beta_dot_dot_dot_h_3rd_case,Model) +
        alpha_dot_dot_dot_1(Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case,Model,S_q, S_q_tilde) +
        4*kroneckerProduct(beta_dot_dot_h_3rd_case, Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case).transpose()*Lambda_0 +
        3*alpha_1(beta_h_3rd_case,Model,S_n, S_q) +
        6*kroneckerProduct((Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case), kroneckerProduct(beta_h_3rd_case, (Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case))).transpose()*M_0;
  
      /// Betas 
      beta_dot_dot_dot_h_3rd_case =  b_1(beta_dot_dot_dot_h_3rd_case,Model) +
        beta_dot_dot_dot_1(Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case,Model,S_q,S_q_tilde) +
        4*Lambda_1.transpose()*kroneckerProduct(beta_dot_dot_h_3rd_case, Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case) +
        3*beta_1(beta_h_3rd_case,Model,S_n, S_q) +
        6*M_1.transpose()*kroneckerProduct(Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case, kroneckerProduct(beta_h_3rd_case, Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case));
      beta_dot_dot_h_3rd_case =  b_1(beta_dot_dot_h_3rd_case,Model) +
        beta_dot_dot_1(Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case,Model,S_q) +
        3*Lambda_1.transpose()*kroneckerProduct(beta_h_3rd_case, Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case);
        beta_h_3rd_case = beta_1(Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case,Model,S_n,S_q) +
        b_1(beta_h_3rd_case,Model);
      b_h_3rd_case = b_1(Gamma_3rd_case.row(max_H-h_ind).transpose() + b_h_3rd_case,Model);
  
      // 4th case
      /// Alphas
      a_h_4th_case = a_h_4th_case + a_1(Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case, Model);
      alpha_h_4th_case = alpha_1(Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case,Model,S_n,S_q) +
        alpha_h_4th_case + a_1(beta_h_4th_case,Model);
      alpha_dot_dot_h_4th_case = alpha_dot_dot_h_4th_case + a_1(beta_dot_dot_h_4th_case,Model) +
        alpha_dot_dot_1(Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case,Model,S_q) +
        3*kroneckerProduct(beta_h_4th_case, Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case).transpose()*Lambda_0;
      alpha_dot_dot_dot_h_4th_case = alpha_dot_dot_dot_h_4th_case + a_1(beta_dot_dot_dot_h_4th_case,Model) +
        alpha_dot_dot_dot_1(Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case,Model,S_q, S_q_tilde) +
        4*kroneckerProduct(beta_dot_dot_h_4th_case, Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case).transpose()*Lambda_0 +
        3*alpha_1(beta_h_4th_case,Model,S_n, S_q) +
        6*kroneckerProduct((Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case), kroneckerProduct(beta_h_4th_case, (Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case))).transpose()*M_0;
  
      /// Betas 
      beta_dot_dot_dot_h_4th_case =  b_1(beta_dot_dot_dot_h_4th_case,Model) +
        beta_dot_dot_dot_1(Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case,Model,S_q,S_q_tilde) +
        4*Lambda_1.transpose()*kroneckerProduct(beta_dot_dot_h_4th_case, Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case) +
        3*beta_1(beta_h_4th_case,Model,S_n, S_q) +
        6*M_1.transpose()*kroneckerProduct(Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case, kroneckerProduct(beta_h_4th_case, Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case));
      beta_dot_dot_h_4th_case =  b_1(beta_dot_dot_h_4th_case,Model) +
        beta_dot_dot_1(Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case,Model,S_q) +
        3*Lambda_1.transpose()*kroneckerProduct(beta_h_4th_case, Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case);
      beta_h_4th_case = beta_1(Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case,Model,S_n,S_q) +
        b_1(beta_h_4th_case,Model);
      b_h_4th_case = b_1(Gamma_4th_case.row(max_H-h_ind).transpose() + b_h_4th_case,Model);

    }
    
    if (std::find(H.data(), H.data() + H.size(), h_ind) != H.data() + H.size()) {
      // h_ind is found in H, print it
      //std::cout << "h_ind = " << h_ind << std::endl;
      
      // in that case, the loadings are expected in the output
      
      
      // Deflator for the first case
      double deflator = (1 - indic_5y_in_Xy) * (1 + static_cast<int>((h_ind - 1) / freq)) +
        indic_5y_in_Xy * std::min(5, 1 + static_cast<int>((h_ind - 1) / freq));
      
      //std::cout << "Deflator: " << deflator << std::endl;
      
      
     // a[count] = matrixXdToNumericVector(a_h_1st_case); //element [,,1][1,1]
     
      // 1st case
      a_Xd(0,count) = a_h_1st_case(0,0) / deflator;
      alpha_Xd(0,count) = alpha_h_1st_case(0,0) / std::pow(deflator, 2);
      alpha_dot_dot_Xd(0,count) = alpha_dot_dot_h_1st_case(0,0) / std::pow(deflator, 3);
      alpha_dot_dot_dot_Xd(0,count) = alpha_dot_dot_dot_h_1st_case(0,0) / std::pow(deflator, 4);
      
      //b_Xd.block(0, count, n+q, 1) = b_h_1st_case / deflator;
      b_Xd.block(count*(n+q), 0, n+q, 1) = b_h_1st_case / deflator;
      beta_Xd.block(count*(n+q), 0, n+q, 1) = beta_h_1st_case / std::pow(deflator, 2);
      beta_dot_dot_Xd.block(count*(n+q), 0, n+q, 1) = beta_dot_dot_h_1st_case / std::pow(deflator, 3);
      beta_dot_dot_dot_Xd.block(count*(n+q), 0, n+q, 1) = beta_dot_dot_dot_h_1st_case / std::pow(deflator, 4);
      
      
      // 2nd case
      a_Xd(1,count) = a_h_2nd_case(0,0);
      alpha_Xd(1,count) = alpha_h_2nd_case(0,0);
      alpha_dot_dot_Xd(1,count) = alpha_dot_dot_h_2nd_case(0,0);
      alpha_dot_dot_dot_Xd(1,count) = alpha_dot_dot_dot_h_2nd_case(0,0);
      
      //b_Xd.block(n+q, count, n+q, 1) = b_h_2nd_case;
      b_Xd.block(count*(n+q), 1, n+q, 1)  = b_h_2nd_case;
      beta_Xd.block(count*(n+q), 1, n+q, 1)  = beta_h_2nd_case;
      beta_dot_dot_Xd.block(count*(n+q), 1, n+q, 1)  = beta_dot_dot_h_2nd_case;
      beta_dot_dot_dot_Xd.block(count*(n+q), 1, n+q, 1)  = beta_dot_dot_dot_h_2nd_case;
      
      // 3rd case
      a_Xd(2,count) = a_h_3rd_case(0,0);
      alpha_Xd(2,count) = alpha_h_3rd_case(0,0);
      alpha_dot_dot_Xd(2,count) = alpha_dot_dot_h_3rd_case(0,0);
      alpha_dot_dot_dot_Xd(2,count) = alpha_dot_dot_dot_h_3rd_case(0,0);
      
      //b_Xd.block((n+q)*2, count, n+q, 1) = b_h_3rd_case;
      b_Xd.block(count*(n+q), 2, n+q, 1) = b_h_3rd_case;
      beta_Xd.block(count*(n+q), 2, n+q, 1) = beta_h_3rd_case;
      beta_dot_dot_Xd.block(count*(n+q), 2, n+q, 1) = beta_dot_dot_h_3rd_case;
      beta_dot_dot_dot_Xd.block(count*(n+q), 2, n+q, 1) = beta_dot_dot_dot_h_3rd_case;
      
    
      // 4th case
      a_Xd(3,count) = a_h_4th_case(0,0);
      alpha_Xd(3,count) = alpha_h_4th_case(0,0);
      alpha_dot_dot_Xd(3,count) = alpha_dot_dot_h_4th_case(0,0);
      alpha_dot_dot_dot_Xd(3,count) = alpha_dot_dot_dot_h_4th_case(0,0);
      
      //b_Xd.block((n+q)*3, count, n+q, 1) = b_h_4th_case;
      b_Xd.block(count*(n+q), 3, n+q, 1) = b_h_4th_case;
      beta_Xd.block(count*(n+q), 3, n+q, 1) = beta_h_4th_case;
      beta_dot_dot_Xd.block(count*(n+q), 3, n+q, 1) = beta_dot_dot_h_4th_case;
      beta_dot_dot_dot_Xd.block(count*(n+q), 3, n+q, 1) = beta_dot_dot_dot_h_4th_case;
      
      count = count + 1;
      
    }
    
  } 
 
  // Create the NumericVector with dimensions 1xhx4 and all elements set to 0
  int h = nb_horizons; // set h to the desired value
  //NumericVector myVec(h*4); // create a NumericVector with length h*4 filled with 0s
  //NumericVector myVec2(h*(n+q)*4); // create a NumericVector with length h*(n+q) filled with 0s
  
  IntegerVector dim = IntegerVector::create(1, h, 4); // create an IntegerVector with the desired dimensions
  IntegerVector dim2 = IntegerVector::create(n+q, h, 4); // create an IntegerVector with the desired dimensions
  
  // Create a, alpha, alpha_dot_dot, alpha_dot_dot_dot, b, beta, beta_dot_dot, beta_dot_dot_dot, 
  NumericVector a = arrayC(matrixXdToNumericVector(a_Xd.transpose()), dim);
  NumericVector alpha = arrayC(matrixXdToNumericVector(alpha_Xd.transpose()), dim);
  NumericVector alpha_dot_dot = arrayC(matrixXdToNumericVector(alpha_dot_dot_Xd.transpose()), dim);
  NumericVector alpha_dot_dot_dot = arrayC(matrixXdToNumericVector(alpha_dot_dot_dot_Xd.transpose()), dim);
  NumericVector b = arrayC(matrixXdToNumericVector(b_Xd), dim2);
  NumericVector beta = arrayC(matrixXdToNumericVector(beta_Xd), dim2);
  NumericVector beta_dot_dot = arrayC(matrixXdToNumericVector(beta_dot_dot_Xd), dim2);
  NumericVector beta_dot_dot_dot = arrayC(matrixXdToNumericVector(beta_dot_dot_dot_Xd), dim2);
  //NumericVector a = arrayC(myVec, dim); // set the dimensions of myVec using the arrayC function
  //NumericVector alpha = arrayC(myVec, dim); // set the dimensions of myVec using the arrayC function
  //NumericVector alpha_dot_dot = arrayC(myVec, dim); // set the dimensions of myVec using the arrayC function
  //NumericVector alpha_dot_dot_dot = arrayC(myVec, dim); // set the dimensions of myVec using the arrayC function
  //NumericVector b = arrayC(myVec2, dim2); // set the dimensions of myVec using the arrayC function
  //NumericVector beta = arrayC(myVec2, dim2); // set the dimensions of myVec using the arrayC function
  //NumericVector beta_dot_dot = arrayC(myVec2, dim2); // set the dimensions of myVec using the arrayC function
  //NumericVector beta_dot_dot_dot = arrayC(myVec2, dim2); // set the dimensions of myVec using the arrayC function
  
  
  
  return List::create(Named("a") = a, Named("b") = b, Named("alpha") = alpha, Named("beta") = beta,
                      Named("alpha.dot.dot") = alpha_dot_dot, Named("beta.dot.dot") = beta_dot_dot,
                      Named("alpha.dot.dot.dot") = alpha_dot_dot_dot, 
                      Named("beta.dot.dot.dot") = beta_dot_dot_dot, 
                      Named("Gamma.1st.case") = Gamma_1st_case, Named("Gamma.2nd.case") = Gamma_2nd_case,
                      Named("Gamma.3rd.case") = Gamma_3rd_case, Named("Gamma.4th.case") = Gamma_4th_case);
 
}

// 
// // [[Rcpp::export]]
// Rcpp::List compute_loadings(MatrixXd Gamma, List Model, VectorXd H, MatrixXd S_n, MatrixXd S_q, MatrixXd S_q_tilde, int indic_5y_in_Xy = 0) {
//   
//   int n = Model["n"];
//   int q = Model["q"];
//   
//   // Extract the length of H
//   int nb_horizons = H.size();
//   
//   // Create arrays of appropriate sizes
//   ArrayXXd a(1, nb_horizons, 4); // create a dynamic array of size 1 x h x 4
//   ArrayXXd alpha(1, nb_horizons, 4); // create a dynamic array of size 1 x h x 4
//   ArrayXXd alpha_dot_dot(1, nb_horizons, 4); // create a dynamic array of size 1 x h x 4
//   ArrayXXd alpha_dot_dot_dot(1, nb_horizons, 4); // create a dynamic array of size 1 x h x 4
//   ArrayXXd b(n+q, nb_horizons, 4); // create a dynamic array of size n+q x h x 4
//   ArrayXXd beta(n+q, nb_horizons, 4); // create a dynamic array of size n+q x h x 4
//   ArrayXXd beta_dot_dot(n+q, nb_horizons, 4); // create a dynamic array of size n+q x h x 4
//   ArrayXXd beta_dot_dot_dot(n+q, nb_horizons, 4); // create a dynamic array of size n+q x h x 4
//   
//   // Fill arrays with zeros
//   a.setZero();
//   alpha.setZero();
//   alpha_dot_dot.setZero();
//   alpha_dot_dot_dot.setZero();
//   b.setZero();
//   beta.setZero();
//   beta_dot_dot.setZero();
//   beta_dot_dot_dot.setZero();
//   
//   //int* a = new int[1 * nb_horizons * 4]; // create a dynamic array of size 1 x h x 4
//   //int* alpha = new int[1 * nb_horizons * 4]; // create a dynamic array of size 1 x h x 4
//   //int* alpha_dot_dot = new int[1 * nb_horizons * 4]; // create a dynamic array of size 1 x h x 4
//   //int* alpha_dot_dot_dot = new int[1 * nb_horizons * 4]; // create a dynamic array of size 1 x h x 4
//   //int* b = new int[(n+q) * nb_horizons * 4]; // create a dynamic array of size 1 x h x 4
//   //int* beta = new int[(n+q) * nb_horizons * 4]; // create a dynamic array of size 1 x h x 4
//   //int* beta_dot_dot = new int[(n+q) * nb_horizons * 4]; // create a dynamic array of size 1 x h x 4
//   //int* beta_dot_dot_dot = new int[(n+q) * nb_horizons * 4]; // create a dynamic array of size 1 x h x 4
//    
//   return List::create(Named("a") = a, Named("b") = b, Named("alpha") = alpha, Named("beta") = beta,
//                       Named("alpha.dot.dot") = alpha_dot_dot, Named("beta.dot.dot") = beta_dot_dot,
//                       Named("alpha.dot.dot.dot") = alpha_dot_dot_dot, 
//                       Named("beta.dot.dot.dot") = beta_dot_dot_dot);
//                                  
// }

