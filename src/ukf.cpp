#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#define M_PI 3.14159265358979323846

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.42;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.38;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Parameter lambda
  lambda_ = 3 - n_aug_;

  // weights to calculate the mean and covariance
  weights_ = VectorXd(2 * n_aug_ + 1); 

  // Prediction of sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  cout << "Processing measurement..." << endl;
  if (!is_initialized_){
    

    // Initializing state vector
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      
      x_(0) = rho*cos(phi);
      x_(1) = rho*sin(phi);
      x_(2) = 0.0;
      x_(3) = 0.0;
      x_(4) = 0.0; 

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
      x_(2) = 0.0;
      x_(3) = 0.0;
      x_(4) = 0.0;
    }

    // Initializing state covariance
    P_ = MatrixXd::Identity(n_x_, n_x_);

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
    std::cout << "Radar update: " << std::endl;
    UpdateRadar(meas_package);
    cout << "End radar update" << endl;
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
    std::cout << "Lidar update: " << std::endl;
    UpdateLidar(meas_package);
    cout << "End laser update" << endl;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  cout << "prediction step" << endl;
  // Initializing variables

  cout << "x_in = " << x_ << endl;
  cout << "P_in = " << P_ << endl;

  MatrixXd Xsig_aug_ = MatrixXd(n_aug_,2*n_aug_+1);
  MatrixXd P_aug_ = MatrixXd(n_aug_,n_aug_);
  VectorXd x_aug_ = VectorXd(n_aug_);

  // Fill augmented state with values
  x_aug_.head(5) = x_;
  x_aug_(5) = 0.0;
  x_aug_(6) = 0.0;

  //cout << "Stage 1" << endl;

  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(n_x_, n_x_) = std_a_*std_a_;
  P_aug_(n_x_ + 1,n_x_ + 1) = std_yawdd_*std_yawdd_;

  // Get the square root of the covariance matrix
  MatrixXd L = P_aug_.llt().matrixL();

  Xsig_aug_.fill(0.0);
  Xsig_aug_.col(0) = x_aug_;
  for(int i = 0; i < n_aug_; i++){
    Xsig_aug_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_)*L.col(i);
    Xsig_aug_.col(n_aug_+i+1) = x_aug_ - sqrt(lambda_ + n_aug_)*L.col(i);
  }

  //cout << "Stage 2" << endl;

  Xsig_pred_.fill(0.0);
  for(int j=0; j < 2*n_aug_ + 1; j++){
    // Identify values of the stae vector
    double p_x = Xsig_aug_(0,j);
    double p_y = Xsig_aug_(1,j);
    double v = Xsig_aug_(2,j);
    double yaw = Xsig_aug_(3,j);
    double yawd = Xsig_aug_(4,j);
    double nu_a = Xsig_aug_(5,j);
    double nu_yawdd = Xsig_aug_(6,j);
    //cout << "Stage 2.1" << endl;
    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.0001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
    //cout << "Stage 2.2" << endl;
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
    //cout << "Stage 2.3" << endl;

    //write predicted sigma point into right column
    
    Xsig_pred_(0,j) = px_p;
    Xsig_pred_(1,j) = py_p;
    Xsig_pred_(2,j) = v_p;
    Xsig_pred_(3,j) = yaw_p;
    Xsig_pred_(4,j) = yawd_p;
    //cout << "Stage 2.4" << endl;
  }

  //cout << "Stage 3" << endl;

  // To predict the state vector we have to set the weights
  weights_(0) = lambda_/(lambda_ + n_aug_);
  weights_.tail(2*n_aug_).setConstant(0.5/(n_aug_ + lambda_));
  cout << "weights = " << weights_ << endl;
  //cout << "Stage 4" << endl;

  x_.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){
    x_ += weights_(i)*Xsig_pred_.col(i);
  }


  //cout << "Stage 5" << endl;

  P_.fill(0.0);

  for(int i = 0; i < 2*n_aug_+1; i++){
    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    if(x_diff(3)>M_PI) x_diff(3)-= 2.*M_PI;
    if(x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

    P_ = P_ + weights_(i)*x_diff*x_diff.transpose();

    cout << "x_out = " << x_ << endl;
    cout << "P_out = " << P_ << endl;
  }

cout << "End" << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.


  You'll also need to calculate the lidar NIS.
  */

  // Variables for the LASER measurement
  //cout << "LASER_1" << endl;
  int n_z_laser = 2;
  MatrixXd H_laser_ = MatrixXd(n_z_laser, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;
  //cout << "LASER_2" << endl;
  MatrixXd R_laser_ = MatrixXd(n_z_laser, n_z_laser);
  R_laser_.fill(0.0);
  R_laser_(0,0) = std_a_*std_a_;
  R_laser_(1,1) = std_yawdd_*std_yawdd_;

  VectorXd z_laser_ = VectorXd(n_z_laser);
  z_laser_(0) = meas_package.raw_measurements_[0];
  z_laser_(1) = meas_package.raw_measurements_[1];
  //cout << "LASER_3" << endl;

  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size,x_size);
  VectorXd y_laser_ = z_laser_ - H_laser_*x_;
  MatrixXd S_laser_ = H_laser_*P_*H_laser_.transpose() + R_laser_;
  MatrixXd K_laser_ = P_*H_laser_.transpose()*S_laser_.inverse();
  //cout << "LASER_4" << endl;

  // Updating state and covariance matrix from LASER measurement
  x_ = x_ + (K_laser_*y_laser_);
  P_ = (I - K_laser_*H_laser_)*P_;
  //cout << "LASER_5" << endl;

  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // Variables for the RADAR measurement
  int n_z_radar = 3;
  MatrixXd Zsig_radar_ = MatrixXd(n_z_radar, 2*n_aug_+1);
  MatrixXd R_radar_ = MatrixXd(n_z_radar, n_z_radar);
  MatrixXd S_radar_ = MatrixXd(n_z_radar, n_z_radar);
  MatrixXd K_radar_ = MatrixXd(n_x_, n_z_radar);
  MatrixXd Tc_radar_ = MatrixXd(n_x_, n_z_radar);

  // transform sigma points into the measurement space
  Zsig_radar_.fill(0.0);
  for(int i=0; i < 2*n_aug_+1; i++){
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v_x = cos(yaw)*v;
    double v_y = sin(yaw)*v;

    // measurement model
    
    Zsig_radar_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_radar_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_radar_(2,i) = (p_x*v_x + p_y*v_y ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot

    // Normalize phi
    if (Zsig_radar_(1,i)> M_PI) Zsig_radar_(1,i) -=2.*M_PI;
    if (Zsig_radar_(1,i)<-M_PI) Zsig_radar_(1,i) +=2.*M_PI;
  }
  cout << "Zsig_radar_ = " << Zsig_radar_ << endl;

  //mean predicted measurement
  VectorXd z_pred_radar_ = VectorXd(3);
  z_pred_radar_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_radar_ = z_pred_radar_ + weights_(i) * Zsig_radar_.col(i);
  }

  cout << "z_pred_radar_ = " << z_pred_radar_ << endl;

  //measurement covariance matrix S
  S_radar_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;

    //angle normalization
    if (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    if (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S_radar_ = S_radar_ + weights_(i) * z_diff * z_diff.transpose();
  }

  cout << "S_radar_ = " << S_radar_ << endl;

  //add measurement noise covariance matrix
  R_radar_ <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S_radar_ = S_radar_ + R_radar_;

  cout << "S_radar_ plus R_radar_ = " << S_radar_ << endl;

  VectorXd z_radar_ = VectorXd(3);
  z_radar_(0) = meas_package.raw_measurements_[0];
  z_radar_(1) = meas_package.raw_measurements_[1];
  z_radar_(2) = meas_package.raw_measurements_[2];

  cout << "z_radar_ = " << z_radar_ << endl;

  Tc_radar_.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;
    //angle normalization
    if (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    if (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    if (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    if (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc_radar_ += weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  K_radar_ = Tc_radar_ * S_radar_.inverse();

  //residual
  //z_diff = z_radar_ - z_pred_radar_;

  //angle normalization
  //if (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  //if (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K_radar_ * (z_radar_ - z_pred_radar_);
  P_ = P_ - K_radar_*S_radar_*K_radar_.transpose();

  if (x_(3)> M_PI) x_(3)-=2.*M_PI;
  if (x_(3)<-M_PI) x_(3)+=2.*M_PI;

  // Radat NIS calculation
  //double NIS_radar_ = z_diff.transpose()*S_radar_.inverse()*z_diff;
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}
