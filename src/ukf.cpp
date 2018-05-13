#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <cfloat>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//set state dimension
const int n_x = 5;

//set augmented dimension
const int n_aug = 7;

//define spreading parameter
const double lambda = 3 - n_aug;

//weight coefficient
const double inv_lambda_na = 1 / (lambda + n_aug);

//the number of sigma point
const int n_sigma_point = 2 * n_aug +1;

//set measurement dimension, radar can measure r, phi, and r_dot
const int n_z = 3;

//set measurement dimension, lidar can measure px and py
const int n_zl = 2;

// previous timestamp
long long previous_timestamp_;

//M_PI*2
const double PI2 = M_PI * 2;

//create augmented mean vector
VectorXd x_aug = VectorXd(7);

//create augmented state covariance
MatrixXd P_aug = MatrixXd(7, 7);

//create sigma point matrix
MatrixXd Xsig_aug = MatrixXd(n_aug, n_sigma_point);

//create matrix with predicted sigma points as columns
MatrixXd Xsig_pred = MatrixXd(n_x, n_sigma_point);

//create matrix for x difference
MatrixXd x_diff = MatrixXd(n_x, n_sigma_point);

//create matrix for x difference with weight
MatrixXd x_diff_weights = MatrixXd(n_x, n_sigma_point);

//create vector for weights
VectorXd weights = VectorXd(n_sigma_point);

//create vector for predicted state
VectorXd x = VectorXd(n_x);

//create covariance matrix for prediction
MatrixXd P = MatrixXd(n_x, n_x);

//create weight matrix
MatrixXd weights_mat = MatrixXd(n_x, n_sigma_point);

//create weight matrix for radar
MatrixXd weights_matr = MatrixXd(n_z, n_sigma_point);

//create weight matrix for lidar
MatrixXd weights_matl = MatrixXd(n_zl, n_sigma_point);

//create matrix for sigma points in measurement space for radar
MatrixXd Zsig = MatrixXd(n_z, n_sigma_point);

//create matrix for sigma points in measurement space for lidar
MatrixXd Zsigl = MatrixXd(n_zl, n_sigma_point);

//mean predicted measurement for radar
VectorXd z_pred = VectorXd(n_z);

//mean predicted measurement for lidar
VectorXd z_predl = VectorXd(n_zl);
  
//measurement covariance matrix S for radar
MatrixXd S = MatrixXd(n_z, n_z);

//measurement covariance matrix S for lidar
MatrixXd Sl = MatrixXd(n_zl, n_zl);

//Radar measurement noise matrix
MatrixXd R = MatrixXd(n_z, n_z);

//Lidar measurement noise matrix
MatrixXd Rl = MatrixXd(n_zl, n_zl);

//create matrix for cross correlation Tc for radar
MatrixXd Tc = MatrixXd(n_x, n_z);

//create matrix for cross correlation Tc for lidar
MatrixXd Tcl = MatrixXd(n_x, n_zl);


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
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
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

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */


 //set weight
 weights.fill(0.5 *inv_lambda_na);
 weights(0) = lambda * inv_lambda_na;

  //weight matrix(repmat)
  for(int i=0; i<n_x; i++){
      weights_mat.row(i) = weights;
  }
  for(int i=0; i<n_z; i++){
      weights_matr.row(i) = weights;
  }
  for(int i=0; i<n_zl; i++){
      weights_matl.row(i) = weights;
  }

 R.fill(0.0);
 R(0,0) = std_radr_   * std_radr_;
 R(1,1) = std_radphi_ * std_radphi_;
 R(2,2) = std_radrd_  * std_radrd_;

 Rl.fill(0.0);
 Rl(0,0) = std_laspx_ * std_laspx_;
 Rl(1,1) = std_laspy_ * std_laspy_;

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

	/*****************************************************************************
	*  Initialization
	****************************************************************************/
	if (!is_initialized_) {

		// first measurement
		cout << "UKF: " << endl;

		P_ = MatrixXd::Identity(n_x, n_x);

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		/**
		Convert radar from polar to cartesian coordinates and initialize state.
		*/
			float px  = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
			float py  = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
			float v   = meas_package.raw_measurements_[2];
			float psi = meas_package.raw_measurements_[1];
			x_ << px, py, v, psi, 0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      		/**
      		Initialize state.
      		*/
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
		}

		previous_timestamp_ = meas_package.timestamp_;
		// done initializing, no need to predict or update
		is_initialized_ = true;
    		return;
	}


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

	//compute the time elapsed between the current and previous measurements
	float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;

	Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  // print the output

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
     if(use_radar_){
    // Radar updates
	UpdateRadar(meas_package);
     }
  } else {
     if(use_laser_){
    // Laser updates
	UpdateLidar(meas_package);
     }
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

/*******************************************************************************
 * Augmantation
 ******************************************************************************/

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  //set first column of sigma point matrix
  Xsig_aug.col(0)  = x_aug;

  //set remaining sigma points
  for (int i = 0; i < n_aug; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug) * A.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * A.col(i);
  }

/*******************************************************************************
 * Sigma Point Prediction
 ******************************************************************************/

  double delta2_0p5 = delta_t * delta_t * 0.5;

  for(int i=0; i<n_sigma_point; i++){
      VectorXd XsigCol      = Xsig_aug.col(i);
      double px     = XsigCol(0);
      double py     = XsigCol(1);
      double v      = XsigCol(2);
      double psi    = XsigCol(3);
      double psidot = XsigCol(4);
      double nua    = XsigCol(5);
      double nupdd  = XsigCol(6);

      double cos_psi        = cos(psi);
      double sin_psi        = sin(psi);
      double delta2_0p5_nua = delta2_0p5 * nua;

      double px_noise       = delta2_0p5_nua * cos_psi;
      double py_noise       = delta2_0p5_nua * sin_psi;
      double v_noise        = delta_t    * nua;
      double psi_noise      = delta2_0p5 * nupdd;
      double psidot_noise   = delta_t    * nupdd;

      if(fabs(psidot) < DBL_EPSILON){
        double px_new     = px + v * cos_psi * delta_t + px_noise;
        double py_new     = py + v * sin_psi * delta_t + py_noise;
        double v_new      = v + v_noise;
        double psi_new    = psi + psi_noise;
        double psidot_new = psidot + psidot_noise; 
        
        Xsig_pred.col(i) << px_new, py_new, v_new, psi_new, psidot_new;
      }
      else{
        double vkpk       = v / psidot;
        double psi_delta  = psidot * delta_t;

        double psi_new    = psi + psi_delta;
        double px_new     = px + vkpk * ( sin(psi_new) - sin_psi ) + px_noise;
        double py_new     = py + vkpk * (-cos(psi_new) + cos_psi ) + py_noise;
        double v_new      = v + v_noise;
        double psidot_new = psidot + psidot_noise;
        psi_new         += psi_noise;

        Xsig_pred.col(i) << px_new, py_new, v_new, psi_new, psidot_new;
      }      
  }


/*******************************************************************************
 * Predict Mean and Covariance
 ******************************************************************************/

  //Predict Mean
  x_      = Xsig_pred * weights;

  x_diff = Xsig_pred.colwise() - x_;

  //angle normalization
  for (int i = 0; i < n_sigma_point; i++) {  //iterate over sigma points
     while (x_diff(3,i) >  M_PI) x_diff(3,i) -= PI2;
     while (x_diff(3,i) < -M_PI) x_diff(3,i) += PI2;
  }

  x_diff_weights = x_diff.array() * weights_mat.array();

  //Predict Covariance 
  P_ = x_diff_weights * x_diff.transpose();

}




/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
cout << "LIDAR " << endl;
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  VectorXd z = VectorXd(n_zl);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

  for(int i=0; i<n_sigma_point; i++){
      VectorXd XsigCol      = Xsig_pred.col(i);
      double px      = XsigCol(0);
      double py      = XsigCol(1);
      Zsigl.col(i) << px, py;
  }

  z_predl  = Zsigl * weights;

  MatrixXd z_diff = Zsigl.colwise() - z_predl;

  MatrixXd z_diff_weights = z_diff.array() * weights_matl.array();

  //calculate innovation covariance matrix S
  Sl = z_diff_weights * z_diff.transpose() + Rl;

  //calculate Cross-correlation Matrix Tc
  Tcl = x_diff_weights * z_diff.transpose();

  //update x and P
  VectorXd y = z - z_predl;
  MatrixXd K = Tcl * Sl.inverse();

  //new estimate
  x_ = x_ + (K * y);
  P_ = P_ - K * Sl * K.transpose();

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
cout << "RADAR " << endl;
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */


/*******************************************************************************
 * Predict Measurement
 ******************************************************************************/

  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

  for(int i=0; i<n_sigma_point; i++){
      VectorXd XsigCol      = Xsig_pred.col(i);
      double px      = XsigCol(0);
      double py      = XsigCol(1);
      double v       = XsigCol(2);
      double psi     = XsigCol(3);

      double cos_psi = cos(psi);
      double sin_psi = sin(psi);

      double rho     = sqrt(px*px + py*py);
      double phi     = atan2(py, px);
      //double phi     = psi;//??
    
      double rhodot;

      if(rho < DBL_EPSILON){
           rhodot = 0.0;
      }
      else{
           rhodot = v * (px * cos_psi + py * sin_psi) / rho;
      }      
      Zsig.col(i) << rho, phi, rhodot;
  }

  z_pred  = Zsig * weights;

  MatrixXd z_diff = Zsig.colwise() - z_pred;

  //angle normalization  
  for (int i = 0; i < n_sigma_point; i++) {  //iterate over sigma points
     while (z_diff(1,i) >  M_PI) z_diff(1,i) -= PI2;
     while (z_diff(1,i) < -M_PI) z_diff(1,i) += PI2;
  }

  MatrixXd z_diff_weights = z_diff.array() * weights_matr.array();

  //calculate innovation covariance matrix S
  S = z_diff_weights * z_diff.transpose() + R;
 
  //calculate Cross-correlation Matrix Tc
  Tc = x_diff_weights * z_diff.transpose();

  //update x and P
  VectorXd y = z - z_pred;
  MatrixXd K = Tc * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  P_ = P_ - K * S * K.transpose();

}
