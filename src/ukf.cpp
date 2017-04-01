#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_package) {
  switch (measurement_package.sensor_type_) {
    case MeasurementPackage::RADAR:
      if (is_initialized_) {
        if (use_radar_) {
          UpdateRadar(measurement_package);
        }
      } else {
        InitializeRadar(measurement_package);
        is_initialized_ = true;
      }
      break;
    case MeasurementPackage::LASER:
      if (is_initialized_) {
        if (use_laser_) {
          UpdateLidar(measurement_package);
        }
      } else {
        InitializeLidar(measurement_package);
        is_initialized_ = true;
      }
      break;
    default:
      cerr << "bad measurement pack sensor type " <<
        measurement_package.sensor_type_ << endl;
      exit(EXIT_FAILURE);
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
}

void UKF::InitializeLidar(MeasurementPackage measurement_package) {
  VectorXd z = measurement_package.raw_measurements_;
  x_ << z(0), z(1), 0, 0, 0;

  // TODO: can probably choose better initial values
  double var_px = std_laspx_ * std_laspx_;
  double var_py = std_laspy_ * std_laspy_;
  double var_v = std_a_ * std_a_;
  double var_yaw = std_yawdd_ * std_yawdd_;
  P_ <<
    var_px,      0,     0,       0,       0,
         0, var_py,     0,       0,       0,
         0,      0, var_v,       0,       0,
         0,      0,     0, var_yaw,       0,
         0,      0,     0,       0, var_yaw;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

void UKF::InitializeRadar(MeasurementPackage measurement_package) {
  VectorXd z = measurement_package.raw_measurements_;
  double rho = z(0);
  double phi = z(1);
  double px = rho * cos(phi);
  double py = rho * sin(phi);

  x_ << px, py, 0, 0, 0;

  // TODO: can definitely choose better initial values
  double var_px = std_laspx_ * std_laspx_;
  double var_py = std_laspy_ * std_laspy_;
  double var_v = std_a_ * std_a_;
  double var_yaw = std_yawdd_ * std_yawdd_;
  P_ <<
    var_px,      0,     0,       0,       0,
         0, var_py,     0,       0,       0,
         0,      0, var_v,       0,       0,
         0,      0,     0, var_yaw,       0,
         0,      0,     0,       0, var_yaw;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
