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

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(1 / (2 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

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
  //
  // Generate (augmented) sigma points.
  //
  int n_sig_aug = 2 * n_aug_ + 1;

  VectorXd x_aug(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  MatrixXd P_aug(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  MatrixXd A_aug = P_aug.llt().matrixL(); // matrix square root
  double d = sqrt(lambda_ + n_x_);

  MatrixXd Xsig_aug(n_aug_, n_sig_aug);
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) =
    (d * A_aug).colwise() + x_aug;
  Xsig_aug.block(0, n_aug_ + 1, n_aug_, n_aug_) =
    (-d * A_aug).colwise() + x_aug;

  //
  // Predict (augmented) sigma points.
  //
  double hdt2 = 0.5 * delta_t * delta_t;
  for (int i = 0; i < n_sig_aug; ++i) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double dyaw = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_dyaw = Xsig_aug(6, i);

    if (fabs(dyaw) > 1e-3) {
      double r = v / dyaw;
      Xsig_pred_(0, i) =
        px + r * (sin(yaw + dyaw * delta_t) - sin(yaw)) +
        hdt2 * cos(yaw) * nu_a;
      Xsig_pred_(1, i) =
        py + r * (-cos(yaw + dyaw * delta_t) + cos(yaw)) +
        hdt2 * sin(yaw) * nu_a;
    } else {
      Xsig_pred_(0, i) =
        px + v * cos(yaw) * delta_t +
        hdt2 * cos(yaw) * nu_a;
      Xsig_pred_(1, i) =
        py + v * sin(yaw) * delta_t +
        hdt2 * sin(yaw) * nu_a;
    }
    Xsig_pred_(2, i) = v + delta_t * nu_a;
    Xsig_pred_(3, i) = yaw + dyaw * delta_t + hdt2 * nu_dyaw;
    Xsig_pred_(4, i) = dyaw + delta_t * nu_dyaw;
  }

  //
  // Update mean and covariance from predicted sigma points.
  //
  x_ = Xsig_pred_ * weights_;

  P_.fill(0.0);
  for (int i = 0; i < n_sig_aug; ++i) {
    VectorXd d = Xsig_pred_.col(i) - x_;
    P_ += weights_(i) * d * d.transpose();
  }
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
