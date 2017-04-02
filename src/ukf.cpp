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
  lambda_ = 6 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(1 / (2 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  time_us_ = 0;

  NIS_radar_ = 0;
  NIS_laser_ = 0;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 6;

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
      cout << "RADAR" << endl;
      if (is_initialized_) {
        if (use_radar_) {
          double delta_t = (measurement_package.timestamp_ - time_us_) / 1e6;
          Prediction(delta_t);
          UpdateRadar(measurement_package);
        }
      } else {
        InitializeRadar(measurement_package);
        is_initialized_ = true;
      }
      break;
    case MeasurementPackage::LASER:
      cout << "LIDAR" << endl;
      if (is_initialized_) {
        if (use_laser_) {
          double delta_t = (measurement_package.timestamp_ - time_us_) / 1e6;
          Prediction(delta_t);
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
  time_us_ = measurement_package.timestamp_;
  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_ << endl;
}

// Based on
// http://commons.apache.org/proper/commons-math/javadocs/api-3.1/org/apache/
//   commons/math3/util/MathUtils.html#normalizeAngle(double,%20double)
double NormalizeAngle(double theta) {
  const double TAU = 2 * M_PI;
  double x = theta;
  double r = theta - TAU * floor((theta + M_PI) / TAU);
  if (fabs(x - r) > 1e-3) {
    cout << "NORM " << x/M_PI << " -> " << r/M_PI << endl;
  }
  return r;
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

  Eigen::LLT<MatrixXd> lltA = P_aug.llt();
  if (lltA.info() != Eigen::ComputationInfo::Success) {
    cerr << "Cholesky failed with ComputationInfo=" << lltA.info() << endl;
    exit(-1);
  }
  MatrixXd A_aug = lltA.matrixL(); // matrix square root
  double d = sqrt(lambda_ + n_aug_);
  cout << "PREDICT A_aug" << endl << A_aug << endl;

  MatrixXd Xsig_aug(n_aug_, n_sig_aug);
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) =
    (d * A_aug).colwise() + x_aug;
  Xsig_aug.block(0, n_aug_ + 1, n_aug_, n_aug_) =
    (-d * A_aug).colwise() + x_aug;

  // cout << "PREDICT delta_t = " << delta_t << endl;
  // cout << "PREDICT Xsig_aug = " << endl << Xsig_aug << endl;

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

    // if (v < 0) {
    //   v *= -1;
    //   yaw = NormalizeAngle(yaw + M_PI);
    // }

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
    Xsig_pred_(3, i) = NormalizeAngle(yaw + dyaw * delta_t + hdt2 * nu_dyaw);
    Xsig_pred_(4, i) = dyaw + delta_t * nu_dyaw;
  }

  // cout << "PREDICT Xsig_pred = " << endl << Xsig_pred_ << endl;

  //
  // Update mean and covariance from predicted sigma points.
  //
  x_ = Xsig_pred_ * weights_;

  P_.fill(0.0);
  for (int i = 0; i < n_sig_aug; ++i) {
    VectorXd d = Xsig_pred_.col(i) - x_;
    d(3) = NormalizeAngle(d(3));
    // cout << "PRED i=" << i << " w=" << weights_(i) << " d=" << d.transpose() << endl << d * d.transpose() << endl;
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
  VectorXd z = measurement_package.raw_measurements_;

  MatrixXd H(2, 5);
  H <<
    1, 0, 0, 0, 0,
    0, 1, 0, 0, 0;

  MatrixXd R(2, 2);
  R <<
    std_laspx_ * std_laspx_, 0,
    0, std_laspy_ * std_laspy_;

  VectorXd z_pred = H * x_;
  VectorXd y = z - z_pred;
  // cout << "LIDAR y=" << y.transpose() << endl;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  // cout << "LIDAR K=" << endl << K << endl;
  // cout << "LIDAR KH=" << endl << K * H << endl;
  // cout << "LIDAR pre-P=" << endl << P_ << endl;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;

  NIS_laser_ = y.transpose() * Si * y;
}

void UKF::InitializeRadar(MeasurementPackage measurement_package) {
  VectorXd z = measurement_package.raw_measurements_;
  double rho = z(0);
  double phi = z(1);
  double px = rho * cos(phi);
  double py = rho * sin(phi);

  x_ << px, py, 0, 0, 0;

  // TODO: can definitely choose better initial values
  double var_px = 1; // std_laspx_ * std_laspx_;
  double var_py = 1; // std_laspy_ * std_laspy_;
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
  VectorXd z = measurement_package.raw_measurements_;
  int n_sig_aug = 2 * n_aug_ + 1;
  int n_z = 3;
  MatrixXd Zsig(n_z, n_sig_aug);

  //
  // Find measurement points for predicted sigma points
  //
  for (int i = 0; i < n_sig_aug; ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double speed = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double dyaw = Xsig_pred_(4, i);

    // transform sigma points into measurement space
    double range = sqrt(px * px + py * py);
    double angle = atan2(py, px);
    double drange = px * cos(yaw) * speed + py * sin(yaw) * speed;
    if (fabs(range) > 1e-3) {
      drange /= range;
    } else {
      drange = 0;
    }

    Zsig(0, i) = range;
    Zsig(1, i) = angle;
    Zsig(2, i) = drange;
  }

  // calculate mean predicted measurement
  VectorXd z_pred = Zsig * weights_;

  // calculate measurement covariance matrix S
  MatrixXd S(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_aug; ++i) {
    VectorXd diff = Zsig.col(i) - z_pred;
    diff(1) = NormalizeAngle(diff(1));
    S += weights_(i) * diff * diff.transpose();
  }
  S(0, 0) += std_radr_ * std_radr_;
  S(1, 1) += std_radphi_ * std_radphi_;
  S(2, 2) += std_radrd_ * std_radrd_;

  // calculate cross correlation matrix
  MatrixXd Tc(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_aug; ++i) {
    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    xdiff(3) = NormalizeAngle(xdiff(3));

    VectorXd zdiff = Zsig.col(i) - z_pred;
    zdiff(1) = NormalizeAngle(zdiff(1));
    Tc += weights_(i) * xdiff * zdiff.transpose();
  }

  // cout << "S = " << endl << S << endl;
  // cout << "Sinv = " << endl << S.inverse() << endl;

  //calculate Kalman gain K
  MatrixXd K = Tc * S.inverse();

  cout << "K = " << endl << K << endl;

  VectorXd dz = z - z_pred;
  dz(1) = NormalizeAngle(dz(1));

  //update state mean and covariance matrix
  x_ = x_ + K * dz;
  x_(3) = NormalizeAngle(x_(3));
  P_ = P_ - K * S * K.transpose();

  NIS_radar_ = dz.transpose() * S.inverse() * dz;
}
