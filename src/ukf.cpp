#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;

const bool UKF::DEFAULT_USE_LASER = true;
const bool UKF::DEFAULT_USE_RADAR = true;
const double UKF::DEFAULT_STD_A = 0.267;
const double UKF::DEFAULT_STD_YAWD = M_PI / 5.9;
const double UKF::DEFAULT_LAMBDA = 1.79 - 7.0;

// Based on
// http://commons.apache.org/proper/commons-math/javadocs/api-3.1/org/apache/
//   commons/math3/util/MathUtils.html#normalizeAngle(double,%20double)
double NormalizeAngle(double theta) {
  const double TAU = 2 * M_PI;
  double x = theta;
  double r = theta - TAU * floor((theta + M_PI) / TAU);
  // if (fabs(x - r) > 1e-3) {
  //   cout << "NORM " << x/M_PI << " -> " << r/M_PI << endl;
  // }
  return r;
}

void UKF::CanonicalizeStateAngle::operator()(
  Eigen::Matrix<double, 5, 1> &vector) const
{
  vector(3) = NormalizeAngle(vector(3)); // 3: phi
}

void UKF::CanonicalizeMeasurementAngle::operator()(
  Eigen::Matrix<double, 3, 1> &vector) const
{
  vector(1) = NormalizeAngle(vector(1)); // 1: phi
}

//
// Radar Sensor
//

UKF::Radar::Radar(UKF::Filter &filter) :
  Filter::Sensor<3, UKF::CanonicalizeMeasurementAngle>(filter)
{
  // Radar measurement noise standard deviation radius in m
  const double std_radr = 0.3;

  // Radar measurement noise standard deviation angle in rad
  const double std_radphi = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  const double std_radrd = 0.3;

  // Measurement covariance matrix:
  R_ <<
    std_radr * std_radr, 0, 0,
    0, std_radphi * std_radphi, 0,
    0, 0, std_radrd * std_radrd;
}

void UKF::Radar::Initialize(const MeasurementVector &z,
  double var_v, double var_yaw)
{
  double rho = z(0);
  double phi = z(1);
  double px = rho * cos(phi);
  double py = rho * sin(phi);

  Filter::StateVector x;
  x << px, py, 0, 0, 0;

  // TODO: can definitely choose better initial values
  double var_px = 1; // std_laspx_ * std_laspx_;
  double var_py = 1; // std_laspy_ * std_laspy_;
  Filter::StateMatrix P;
  P <<
    var_px,      0,     0,       0,       0,
         0, var_py,     0,       0,       0,
         0,      0, var_v,       0,       0,
         0,      0,     0, var_yaw,       0,
         0,      0,     0,       0, var_yaw;

  filter_.Initialize(x, P);
}

double UKF::Radar::Update(const MeasurementVector &z) {
  //
  // Find measurement points for predicted sigma points
  //
  MeasurementSigmaMatrix Z_sigma;
  for (size_t i = 0; i < Z_sigma.cols(); ++i) {
    double px = filter_.sigma_matrix()(0, i);
    double py = filter_.sigma_matrix()(1, i);
    double speed = filter_.sigma_matrix()(2, i);
    double yaw = filter_.sigma_matrix()(3, i);
    double dyaw = filter_.sigma_matrix()(4, i);

    // transform sigma points into measurement space
    double range = sqrt(px * px + py * py);
    double angle = atan2(py, px);
    double drange = px * cos(yaw) * speed + py * sin(yaw) * speed;
    if (fabs(range) > 1e-3) {
      drange /= range;
    } else {
      drange = 0;
    }

    Z_sigma(0, i) = range;
    Z_sigma(1, i) = angle;
    Z_sigma(2, i) = drange;
  }

  return UpdateUKF(z, Z_sigma, R_);
}

//
// Laser Sensor
//

UKF::Laser::Laser(UKF::Filter &filter) : Filter::Sensor<2>(filter) {
  // Laser measurement noise standard deviation position1 in m
  double std_laspx = 0.15;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy = 0.15;

  // Measurement covariance matrix:
  R_ <<
    std_laspx * std_laspx, 0,
    0, std_laspy * std_laspy;

  // Measurement matrix:
  H_ <<
    1, 0, 0, 0, 0,
    0, 1, 0, 0, 0;
}

void UKF::Laser::Initialize(const MeasurementVector &z, double var_v,
  double var_yaw)
{
  Filter::StateVector x;
  x << z(0), z(1), 0, 0, 0;

  Filter::StateMatrix P;
  P <<
    R_(0),      0,     0,       0,       0,
         0, R_(3),     0,       0,       0,
         0,      0, var_v,       0,       0,
         0,      0,     0, var_yaw,       0,
         0,      0,     0,       0, var_yaw;

  filter_.Initialize(x, P);
}

double UKF::Laser::Update(const MeasurementVector &z) {
  return Filter::Sensor<2>::Update(z, H_, R_);
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF(bool use_laser, bool use_radar, double std_a, double std_yawdd,
    double lambda) :
  NIS_radar_(0),
  NIS_laser_(0),
  use_laser_(use_laser),
  use_radar_(use_radar),
  std_a_(std_a),
  std_yawdd_(std_yawdd),
  is_initialized_(false),
  time_us_(0),
  filter_(lambda),
  radar_(filter_),
  laser_(filter_) {}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_package) {
  double delta_t = (measurement_package.timestamp_ - time_us_) / 1e6;
  switch (measurement_package.sensor_type_) {
    case MeasurementPackage::RADAR:
      // cout << "RADAR" << endl;
      if (is_initialized_) {
        if (use_radar_) {
          Prediction(delta_t);
          NIS_radar_ = radar_.Update(measurement_package.raw_measurements_);
        }
      } else {
        radar_.Initialize(measurement_package.raw_measurements_,
          std_a_ * std_a_, std_yawdd_ * std_yawdd_);
        is_initialized_ = true;
      }
      break;
    case MeasurementPackage::LASER:
      // cout << "LIDAR" << endl;
      if (is_initialized_) {
        if (use_laser_) {
          Prediction(delta_t);
          NIS_laser_ = laser_.Update(measurement_package.raw_measurements_);
        }
      } else {
        laser_.Initialize(measurement_package.raw_measurements_,
          std_a_ * std_a_, std_yawdd_ * std_yawdd_);
        is_initialized_ = true;
      }
      break;
    default:
      cerr << "bad measurement pack sensor type " <<
        measurement_package.sensor_type_ << endl;
      exit(EXIT_FAILURE);
  }
  time_us_ = measurement_package.timestamp_;
  // cout << "x_ = " << endl << filter_.state() << endl;
  // cout << "P_ = " << endl << filter_.covariance() << endl;
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
  Filter::AugmentedNoiseMatrix Q;
  Q <<
    std_a_ * std_a_, 0,
    0, std_yawdd_ * std_yawdd_;

  Filter::AugmentedStateSigmaMatrix Xsig_aug =
    filter_.GenerateSigmaPoints(Q);

  //
  // Predict (augmented) sigma points.
  //
  double hdt2 = 0.5 * delta_t * delta_t;
  Filter::StateSigmaMatrix Xsig_pred_; // TODO rename
  for (size_t i = 0; i < Xsig_aug.cols(); ++i) {
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
    Xsig_pred_(3, i) = NormalizeAngle(yaw + dyaw * delta_t + hdt2 * nu_dyaw);
    Xsig_pred_(4, i) = dyaw + delta_t * nu_dyaw;
  }

  filter_.PredictUKF(Xsig_pred_);
}
