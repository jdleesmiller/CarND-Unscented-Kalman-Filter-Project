#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;

// Default parameters for the constructor.
const bool UKF::DEFAULT_USE_LASER = true;
const bool UKF::DEFAULT_USE_RADAR = true;
const double UKF::DEFAULT_STD_A = 0.484;
const double UKF::DEFAULT_STD_YAWD = M_PI / 5.25;
const double UKF::DEFAULT_LAMBDA = 2.05 - 7.0;

// Initial variances for components of the state that we don't observe.
const double INITIAL_SPEED_VARIANCE = 25;
const double INITIAL_YAW_VARIANCE = M_PI * M_PI / 16;
const double INITIAL_YAWD_VARIANCE = M_PI * M_PI / 64;

// For avoiding division by zero: treat numbers smaller than this as zero.
const double EPSILON = 1e-3;

//
// Angle Canonicalization
//

// Map an arbitrary angle to [-pi, pi). Based on
// http://commons.apache.org/proper/commons-math/javadocs/api-3.1/org/apache/
//   commons/math3/util/MathUtils.html#normalizeAngle(double,%20double)
double CanonicalizeAngle(double theta) {
  const double TAU = 2 * M_PI;
  return theta - TAU * floor((theta + M_PI) / TAU);
}

void UKF::CanonicalizeStateAngle::operator()(
  Eigen::Matrix<double, 5, 1> &vector) const
{
  vector(3) = CanonicalizeAngle(vector(3)); // 3: phi
}

void UKF::CanonicalizeMeasurementAngle::operator()(
  Eigen::Matrix<double, 3, 1> &vector) const
{
  vector(1) = CanonicalizeAngle(vector(1)); // 1: phi
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

void UKF::Radar::Initialize(const MeasurementVector &z)
{
  double rho = z(0);
  double phi = z(1);
  double var_rho = R_(0, 0);
  double var_phi = R_(1, 1);

  double px = rho * cos(phi);
  double py = rho * sin(phi);

  // In addition to the position, we can also obtain an estimate of the variance
  // due to measurement error in the radar sensor. Let's take the variance of px
  // as an example; the variance for py is computed similarly.
  //
  // We need to propagate error through a product, rho * cos(phi), and a cosine.
  // Fortunately, Wikipedia tells us exactly how to do this:
  // https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas
  // (But note that all of these rules are approximate; they hold only when
  // the measurement error is small relative to the value being measured.)
  //
  // The radar R_ matrix is diagonal, which tells us that measurements of rho
  // and phi are uncorrelated. The formula for propagating error through the
  // product AB for normal random variables A and B therefore simplifies to:
  //
  // Var(AB) = Var(A) B^2 + Var(B) A^2
  //
  // The rule for cos(C) for normal random variable C is:
  //
  // Var(cos(C)) = Var(C) sin^2(C)
  //
  // Combining the above, we obtain
  //
  // Var(rho cos(phi))
  //   = Var(rho) cos^2(phi) + Var(cos(phi)) rho^2
  //   =                 ... + Var(phi) sin^2(phi) rho^2
  //   =                 ... + Var(phi) (py)^2
  //
  // Note: it turns out is possible to obtain an exact expression for
  // Var(cos(phi)): http://nbviewer.jupyter.org/gist/dougalsutherland/8513749
  // so we could try that if this approximation causes problems.
  double cos2phi = cos(phi) * cos(phi);
  double sin2phi = sin(phi) * sin(phi);
  double var_px = var_rho * cos2phi + var_phi * py * py;
  double var_py = var_rho * sin2phi + var_phi * px * px;

  Filter::StateVector x;
  x << px, py, 0, 0, 0;

  Filter::StateMatrix P;
  double var_v = INITIAL_SPEED_VARIANCE;
  double var_yaw = INITIAL_YAW_VARIANCE;
  double var_yawd = INITIAL_YAWD_VARIANCE;
  P <<
        var_px,      0,     0,       0,        0,
             0, var_py,     0,       0,        0,
             0,      0, var_v,       0,        0,
             0,      0,     0, var_yaw,        0,
             0,      0,     0,       0, var_yawd;

  filter_.Initialize(x, P);
}

double UKF::Radar::Update(const MeasurementVector &z) {
  // Find measurement points for predicted sigma points
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
    if (fabs(range) > EPSILON) {
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

void UKF::Laser::Initialize(const MeasurementVector &z)
{
  Filter::StateVector x;
  x << z(0), z(1), 0, 0, 0;

  Filter::StateMatrix P;
  double var_v = INITIAL_SPEED_VARIANCE;
  double var_yaw = INITIAL_YAW_VARIANCE;
  double var_yawd = INITIAL_YAWD_VARIANCE;
  P <<
    R_(0, 0),        0,     0,       0,        0,
           0, R_(1, 1),     0,       0,        0,
           0,        0, var_v,       0,        0,
           0,        0,     0, var_yaw,        0,
           0,        0,     0,       0, var_yawd;

  filter_.Initialize(x, P);
}

double UKF::Laser::Update(const MeasurementVector &z) {
  // This is just a standard linear update.
  return Filter::Sensor<2>::Update(z, H_, R_);
}

//
// Filter
//

/**
 * Make the 'Q' matrix to add to the lower right block of the augmented
 * covariance (P_aug) matrix.
 */
UKF::Filter::AugmentationMatrix MakeQMatrix(double std_a, double std_yawdd)
{
  UKF::Filter::AugmentationMatrix Q;
  Q <<
    std_a * std_a,                     0,
                0, std_yawdd * std_yawdd;
  return Q;
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
  is_initialized_(false),
  time_us_(0),
  filter_(lambda, MakeQMatrix(std_a, std_yawdd)),
  radar_(filter_),
  laser_(filter_) {}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage measurement_package) {
  double delta_t = (measurement_package.timestamp_ - time_us_) / 1e6;
  switch (measurement_package.sensor_type_) {
    case MeasurementPackage::RADAR:
      // cout << "RADAR" << endl;
      if (is_initialized_) {
        if (use_radar_) {
          Predict(delta_t);
          NIS_radar_ = radar_.Update(measurement_package.raw_measurements_);
        }
      } else {
        radar_.Initialize(measurement_package.raw_measurements_);
        is_initialized_ = true;
      }
      break;
    case MeasurementPackage::LASER:
      // cout << "LIDAR" << endl;
      if (is_initialized_) {
        if (use_laser_) {
          Predict(delta_t);
          NIS_laser_ = laser_.Update(measurement_package.raw_measurements_);
        }
      } else {
        laser_.Initialize(measurement_package.raw_measurements_);
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

void UKF::Predict(double delta_t) {

  // Predict (augmented) sigma points.
  Filter::AugmentedStateSigmaMatrix Xsig_aug =
    filter_.GenerateAugmentedSigmaPoints();
  Filter::StateSigmaMatrix Xsig_pred;
  double hdt2 = 0.5 * delta_t * delta_t;
  for (size_t i = 0; i < Xsig_aug.cols(); ++i) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double dyaw = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_dyaw = Xsig_aug(6, i);

    if (fabs(dyaw) > EPSILON) {
      double r = v / dyaw;
      Xsig_pred(0, i) =
        px + r * (sin(yaw + dyaw * delta_t) - sin(yaw)) +
        hdt2 * cos(yaw) * nu_a;
      Xsig_pred(1, i) =
        py + r * (-cos(yaw + dyaw * delta_t) + cos(yaw)) +
        hdt2 * sin(yaw) * nu_a;
    } else {
      Xsig_pred(0, i) =
        px + v * cos(yaw) * delta_t +
        hdt2 * cos(yaw) * nu_a;
      Xsig_pred(1, i) =
        py + v * sin(yaw) * delta_t +
        hdt2 * sin(yaw) * nu_a;
    }
    Xsig_pred(2, i) = v + delta_t * nu_a;
    Xsig_pred(3, i) = CanonicalizeAngle(yaw + dyaw * delta_t + hdt2 * nu_dyaw);
    Xsig_pred(4, i) = dyaw + delta_t * nu_dyaw;
  }

  filter_.PredictUKF(Xsig_pred);
}
