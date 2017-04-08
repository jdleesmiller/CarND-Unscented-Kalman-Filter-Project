#ifndef UKF_H
#define UKF_H

#include <vector>
#include <string>
#include <iostream>

#include "Eigen/Dense"

#include "measurement_package.h"
#include "kalman_filter.h"

class UKF {
  struct CanonicalizeStateAngle {
    void operator()(Eigen::Matrix<double, 5, 1> &vector) const;
  };

  struct CanonicalizeMeasurementAngle {
    void operator()(Eigen::Matrix<double, 3, 1> &vector) const;
  };

  typedef UnscentedKalmanFilter<5, 7, CanonicalizeStateAngle> Filter;

public:

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  const static bool DEFAULT_USE_LASER;
  const static bool DEFAULT_USE_RADAR;
  const static double DEFAULT_STD_A;
  const static double DEFAULT_STD_YAWD;
  const static double DEFAULT_LAMBDA;

  /**
   * Constructor.
   *
   * @param use_laser if false, lidar measurements will be ignored except during
   *                  init
   * @param use_radar if false, radar measurements will be ignored except during
   *                  init
   * @param std_a process noise standard deviation: longitudinal acceleration
   *              in m/s^2
   * @param std_yawdd process noise standard deviation yaw acceleration in
   *                  rad/s^2
   * @param lambda scaling parameter for UKF sigma points
   */
  UKF(
    bool use_laser = DEFAULT_USE_LASER,
    bool use_radar = DEFAULT_USE_RADAR,
    double std_a = DEFAULT_STD_A,
    double std_yawdd = DEFAULT_STD_YAWD,
    double lambda = DEFAULT_LAMBDA);

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage measurement_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  const Filter::StateVector &state() const { return filter_.state(); }

private:
  /**
   * The radar sensor: handles the nonlinear measurement function and Jacobian
   * calculation for updating the filter. If the first measurement is a radar
   * measurement, this class also handles initializing the filter using that
   * first radar measurement.
   *
   * The radar returns a three-dimensional measurement vector: range (rho),
   * angle (phi) and radial speed (rho_dot).
   */
  struct Radar : public Filter::Sensor<3, CanonicalizeMeasurementAngle>
  {
    explicit Radar(Filter &filter);

    void Initialize(const MeasurementVector &z);

    double Update(const MeasurementVector &z);

  private:
    // Measurement covariance matrix for radar
    MeasurementMatrix R_;
  };

  /**
   * The laser sensor: much simpler than the Radar sensor, because the
   * measurement function is linear. If the first measurement is a laser
   * measurement, this class also handles initializing the filter using that
   * first laser measurement.
   *
   * The laser returns a two-dimensional measurement vector: x and y position.
   */
  struct Laser : public Filter::Sensor<2>
  {
    explicit Laser(Filter &filter);

    void Initialize(const MeasurementVector &z);

    double Update(const MeasurementVector &z);

  private:
    // Measurement covariance matrix for laser
    MeasurementMatrix R_;

    // Measurement matrix for laser
    MeasurementStateMatrix H_;
  };

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  double std_a_;

  double std_yawdd_;

  // check whether the tracking toolbox was initiallized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long time_us_;

  // We'll have a single Filter object, and it will be updated by the two
  // Sensor objects, Radar and Laser, which represent our two sensors to be
  // fused.
  Filter filter_;
  Radar radar_;
  Laser laser_;
};

#endif /* UKF_H */
