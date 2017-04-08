#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include "Eigen/Dense"
#include "Eigen/LU"

/**
 * Functor that does nothing to a state, in case the state does not require
 * any canonicalization to deal with e.g. angles.
 */
template <int VectorSize>
struct CanonicalizeIdentity {
  typedef Eigen::Matrix<double, VectorSize, 1> Vector;
  void operator()(Vector&) const { }
};

/**
 * Generic Kalman Filter class.
 *
 * The filter maintains the system state and the variance as its state. The
 * updates to the state are managed by Sensor objects. Multiple Sensor objects
 * can reference the same KalmanFilter.
 *
 * The template class is parameterised by the dimension of the state vector.
 * Each Sensor template class is also parameterised by the dimension of the
 * measurement vector. This makes it possible to use fixed-size matrix math
 * for improved compile-time correctness checking and (sometimes) better
 * performance.
 */
template <int StateSize, int AugmentedStateSize,
  typename CanonicalizeState = CanonicalizeIdentity<StateSize>>
struct UnscentedKalmanFilter {
  static const int SigmaSize = 2 * AugmentedStateSize + 1;

  typedef Eigen::Matrix<double, StateSize, 1> StateVector;
  typedef Eigen::Matrix<double, StateSize, StateSize> StateMatrix;

  typedef Eigen::Matrix<double, SigmaSize, 1> SigmaVector;
  typedef Eigen::Matrix<double, SigmaSize, SigmaSize> SigmaMatrix;
  typedef Eigen::Matrix<double, StateSize, SigmaSize> StateSigmaMatrix;

  typedef Eigen::Matrix<double, AugmentedStateSize, 1> AugmentedStateVector;
  typedef Eigen::Matrix<double, AugmentedStateSize, AugmentedStateSize>
    AugmentedStateMatrix;
  typedef Eigen::Matrix<double, AugmentedStateSize, SigmaSize>
    AugmentedStateSigmaMatrix;
  typedef Eigen::Matrix<double,
    AugmentedStateSize - StateSize,
    AugmentedStateSize - StateSize> AugmentationMatrix;

  /**
   * A generic sensor for a Kalman Filter. The Sensor is responsible for
   * updating the state of the filter based on its measurements.
   */
  template <int MeasurementSize,
    typename CanonicalizeMeasurement = CanonicalizeIdentity<MeasurementSize> >
  struct Sensor {
    typedef UnscentedKalmanFilter<StateSize, AugmentedStateSize,
      CanonicalizeState> Filter;
    typedef Eigen::Matrix<double, MeasurementSize, 1> MeasurementVector;
    typedef Eigen::Matrix<double, MeasurementSize, MeasurementSize>
      MeasurementMatrix;
    typedef Eigen::Matrix<double, MeasurementSize, SigmaSize>
      MeasurementSigmaMatrix;
    typedef Eigen::Matrix<double, MeasurementSize, StateSize>
      MeasurementStateMatrix;
    typedef Eigen::Matrix<double, StateSize, MeasurementSize>
      StateMeasurementMatrix;

    /**
     * Constructor: provide a reference to the KalmanFilter that this Sensor
     * updates.
     */
    Sensor(Filter &filter) : filter_(filter) { }

    /**
     * Update a linear Kalman Filter with a new measurement.
     *
     * @param z measurement
     * @param H measurement matrix
     * @param R measurement covariance matrix
     */
    double Update(
      const MeasurementVector &z,
      const MeasurementStateMatrix &H,
      const MeasurementMatrix &R) {
      return UpdateEKF(z, H * filter_.state(), H, R);
    }

    /**
     * Update an Extended Kalman Filter with a new measurement.
     *
     * @param z measurement
     * @param h measurement predicted (possibly non-linearly) from the state
     * @param H measurement matrix (possibly linearized)
     * @param R measurement covariance matrix
     */
    double UpdateEKF(
      const MeasurementVector &z,
      const MeasurementVector &h,
      const MeasurementStateMatrix &H,
      const MeasurementMatrix &R)
    {
      const StateMatrix &P = filter_.covariance();
      MeasurementVector y = z - h;
      MeasurementMatrix S = H * P * H.transpose() + R;

      // Avoid explicit matrix inversion. Starting from
      // K = P H^T S^-1
      // postmultiply by S to get
      // K S = P H^T
      // which is equivalent to
      // S^T K^T = H P^T
      // which is in standard form (Ax=B)
      Eigen::FullPivLU<MeasurementMatrix> lu(S.transpose());
      StateMeasurementMatrix K = lu.solve(H * P.transpose()).transpose();

      filter_.UpdateOptimal(K * y, K * H);

      // Compute the normalized innovation score y^T S^{-1} y. We already have
      // the LU factorization of S^T, so we would like to reuse it. To do this,
      // if we start with
      // N = y^T S^{-1} y
      // then
      // N = (S^{-1}^T y)^T y
      // so
      // N = (S^T^{-1} y)^T y
      // since the transpose of an inverse equals the inverse of a transpose.
      // Let c = S^T^{-1} y and premultiply by S^T to obtain
      // S^T c = y
      // which is in standard form. If we solve for c, we can then use
      // N = c^T y
      // to obtain the NIS score without the explicit matrix inverse.
      MeasurementVector c = lu.solve(y);
      return c.transpose() * y;
    }

    /**
     * TODO
     * @param  z       [description]
     * @param  Z_sigma [description]
     * @param  R       [description]
     * @return the normalized innovations squared (NIS) score
     */
    double UpdateUKF(
      const MeasurementVector &z,
      const MeasurementSigmaMatrix &Z_sigma,
      const MeasurementMatrix &R)
    {
      // calculate mean predicted measurement
      MeasurementVector z_pred = Z_sigma * filter_.sigma_weights();

      // calculate measurement covariance matrix S
      // calculate cross correlation matrix Tc
      MeasurementMatrix S(R);
      StateMeasurementMatrix Tc;
      Tc.fill(0.0);

      for (size_t i = 0; i < SigmaSize; ++i) {
        MeasurementVector zdiff = Z_sigma.col(i) - z_pred;
        canonicalize_measurement_(zdiff);

        StateVector xdiff = filter_.sigma_matrix().col(i) - filter_.state();
        canonicalize_state_(xdiff);

        S += filter_.sigma_weights()(i) * zdiff * zdiff.transpose();
        Tc += filter_.sigma_weights()(i) * xdiff * zdiff.transpose();
      }

      // Calculate Kalman gain K = T_c S^{-1}. To avoid the explicit inverse,
      // we can rewrite this as K S = T_c, and then
      // S^T K^T = T_c^T
      // which is in standard form.
      Eigen::FullPivLU<MeasurementMatrix> lu(S.transpose());
      StateMeasurementMatrix K = lu.solve(Tc.transpose()).transpose();

      MeasurementVector y = z - z_pred;
      canonicalize_measurement_(y);

      // update state mean and covariance matrix
      filter_.Update(K * y, K * S * K.transpose());

      // Compute the normalized innovation score y^T S^{-1} y. See notes in
      // UpdateEKF about how to reuse the LU factorization of S^T.
      MeasurementVector c = lu.solve(y);
      return c.transpose() * y;
    }

  protected:
    CanonicalizeState canonicalize_state_;
    CanonicalizeMeasurement canonicalize_measurement_;
    Filter &filter_;
  };

  /**
   * Constructor: the state and covariance matrices are uninitialized (may
   * contain values from uninitialized memory).
   */
  UnscentedKalmanFilter(double lambda, const AugmentationMatrix &Q) :
    UnscentedKalmanFilter(StateVector(), StateMatrix(), lambda, Q) {}

  /**
   * Constructor.
   *
   * @param x initial state
   * @param P initial covariance
   * @param lambda for calculation of weights
   */
  UnscentedKalmanFilter(
    const StateVector &x, const StateMatrix &P,
    double lambda, const AugmentationMatrix &Q)
    : canonicalize_(), I_(StateMatrix::Identity()),
      x_(x), P_(P), lambda_(lambda), Q_(Q),
      sigma_weights_(MakeSigmaWeights(lambda)),
      sigma_matrix_(StateSigmaMatrix::Zero()) // initialize for safety
  { }

  const StateVector &state() const { return x_; }

  const StateMatrix &covariance() const { return P_; }

  const double &lambda() const { return lambda_; }

  const SigmaVector &sigma_weights() const { return sigma_weights_; }

  const StateSigmaMatrix &sigma_matrix() const { return sigma_matrix_; }

  /**
   * TODO
   *
   * @return   [description]
   */
  AugmentedStateSigmaMatrix GenerateSigmaPoints() {
    AugmentedStateVector x_aug;
    x_aug.fill(0.0);
    x_aug.head(StateSize) = x_;

    AugmentedStateMatrix P_aug;
    P_aug.fill(0.0);
    P_aug.topLeftCorner(StateSize, StateSize) = P_;
    P_aug.bottomRightCorner(Q_.rows(), Q_.cols()) = Q_;

    Eigen::LLT<AugmentedStateMatrix> lltA = P_aug.llt();
    if (lltA.info() != Eigen::ComputationInfo::Success) {
      std::cerr << "Cholesky failed with ComputationInfo=" << lltA.info() <<
        std::endl;
      exit(-1);
    }
    AugmentedStateMatrix A_aug = lltA.matrixL(); // matrix square root
    double d = sqrt(lambda_ + AugmentedStateSize);

    AugmentedStateSigmaMatrix X_sigma_aug;
    X_sigma_aug.col(0) = x_aug;
    X_sigma_aug.block(0, 1,
      AugmentedStateSize, AugmentedStateSize) = (d * A_aug).colwise() + x_aug;
    X_sigma_aug.block(0, AugmentedStateSize + 1,
      AugmentedStateSize, AugmentedStateSize) = (-d * A_aug).colwise() + x_aug;

    return X_sigma_aug;
  }

  /**
   * Predict TODO
   *
   * @param X_sigma the predicted sigma matrix
   */
  void PredictUKF(const StateSigmaMatrix &X_sigma)
  {
    // Save the predicted sigma matrix for use in sensor updates.
    sigma_matrix_ = X_sigma;

    x_ = X_sigma * sigma_weights_;
    canonicalize_(x_);

    P_.fill(0.0);
    for (size_t i = 0; i < SigmaSize; ++i) {
      StateVector diff = X_sigma.col(i) - x_;
      canonicalize_(diff);
      P_ += sigma_weights_(i) * diff * diff.transpose();
    }
  }

  /**
   * Initialize (or reinitialize) the state of the filter.
   *
   * @param x state
   * @param P covariance matrix
   */
  void Initialize(const StateVector &x, const StateMatrix &P) {
    x_ = x;
    P_ = P;
  }

  /**
   * Update the state and covariance estimates using the update formula for
   * the optimal gain. The caller is responsible for computing the Kalman gain
   * and using it to scale the residual and the measurement matrix. This allows
   * this method to be independent of the dimension of the measurement vector.
   *
   * @param Ky residual, premultiplied by the Kalman gain, K
   * @param KH measurement matrix, premultiplied by the Kalman gain, K
   */
  void UpdateOptimal(const StateVector &Ky, const StateMatrix &KH) {
    x_ = x_ + Ky;
    canonicalize_(x_);
    P_ = (I_ - KH) * P_;
  }

  /**
   * Update the state and covariance estimates using the given scaled state
   * and covariance updates.
   *
   * @param Ky the residual, y, premultiplied by the Kalman gain, K
   * @param KP the covariance matrix update, after scaling by the Kalman gain
   */
  void Update(const StateVector &Ky, const StateMatrix &KP) {
    x_ = x_ + Ky;
    canonicalize_(x_);
    P_ = P_ - KP;
  }

private:
  // Create the weights for the sigma points. The weights do not change, so we
  // do it once when the filter is initialized.
  static SigmaVector MakeSigmaWeights(double lambda) {
    SigmaVector weights;
    weights.fill(1 / (2 * (lambda + AugmentedStateSize)));
    weights(0) = lambda / (lambda + AugmentedStateSize);
    return weights;
  }

  // Functor to canonicalize the state.
  const CanonicalizeState canonicalize_;

  // Identity matrix (for UpdateOptimal)
  const StateMatrix I_;

  // State vector
  StateVector x_;

  // State covariance matrix
  StateMatrix P_;

  // State covariance augmentation matrix
  const AugmentationMatrix Q_;

  // Sigma point spreading parameter
  const double lambda_;

  // Weight vector for sigma point matrix
  const SigmaVector sigma_weights_;

  // Predicted state sigma point matrix
  StateSigmaMatrix sigma_matrix_;
};

#endif /* KALMAN_FILTER_H_ */
