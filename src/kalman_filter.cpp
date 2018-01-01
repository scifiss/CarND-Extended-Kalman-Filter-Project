#include "kalman_filter.h"

using namespace std;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	updateState(y);

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd z_pred = ComputeRadarFromState();
  VectorXd y = z - z_pred;  // rho, phi, drho
  float pi =  3.14159265359;
  float pi2 = 6.28318530718;
  while (y(1) > pi )
  {
      cout << " phi > pi "  << endl;
      y(1) = y(1) - pi2;
  }
  while (y(2)< -pi)
  {
      cout << " phi < -pi" << endl;
      y(1) = y(1) + pi2;
  }

  updateState(y);

}

VectorXd KalmanFilter::ComputeRadarFromState()
{
    float px = x_[0];
    float py = x_[1];
    float vx = x_[2];
    float vy = x_[3];

    VectorXd z_radar = VectorXd(3);
    float rho = sqrt(px*px + py*py);
    float phi = atan2(py,px);
    float drho;
    if (rho < 1e-4)
    {
        phi = 0;
        drho = 0;

    }
    else
    {
        phi = atan2(py,px);
        drho = (px * vx + py * vy)/rho;

    }

    z_radar << rho, phi, drho;
    return z_radar;

}

void KalmanFilter::updateState(const VectorXd &y)
{

    MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;
    //new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}
