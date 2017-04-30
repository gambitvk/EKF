#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
 //output 
  x_ = x_in;
  P_ = P_in;

  //change depends on snsoir
  H_ = H_in;
  R_ = R_in;

  //static after initi
  F_ = F_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() 
{
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) 
{
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    updateCommon(y);
}

Eigen::VectorXd KalmanFilter::getPolar()
{
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);
    Eigen::VectorXd  result(3);
    result(0) = sqrt(px*px + py*py);
    result(1) = atan2(py,px);
    result(2) = (px*vx + py*vy) / result(0);
    return result;
}

void KalmanFilter::updateCommon(Eigen::VectorXd &y)
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

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    VectorXd y = z - getPolar();
    updateCommon(y);
}
