#include "kalman_filter.h"
#include "iostream"

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
    float c1 = px*px + py*py;
    float c2 =  (px*vx + py*vy);
    Eigen::VectorXd  result(3);
    result(0) = sqrt(c1);
    if(fabs(result(0) < 0.0001))
    {
        result(0) = 0.0001;
    }
        
    result(1) = atan2(py,px) ;
    result(2) =  (c2 / result(0));

    return result;
}

float KalmanFilter::normalize( float y)
{
    return atan2(sin(y),cos(y));
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
    y(1) = normalize(y(1));
    updateCommon(y);
}

void KalmanFilter::setX(const double &r,const  double &p ,const  double &rr)
{
    x_ << r*cos(p), r* sin(p) , rr * cos(p) , rr * sin(p);
}


void KalmanFilter::setX(const double &x,const  double &y)
{
    x_ << x,y,0.0,0.0;
}
