#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    if(estimations.size() == 0 ||
       estimations.size() != ground_truth.size())
    {
        std::cout << "Invalid estimation or ground_truth data" << std::endl;
        return rmse;
    }

    for(unsigned int i=0; i < estimations.size() ; ++i)
    {
        VectorXd res = estimations[i] - ground_truth[i];

        res = res.array() * res.array();
        rmse += res;
    }

    rmse = rmse/estimations.size();
    return rmse.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //check division by zero

    if(px*py*vx*vy == 0)
    {
        std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl; 
        return Hj;
    }
    //compute the Jacobian matrix
    Hj(2,2) = Hj(0,0) = px / sqrt(cal1(px,py));
    Hj(2,3) = Hj(0,1) = py / sqrt(cal1(px,py));

    Hj(1,0) = -py / cal1(px,py);
    Hj(1,1) = px / cal1(px,py);

    Hj(2,0) = cal2(py,px,vx,vy) / pow(cal1(px,py),3/2);
    Hj(2,1) = cal2(px,py,vy,vx) / pow(cal1(px,py),3/2);

    return Hj;
}
