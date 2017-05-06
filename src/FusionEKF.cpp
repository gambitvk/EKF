#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  G = MatrixXd(4,2);
  Qv = MatrixXd(2,2);
  
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ << 1,0,0,0,
              0,1,0,0;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  noise_ax = 9;
  noise_ay = 9;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) 
    {
        /**
TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

         //state covariance matrix P
         ekf_.P_ = MatrixXd(4, 4);
         ekf_.P_ << 1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1000, 0,
                    0, 0, 0, 1000;
        
         //measurement covariance
         ekf_.R_ = MatrixXd(2, 2);
         ekf_.R_ << 0.0225, 0,
                    0, 0.0225;
        
         //measurement matrix
         ekf_.H_ = MatrixXd(2, 4);
         ekf_.H_ << 1, 0, 0, 0,
                    0, 1, 0, 0;
        
         //the initial transition matrix F_
         ekf_.F_ = MatrixXd(4, 4);
         ekf_.F_ << 1, 0, 1, 0,
                    0, 1, 0, 1,
                    0, 0, 1, 0,
                    0, 0, 0, 1;
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
        {
            /**
              Convert radar from polar to cartesian coordinates and initialize state.
              */
            ekf_.setX(
                    measurement_pack.raw_measurements_(0),
                    measurement_pack.raw_measurements_(1),
                    measurement_pack.raw_measurements_(2)
                    );
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
        {
            /**
              Initialize state.
              */
            ekf_.setX(measurement_pack.raw_measurements_(0),
                      measurement_pack.raw_measurements_(1));
        }

        previous_timestamp_ = measurement_pack.timestamp_;
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
TODO:
     * Update the state transition matrix F according to the new elapsed time.
     - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */

    //compute the time elapsed between the current and previous measurements
    float dt =  processTime(measurement_pack.timestamp_);
    //1. Modify the F matrix so that the time is integrated
    setF(dt);
    //2. Set the process covariance matrix Q
    setQ(dt);

    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
    {
        // Radar updates
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
    else 
    {
        // Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}

float FusionEKF::processTime(const long long &t)
{
    long long pt =  previous_timestamp_;
    previous_timestamp_ = t;
    return float((t - pt) / 1000000.0); //dt - expressed in seconds
}

void FusionEKF::setF(float &dt)
{
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;
}

void FusionEKF::setQ(float &dt)
{
    float dt2 =  float(dt*dt) /2.0;
    G << dt2,0,
         0,dt2,
         dt,0,
         0,dt;
    Qv << noise_ax ,0,
          0, noise_ay;
    ekf_.Q_ = G * Qv * G.transpose();
}
