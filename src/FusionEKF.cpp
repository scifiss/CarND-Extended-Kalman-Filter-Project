
#include <iostream>
#include "Eigen/Dense"

#include "FusionEKF.h"
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
Tools tools;


/*
 * Constructor.
 */
FusionEKF::FusionEKF()
{
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    H_radar_ = MatrixXd(3, 4);  // to be filled
//    cout << "FusionEKF" <<endl;

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
             0, 0.0225;


    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
             0, 0.0009, 0,
             0, 0, 0.09;

    /**
    TODO:
      * Finish initializing the FusionEKF.
      * Set the process and measurement noises
    */
    H_laser_ << 1,0,0,0,
             0,1,0,0;
//              cout << "H_laser_" <<endl;
    // state transition matrix
    ekf_.F_ = MatrixXd(4, 4);  // to be filled
//     cout << "ekf_.F_ " <<endl;
    // State covariance matrix P
    ekf_.P_ = MatrixXd(4, 4);
//     cout << "ekf_.P_ " <<endl;
    ekf_.P_ << 1,0,0,0,
               0,1,0,0,
               0,0,1000,0,
               0,0,0,1000;
    // Process covariance matrix
    ekf_.Q_ = MatrixXd(4, 4);  // to be filled
//cout << "ekf_.Q_ " <<endl;
    noise_ax = 9.0;
    noise_ay = 9.0;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{



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
//        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        float px,py,vx,vy;
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
        {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            float rho = measurement_pack.raw_measurements_(0);
            float phi = measurement_pack.raw_measurements_(1);
            float drho = measurement_pack.raw_measurements_(2);
            px = rho * cos(phi);
            py = rho * sin(phi);
            // assume direction of drho has small deviation from v
            vx = drho * cos(phi);
            vy = drho * sin(phi);
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
        {
            /**
            Initialize state.
            */

            px = measurement_pack.raw_measurements_(0);
            py = measurement_pack.raw_measurements_(1);
            // assume 0 since nothing is known about the velocity
            vx = 0;
            vy = 0;

        }

        ekf_.x_ << px,py,vx,vy;
        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        cout << "EKF initialized" << endl;

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
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) * 1e-6;
    previous_timestamp_ = measurement_pack.timestamp_;

    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;
    // state transition matrix
     ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, dt, 0,
          0, 1, 0, dt,
          0, 0, 1, 0,
          0, 0, 0, 1;
    // process covariance
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
               0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
               dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
               0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

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
        H_radar_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = H_radar_;
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
