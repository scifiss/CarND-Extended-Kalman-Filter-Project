#include <iostream>

#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
//  cout << "Tools::CalculateRMSE" <<endl;

  VectorXd rmse(4);
  rmse << 0,0,0,0;  // for px, py, vx, vy components

  //each component c has 1...n estimations,
  // a(i) = est(i) - truth(i)
  // rmse(c) = sqrt(  sum( a(1)^2+a(2)^2+...+a(n)^2 )/n  )

  // check the validation of inputs
  if (estimations.size() == 0
      || estimations.size() != ground_truth.size() )
  {

      cout << "Sizes of estimation and ground truth don't match or is zero." << endl;
      return rmse;
  }
  // squared sum of residuals for each of the n estimations
  for (unsigned int i=0; i< estimations.size(); i++)
  {
      // the ith data points, each has px, py, vx, vy components
      VectorXd residual = estimations[i] - ground_truth[i];
      // each estimation, ground truth, and residual
      residual = residual.array() * residual.array();
      rmse += residual;

  }
  // mean squared sum
  rmse = rmse/estimations.size();
  // root of each component
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
//cout << "Tools::CalculateJacobian" <<endl;
  // the following code is from the course
  	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "Error: CalculateJacobian (). Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}
