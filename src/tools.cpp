#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	// ... your code here
	if(estimations.size() == 0){
	    cout << "Estimation vector size is zero." << endl;
	    return rmse;
	}
	if(estimations.size() != ground_truth.size()){
	    cout << "Estimtion vector size isn't match ground truth vector size." << endl;
	    return rmse;
	}
	

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
    		VectorXd tmp(4);
        	tmp =  (ground_truth[i] - estimations[i]);
		tmp = tmp.array() * tmp.array();
       		rmse += tmp;
	}

	//calculate the mean
    	rmse /= estimations.size();

	//calculate the squared root
    	rmse = rmse.array().sqrt();

	return rmse;
}