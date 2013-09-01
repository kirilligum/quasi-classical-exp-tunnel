#pragma once
#include <vector>
#include <functional>
#include <map>
#include <string>
//#include <iostream>
//#include <fstream>
//#include <boost/range/algorithm.hpp>
//#include <boost/range/algorithm_ext.hpp>
//#include <boost/numeric/odeint.hpp>

#include "tunnel.hpp"

using namespace std;


struct ode_step {
  const double dt;
  tunnel sys;
  //function<void(const vector<double>&,vector<double>&,const double)> eom;
  //function<double(const vector<double>&)> ham;
  map<string,double> param;
  double old_energy;
  double energy_threshold;
  ode_step(
    tunnel sys_,
    map<string,double> param_,
    double dt_,
    double old_energy_,
    double energy_threshold_
    //function<void(const vector<double>&,vector<double>&,const double)> eom_, 
    //function<double(const vector<double>&)> hamiltonian_energy
    );
  vector<double> operator()(vector<double> a, vector<double> b);
};

