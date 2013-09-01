#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <functional>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/numeric.hpp>
#include <boost/numeric/odeint.hpp>

#include "integrate_traj.hpp"
//#include "ran_pi.hpp"
//#include "p_from_n.hpp"
//#include "x_from_n.hpp"
#include "tunnel.hpp"
//#include "pcet.hpp"
//#include "spin_boson.hpp"
#include "ode_step.hpp"
//#include "j.hpp"
//#include "c.hpp"
#include "o.hpp"

using namespace std;
using namespace std::placeholders;

vector<vector<double>> integrate_traj(vector<double> state, int n_times, double end_time, map<string,double> sys_param) {
  vector<vector<double>> traj(n_times,vector<double>(state));
  tunnel sys{sys_param};
  auto energy = bind(&tunnel::energy,&sys,_1);
  //auto eom = bind(&tunnel::eom,&sys,_1,_2,_3);
  double old_energy = abs(energy(state));
  double energy_threshold = sys_param.at("energy_threshold_factor_absolute");
  boost::partial_sum(traj, traj.begin(),ode_step(sys,sys_param,end_time/n_times,old_energy,energy_threshold));
  //boost::partial_sum(traj, traj.begin(),ode_step(sys_param,end_time/n_times,old_energy,energy_threshold,eom,energy));
  traj.front().at(9)= energy(state);
  return traj;
}

