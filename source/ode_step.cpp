#include <vector>
#include <functional>
//#include <iostream>
//#include <fstream>
//#include <boost/range/algorithm.hpp>
//#include <boost/range/algorithm_ext.hpp>
#include <boost/numeric/odeint.hpp>

#include "ode_step.hpp"
#include "tunnel.hpp"
#include "p_from_n.hpp"
#include "x_from_n.hpp"
#include "state.hpp"

using namespace std;
using namespace std::placeholders;


ode_step::ode_step(
    tunnel sys_,
    map<string,double> param_,
    double dt_,
    double old_energy_,
    double energy_threshold_
    //function<void(const vector<double>&,vector<double>&,const double)> eom_, 
    //function<double(const vector<double>&)> hamiltonian_energy
    )
  : param(param_),dt(dt_)
  , sys(sys_)
  , old_energy(old_energy_)
  , energy_threshold(energy_threshold_)
  //, eom(eom_)
  //, ham(hamiltonian_energy){}
  {}

vector<double> ode_step::operator()(vector<double> a, vector<double> b) {
  using namespace boost::numeric::odeint;
  vector<double> current(a);
  double n1_shift = param.at("n1_shift");
  double n2_shift = param.at("n2_shift");
  double current_energy = sys.energy(current);
  double energy_step = abs(old_energy-current_energy);
  //if(max(sys.v1(current),sys.v2(current))>=max(current[11],current[12]) ) {
  //if(abs(sys.v1(current)*get_n{n1_shift}(current,1))>=max(current[11],current[12]) ) {
  if(abs(sys.v2(current)*1)>=max(current[11],current[12]) ) {
    current[13] = 1.0;
    double q1=atan2(current[3],current[4]);
    current[3]=p_from_n{n1_shift}(1.0,q1);
    current[4]=x_from_n{n1_shift}(1.0,q1);
    double q2=atan2(current[5],current[6]);
    current[5]=p_from_n{n2_shift}(0.0,q2);
    current[6]=x_from_n{n2_shift}(0.0,q2);
  }
  if(abs(sys.v1(current)*1)>=max(current[11],current[12]) ) {
    current[13] = 1.0;
    double q1=atan2(current[3],current[4]);
    current[3]=p_from_n{n1_shift}(0.0,q1);
    current[4]=x_from_n{n1_shift}(0.0,q1);
    double q2=atan2(current[5],current[6]);
    current[5]=p_from_n{n2_shift}(1.0,q2);
    current[6]=x_from_n{n2_shift}(1.0,q2);
  }
  //if(energy_step>energy_threshold) current[10]=1.0;
  //if(current[10] ==0){
    auto eom = bind(&tunnel::eom,&sys,_1,_2,_3);
    bulirsch_stoer<std::vector<double>> controlled_stepper(1.0e-9,1.0e-9,0.0,0.0);
    //typedef runge_kutta_cash_karp54<std::vector<double>> error_stopper_type;
    //typedef bulirsch_stoer<error_stopper_type> controlled_stepper_type;
    //typedef controlled_runge_kutta<error_stopper_type> controlled_stepper_type;
    //controlled_stepper_type controlled_stepper(default_error_checker<double>(1.0e-6,1.0e-6,1.0,1.0));
    //controlled_stepper_type controlled_stepper(default_error_checker<double>(1.0e-9,1.0e-9,1.0,1.0));
    integrate_const(controlled_stepper,eom, current, a[0],a[0]+dt,dt*1e-3);
    current[7] = get_n{n1_shift}(current,1);
    current[8] = get_n{n2_shift}(current,2);
    current.at(9)=sys.energy(current);
    current[14]=sys.v1(current);
    current[15]=sys.v2(current);
    old_energy = current_energy;
    energy_threshold = param.at("energy_threshold_factor_absolute")+ param.at("energy_threshold_factor_relative")*energy_step;
  //}
  current[0]+=dt;
  return current;
}
