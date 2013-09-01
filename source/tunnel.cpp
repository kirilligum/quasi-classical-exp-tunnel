#include <iostream>
#include "tunnel.hpp"
#include "p_from_n.hpp"
#include "x_from_n.hpp"
#include "state.hpp"

using namespace std;

vector<double> tunnel::make_initial_state(const map<string,double> param, int seed){
  std::mt19937 g(++seed);
  std::uniform_real_distribution<double> ran_pi(0.0,2*M_PI);
  double n1_bin_start = param.at("n1_bin_start");
  double n2_bin_start = param.at("n2_bin_start");
  double n1_shift = param.at("n1_shift");
  double n2_shift = param.at("n2_shift");
  double n1 = param.at("n1");
  double n2 = param.at("n2");
  double pp = param.at("pp");
  double qp = param.at("qp");
  uniform_real_distribution<double> ran_bin1(-n1_bin_start,n1_bin_start);
  uniform_real_distribution<double> ran_bin2(-n2_bin_start,n2_bin_start);
  const double q1=ran_pi(g),q2=ran_pi(g);
  vector<double> state (1+2+4+2+1+1+2+1+2,0.0);// time+pxpx+n1+n2+energy+traj_state+v1i+v2i+nn_dyn+v1+v2
  state[0]=0.0;
  normal_distribution<double> normal(0.0,1.0);
  state[1]=pp+normal(g);
  state[2]=qp+normal(g);
  state[3]=p_from_n{n1_shift}(n1,q1);
  state[4]=x_from_n{n1_shift}(n1,q1);
  state[5]=p_from_n{n2_shift}(n2,q2);
  state[6]=x_from_n{n2_shift}(n2,q2);
  state[7] = get_n{n1_shift}(state,1);
  state[8] = get_n{n2_shift}(state,2);
  double v1e=v1(state);
  double v2e=v2(state);
  cout << "   v1  " << v1e << "  v2e   " << v2e << endl;
  double vvnorm = (v1e/v2e+v2e/v1e);
  n1+=ran_bin1(g)*v2e/v1e/vvnorm;
  n2+=ran_bin2(g)*v1e/v2e/vvnorm;
  state[3]=p_from_n{n1_shift}(n1,q1);
  state[4]=x_from_n{n1_shift}(n1,q1);
  state[5]=p_from_n{n2_shift}(n2,q2);
  state[6]=x_from_n{n2_shift}(n2,q2);
  state[7] = get_n{n1_shift}(state,1);
  state[8] = get_n{n2_shift}(state,2);
  state.at(9)=0.0;
  state.at(10)=0.0;// 0 if traj is good, 1 is traj is bad
  state[11]=v1(state);// needed to see when the n dynamics should be truncated
  state[12]=v2(state);// needed to see when the n dynamics should be truncated
  state[13]=0.0;// 1 truncates the dynamics
  return state;
}

double tunnel::v1(vector<double> s) {
  double qp = s[2];
  return 0.1*exp(qp);
}

double tunnel::v2(vector<double> s) {
  double qp = s[2];
  return 0.1*exp(-qp);
}

double tunnel::energy(vector<double> s) {
  double m=param.at("m");
  double nn_coup=param.at("nn_coup");
  double empty
    , pp = s[1]
    , qp = s[2]
    , n1_shift = param.at("n1_shift")
    , n2_shift = param.at("n2_shift")
    , n1 = get_n{n1_shift}(s,1)
    , n2 = get_n{n2_shift}(s,2)
    , nn = get_nn(s)
    , kin = pp*pp*0.5/m
    , v12 = nn_coup*exp(-qp*qp)
    //, v12 = 0.0
    ;
  //return kin + v1*n1 + v2*n2 +v12*nn;
  //if(n2!=0.0) cout <<"error n2 is " << n2 << " \n";
  return kin + v1(s)*n1 + v2(s)*n2 +v12*nn;
  //return kin;
}

void tunnel::eom(const vector<double> &s, vector<double> &ds, const double) {
  double m=param.at("m");
  double nn_coup=param.at("nn_coup");
  double empty
    , pp=s[1]
    , qp=s[2]
    , n1_shift = param.at("n1_shift")
    , n2_shift = param.at("n2_shift")
    , n1 = get_n{n1_shift}(s,1)
    , n2 = get_n{n2_shift}(s,2)
    , nn = get_nn(s)
    , kin = pp*pp*0.5
    , v12 = nn_coup*exp(-qp*qp)
    , dv1_dr = 0.1*exp(qp)
    , dv2_dr = -0.1*exp(-qp)
    , dv12_dr = -2*nn_coup*qp*exp(-qp*qp)
    //, v12 = 0.0
    //, dv1_dr = 1.0
    //, dv2_dr = 1.0
    //, dv12_dr = 0.0
    , dn1_dp1 = get_dn_dp(s,1)
    , dn1_dx1 = get_dn_dx(s,1)
    , dn2_dp2 = get_dn_dp(s,2)
    , dn2_dx2 = get_dn_dx(s,2)
    , dnn_dp1 = get_dnn_dp(s,1)
    , dnn_dx1 = get_dnn_dx(s,1)
    , dnn_dp2 = get_dnn_dp(s,2)
    , dnn_dx2 = get_dnn_dx(s,2)
    ;
  //for(auto &i:ds) i=0.0;
  ds[1] =
    dv1_dr*n1+dv2_dr*n2+dv12_dr*nn;
  ds[1] *=-1; // - d_sys_ham_dr;
  ds[2] = // d_sys_ham_dp;
    pp/m;
  //if(s[13] == 0) {
  //if(max(v1(s),v2(s))<=max(s[11],s[12]) ){
    ds[3] = // d_sys_ham_dx1;
      v1(s)*dn1_dx1+v12*dnn_dx1;
    ds[3] *=-1;
    ds[4] = // d_sys_ham_dp1;
      v1(s)*dn1_dp1+v12*dnn_dp1;
    ds[5] = // d_sys_ham_dx2;
      v2(s)*dn2_dx2+v12*dnn_dx2;
    ds[5] *=-1;
    ds[6] = // d_sys_ham_dp2;
      v2(s)*dn2_dp2+v12*dnn_dp2;
  //} else {
    //ds[3]=0.0;
    //ds[4]=0.0;
    //ds[5]=0.0;
    //ds[6]=0.0;
  //}

  ds[0] = 0.0;
  ds[7] = 0.0;
  ds[8] = 0.0;
  ds[9] = 0.0;
  ds[10] = 0.0;
  ds[0] = 0.0;
}
