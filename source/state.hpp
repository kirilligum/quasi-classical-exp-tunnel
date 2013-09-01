#pragma once
#include <vector>
using namespace std;

double get_nn(vector<double> s) ;
double get_dn_dp(vector<double> s, int i) ;
double get_dn_dx(vector<double> s, int i) ;
double get_dnn_dp(vector<double> s, int i) ;
double get_dnn_dx(vector<double> s, int i) ;

struct get_n {
  double n_shift;
  double operator()(vector<double> s,int i) ;
};

//struct get_n1 {
  //double n_shift;
  //double operator()(vector<double> s) ;
//};

//struct get_n2 {
  //double n_shift;
  //double operator()(vector<double> s) ;
//};
