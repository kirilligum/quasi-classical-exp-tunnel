#pragma once
#include <map>
#include <vector>
#include <random>

using namespace std;

struct tunnel {
  map<string,double> param;
  vector<double> make_initial_state(const map<string,double> param, int seed);
  double v1(vector<double> s) ;
  double v2(vector<double> s) ;
  double energy(vector<double> state);
  void eom(const vector<double> &state, vector<double> &dstates, const double);
};
