#pragma once
#include <vector>
#include <map>

using namespace std;

vector<vector<double>> integrate_traj(vector<double> state, int n_times, double end_time, map<string,double> sys_param) ;
