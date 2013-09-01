#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <numeric>
#include <iomanip>
//#include <boost/range/algorithm.hpp>
//#include <boost/range/algorithm_ext.hpp>
//#include <boost/range/numeric.hpp>
//#include <boost/math/constants/constants.hpp>
//#include <boost/numeric/odeint.hpp>

#include "o.hpp"
//#include "make_wb.hpp"
//#include "make_cb.hpp"
//#include "make_bath.hpp"
//#include "make_initial_state.hpp"
#include "state.hpp"
#include "integrate_traj.hpp"
#include "box.hpp"
#include "electronic_population.hpp"
#include "tunnel.hpp"
#include "lua_read.hpp"

using namespace std;
using namespace std::placeholders;

int main(int argc, const char *argv[]) {
  using namespace std; 
  typedef vector<double> vd; typedef vector<vector<double>> vvd; typedef vector<vector<vector<double>>> vvvd;
  cout << "Reading parameters from " << argv[1] << " ... \n";

  ofstream ofp("param.out");
  auto sp = Lua_read(argv[1]).read_config(vector<string>{"m","n1_shift","n2_shift","nn_coup","energy_threshold_factor_relative","energy_threshold_factor_absolute"});
  ofp << "map of system's parameters: \n"; for(auto i: sp) ofp<< setw(20) << i.first << " = " << i.second << "\n";
  
  auto ip = Lua_read(argv[1]).read_config(vector<string>{"pp","qp","n1","n2","n1_bin_start","n2_bin_start","n1_shift","n2_shift"});
  ofp << "map of initial parameters: \n"; for(auto i: ip) ofp<< setw(20) << i.first << " = " << i.second << "\n";
  
  auto tp = Lua_read(argv[1]).read_config(vector<string>{"n_times","end_time","time_step"});
  ofp << "map of integration parameters: \n"; for(auto i: tp) ofp<< setw(20) << i.first << " = " << i.second << "\n";
  
  auto mp = Lua_read(argv[1]).read_config(vector<string>{"seed","n1_shift","n2_shift","n1_bin_start","n2_bin_start","n1_bin_end","n2_bin_end","n_trajs"});
  ofp << "map of method parameters: \n"; for(auto i: mp) ofp<< setw(20) << i.first << " = " << i.second << "\n";
  
  auto seed = mp.at("seed");


  //cout << "Dealing with bath ... \n";
  //auto wb = make_wb(sp["n_modes"], sp["dw"]);  o(wb,"wb.dat"); 
  //auto cb = make_cb(wb,sp["eta"],sp["w_c"],sp["dw"]); o(cb,"cb.dat"); 

  cout << "Making initial states ... \n";
  tunnel sys{ip};
  vector<vector<double>> states(mp.at("n_trajs"));
  for(auto &i:states) {
    i= sys.make_initial_state(ip,++seed); 
  }
  o(states,"states.dat");

  
  //vector<double> dich;
  //transform(states.begin(),states.end(),back_inserter(dich),[](vector<double>v){return (v[1]);});
  //sort(dich.begin(),dich.end());
  //ofstream of("dich.dat");
  //for(auto i:dich) of<<i<<"\n";

  cout << "Integrating trajectories \n";
  vector<vector<vector<double>>> trajs;
  transform(states.begin(),states.end(),back_inserter(trajs), bind(integrate_traj,_1,tp["n_times"],tp["end_time"],sp));

  double an1=0, an2=0;
  for(auto is: trajs) {
    an1+=is.front().at(7);
    an2+=is.front().at(8);
  }
  cout << "tot n = " << an1 << "  " << an2 << "\n";
  cout << "ave n = " << an1/states.size() << "  " << an2/states.size() << "\n";
  
  auto traj = trajs.front();
  for(auto &i:traj) {
    i.push_back(get_n{mp.at("n1_shift")}(i,1));
    i.push_back(get_n{mp.at("n2_shift")}(i,2));
  }
  o(traj,"traj.dat");

  cout << "binning ... \n";
  auto el_pop = electronic_population(trajs,ip.at("n1"),ip.at("n2"),mp.at("n1_bin_end"),mp.at("n2_bin_end"),mp.at("n1_shift"),mp.at("n2_shift"));
  o(el_pop ,"el_pop.dat");
  
  return 0;
}
