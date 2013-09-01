#include "state.hpp"

#include <vector>
#include <iostream>
using namespace std;

double get_dnn_dp(vector<double> s, int i) {
  if(i==1) {
    return s[5];
  }else if(i==2) {
    return s[3];
  } else {
    cout << " error wrong N in nN \n";
    return 1337.35505;
  }
}

double get_dnn_dx(vector<double> s, int i) {
  if(i==1) {
    return s[6];
  }else if(i==2) {
    return s[4];
  } else {
    cout << " error wrong N in nN \n";
    return 1337.35505;
  }
}

double get_nn(vector<double> s) {
  return (s[3]*s[5]+s[4]*s[6]);// p1*p2+x1*x2
}

double get_dn_dp(vector<double> s, int i) {
  if(i==1 || i==2) {
    auto j=i;
    j*=2;
    j+=1;
    return s[j];
  } else {
    cout << " error wrong N in nN \n";
    return 1337.35505;
  }
}

double get_dn_dx(vector<double> s, int i) {
  if(i==1 || i==2) {
    auto j=i;
    j*=2;
    j+=1;
    return s[j+1];
  } else {
    cout << " error wrong N in nN \n";
    return 1337.35505;
  }
}

double get_n::operator()(vector<double> s,int i) {
  //cout << "s5 " << s[5] << "  s6 " << s[6] << "\n";
  if(i==1 || i==2) {
    auto j=i;
    j*=2;
    j+=1;
    return 0.5*(s[j]*s[j]+s[j+1]*s[j+1])-n_shift;// (p*p+x*x)/2-n_shift
  } else {
    cout << " error wrong N in nN \n";
    return 1337.35505;
  }
}

