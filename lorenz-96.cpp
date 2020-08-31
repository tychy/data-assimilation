#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#define ll long long
#define rep(i,n) for(int (i)=0;(i)<(n);(i)++)
#define Graph vector<vector<int>>;
#define iterG(next_v, G, v) for(auto next_v : G[v]


std::vector<double> adder(std::vector<double> v1, std::vector<double> v2){
    //vector全体の加算器
    std::vector<double> v3;
    std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(v3), std::plus<double>());
    return v3;
}
std::vector<double> muler(std::vector<double> v1, double s){
    //vector全体の加算器
    std::vector<double> v3(v1.size());
    rep(i, v1.size()){
        v3[i] = v1[i] * s;
    }
    return v3;
}

std::vector<double> lorenz96(std::vector<double> v,std::vector<double> k,double t, double F){
    //引数　x,kとタイムステップ
    //返り値 f(v, t)
    int N = v.size();
    int ia, ib, ic;
    std::vector<double> v_in(N, 0.0);
    std::vector<double> v_out(N, 0.0);
    //v_in=v + k * t
    v_in = adder(v, muler(k,t));
    //微分方程式の右辺を計算
    rep(i, N){
        if(i = N - 1){
            ia = 0;
        }else{
            ia = i + 1;
        }
        if(i == 0){
            ib = N - 2;
            ic = N - 1;
        }else if(i == 1){
            ib = N - 1;
            ic = 0;
        }else{
        ib = i - 2;
        ic = i - 1;
        }
        v_out[i] = (v_in[ia] - v_in[ib]) * v_in[ic] - v_in[i] + F;
    }
    return v_out;
}
void printer(std::vector<double> v){
    rep(i, v.size()){
        std::cout << v[i] << ",";
    }
    std::cout << std::endl;
}
int main(){
    //定義
    int N = 40;
    double F = 8;
    double t_max = 1000;
    double dt = 0.2;
    std::vector<double> x(N,1.0);
    std::vector<double> k1(N,0.0);
    std::vector<double> k2(N,0.0);
    std::vector<double> k3(N,0.0);
    std::vector<double> k4(N,0.0);
    // loop
    for(double t=0.0;t+=dt;t<=t_max){
        k1 = lorenz96(x, x, 0, F);//２つ目のxはダミーで0によって消える
        k2 = lorenz96(x, k1, dt/2, F);
        k3 = lorenz96(x, k2, dt/2, F);
        k4 = lorenz96(x, k3, dt, F);
        x = adder(x, muler(adder(adder(k1, k4),muler(adder(k2, k3), 2.0)), dt/6));
        printer(x);
        printer(k1);

    }
    return 0;
}

