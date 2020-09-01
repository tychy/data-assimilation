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
void printer(std::vector<double> v){
    rep(i, v.size()-1){
        std::cout << v[i] << " ";
    }
    std::cout << v[v.size()-1]<< std::endl;
}
std::vector<double> lorenz96(std::vector<double> v,std::vector<double> k,double t, double F){
    //引数　x,kとタイムステップ
    //返り値 f(v, t)
    int N = v.size();
    std::vector<double> v_in(N, 0.0);
    std::vector<double> v_out(N, 0.0);
    v_in = adder(v, muler(k,t));
    //std::cout << "v_in time step" << t << std::endl;
    //printer(v_in);
    //微分方程式の右辺を計算
    rep(i, N){
        v_out[i] = (v_in[(i+1)%N] - v_in[(i - 2 + N) % N]) * v_in[(i - 1 + N)%N] - v_in[i] + F;
    }
    //std::cout << "v_out time step" << t << std::endl;
    //printer(v_out);
    return v_out;
}

int main(){
    //定義
    int N = 40;
    double F = 8.0;
    double t_max = 2010;
    double t_th = 2000;
    double dt = 0.01;
    std::vector<double> x(N,F);//平衡
    //誤差
    x[0] += 0.01;
    std::vector<double> k1(N,0.0);
    std::vector<double> k2(N,0.0);
    std::vector<double> k3(N,0.0);
    std::vector<double> k4(N,0.0);
    std::vector<double> kmid(N,0.0);

    // loop
    for(double t=0.0;t<=t_max;t+=dt){
        //std::cout << "-------------STEP" << t/dt<< "start-------------"<< std::endl; 
        k1 = lorenz96(x, x, 0, F);//２つ目のxはダミーで0によって消える
        k2 = lorenz96(x, k1, dt/2, F);
        k3 = lorenz96(x, k2, dt/2, F);
        k4 = lorenz96(x, k3, dt, F);
        kmid = muler(adder(adder(k1, k4),muler(adder(k2, k3), 2.0)), dt/6);
        x = adder(x, kmid);
        /*
        std::cout << "k1" << std::endl; 
        printer(k1);
        std::cout << "k2" << std::endl; 
        printer(k2);
        std::cout << "k3" << std::endl; 
        printer(k3);
        std::cout << "k4" << std::endl; 
        printer(k4);
        std::cout << "kmid" << std::endl; 
        printer(kmid);
        std::cout << "x" << std::endl; 
        */
        if(t >= t_th){
            printer(x);
        }
    }
    return 0;
}

