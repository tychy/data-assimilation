#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>

#define ll long long
#define rep(i,n) for(int (i)=0;(i)<(n);(i)++)
#define Graph vector<vector<int>>;
#define iterG(next_v, G, v) for(auto next_v : G[v]
void printer(std::vector<double> v){
    rep(i, v.size()-1){
        std::cout << v[i] << " ";
    }
    std::cout << v[v.size()-1]<< std::endl;
}
double gen_mt(){
    //[0,1]の一様分布を作成
    std::mt19937 mt{ std::random_device{}() };
    std::uniform_real_distribution<double> dist(0, 1);
    return dist(mt);
}
double gen_norm(){
    constexpr double pi = std::acos(-1);
    return std::sqrt(-2*std::log(gen_mt())) * std::cos(2 * pi * gen_mt());
}
std::vector<double> mtv(int s){
    //sサイトのそれぞれのサイトに置いて独立な正規分布を作る
    std::vector<double> v(s);
    rep(i, s){
        v[i] = gen_norm();
    }
    return v;
}

int main(){
    int N = 40;
    for ( int i = 0 ; i < 1460 ; ++i )
        printer(mtv(N));
    return 0;
}

