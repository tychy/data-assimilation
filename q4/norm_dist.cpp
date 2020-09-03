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
int main(){
    for ( int i = 0 ; i != 100000 ; ++i )
        std::cout << gen_norm() << std::endl;
    return 0;
}

