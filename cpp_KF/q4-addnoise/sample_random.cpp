#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
using namespace std;
#define ll long long
#define rep(i,n) for(int (i)=0;(i)<(n);(i)++)
#define Graph vector<vector<int>>;
#define iterG(next_v, G, v) for(auto next_v : G[v]
int main(){
    std::mt19937 mt{ std::random_device{}() };
    std::uniform_int_distribution<int> dist(1, 6);
    for ( int i = 0 ; i != 10 ; ++i )
      std::cout << dist(mt) << std::endl;
    return 0;
}

