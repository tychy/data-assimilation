#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;
#define ll long long
#define rep(i,n) for(int (i)=0;(i)<(n);(i)++)
#define Graph vector<vector<int>>;
#define iterG(next_v, G, v) for(auto next_v : G[v]
//ベクトル演算
std::vector<std::vector <double>> transpose(vector<vector<double> > b){
    if (b.size() == 0)
        return b;
    vector<vector<double> > trans_vec(b[0].size(), vector<double>());
    for (int i = 0; i < b.size(); i++){
        for (int j = 0; j < b[i].size(); j++){
            trans_vec[j].push_back(b[i][j]);
        }
    }
    return trans_vec;
}

std::vector<std::vector <double>> genI(int N){
    std::vector<std::vector <double>> I(N, std::vector<double>(N,0.0));
    rep(i, N){
        I.at(i).at(i) = 1.0;
    }
    return I;
}
std::vector<std::vector <double>> genZ(int N){
    std::vector<std::vector <double>> I(N, std::vector<double>(N,0.0));
    return I;
}
//operator
std::vector<double> operator+(const std::vector<double> &v1,const std::vector<double> &v2){
    //1d
    std::vector<double> ans = v1;
    rep(i, v1.size())
        ans[i] += v2[i];
    return ans;
}
std::vector<std::vector <double>> operator+(const std::vector<std::vector <double>> &v1,const std::vector<std::vector <double>> &v2){
    //2d
    std::vector<std::vector <double>> ans = genZ(v1.size());   
    rep(i, v1.size()){
        ans[i] = v1[i] + v2[i];
    }
    return ans;
}
double operator*(const std::vector<double> &v1,const std::vector<double> &v2){
    //1d ちょっと特殊 n * n = scalar
    double ans = 0.0;
    rep(i, v1.size())
        ans += v1[i] * v2[i];
    return ans;
}
std::vector<std::vector <double>> operator*(const std::vector<std::vector <double>> &v1,const std::vector<std::vector <double>> &v2){
    //2d
    std::vector<std::vector <double>> ans = genZ(v1.size());   
    std::vector<std::vector <double>> mid = transpose(v2);  

    rep(i, v1.size()){
        rep(j, v1[i].size()){
            ans[i][j] = v1[i] * mid[j];
        }
    }
    return ans;
}
void printer(std::vector<double> v){
    rep(i, v.size()-1){
        std::cout << v[i] << " ";
    }
    std::cout << v[v.size()-1]<< std::endl;
}
void printernn(std::vector<std::vector <double>> v){
    // n n行列のプリント
    rep(i, v.size())
        printer(v.at(i));
}

int main(){
    int N = 2;
    std::vector<std::vector <double>> pf = genI(N);
    std::vector<std::vector <double>> pa = genI(N);
    pf[0][0] += 1;
    pf[1][0] += 1;

    pa[0][0] += 1;
    std::cout << "A" << std::endl;
    printernn(pf);
    std::cout << "B" << std::endl;
    printernn(pa);
    std::cout << "A^T" << std::endl;
    printernn(transpose(pf));
    std::cout << "A+B" << std::endl;
    printernn(pf*pa);
    return 0;
}

