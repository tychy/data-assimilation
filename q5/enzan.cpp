#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include "enzan.h"
#define rep(i,n) for(int (i)=0;(i)<(n);(i)++)
//ベクトル演算
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
std::vector<std::vector <double>> transpose(std::vector<std::vector<double> > b){
    if (b.size() == 0)
        return b;
    std::vector<std::vector<double> > trans_vec(b[0].size(), std::vector<double>());
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
std::vector<double> operator-(const std::vector<double> &v1,const std::vector<double> &v2){
    //1d
    std::vector<double> ans = v1;
    rep(i, v1.size())
        ans[i] -= v2[i];
    return ans;
}
std::vector<std::vector <double>> operator-(const std::vector<std::vector <double>> &v1,const std::vector<std::vector <double>> &v2){
    //2d
    std::vector<std::vector <double>> ans = genZ(v1.size());   
    rep(i, v1.size()){
        ans[i] = v1[i] - v2[i];
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
    //v1[0].size != v2.sizeのときアラートを出したい
    std::vector<std::vector <double>> ans(v1.size());   
    std::vector<std::vector <double>> mid = transpose(v2);
    rep(i, mid.size()){
            rep(j, v1.size()){
                ans[j].push_back(v1[j] * mid[i]);
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
std::vector<double> Axeqb(const std::vector<std::vector<double>> &v, std::vector<double> &b){
    //Ax = bのxを返す
    int N = b.size();
    std::vector<std::vector<double>> A_c = v;
    std::vector<double> x(N, 0.0);
    std::vector<double> buff(N, 0.0);//枢軸の入れ替え用のバッファ
    std::vector<double> b_c = b;

    //std::vector<int> p(N, 0);//枢軸の操作
    //rep(i, N)
    //    p[i] = i;
    int p_cur;
    double p_max, ans, s;
    rep(i, N-1){
        p_max = A_c[i][i];
        p_cur = i;
        for(int j=i;j<N;j++){
            if (A_c[j][i] > p_max){
                p_max = A_c[j][i];
                p_cur = j;
            }
        }
        //もしp_maxが0ならば正則出ないと言えるのでエラーを出したい
        //入れ替え
        if(p_cur != i){
            //A_c b_c 
            buff = A_c[i];
            A_c[i] = A_c[p_cur];
            A_c[p_cur] = buff;

            s = b_c[i];
            b_c[i] = b_c[p_cur];
            b_c[p_cur] = s;
        }

        //前進消去
        for(int k=i+1;k<N;k++){
            b_c[k] = b_c[k] - b_c[i] * A_c[k][i]/A_c[i][i];
            A_c[k] = A_c[k] - muler(A_c[i], A_c[k][i]/A_c[i][i]);

        }
    }
    //後退代入
    for(int i = N-1;i >= 0;i--){
        ans = b_c[i];
        for(int k = i + 1;k < N;k++){
            ans -= A_c[i][k] * x[k];
        }
        ans = ans / A_c[i][i];
        x[i] = ans;
    }
    return x;
}


