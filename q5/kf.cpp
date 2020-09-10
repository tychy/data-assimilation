#include <iostream>
#include <complex>
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
    //1dvector全体の加算器
    std::vector<double> v3;
    std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(v3), std::plus<double>());
    return v3;
}
std::vector<double> minuser(std::vector<double> v1, std::vector<double> v2){
    //1dvector全体
    std::vector<double> v3;
    std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(v3), std::minus<double>());
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
//ベクトル演算
std::vector<std::vector <double>> transpose(vector<vector<double> > b){
    if (b.size() == 0)
        return;
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
double norm(std::vector<double> v){
    double accum = 0.;
    for (int i = 0; i < v.size(); ++i) {
        accum += v[i] * v[i];
    }
    return sqrt(accum);
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

std::vector<double> M_6h(std::vector<double> x){
    // 6h進める
    //定義
    int N = 40;
    double F = 8.0;
    double t_max = 0.05;//
    double dt = 0.01;
    std::vector<double> k1(N,0.0);
    std::vector<double> k2(N,0.0);
    std::vector<double> k3(N,0.0);
    std::vector<double> k4(N,0.0);
    std::vector<double> kmid(N,0.0);

    // loop
    for(double t=0.01;t<=t_max;t+=dt){
        //std::cout << "-------------STEP" << t/dt<< "start-------------"<< std::endl; 
        k1 = lorenz96(x, x, 0, F);//２つ目のxはダミーで0によって消える
        k2 = lorenz96(x, k1, dt/2, F);
        k3 = lorenz96(x, k2, dt/2, F);
        k4 = lorenz96(x, k3, dt, F);
        kmid = muler(adder(adder(k1, k4),muler(adder(k2, k3), 2.0)), dt/6);
        x = adder(x, kmid);
    }
    return x;
}

int main(){
    int N = 40;
    std::vector<double> xf(N,0.0);//平衡
    std::vector<double> xa(N,0.0);
    std::vector<std::vector <double>> pf = genI(N);
    std::vector<std::vector <double>> pa = genI(N);
    std::vector<std::vector <double>> R = genI(N);
    std::vector<std::vector <double>> H = genI(N);
    std::vector<std::vector <double>> K = genZ(N);
    //printernn(pf);    
    for(int i=0;i<1460;i++){
        K = pf * transpose(H) * inverse(H * pf * transpose(H) + R);
        xa = xf + K * (y - H * xf);
        printer(xf);
        xf = M_6h(xf);
    }
    return 0;
}

