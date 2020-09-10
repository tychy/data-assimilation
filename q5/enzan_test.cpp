#include <iostream>
#include <vector>
#include <cmath>
#define rep(i,n) for(int (i)=0;(i)<(n);(i)++)
#include "enzan.h"

int main(){
    int N = 3;
    std::vector<std::vector <double>> pf = genZ(N);
    std::vector<std::vector <double>> pa = genI(N);
    std::vector<double> x(N, 1.0);
    pf[0][0] += 2;
    pf[0][1] += -2;
    pf[0][2] += 3;
    pf[1][0] += 3;
    pf[1][1] += 2;
    pf[1][2] += -4;
    pf[2][0] += 4;
    pf[2][1] += -3;
    pf[2][2] += 2;
    
    x[0] = 7;
    x[1] = -5;
    x[2] =  4;
    pa[0][0] += 1;
    std::cout << "A" << std::endl;
    printernn(pf);
    std::cout << "B" << std::endl;
    printernn(pa);
    std::cout << "A^T" << std::endl;
    printernn(transpose(pf));
    std::cout << "A+B" << std::endl;
    printernn(pf*pa);
    std::cout << "x" << std::endl;
    printer(x);

    std::cout << "Ax = b then x?" << std::endl;
    printer(Axeqb(pf, x));

 
    return 0;
}

