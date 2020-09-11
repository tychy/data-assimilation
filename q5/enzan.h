#ifndef ENZAN_HPP
#define ENZAN_HPP

#include <vector>

std::vector<double> adder(std::vector<double> v1, std::vector<double> v2);
std::vector<double> muler(std::vector<double> v1, double s);
std::vector<std::vector <double>> transpose(std::vector<std::vector<double> > b);
std::vector<std::vector <double>> transpose1d(std::vector<double> b);
std::vector<double> retranspose1d(std::vector<std::vector<double>> b);
std::vector<std::vector <double>> genI(int N);
std::vector<std::vector <double>> genZ(int N);
std::vector<double> operator+(const std::vector<double> &v1,const std::vector<double> &v2);
std::vector<std::vector <double>> operator+(const std::vector<std::vector <double>> &v1,const std::vector<std::vector <double>> &v2);
std::vector<double> operator-(const std::vector<double> &v1,const std::vector<double> &v2);
std::vector<std::vector <double>> operator-(const std::vector<std::vector <double>> &v1,const std::vector<std::vector <double>> &v2);
double operator*(const std::vector<double> &v1,const std::vector<double> &v2);
std::vector<std::vector <double>> operator*(const std::vector<std::vector <double>> &v1,const std::vector<std::vector <double>> &v2);
void printer(std::vector<double> v);
void printernn(std::vector<std::vector <double>> v);
std::vector<double> Axeqb(std::vector<std::vector<double>> v, std::vector<double> b);

#endif
