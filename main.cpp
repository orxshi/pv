#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

//const int n = 120;
const int n = 4320;
const double GM = 1.4;
const double cp = 1.004;
const double R = 0.287;
const double TL = 300;
const double TH = 600;
const double VL = 0.1;
const double cr = 6; // compression ratio
const double m = 1;
        
double P[n];
double V[n];
double T[n];

int nquad = n / 4;
double VH = 0.1 * cr;

double dw(int j)
{
    if (j == 0)
    {
        return 0.;
    }

    return P[j] * (V[j] - V[j-1]);
}

double work = 0.;

void isothermal_heat_rejection(int k, double VI)
{
    double inc = (VH - VI) / (nquad+1);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VH - i * inc;
        T[j] = TL;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }

    //work += m * R * TL * std::log(VI/VH);
}

void isochoric_compression(int k)
{
    double inc = (TH - TL) / (nquad);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VL;
        T[j] = TL + i * inc;
        P[j] = m * R * T[j] / V[j];
    }
}

void isentropic_compression(int k)
{
    double inct = (TH - TL) / (nquad);
    double VI = VL * std::pow(TH/TL,1/(GM-1));
    double incv = (VI - VL) / (nquad);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VI - i * incv;
        T[j] = TL + i * inct;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }
}

void isentropic_expansion(int k)
{
    double inct = (TH - TL) / (nquad);
    double VI = VH * std::pow(TL/TH,1/(GM-1));
    double incv = (VH - VI) / (nquad);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VI + i * incv;
        T[j] = TH - i * inct;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }
}

void isothermal_heat_addition(int k, double VI)
{
    double inc = (VI - VL) / (nquad);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VL + i * inc;
        T[j] = TH;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }
}

void isochoric_expansion(int k)
{
    double inc = (TH - TL) / (nquad);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VH;
        T[j] = TH - i * inc;
        P[j] = m * R * T[j] / V[j];
    }
}

void print(std::string filename)
{
    std::ofstream out;
    out.open(filename);

    for (int i=0; i<n; ++i)
    {
        out << i;
        out << " ";
        out << P[i];
        out << " ";
        out << V[i];
        out << " ";
        out << T[i];
        out << std::endl;
    }

    out.close();
}

void stirling()
{
    work = 0.;

    isothermal_heat_rejection(0, VL);
    isochoric_compression(1);
    isothermal_heat_addition(2, VH);
    isochoric_expansion(3);

    print("stirling");

    std::cout << "net work from stirling: " << work << std::endl;
}

void carnot()
{
    work = 0.;

    isentropic_compression(0);
    isentropic_expansion(2);
    isothermal_heat_addition(1, V[nquad*2]);
    isothermal_heat_rejection(3, V[0]);

    print("carnot");

    std::cout << "net work from carnot: " << work << std::endl;
}


int main()
{
    stirling();
    carnot();
        

    return 0;
}
