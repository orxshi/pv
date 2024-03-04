#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cassert>

//const int n = 120;
const int n = 4320;
const double GM = 1.4;
const double cp = 1.004;
const double R = 0.287;
const double TL = 300;
const double TH = 600;
const double VL = 0.1;
const double cr = 6; // compression ratio
const double VH = cr * VL;
const double m = 1;
        
double P[n];
double V[n];
double T[n];

int nquad = n / 4;

double dw(int j)
{
    if (j == 0)
    {
        return 0.;
    }

    return P[j] * (V[j] - V[j-1]);
}

//double dq(int j)
//{
//    if (j == 0)
//    {
//        return 0.;
//    }
//
//    return T[j] * (S[j] - S[j-1]);
//}

double work = 0.;

void isobaric_heat_rejection(int k, double PC)
{
    double inc = (TH - TL) / (nquad+1);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        T[j] = TH - i * inc;
        P[j] = PC;
        V[j] = m * R * T[j] / P[j];

        work += dw(j);
    }
}

void isobaric_heat_addition(int k, double PC)
{
    double inc = (TH - TL) / (nquad+1);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        T[j] = TL + i * inc;
        P[j] = PC;
        V[j] = m * R * T[j] / P[j];

        work += dw(j);
    }
}

void isothermal_heat_rejection(int k, double Vmin, double Vmax)
{
    double inc = (Vmax - Vmin) / (nquad+1);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        //V[j] = VH - i * inc;
        V[j] = Vmax - i * inc;
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

double isentropic_full_compression(int k)
{
    double incv = (VH - VL) / (nquad);
    double TI = TL * std::pow(VH/VL,(GM-1));
    double inct = (TI - TL) / (nquad);

    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VH - i * incv;
        T[j] = TL + i * inct;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }

    return TI;
}

void isentropic_full_expansion(int k, double TH)
{
    double incv = (VH - VL) / (nquad);
    double TI = TH * std::pow(VH/VL,(GM-1));
    double inct = (TH - TI) / (nquad);

    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VL + i * incv;
        T[j] = TH - i * inct;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
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

void isentropic_compression(int k, double TI)
{
    double inct = (TI - TL) / (nquad);
    double VI = VL * std::pow(TI/TL,1/(GM-1));
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
    assert(VI < VH);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VI + i * incv;
        T[j] = TH - i * inct;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }
}

void isentropic_expansion(int k, double TI)
{
    double inct = (TI - TL) / (nquad);
    double VI = VH * std::pow(TL/TH,1/(GM-1));
    double incv = (VH - VI) / (nquad);
    assert(VI < VH);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VI + i * incv;
        T[j] = TI - i * inct;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }
}

void isothermal_heat_addition_Q(int k, double Vmin, double QH)
{
    double VI = VL * std::exp(QH/m/R/TH); 
    double inc = (VI - Vmin) / (nquad);
    assert(VI > Vmin);

    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = Vmin + i * inc;
        T[j] = TH;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }
}

void isothermal_heat_addition(int k, double Vmin, double Vmax)
{
    double inc = (Vmax - Vmin) / (nquad);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        //V[j] = VL + i * inc;
        V[j] = Vmin + i * inc;
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

void isochoric_heat_addition(int k, double TI, double QH)
{
    double TH_real = TI + QH / m / cp; 
    double inc = (TH_real - TI) / (nquad);
    assert(inc > 0);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VL;
        T[j] = TI + i * inc;
        P[j] = m * R * T[j] / V[j];
    }
}

void isochoric_heat_addition(int k, double TI)
{
    double inc = (TH - TI) / (nquad);
    assert(inc > 0);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VL;
        T[j] = TI + i * inc;
        P[j] = m * R * T[j] / V[j];
    }
}

void isochoric_heat_rejection(int k, double TI)
{
    double inc = (TI - TL) / (nquad);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VH;
        T[j] = TI - i * inc;
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

    isothermal_heat_rejection(0, VL, VH);
    isochoric_compression(1);
    isothermal_heat_addition(2, VL, VH);
    isochoric_expansion(3);

    print("stirling");

    std::cout << "net work from stirling: " << work << std::endl;
}

void carnot1()
{
    work = 0.;

    isentropic_compression(0);
    isentropic_expansion(2);
    isothermal_heat_addition(1, VL, V[nquad*2]);
    isothermal_heat_rejection(3, V[0], VH);

    print("carnot");

    std::cout << "net work from carnot: " << work << std::endl;
}

void carnot2()
{
    // does not seem possible.
    work = 0.;

    isothermal_heat_addition_Q(1, VL, 100);
    isentropic_expansion(2, T[2*nquad-1]);
    isentropic_compression(0, T[2*nquad-1]);
    isothermal_heat_rejection(3, V[0], VH);

    print("carnot");

    std::cout << "net work from carnot: " << work << std::endl;
}

void ericsson()
{
    work = 0.;

    double PL = m * R * TH / VH; // 4
    double PH = m * R * TL / VL; // 2

    isobaric_heat_rejection(0, PL);
    isothermal_heat_rejection(1, VL, V[nquad*1-1]);
    isobaric_heat_addition(2, PH);
    isothermal_heat_addition(3, V[nquad*3-1], VH);

    print("ericsson");

    std::cout << "net work from ericsson: " << work << std::endl;
}

void otto()
{
    // setting compression ratio and temperature range causes contradiction sometimes. If suitable parameters are found, in this case, Carnot fails. Failure happens this way: Max volume may be exceeded or max temperature may be exceeded.
    work = 0.;

    double TI = isentropic_full_compression(0);
    isochoric_heat_addition(1, T[nquad-1]);
    isentropic_full_expansion(2, T[2*nquad-1]);
    isochoric_heat_rejection(3, T[3*nquad-1]);

    print("otto");

    std::cout << "net work from otto: " << work << std::endl;
}


int main()
{
    stirling();
    carnot1();
    ericsson();
    //otto(); // plot otto. you will see that it gives wrong plot.

    //https://www.quora.com/Why-do-we-add-heat-at-a-constant-pressure-in-a-gas-turbine-plant
    // https://en.wikipedia.org/wiki/Brayton_cycle
    // https://dieselnet.com/tech/diesel_engines.php
    // https://engineering.stackexchange.com/questions/32204/why-does-combustion-happen-at-constant-pressure-in-jet-engines
    // https://en.wikipedia.org/wiki/Pulsejet
    // Ramjet
    // mixed pdv and vdp
    // https://www.physicsforums.com/threads/why-does-heat-addition-happen-at-constant-pressure-in-diesel-cycle.746745/
    // https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node40.html
        

    return 0;
}
