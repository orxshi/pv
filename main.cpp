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
const double cv = cp - R;
const double TL = 300;
//double TH = 600;
const double VL = 0.1;
//const double cr = 6; // compression ratio
//const double VH = cr * VL;
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

double isentropic_V(double V1, double T1, double T2) 
{
    return V1 * std::pow(T1 / T2, 1 / (GM-1));
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

void isobaric_heat_rejection(int k, double PC, double TH)
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

void isobaric_heat_addition(int k, double PC, double TH)
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
    assert(inc > 0);

    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = Vmax - i * inc;
        T[j] = TL;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }

    //work += m * R * TL * std::log(VI/VH);
}

void isochoric_compression(int k, double TH)
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

double isentropic_full_compression(int k, double VH)
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

void isentropic_full_expansion(int k, double TH, double VH)
{
    double incv = (VH - VL) / (nquad);
    assert(incv > 0);
    //double TI = TH * std::pow(VH/VL,(GM-1));
    double TI = TH * std::pow(VL/VH,(GM-1));
    std::cout << TI << std::endl;
    double inct = (TH - TI) / (nquad);
    assert(inct > 0);

    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = VL + i * incv;
        T[j] = TH - i * inct;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }
}

//void isentropic_compression(int k, double TL, double TH)
//{
//    double inct = (TH - TL) / (nquad);
//    double VH = VL * std::pow(TL / TH, 1 / (GM-1));
//    double incv = (VI - VL) / (nquad);
//    for (int i=0; i<nquad; ++i)
//    {
//        int j = i + k * nquad;
//        V[j] = VI - i * incv;
//        T[j] = TL + i * inct;
//        P[j] = m * R * T[j] / V[j];
//
//        work += dw(j);
//    }
//}

double isentropic_compression(int k, double T1, double T2, double V2)
{
    // compression from (T1, V1) to (T2, V2) where V1 > V2 and T1 < T2.
    
    double inct = (T2 - T1) / (nquad);
    assert(inct > 0);
    double V1 = isentropic_V(V2, T2, T1);
    double incv = (V1 - V2) / (nquad);
    assert(incv > 0);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = V1 - i * incv;
        T[j] = T1 + i * inct;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }

    return V1;
}

double isentropic_expansion_low(int k, double T1, double T2, double V1)
{
    // expansion from (T1, V1) to (T2, V2) where V2 > V1 and T1 > T2.

    double inct = (T1 - T2) / (nquad);
    assert(inct > 0);
    double V2 = isentropic_V(V1, T1, T2);
    double incv = (V2 - V1) / (nquad);
    assert(incv > 0);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = V1 + i * incv;
        T[j] = T1 - i * inct;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }

    return V2;
}

double isentropic_expansion_high(int k, double T1, double T2, double V2)
{
    // expansion from (T1, V1) to (T2, V2) where V2 > V1 and T1 > T2.

    double inct = (T1 - T2) / (nquad);
    assert(inct > 0);
    double V1 = isentropic_V(V2, T2, T1);
    double incv = (V2 - V1) / (nquad);
    assert(incv > 0);
    for (int i=0; i<nquad; ++i)
    {
        int j = i + k * nquad;
        V[j] = V1 + i * incv;
        T[j] = T1 - i * inct;
        P[j] = m * R * T[j] / V[j];

        work += dw(j);
    }

    return V1;
}

//void isentropic_expansion(int k, double TI, double VH)
//{
//    // expansion from (T1, V1) to (T2, V2) where V2 > V1 and T2 < T1.
//
//    double inct = (T2 - T1) / (nquad);
//    assert(inct > 0);
//    double VI = VH * std::pow(TL/TI,1/(GM-1)); // mistake?
//    double incv = (VH - VI) / (nquad);
//    assert(VI < VH);
//    for (int i=0; i<nquad; ++i)
//    {
//        int j = i + k * nquad;
//        V[j] = VI + i * incv;
//        T[j] = TI - i * inct;
//        P[j] = m * R * T[j] / V[j];
//
//        work += dw(j);
//    }
//}

//void isothermal_heat_addition_Q(int k, double Vmin, double QH)
//{
//    double VI = VL * std::exp(QH/m/R/TH); 
//    double inc = (VI - Vmin) / (nquad);
//    assert(VI > Vmin);
//
//    for (int i=0; i<nquad; ++i)
//    {
//        int j = i + k * nquad;
//        V[j] = Vmin + i * inc;
//        T[j] = TH;
//        P[j] = m * R * T[j] / V[j];
//
//        work += dw(j);
//    }
//}

void isothermal_heat_addition(int k, double Vmin, double Vmax, double TH)
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

void isochoric_expansion(int k, double TH, double VH)
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

void isochoric_heat_addition_QH(int k, double TI, double QH)
{
    double TH_real = TI + QH / m / cv;
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

void isochoric_heat_addition_TH(int k, double TI, double TH)
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

void isochoric_heat_rejection(int k, double TI, double VH)
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

void stirling(double cr, double TH)
{
    work = 0.;

    double VH = cr * VL;

    isothermal_heat_rejection(0, VL, VH);
    isochoric_compression(1, TH);
    isothermal_heat_addition(2, VL, VH, TH);
    isochoric_expansion(3, TH, VH);

    print("stirling");

    std::cout << "net work from stirling: " << work << std::endl;
}

void carnot1(double cr, double TH)
{
    work = 0.;

    double VH = cr * VL;

    double V2 = isentropic_compression(0, TL, TH, VL);
    isothermal_heat_rejection(3, V2, VH);
    double V4 = isentropic_expansion_high(2, TH, TL, VH);
    isothermal_heat_addition(1, VL, V4, TH);

    print("carnot");

    std::cout << "net work from carnot: " << work << std::endl;
}

//void carnot2()
//{
//    // does not seem possible.
//    work = 0.;
//
//    isothermal_heat_addition_Q(1, VL, 100);
//    isentropic_expansion(2, T[2*nquad-1]);
//    isentropic_compression(0, T[2*nquad-1]);
//    isothermal_heat_rejection(3, V[0], VH);
//
//    print("carnot");
//
//    std::cout << "net work from carnot: " << work << std::endl;
//}

void ericsson(double cr, double TH)
{
    work = 0.;

    double VH = cr * VL;

    double PL = m * R * TH / VH; // 4
    double PH = m * R * TL / VL; // 2

    isobaric_heat_rejection(0, PL, TH);
    isothermal_heat_rejection(1, VL, V[nquad*1-1]);
    isobaric_heat_addition(2, PH, TH);
    isothermal_heat_addition(3, V[nquad*3-1], VH, TH);

    print("ericsson");

    std::cout << "net work from ericsson: " << work << std::endl;
}

//void diesel()
//{
//    work = 0.;
//
//    double PL = m * R * TH / VH; // 4
//    double PH = m * R * TL / VL; // 2
//
//    double TI = isentropic_full_compression(0);
//    isobaric_heat_addition(2, PH);
//    isentropic_full_expansion(2, T[2*nquad-1]);
//    isochoric_heat_rejection(3, T[3*nquad-1]);
//
//    isobaric_heat_rejection(0, PL);
//    isothermal_heat_rejection(1, VL, V[nquad*1-1]);
//
//    print("otto");
//
//    std::cout << "net work from otto: " << work << std::endl;
//}

void otto(double cr, double TH)
{
    // setting compression ratio and temperature range causes contradiction sometimes. If suitable parameters are found, in this case, Carnot fails. Failure happens this way: Max volume may be exceeded or max temperature may be exceeded.
    work = 0.;

    double VH = cr * VL;

    double TI = isentropic_full_compression(0, VH);
    isochoric_heat_addition_TH(1, T[nquad-1], TH);
    isentropic_full_expansion(2, T[2*nquad-1], VH);
    isochoric_heat_rejection(3, T[3*nquad-1], VH);

    print("otto");

    std::cout << "net work from otto: " << work << std::endl;
}

void otto_TH_QR(double TH, double QR)
{
    // otto cycle with given max temperature and heat rejection.
    
    work = 0.;

    double T4 = TL + QR / m / cv;
    double VH = VL * std::pow(TH / T4, 1. / (GM-1));

    isentropic_full_expansion(0, TH, VH);
    isochoric_heat_rejection(1, T4, VH);
    double T2 = isentropic_full_compression(2, VH);
    isochoric_heat_addition_TH(3, T2, TH);

    print("otto");

    std::cout << "net work from otto: " << work << std::endl;
}

void diesel_TH_QR(double TH, double QR)
{
    // diesel cycle with given max temperature and heat rejection.
    // Just Tmax and QR seems like not enough.
    
    work = 0.;

    double T4 = TL + QR / m / cv;
    std::cout << "T4: " << T4 << std::endl;
    double V3 = VL * std::pow(TH / T4, 1. / (GM-1)); // wrong. not VL bu VH and I dont know VH!
    std::cout << "V3: " << V3 << std::endl;
    double P3 = m * R * TH / V3;
    std::cout << "P3: " << P3 << std::endl;
    double P2 = P3;
    double T2 = P2 * VL / m / R;
    std::cout << "T2: " << T2 << std::endl;

    isentropic_compression(0, TL, T2, VL);
    isobaric_heat_addition(1, P2, TH);
    double V4 = isentropic_expansion_low(2, TH, T4, V3);
    isochoric_heat_rejection(4, T4, V4);

    print("diesel");

    std::cout << "net work from diesel: " << work << std::endl;
}


int main()
{
    //stirling(6, 600);
    //carnot1(6, 600);
    //ericsson(6, 600);
    //otto_TH_QR(1000, 300);
    diesel_TH_QR(600, 10);

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
