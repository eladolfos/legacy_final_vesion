#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <cmath>

extern "C" {
    void lhapdfstart_(char* name, int* set, int length);
    double alphaspdf6_(const double &Q);
    void orderas6_(int& order);
    double xfx6_(const int &id, const double &x1, const double &Q);
    void nfl6_(const double &Q, int &nf);
    double q0pdf_();
}

#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
LHAPDF::PDF* pdf;
#endif

void lhapdfstart_(char* name, int* set, int length){
    int Iset = *set;
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
    pdf = LHAPDF::mkPDF(name, Iset);
#else
    LHAPDF::initPDFSet(name, LHAPDF::LHGRID, Iset);
#endif
}

double alphaspdf6_(const double &Q){
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
    double alphas = pdf -> alphasQ2(Q*Q);
#else
    double alphas = LHAPDF::alphasPDF(Q);
#endif
    return alphas;
}

void orderas6_(int &order) {
    order = pdf -> orderQCD();
}

double xfx6_(const int &id, const double &x1, const double &Q){
    int pid = id;
    int sign = pid/std::abs(pid);
    if(std::abs(pid) == 1) pid = sign*2;
    else if(std::abs(pid) == 2) pid = sign*1;
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
    return pdf -> xfxQ2(pid, x1, Q*Q)/x1;
#else
    if(id == 21)
    {
        return LHAPDF::xfx(x1,Q,0);
    }
    else
    {
        return LHAPDF::xfx(x1,Q,pid);
    }
#endif
}

void nfl6_(const double &Q, int &neff) {
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
    neff = pdf -> alphaS().numFlavorsQ(Q);
    if(neff > 5) neff = 5;
#else
    static_assert(false, "Need LHAPDF6");
#endif
}

double q0pdf_() {
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
    return pdf -> qMin();
#else
    return sqrt(LHAPDF::getQ2min(0));
#endif
}
