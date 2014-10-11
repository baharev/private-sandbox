//=============================================================================
//
//  This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
//  Copyright (C) 2014  Ali Baharev
//  All rights reserved. E-mail: <my_first_name.my_last_name@gmail.com>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//=============================================================================

#include <iostream>
#include <assert.h>
#include "affine.hpp"
#include "rmatrix.hpp"
#include "ivector.hpp"
#include "matinv.hpp"

using namespace std;
using namespace cxsc;
using namespace aa_lib;

// x_new := x intersect F(x); iteration with interval Newton method and with
// affine linearization and -A^(-1)*b contraction. (LV Kolev, An improved
// interval linearization for solving nonlinear problems NUMA 2004)
//
// Newton iteration stalls from DELTA = 0.46772 but succeeds with 0.46771;
// starts to deteriorate at 0.253 (3 vs. 4 iterations).
// Affine still succeeds with DELTA = 0.8938 in 99 iterations.

const double DELTA = 0.4; // box = [x-delta, x+delta]
const double TOLER = 1.0e-3;

double solution[] = {
        -0.477446285098725,
        -0.520433482442136,
        -0.558370943009135,
        -0.592117117250195,
        -0.622301248771590,
        -0.650541249094133,
        -0.648055016625447,
        -0.645614708220959,
        -0.643573559963035,
        -0.642165324775287,
        -0.641534421218548,
        -0.641838340629904,
        -0.642048212981229,
        -0.642138128758668,
        -0.643065534204966,
        -0.614079591979722
};

// Broyden banded
// AMPL code:
// param N := 16;
// var x{1..N} >= -100, <= 101;
// con{i in 1..N}:
//   x[i]*(5*x[i]^2+2)+1-sum{j in 1..N: j>=max(1,i-5) and j<=min(N,i+1)} (x[j]^2+x[j]) = 0;
void f_aa_broyden(const affine x[], affine r[]) {
    for (int i=0; i<16; ++i) {
        affine sum(0.0);
        for (int j=0; j<16; ++j) {
            if (j>=std::max(0, i-5) && j<=std::min(15, i+1))
                sum += (aa_lib::sqr(x[j])+x[j]);
        }
        r[i] = x[i]*(mult_const(5.0, aa_lib::sqr(x[i]))+affine(2.0))-sum+affine(1.0);
    }
}

// returns x_new: x intersect -A^(-1)*b
ivector aa_iterate(ivector x) {
    const int nvars = VecLen(x);

    affine::set_max_used_index(0);

    affine* aa_var = new affine[nvars];
    for (int i=0; i<nvars; ++i) {
        aa_var[i] = affine(_double(Inf(x[i+1])), _double(Sup(x[i+1])));
    }

    affine* linearized_constraint = new affine[nvars];

    f_aa_broyden(aa_var, linearized_constraint);

//    cout << endl << "The affine linearization\n" << endl;
//    for (int i=0; i<nvars; ++i) {
//        cout << "Eq #" << i << endl;
//        cout << linearized_constraint[i] << flush;
//    }

    rmatrix A(nvars, nvars);
    for (int i=1; i<=nvars; ++i)
        for (int j=1; j<=nvars; ++j)
            A[i][j] = 0.0;

    ivector b(nvars);

    for (int i=0; i<nvars; ++i) {
        affine& r     = linearized_constraint[i];
        for (int k=0; k< r.n; ++k) {
            int j      = r.index[k];
            double val = r.value[k];
            if (j==0) {
                b[i+1] = val;
            }
            else if (j<=nvars) {
                A[i+1][j] = val;
            }
            else {
                double fval = std::fabs(val);
                b[i+1] += interval(-fval, fval);
            }
        }
    }

    cout << "x:\n" << x << endl;
    cout << "A:\n" << A << endl;
    cout << "b:\n" << b << endl;

    rmatrix R(nvars, nvars);
    int err = 0;
    MatInv(A, R, err);
    if (err) {
        cerr << "Inversion error" << endl;
        exit(1);
    }

    ivector eps_new = R*(-b);
    cout << "eps_new:\n" << eps_new << endl;

    ivector xnew(nvars);

    for (int i=0; i<nvars; ++i) {
        double* a = aa_var[i].value;
        // HACK: Leave too narrow intervals intact
        interval eps =(a[1]<1.0e-5) ? interval(-1,1) : (eps_new[i+1] & interval(-1,1));
        xnew[i+1] = interval(a[0] + a[1]*eps);
    }

    delete[] linearized_constraint;
    delete[] aa_var;

    return xnew;
}

//==============================================================================
#include <nlinsys.hpp>

GTvector f_broyden(const GTvector& x) {
    int n = x.Dim();
    GTvector Result(n);
    for (int i=1; i<=n; ++i) {
        GradType sum(n);
        sum = 0.0;
        for (int j=1; j<=n; ++j) {
            if (j>=std::max(1, i-5) && j<=std::min(n, i+1))
                sum = sum + (sqr(x[j])+x[j]);
        }
        Result[i] = x[i]*(5.0*sqr(x[i])+2.0)-sum+1.0;
    }
    return Result;
}

// As taken from C-XSC but dropped the extended Newton part
static ivector NewtonStep(GTvector_FctPtr f,
                          ivector         Y,
                          imatrix&        JfY,
                          real            Epsilon)
{
    const int     n = Ub(Y);
    rvector       c(n);
    ivector       fC(n), b(n), Y_minus_c(n);
    rmatrix       R(n,n);
    imatrix       A(n,n);
    int           i, i0, j, InvErr;
    idotprecision Accu;
    interval      h;

    c = mid(Y);  fEvalJ(f,_ivector(c),fC);         // Midpoint evaluation of 'f'
    //---------------------------
    MatInv(mid(JfY),R,InvErr);                     // Invert the midpoint matrix
    if (InvErr) R = Id(R);                         //---------------------------

    A = R * JfY;  b = R * fC;         // Compute data for Gauss-Seidel step
    Y_minus_c = Y - c;                //-----------------------------------

    i0 = 0;                           // Initializations, A[i0,i0] contains zero
    //----------------------------------------
    for (i = 1; i <= n; i++ ) {                // Interval-Gauss-Seidel-step for
        // non-zero A[i,i] elements
        //-------------------------------
        if ( in(0.0,A[i][i]) || (RelDiam(Y[i]) <= Epsilon) ) {
            i0 = i;       // Largest i with 0 in A[i,i]
            continue;     // Next i
        }

        Accu = b[i];
        for (j = 1; j < i; j++) accumulate(Accu,A[i][j],Y_minus_c[j]);
        for (j = i+1; j <= n; j++) accumulate(Accu,A[i][j],Y_minus_c[j]);
        rnd(Accu,h);
        h = c[i] - h / A[i][i];

        if ( Disjoint(Y[i],h) ) { exit(0); }                // No solution possible
        //---------------------
        Y[i] = Y[i] & h;  Y_minus_c[i] = Y[i] - c[i];
    } // for (i...

    return Y;
}

static ivector XINewton(GTvector_FctPtr f,
                        ivector         y,
                        real&           Epsilon)
{
    const int  n = Ub(y);
    ivector    fy(n);
    imatrix    Jfy(n,n);

    fJEvalJ(f, y, fy, Jfy);               // Compute f(y) and Jf(y)
    if ( !in(0,fy) )
        exit(0);
    //------------------------------------
    return NewtonStep(f,y,Jfy,Epsilon);
}

ivector newton_iterate(ivector x) {
    real Epsilon = TOLER;
    return XINewton(f_broyden, x, Epsilon);
}

//==============================================================================

void linearization(ivector (*iterate)(ivector)) {

    const int nvars = 16;
    ivector x(nvars);

    for (int i=0; i<nvars; ++i) {
        // Pseudo random perturbation for debugging purposes
        //double lb = solution[i]-DELTA*(    i%2? (0.45+i/20.0): (1.0-i/20.0));
        //double ub = solution[i]+DELTA*((i+1)%2? (0.30+i/22.0): (1.0-i/22.0));
        double lb = solution[i]-DELTA;
        double ub = solution[i]+DELTA;
        x[i+1] = interval(lb, ub);
    }

    real maxreldiam = MaxRelDiam(x);
    int counter = 0;

    while (maxreldiam >= TOLER) {
        cout << "max rel. diam. before iteration: " << maxreldiam << endl;
        x = iterate(x);
        maxreldiam = MaxRelDiam(x);
        cout << "x new:\n" << x << endl;
        ++counter;
    }

    real maxdev = 0;

    for (int i=0; i<nvars; ++i) {
        maxdev = Max(Mid(x[i+1])-solution[i], maxdev);
    }

    cout << "max rel. diam. at exit: " << maxreldiam << endl;
    cout << "max deviation from sol: " << maxdev << endl;
    cout << "No. iter: " << counter << endl;
}

void newton_linearization() {
    linearization(newton_iterate);
}

void affine_linearization() {
    linearization(aa_iterate);
}

int main(int argc, char* argv[]) {
    cout << "Delta = " << DELTA << endl;
    cout << "=========================================================="<< endl;
    cout << "Affine" << endl;
    affine_linearization();
    cout << "=========================================================="<< endl;
    cout << "Interval Newton" << endl;
    newton_linearization();
}
