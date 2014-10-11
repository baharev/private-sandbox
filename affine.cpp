//=============================================================================
//
//  This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
//  Copyright (C) 2010-2014  Ali Baharev
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

#include "affine.hpp"
#include <cmath>
#include <iostream>
#include <cstdlib>
//#include <iomanip>
#include <limits>
#include <assert.h>

//=============================================================================
//
// If you wish to understand the logic of the code below then you will need to
// have a deep understanding of the mixed affine arithmetic / interval
// arithmetic model. The brilliant monograph of Stolfi and Figueiredo entitled
// Self-Validated Numerical Methods and Applications (1997) will give you the
// necessary information.
//
//=============================================================================

namespace aa_lib {

bool affine::isValid(true);

int affine::max_used_index(0);

//-----------------------------------------------------------------------------
//
// An excellent guideline for implementing numerical types is Item 20 in
// Sutter: Exceptional C++. Some of the member functions are now written
// following his advice; e.g. print below, and print is called by operator<<.
//
//-----------------------------------------------------------------------------
//
// One may wish to wrap up the index value pairs in a struct (or in an 
// std::pair from the STL). Well, the previous implementation used 
// std::map<int, cxsc::real> and I did not find any benefit in having a 
// pair<int, real> compared to the current approach (now index and value are 
// treated separately). Please feel free to differ.
//
//-----------------------------------------------------------------------------

std::ostream& affine::print(std::ostream& os) const {

	using std::endl;

  if (!isValid)
	os << "empty_interval" << endl << endl;
  else {
	for (int i=0; i<n; ++i)
		os << index[i] << ": " << value[i] << endl;

	os << "IA: " << "[ " << lb << ", " << ub << "]" << endl << endl;
  }

  return os;
}

std::ostream& operator<<(std::ostream& os, const affine & a) { return a.print(os); }

void error(const char* msg) {
	std::cout << "Error: " << msg << "!" << std::endl;
	std::exit(EXIT_FAILURE);
}

affine::affine(const double val) :
	index(new int[1]), value(new double[1]), lb(val), ub(val), n(1), nmax(1)
{
	index[0] = 0;
	value[0] = val;
}

affine::affine(const double inf, const double sup) : 
	index(new int[2]), value(new double[2]), lb(inf), ub(sup), n(2), nmax(2)
{
	if (lb >= ub)
		error("lb >= ub in explicit constructor");

	index[0] = 0;
	index[1] = ++max_used_index;
	value[0] = (inf+sup)/2.0;
	value[1] = (sup-inf)/2.0;
}

affine::affine(const affine& aa) : index(new int[aa.nmax]), value(new double[aa.nmax]),
								   lb(aa.lb), ub(aa.ub), n(aa.n), nmax(aa.nmax)
{
	int*    const idx = index;
	double* const val = value;

	const int*    const i = aa.index;
	const double* const d = aa.value;

	const int n_ = aa.n;

    for (int j=0; j<n_; ++j) {
        idx[j] = i[j];
		val[j] = d[j];
	}
}

//-----------------------------------------------------------------------------
// An elegant way to implement operator= is presented at Item 13 in Sutter:
// Exceptional C++ (create-a-temporary-and-swap idiom). Self-assignment is also
// of concern, see Item 11 in Meyers: Effective C++, Third Edition. These
// guidelines would be followed on a code rewrite. Now, the code is ugly but
// seems to be bug free.
affine& affine::operator=(const affine& rhs)
{
	const int n_ = rhs.n;

	if (n_ > nmax) {

		delete[] index;
		delete[] value;

		index = new int[n_];
		value = new double[n_];
	}

	lb = rhs.lb;
	ub = rhs.ub;
	n = n_;
	nmax = n_;

	int*    const idx = index;
	double* const val = value;

	const int*    const i = rhs.index;
	const double* const d = rhs.value;

    for (int j=0; j<n_; ++j) {
        idx[j] = i[j];
		val[j] = d[j];
	}

	return *this;
}

// Probing needs this
void affine::set_bounds(const double l, const double u) {

	assert(n == 2);
	assert(nmax >= n);

	if (!isValid)
		error("set_bounds() called when state is invalid");

    if (l >= u )
		error("empty or degenerate interval in set_bounds()");

	lb = l;

	ub = u;

	value[0] = (ub + lb)/2.0;
	value[1] = (ub - lb)/2.0;
}

bool affine::set_range(const double l, const double u) {

	assert(n > 2);
	assert(nmax >= n);
	assert(lb <= ub);
	assert(index[0] == 0);

	if (!isValid)
		error("set_range() called when state is invalid");

    if (l >= u )
		error("empty or degenerate interval in set_range()");

    if (l > lb)
    	lb = l;

    if (u < ub)
    	ub = u;

    // TODO Perhaps propagate?

    isValid = (lb <= ub);

	return (isValid);
}

bool affine::intersect_domain(const double l, const double u) {

	assert(n == 2);
	assert(nmax >= n);
	assert(lb <= ub);
	assert(index[0] == 0);

	if (!isValid)
		error("set_range() called when state is invalid");

    if (l >= u )
		error("empty or degenerate interval in set_range()");

    bool has_changed = false;

    if (l > lb) {
    	lb = l;
    	has_changed = true;
    }

    if (u < ub) {
    	ub = u;
    	has_changed = true;
    }

    isValid = (lb <= ub);

    if (!has_changed || !isValid)
    	return isValid;

	value[0] = (ub + lb)/2.0;
	value[1] = (ub - lb)/2.0;

	assert(isValid);

	return isValid;
}

//=============================================================================
//
// "Duplication is Evil"
//
//=============================================================================
//
// I am not satisfied with the code mainly because approximately 25% of the
// code can be eliminated by eliminating the repeating code fragments, more on
// this in a second. Another cause for my dissatisfaction is the weird looking
// and unnatural approach found at many parts. Reasons for the latter: it was a
// research code, I did not know exactly what and how I should implement. The
// code went through several find-and-replace "reviews" as a result (replacing
// std::map<int, cxsc::real> with built-in types, then replacing the memory
// pool, etc). Last but not least the core of the code was my first C++
// project.
//
// Here is how the repeating code fragments can (and should) be eliminated by
// applying the strategy pattern (Gamma, Helm, Johnson, Vlissides: Design
// Patterns). The major part of the code of operator+, operator-, operator* and
// operator/ (addition, subtraction, multiplication, division) is the same but
// there is some variation. The common part should go into a function and the 
// varying parts would be calls to a function object. The function object
// (functor; std::binary_function from the STL) should perform the concrete
// operation.
//
// A very similar approach is followed by the std::sort algorithm, and
// interestingly enough there is no abstraction penalty (Item 46 in Meyers:
// Effective STL): the corresponding function call can be inlined.
//
// The unary operations are even more easier to fix: just a const double has to
// be passed to the function containing the common code fragment.
//
//-----------------------------------------------------------------------------
//
// Usually, operator+ should be implemented in terms of operator+=, see Item 20
// in Sutter: Exceptional C++. But think about it for a moment in case of the 
// affine data type. Roughly speaking, the affine data type is a sparse
// vector for which += is tricky as opposed to complex numbers for example (as
// in Item 20 referenced above). Actually, here the situation is reversed:
// operator+= is implemented in terms of operator+.
//
//-----------------------------------------------------------------------------

const affine operator+(const affine& lhs, const affine& rhs) {

	if (!affine::isValid)
		return affine(false);

	assert((lhs.n > 0) && (rhs.n > 0));
	assert((lhs.lb <= lhs.ub)&&(rhs.lb <= rhs.ub));
	assert((lhs.index[0] == 0) && (rhs.index[0] == 0));
	assert(lhs.nmax >= lhs.n);
	assert(rhs.nmax >= rhs.n);

	const int size = lhs.nmax + rhs.nmax - 1;

	affine result(false, size);

	// If the code below looks very weird to you then you are right!
	// This code was originally written with std::map iterators
	// and went through several find-and-replace "reviews" later when
	// I replaced the map container with built-in types.
	// This issue will be resolved on the next rewrite as discussed above.

	int*    ri = result.index;
	double* rd = result.value;

	const int* p = lhs.index;
	const int* q = rhs.index;

	const int* const p_end = lhs.index+lhs.n;
	const int* const q_end = rhs.index+rhs.n;

	const double* r = lhs.value;
	const double* s = rhs.value;

	const int int_max(std::numeric_limits<int>::max());

	int k = 0;

	double mid = 0.0;
	double rad = 0.0;

	do {
		assert (k < size);

		int i, j, idx;

		if (p != p_end)
			i = *p;
		else
			i = int_max;

		if (q != q_end)
			j = *q;
		else
			j = int_max;

		double x, y;

		if      (i < j) {
			idx = i;
			x = *r;
			y = 0.0;
			++p;
			++r;
		}
		else if (i > j) {
			idx = j;
			x = 0.0;
			y = *s;
			++q;
			++s;
		}
		else {
			idx = i;
			x = *r;
			y = *s;
			++p;
			++r;
			++q;
			++s;
		}

		ri[k] = idx;
		double tmp = x+y;    
		rd[k] = tmp;

		if (k != 0)
			rad += fabs(tmp);	 
		else
			mid = tmp;

		++k;
	}
	while (!((p == p_end) && (q == q_end)));

	double lb = std::max( mid - rad, lhs.lb + rhs.lb);
	double ub = std::min( mid + rad, lhs.ub + rhs.ub);

	if (lb > ub)
		affine::isValid = false;

	result.lb = lb;
	result.ub = ub;

	result.n = k;

	return result;

}

const affine operator-(const affine& aa) {

	if (!affine::isValid)
		return affine(false);

	assert((aa.lb <= aa.ub));
    assert((aa.index[0] == 0));
	assert(aa.nmax >= aa.n);

	const int size = aa.nmax;
	affine result(false, size);

	result.lb = -aa.ub;
	result.ub = -aa.lb;
	const int n_ = aa.n;

	result.n = n_;

	int*    const idx = result.index;
	double* const val = result.value;

	const int*    const i = aa.index;
	const double* const d = aa.value;

    for (int j=0; j<n_; ++j) {
        idx[j] =  i[j];
		val[j] = -d[j];
	}

	return result;
}

const affine operator-(const affine& lhs, const affine& rhs) {

	if (!affine::isValid)
		return affine(false);

	assert((lhs.n > 0) && (rhs.n > 0));
	assert((lhs.lb <= lhs.ub)&&(rhs.lb <= rhs.ub));
	assert((lhs.index[0] == 0) && (rhs.index[0] == 0));
	assert(lhs.nmax >= lhs.n);
	assert(rhs.nmax >= rhs.n);

	const int size = lhs.nmax + rhs.nmax - 1;

	affine result(false, size);

	int*    ri = result.index;
	double* rd = result.value;

	const int* p = lhs.index;
	const int* q = rhs.index;

	const int* const p_end = lhs.index+lhs.n;
	const int* const q_end = rhs.index+rhs.n;

	const double* r = lhs.value;
	const double* s = rhs.value;

	const int int_max(std::numeric_limits<int>::max());

	int k = 0;

	double mid = 0.0;
	double rad = 0.0;

	do {
		assert (k < size);

		int i, j, idx;

		if (p != p_end)
			i = *p;
		else
			i = int_max;

		if (q != q_end)
			j = *q;
		else
			j = int_max;

		double x, y;

		if      (i < j) {
			idx = i;
			x = *r;
			y = 0.0;
			++p;
			++r;
		}
		else if (i > j) {
			idx = j;
			x = 0.0;
			y = *s;
			++q;
			++s;
		}
		else {
			idx = i;
			x = *r;
			y = *s;
			++p;
			++r;
			++q;
			++s;
		}

		ri[k] = idx;
		double tmp = x-y;    
		rd[k] = tmp;

		if (k != 0)
			rad += fabs(tmp);	 
		else
			mid = tmp;

		++k;
	}
	while (!((p == p_end) && (q == q_end)));

	double lb = std::max( mid - rad, lhs.lb - rhs.ub);
	double ub = std::min( mid + rad, lhs.ub - rhs.lb);

	if (lb > ub)
		affine::isValid = false;

	result.lb = lb;
	result.ub = ub;

	result.n = k;

	return result;
}

const affine mult_const(const double c, const affine& x) {

	if (!affine::isValid)
		return affine(false);

	if (x.n==1)
		return affine(c*(x.lb));

	assert(x.lb <= x.ub);
	assert(x.index[0] == 0);
	assert(x.nmax >= x.n);

	const int size = x.nmax;

	affine result(false, size);

	double lb = c*(x.lb);
	double ub = c*(x.ub);

	if (lb > ub) {
		const double temp = lb;
		lb = ub;
		ub = temp;
	}

	result.lb = lb;
	result.ub = ub;

	const int n_ = x.n;
	result.n  = n_;

	int*    const idx = result.index;
	double* const val = result.value;

	const int*    const i = x.index;
	const double* const d = x.value;

    for (int j=0; j<n_; ++j) {
        idx[j] =    i[j] ;
		val[j] = c*(d[j]);
	}

	return result;
}

// Based on the C-XSC source code
void mult(const double xl, const double xu, const double yl, const double yu, double& zl, double& zu) {

  if (xl >=0) {                        /*  0 <= [x]                 */

    if (yl >=0)                        /*  0 <= [y]                 */
      zl=xl*yl;
	else                               /*  [y] <= 0  or  0 \in [y]  */
      zl=xu*yl;

    if (yu <=0)                        /*  [y] <= 0                 */
      zu=xl*yu;  
	else                               /*  0 <= [y]  or  0 \in [y]  */
      zu=xu*yu;

  } else if (xu<=0) {                  /*  [x] <= 0                 */

    if (yu<=0)                         /*  [y] <= 0                 */
      zl=xu*yu;
    else                               /*  0 <= [y]  or  0 \in [y]  */
      zl=xl*yu; 

    if (yl>=0)                         /*  0 <= [y]                 */
      zu=xu*yl;
    else                               /*  [y] <= 0  or  0 \in [y]  */
      zu=xl*yl;

  } else {                             /*  0 \in [x]                */

    if (yl>=0) {                       /*  0 <= [y]                 */
      zl=xl*yu;
      zu=xu*yu;
    } else if (yu<=0) {                /*  [y] <= 0                 */
      zl=xu*yl;
      zu=xl*yl;
    } else {                           /*  0 \in [x], 0 \in [y]     */
		zl=std::min(xl*yu, xu*yl);
		zu=std::max(xl*yl, xu*yu);
    }

  }

}

const affine operator*(const affine& lhs, const affine& rhs) {

	if (!affine::isValid)
		return affine(false);

	assert((lhs.n > 0) && (rhs.n > 0));
	assert((lhs.lb <= lhs.ub)&&(rhs.lb <= rhs.ub));
	assert((lhs.index[0] == 0) && (rhs.index[0] == 0));
	assert(lhs.nmax >= lhs.n);
	assert(rhs.nmax >= rhs.n);

	// The hideous tagged class approach, ouch!
	// On the next rewrite, replace this 
	// with appropriate operator* overloading
	if      (lhs.n == 1)
		return mult_const(lhs.lb, rhs);
	else if (rhs.n == 1)
		return mult_const(rhs.lb, lhs);

	const int size = lhs.nmax + rhs.nmax;
	affine result(false, size);

  int*    const ri = result.index;
  double* const rd = result.value;

  const int* p = lhs.index+1;
  const int* q = rhs.index+1;

  const int* const p_end = lhs.index+lhs.n;
  const int* const q_end = rhs.index+rhs.n;

  const double* r = lhs.value;
  const double* s = rhs.value;

  const int int_max(std::numeric_limits<int>::max());

  //
  //------------------------------

  double x0 = *r++;
  double y0 = *s++;

  double rad_x, rad_y, c, rad;
  rad_x = rad_y = c = rad = 0.0;

  int k = 1;

  while (!((p == p_end) && (q == q_end))) {

	assert (k < size);

    int i, j, idx;

    if (p != p_end)
      i = *p;
    else
      i = int_max;

    if (q != q_end)
      j = *q;
    else
      j = int_max;

    double x_i, y_i;

    if      (i < j) {
      idx = i;
      x_i = *r;
      y_i = 0.0;
      ++p;
      ++r;
    }
    else if (i > j) {
      idx = j;
      x_i = 0.0;
      y_i = *s;
      ++q;
      ++s;
    }
    else {
      idx = i;
      x_i = *r;
      y_i = *s;
      ++p;
      ++r;
      ++q;
      ++s;
    }

	rad_x += fabs(x_i);
	rad_y += fabs(y_i);

	c += x_i*y_i;

    ri[k] = idx;

	double tmp = x0*y_i+y0*x_i;    
	rd[k] = tmp;
    rad += fabs(tmp);

    ++k;
  }

  assert (k < size);

  c /= 2.0;

  double mid = x0*y0 + c;
  result.value[0] = mid;

  double delta = rad_x*rad_y;
  delta -= fabs(c);

  assert(delta >= 0.0);
  // ??? if (delta < 0.0) ???
  if (delta < 0.0)
	  delta = 0.0;

  rad += fabs(delta);

  double lb, ub;

  mult(lhs.lb, lhs.ub, rhs.lb, rhs.ub, lb, ub);

  //cout << "mid+/-rad:  [" << mid-rad << ", " << mid+rad << "]" << endl;
  //cout << "X*Y      :  [" << lb      << ", " << ub      << "]" << endl;

  lb = std::max( lb, mid - rad);
  ub = std::min( ub, mid + rad);

  if (lb > ub)
	  affine::isValid = false;

  result.lb = lb;
  result.ub = ub;

  result.index[0] = 0;

  result.index[k] = ++affine::max_used_index;
  result.value[k] = delta;

  ++k;

  result.n = k;

  return result;
}

//----------------------------------------------------------------------
//  Reciprocal function using minimax / Chebysev approximation
//  The implementation is based on the paper of L. V. Kolev
//  "An improved interval linearization for solving non-linear problems"
//  Numerical Algorithms, 37, pp.213-224 (2004)

const affine reciprocal(const affine& y) {
	
	if (!affine::isValid)
		return affine(false);

	if (y.n==1) {

		return affine(1/(y.lb));
	}

	assert(y.lb <= y.ub);
	assert(y.index[0] == 0);
	assert(y.nmax >= y.n);

	const int size = y.nmax + 1;
	affine result(false, size);

	const double y_inf = y.lb;

	const double y_sup = y.ub;

	if ( (y_inf <= 0.0) && (0.0 <= y_sup) )
		error("division by zero");

	const double s  = -1.0/(y_inf*y_sup);

	const double y1 = -std::sqrt(-1.0/s);

	const double y2 = -y1;

	const double ys = (y_inf>0.0)?y2:y1;

	result.lb = 1.0/y_sup;
	result.ub = 1.0/y_inf;

	const double f_inf = 1.0/ys    - s*ys   ;
	const double f_sup = 1.0/y_inf - s*y_inf;

	const double f0 = (f_inf + f_sup) / 2.0;

	const double rf = f_sup-f0;

	const int*    const i = y.index;
	const double* const d = y.value;

	int*    const ri = result.index;
	double* const rd = result.value;

	const int n = y.n;

	for (int j=0; j<n; ++j) {
		ri[j] =    i[j] ;
		rd[j] = s*(d[j]);
	}

	rd[0] += f0;

	ri[n] = ++affine::max_used_index;
	rd[n] = rf;

	result.n = n+1;

	return result;
}

// Based on the C-XSC source code
void division(const double xl, const double xu, const double yl, const double yu, double& zl, double& zu) {

  if (yl>0) { 

    if (xl>=0)
      zl=xl/yu;
	else 
      zl=xl/yl;

    if (xu<=0)
      zu=xu/yu;
	else 
      zu=xu/yl;

  } else if (yu<0) {

    if (xu<=0)
      zl=xu/yl;
	else 
      zl=xu/yu;

    if (xl>=0)     
      zu=xl/yl;
	else 
      zu=xl/yu;
  } else {
    error("division by zero in IA");
  }

}

//  FIXME Should use operator*(double, affine) for affine P(X - c*Y) and
//  eliminate loop
//  FIXME Handle trivial cases -- constants
const affine operator/(const affine& lhs, const affine& rhs) {

	if (!affine::isValid)
		return affine(false);

	assert((lhs.n > 0) && (rhs.n > 0));
	assert((lhs.lb <= lhs.ub)&&(rhs.lb <= rhs.ub));
	assert((lhs.index[0] == 0) && (rhs.index[0] == 0));
	assert(lhs.nmax >= lhs.n);
	assert(rhs.nmax >= rhs.n);

	const int start_index = affine::max_used_index;

	const int size = rhs.nmax + lhs.nmax - 1;

	affine P(false, size);

	affine Q(reciprocal(rhs));

	//cout << "1/y: " << endl << Q << endl;

	int*    pi = P.index;
	double* pd = P.value;

	const int*    xi = lhs.index;
	const double* xd = lhs.value;

	const int*    yi = rhs.index;
	const double* yd = rhs.value;

	const int* const x_end = xi+lhs.n;
	const int* const y_end = yi+rhs.n;

	const int int_max(std::numeric_limits<int>::max());

	// TODO c = +/- INF or NaN ???
	double c = (*xd)/(*yd);

	++xi; ++xd;
	++yi; ++yd;


	*pi++ =   0;
	*pd++ = 0.0;

	int k = 1;

	double rad_p = 0.0;

	while (!((xi == x_end) && (yi == y_end))) {

		assert( k < size);

		int i, j, index;

		if (xi != x_end)
			i = *xi;
		else
			i = int_max;

		if (yi != y_end)
			j = *yi;
		else
			j = int_max;

		double x_i, y_i;

		if      (i < j) {
			index = i;
			x_i = *xd;
			y_i = 0.0;
			++xi; ++xd;
		}
		else if (i > j) {
			index = j;
			x_i = 0.0;
			y_i = *yd;
			++yi; ++yd;
		}
		else {
			index = i;
			x_i = *xd;
			y_i = *yd;
			++xi; ++xd;
			++yi; ++yd;
		}

		double p_i = x_i-c*y_i;

		*pi++ = index;
		*pd++ = p_i;

		rad_p += fabs(p_i);

		++k;
	}

	assert( k <= size);

	P.n = k;


	double cY_l = c*rhs.lb;
	double cY_u = c*rhs.ub;

	if (cY_l > cY_u) std::swap(cY_l, cY_u);

	// TODO Should check at other parts ???
	double lb_p = std::max(-rad_p, lhs.lb - cY_u);
	double ub_p = std::min( rad_p, lhs.ub - cY_l);

	if (lb_p > ub_p) {
		affine::isValid = false;
		return affine(false);
	}

	P.lb = lb_p;
	P.ub = ub_p;

	affine V(P*Q);

	*(V.value) += c;

	V.lb += c;

	V.ub += c;

	double lb, ub;

	division(lhs.lb, lhs.ub, rhs.lb, rhs.ub, lb, ub);

	if (V.lb < lb)
		V.lb = lb;

	if (V.ub > ub)
		V.ub = ub;

	if (V.lb > V.ub)
		affine::isValid = false;

	const int new_indices = affine::max_used_index - start_index;

	if (new_indices == 2) { // Condense them!

		--affine::max_used_index;

		int n = V.n;

		V.index[n-2] = affine::max_used_index;
		const double before_last = V.value[n-2];
		assert(fabs(before_last) < 1.0e-12);
		V.value[n-2] = before_last + V.value[n-1];

		V.n -= 1;
	}

	return V;
}

const affine ln(const affine& x) {

	if (!affine::isValid)
		return affine(false);

	assert(x.n > 1);
	assert(x.lb <= x.ub);
	assert(x.index[0] == 0);
	assert(x.nmax >= x.n);

	const int size = x.nmax + 1;

	affine result(false, size);

	double a = x.lb;

	if (a <= 0.0)
		error("negative argument in ln()");

	const double b = x.ub;

	const double log_b = log(b);

	const double alpha =  log(a/b)/(a-b);

	const double fu    = -log(alpha);

	const double ru    = -b*alpha+log_b+1.0;

	const double zeta  =  (fu+ru)/2.0-1.0;

	const double delta =  (fu-ru)/2.0;

	const int*    const xi = x.index;
	const double* const xd = x.value;

	int*    const yi = result.index;
	double* const yd = result.value;

	const int n_ = x.n;
	
	for (int j=0; j<n_; ++j) {

		yi[j] =        xi[j] ;
		yd[j] = alpha*(xd[j]);	
	}

	yi[n_] = ++affine::max_used_index;
	
	yd[n_] = delta;

	yd[0] += zeta;

	result.n = x.n + 1;

	result.lb = log(a);

	result.ub = log_b;

	return result;
}

const affine exp(const affine& x) {

	if (!affine::isValid)
		return affine(false);

	assert(x.n > 1);
	assert(x.lb <= x.ub);
	assert(x.index[0] == 0);
	assert(x.nmax >= x.n);

	const double a = x.lb;

	const double b = x.ub;

	const double e_a = std::exp(a);

	// FIXME Find a better way, deal with degenerate cases (b-a) being tiny
	if (a==b)
		return affine(e_a);

	const double e_b = std::exp(b);

	const double alpha = (e_b-e_a)/(b-a);

	const double ln_alpha = log(alpha);

	const double zeta  = (-alpha*ln_alpha+e_b-alpha*b+alpha)/2.0;

	const double delta = ( alpha*ln_alpha+e_b-alpha*b-alpha)/2.0;

	const int*    const xi = x.index;
	const double* const xd = x.value;

	const int size = x.nmax + 1;

	affine result(false, size);

	int*    const yi = result.index;
	double* const yd = result.value;

	const int n = x.n;

	for (int j=0; j<n; ++j) {

		yi[j] =        xi[j] ;
		yd[j] = alpha*(xd[j]);
	}
	
	yi[n] = ++affine::max_used_index;
	
	yd[n] = delta;

	yd[0] += zeta;

	result.n = n + 1;

	result.lb = e_a;

	result.ub = e_b;

	return result;
}

const affine sqr(const affine& x) {
	
	if (!affine::isValid)
		return affine(false);

	assert(x.n > 1);
	assert(x.lb <= x.ub);
	assert(x.index[0] == 0);
	assert(x.nmax >= x.n);

	const int size = x.nmax + 1;

	affine result(false, size);

	const double a = x.lb;

	const double b = x.ub;

	const double a_2_p_b_2 = a*a+b*b;

	const double a_b = a*b;

	const double alpha = a + b;

	const double zeta  =  -(a_2_p_b_2+6.0*a_b)/8.0;

	const double delta =   (a_2_p_b_2-2.0*a_b)/8.0;

	const int*    const xi = x.index;
	const double* const xd = x.value;

	int*    const yi = result.index;
	double* const yd = result.value;

	const int n = x.n;
	
	for (int j=0; j<n; ++j) {

		yi[j] =        xi[j] ;
		yd[j] = alpha*(xd[j]);
	}
	
	yi[n] = ++affine::max_used_index;
	
	yd[n] = delta;

	yd[0] += zeta;

	result.n = n + 1;

	//----------------------------

	double c = fabs(a);

	double d = fabs(b);

	if (c > d) {
		double temp = c;
		c = d;
		d = temp;
	}

	double lb;

	if ((a <= 0.0) && (0.0 <= b))
		lb = 0.0;
	else
		lb = c*c;

	result.lb = lb;

	result.ub = d*d;

	return result;
}

// void bisect(affine x[], affine y[], const int index, const int n) {
//   
//   const int k = index - 1;
//   
//   for (int i=0; i<n; ++i)    
//     y[i] = x[i];
// 
//   assert(isValid);
//   assert((1<=index)&&(index<=N_VARS));
//   assert(x[k].n == 2);
//   assert(x[k].index[0] == 0);
//   assert(x[k].index[1] == index);
//   
//   const double lb = x[k].lb;
// 
//   const double ub = x[k].ub;
//   
//   const double midpoint_l = ((3.0*lb+ub)/4.0);
// 
//   const double midpoint_u = ((3.0*ub+lb)/4.0);
// 
//   const double midpoint   = (lb+ub)/2.0;
//   
//   const double rad = (ub-lb)/4.0;
// 
//   x[k].value[0] = midpoint_l;
//   
//   x[k].value[1] = rad;
// 
//   x[k].lb = lb;
// 
//   x[k].ub = midpoint;
// 
//   y[k].value[0] = midpoint_u;
// 
//   y[k].value[1] = rad;
// 
//   y[k].lb = midpoint;
// 
//   y[k].ub = ub;  
//   
// }

void bisect(affine x[], affine y[], const int index, const int n) {

  assert(affine::isValid);

  const int k = index - 1;

  for (int i=0; i<n; ++i)    
    y[i] = x[i];

  const double lb = x[k].lb;

  const double ub = x[k].ub;

  double midpoint_l, midpoint_u, midpoint, rad_l, rad_u;

//--------------------------------------------

  // !!! <= changed to <
  if ( lb <= 0.0 && 0.0 <= ub) {

    midpoint = (fabs(lb) > fabs(ub)) ? (lb/10.0) : (ub/10.0);

    midpoint_l = (lb + midpoint)/2.0;

    midpoint_u = (midpoint + ub)/2.0;

    rad_l = (midpoint - lb)/2.0;

    rad_u = (ub - midpoint)/2.0;

  }
  else {


    midpoint_l = ((3.0*lb+ub)/4.0);

    midpoint_u = ((3.0*ub+lb)/4.0);

    midpoint   = (lb+ub)/2.0;

    rad_l = rad_u = (ub-lb)/4.0;
  }

//---------------------------------------------

  x[k].value[0] = midpoint_l;
  
  x[k].value[1] = rad_l;

  x[k].lb = lb;

  x[k].ub = midpoint;

  y[k].value[0] = midpoint_u;

  y[k].value[1] = rad_u;

  y[k].lb = midpoint;

  y[k].ub = ub;  
  
}

bool is_narrow(const double l, const double u) {

	assert( l < u );

	bool isNarrow = false;

	const double diam = u - l;
	
	const double fl = fabs(l);
	const double fu = fabs(u);

	const double mv = (fl > fu) ? fl : fu;

	const double rdiam = ( mv==0.0 ) ? 0.0 : (diam/mv) ;

	if (rdiam < tol_solved_lp || diam < tol_solved_lp) {

		isNarrow = true;
	}

	return isNarrow;
}

// DAG needs these
bool operator!=(const double val,  const affine& aa) {

	return !(aa.inf()==val && aa.sup()==val && aa.n_elem() == 1);
}

// Sparse vector! a += b ==> a = a + b
// See also the comment block above operator+
affine& affine::operator+=(const affine& aa){
	
	(*this) = (*this) + aa;
	return (*this);
}

bool in(const double val, const affine& x) {
	return (affine::is_valid() && ((x.inf() <= val) && (val <= x.sup())));
}

//-----------------------------------------------------------------------------
// Dummy functions
bool Disjoint(const affine &, const affine &) {
	error("affine forms are not available for propagation");
	return false;
}

const affine operator&(const affine &, const affine &) {
	error("affine forms are not available for propagation");
	return affine();
}

bool ext_div( affine&, const affine&, const affine& ) {
	error("affine forms are not available for propagation");	
	return false;
}

const affine sqrt(const affine& ) {
	error("affine forms are not available for propagation");	
	return affine();
}

double Inf(const affine& ) {
	error("affine forms are not available for propagation");
	return 0.0;
}

double Sup(const affine& ) {
	error("affine forms are not available for propagation");
	return 0.0;
}

}
