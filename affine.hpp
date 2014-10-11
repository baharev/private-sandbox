//=============================================================================
//
//  This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
//  Copyright (C) 2010  Ali Baharev
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

#ifndef __AFFINE_HPP
#define __AFFINE_HPP

#include <iosfwd>

namespace aa_lib {

class affine;

const affine sqr(const affine& x);

class affine {

	public:

		affine() : index(new int[2]), value(new double[2]), n(0), nmax(2) { }

		explicit affine(const double val);

		affine(const double inf, const double sup);

		affine(const affine& other);

		affine& operator=(const affine& rhs);

		friend const affine operator+(const affine& lhs, const affine& rhs);

		friend const affine operator-(const affine& x);
		
		friend const affine operator-(const affine& lhs, const affine& rhs);

		friend const affine operator*(const affine& lhs, const affine& rhs);

		friend const affine operator/(const affine& lhs, const affine& rhs);

		friend const affine ln(const affine& x);

		friend const affine exp(const affine& x);

		friend const affine sqr(const affine& x);

		double inf() const { return lb; }

		double sup() const { return ub; }

		int n_elem() const { return n; }

		static void set_max_used_index(const int idx) { max_used_index = idx; }

		static int get_max_used_index() { return max_used_index; }

		static void set_valid() { isValid = true; }

		static bool  is_valid() { return isValid; }

		friend class LP_op;

		friend class affine_propagator;

		//======================================================================
		//
		// YOU ARE NOT SUPPOSED TO USE ANY OF THE FUNCTIONS BELOW
		//
		//======================================================================
		
		affine& operator+=(const affine& );

		void set_bounds(const double lb, const double ub);

		bool set_range(const double lb, const double ub);

		bool intersect_domain(const double lb, const double ub);

		friend const affine reciprocal(const affine& x);

		friend const affine mult_const(const double c, const affine& x);

		friend void bisect(affine [], affine [], const int , const int );

        std::ostream& print(std::ostream& os) const;

		~affine() { delete[] index; delete[] value; }

		int*    index;

		double* value;

		int n;

	private:

		explicit affine(bool) : index(0), value(0), n(0), nmax(0) { }

		affine(bool, const int size) : index(new int[size]), value(new double[size]), n(0), nmax(size) { }

		double lb;

		double ub;

		int nmax;

		static bool isValid;

		static int max_used_index;

};

std::ostream& operator<<(std::ostream& , const affine &);

void mult(const double, const double, const double, const double, double&, double&);

void division(const double , const double , const double , const double , double& , double& );

const double tol_solved_lp(1.0e-4);

bool is_narrow(const double lo, const double up);

bool operator!=(const double ,  const affine& );

bool in(const double , const affine& );

// All of these are needed by the DAG
// All of them are dummy functions
bool Disjoint(const affine& , const affine& );

const affine operator&(const affine&,  const affine& );

bool ext_div( affine&, const affine&, const affine& );

const affine sqrt(const affine& );

double Inf(const affine& );

double Sup(const affine& );

}

#endif


