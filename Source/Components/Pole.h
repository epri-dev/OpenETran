/*
  Copyright (c) 1992, 1994, 1998, 2002, 2011, 2012,  
  Electric Power Research Institute, Inc.
  All rights reserved.
  
  This file is part of OpenETran.

  OpenETran is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, using only version 3 of the License.

  OpenETran is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenETran.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef pole_included
#define pole_included

// extern char pair_token[];
// extern char pole_token[];
extern char pole_label_token[];
extern char phase_label_token[];

/* pole and node numbers for branch connections, used in parsing input */
extern int assign_i;
extern int assign_j;
extern int assign_k;

/* as used by LPDW, some of the poles have no insulators, arresters, grounds,
or other compenents.  In that case, we don't need to solve for phase voltages
at each time step, and the traveling waves in modal coordinates simply pass
through */

struct pole {
	int location; /* pole number, 1 to number_of_poles */
	int dirty; /* TRUE if the Ybus matrix has been modified - retriangulate */
	int solve; /* TRUE if we need to solve for phase voltages at this pole */
	int num_nonlinear;
	struct arrbez **backptr;
	gsl_vector *voltage; /* node voltages - dimensioned n+1 so voltage[0] is ground */
	gsl_vector *injection; /* vector of parallel current injections - also dimensiond n+1 */
	/* vmode and imode dimensioned number_of_nodes, don't know #conductors when allocated */
	gsl_vector *vmode; /* node voltages in modal coordinates */
	gsl_vector *imode; /* current injections in modal coordinates */
	gsl_permutation *perm; /* stores row operations for triangularizing Ybus */
	gsl_matrix *Ybus; /* nodal admittance matrix */
	gsl_matrix *y; /* triangularized Ybus */
	gsl_matrix *rcols;
	gsl_matrix *Rthev;
	gsl_permutation *jperm;
	gsl_vector *vnew;
	gsl_vector *inew;
	gsl_vector *f;
	gsl_matrix *jacobian;
	struct pole *next;
};

extern struct pole *pole_head, *pole_ptr;

int init_pole_list (void);
int reset_assignments (void);
int update_assignments (int i, int j, int k);
int next_assignment (int *i, int *j, int *k);
struct pole *new_pole (int location);
void terminate_pole (struct pole *ptr, struct span *defn);

void do_all_poles (void (*verb) (struct pole *));

struct pole *find_pole (int location);
struct span *find_pole_defn (struct pole *ptr);
void triang_pole (struct pole *ptr);
void solve_pole (struct pole *ptr);
void build_rthev (struct pole *ptr);
void zero_pole_injection (struct pole *ptr);
void calc_pole_vmode (struct pole *ptr); /* only for non-network systems */
void inject_pole_imode (struct pole *ptr); /* only for non-network systems */
void add_y (struct pole *ptr, int j, int k, double y);
void print_pole_data (struct pole *ptr);

#endif