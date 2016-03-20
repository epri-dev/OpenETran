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

/* This module contains the time-step loop and supporting functions for
LPDW's transient simulation engine, based on H. W. Dommel's IEEE papers.
Called from the driver and xfmr modules. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "../OETypes.h"
#include "../Parser.h"
#include "../ReadUtils.h"
#include "Meter.h"
#include "../WritePlotFile.h"
#include "Line.h"
#include "ArrBez.h"
#include "Source.h"
#include "Pole.h"

#undef LOG_POLES_AND_LINES
#undef LOG_ARRBEZ

struct pole *pole_head, *pole_ptr;

/* &&&&  pole functions   */

struct pole *find_pole (int location)  /* return pointer to pole # location */
{
	pole_ptr = pole_head;
	while (((pole_ptr = pole_ptr->next) != NULL)) {
		if (pole_ptr->location == location) {
			return (pole_ptr);
		}
	}
	return (NULL);
}

struct span *find_pole_defn (struct pole *ptr)
{
	struct line *ln;
	if (using_network) {
		ln = line_head;
		while (((ln = ln->next) != NULL)) {
			if (ln->left == ptr || ln->right == ptr) {
				return ln->defn;
			}
		}
	}
	return span_head;
}

/* add a branch to the pole's y matrix */

void add_y (struct pole *ptr, int j, int k, double y)
{
	if (j != 0) {
		*gsl_matrix_ptr (ptr->Ybus, j-1, j-1) += y;
	}
	if (k != 0) {
		*gsl_matrix_ptr (ptr->Ybus, k-1, k-1) += y;
	}
	if (j && k) {
		*gsl_matrix_ptr (ptr->Ybus, j-1, k-1) -= y;
		*gsl_matrix_ptr (ptr->Ybus, k-1, j-1) -= y;
	}
	ptr->dirty = TRUE;
}

void print_pole_data (struct pole *ptr)
{
	int i, j;
	if (op) {
		fprintf (op, "pole %d, step %d, t = %g\n", ptr->location, step, t);
		fprintf (op, "\t         Vmode         Imode       Voltage    Injection\n");
		fprintf (op, "\t                             %13.5e %13.5e\n", 
			gsl_vector_get (ptr->voltage, 0),
			gsl_vector_get (ptr->injection, 0));
		for (i = 0; i < number_of_nodes; i++) {
			fprintf (op, "\t %13.5e %13.5e %13.5e %13.5e\n", 
				gsl_vector_get (ptr->vmode, i),
				gsl_vector_get (ptr->imode, i),
				gsl_vector_get (ptr->voltage, i+1),
				gsl_vector_get (ptr->injection, i+1));
		}
		fprintf (op, "\tYbus\n");
		for (i = 0; i < number_of_nodes; i++) {
			fprintf (op, "\t");
			for (j = 0; j < number_of_nodes; j++) {
				fprintf (op, " %14.5e", gsl_matrix_get (ptr->Ybus, i, j));
			}
			fprintf (op, "\n");
		}
		fprintf (op, "\ty\n");
		for (i = 0; i < number_of_nodes; i++) {
			fprintf (op, "\t");
			for (j = 0; j < number_of_nodes; j++) {
				fprintf (op, " %14.5e", gsl_matrix_get (ptr->y, i, j));
			}
			fprintf (op, "\n");
		}
	}
}

static struct arrbez *match_arrbez (struct pole *ptr, int match_num_nonlinear)
{
	struct arrbez *aptr = arrbez_head;
	while (aptr) {
		if (aptr->parent == ptr) {
			if (aptr->pole_num_nonlinear == match_num_nonlinear) {
				return aptr;
			}
		}
		aptr = aptr->next;
	}
	return NULL;
}

void build_rthev (struct pole *ptr)  // must call after ludcmp on ptr->y
{
	int i, j, k, m;
	struct arrbez *aptr;
	gsl_vector_view rhs;

	gsl_matrix_set_zero (ptr->rcols);
	gsl_matrix_set_zero (ptr->Rthev);

	for (i = 0; i < ptr->num_nonlinear; i++) {
		aptr = ptr->backptr[i];
		if (!aptr) {
			if (logfp) fprintf (logfp, "can't find matching arrbez for Thevenin reduction\n");
			oe_exit (ERR_LT_STOPPED);
		}
		gsl_vector_set_zero (ptr->voltage);
		k = aptr->from;
		m = aptr->to;
		if (k > 0) {
			gsl_vector_set (ptr->voltage, k, 1.0);
		}
		if (m > 0) {
			gsl_vector_set (ptr->voltage, m, -1.0);
		}
		rhs = gsl_vector_subvector (ptr->voltage, 1, number_of_nodes);
		gsl_linalg_LU_svx (ptr->y, ptr->perm, &rhs.vector);
		for (j = 0; j < number_of_nodes; j++) {
			gsl_matrix_set (ptr->rcols, i, j, gsl_vector_get (ptr->voltage, j+1));
		}
	}
	for (i = 0; i < ptr->num_nonlinear; i++) {
		for (j = 0; j < ptr->num_nonlinear; j++) {
			aptr = ptr->backptr[j];
			k = aptr->from;
			m = aptr->to;
			if (k > 0) {
				*gsl_matrix_ptr (ptr->Rthev, i, j) += *gsl_matrix_ptr (ptr->rcols, i, k-1);
			}
			if (m > 0) {
				*gsl_matrix_ptr (ptr->Rthev, i, j) -= *gsl_matrix_ptr (ptr->rcols, i, m-1);
			}
		}
	}
}

/* factor the pole Ybus matrix for time-step solutions */
#define Y_OPEN        1.0e-9      /* admittance for an "open circuit" */

void triang_pole (struct pole *ptr)
{
	int i, j;
	int signum;
	struct arrbez *aptr; 

	if (ptr->num_nonlinear > 0 && !ptr->Rthev) {
		ptr->rcols = gsl_matrix_calloc (ptr->num_nonlinear, number_of_nodes);
		ptr->Rthev = gsl_matrix_calloc (ptr->num_nonlinear, ptr->num_nonlinear);
		ptr->backptr = (struct arrbez **) malloc (ptr->num_nonlinear * sizeof (struct arrbez *));
		ptr->inew = gsl_vector_calloc (ptr->num_nonlinear);
		ptr->vnew = gsl_vector_calloc (ptr->num_nonlinear);
		ptr->f = gsl_vector_calloc (ptr->num_nonlinear);
		ptr->jperm = gsl_permutation_alloc (ptr->num_nonlinear);
		ptr->jacobian = gsl_matrix_calloc (ptr->num_nonlinear, ptr->num_nonlinear);
		for (j = 0; j < ptr->num_nonlinear; j++) {
			aptr = match_arrbez (ptr, j+1);
			if (!aptr) {
				if (logfp) fprintf (logfp, "can't find matching arrbez for Thevenin reduction\n");
				oe_exit (ERR_LT_STOPPED);
			}
			ptr->backptr[j] = aptr;
		}
	}	
	if (ptr->dirty && ptr->solve) { /* factor only if we need to */
		gsl_matrix_memcpy (ptr->y, ptr->Ybus);
		for (i = 0; i < number_of_nodes; i++) {
			if (*gsl_matrix_ptr (ptr->y, i, i) <= 0.0) {
				gsl_matrix_set (ptr->y, i, i, Y_OPEN);
			}
		}
		gsl_linalg_LU_decomp (ptr->y, ptr->perm, &signum);
		if (ptr->num_nonlinear > 0) {
			build_rthev (ptr);
		}
		ptr->dirty = FALSE;
	}
}

/* solve for node voltages by back substitution */

#define MAX_NR_ITER   100
#define NR_TOLX      1e-8
#define NR_TOLF      1e-8

/* this version used when solving for voltages */

static void form_newton_raphson (struct pole *ptr, gsl_matrix *jacobian, gsl_vector *vnew, gsl_vector *f,
								 double *bezval, double *bezd1, double *voc)
{
	int i, j;

	gsl_matrix_memcpy (jacobian, ptr->Rthev);
	for (i = 0; i < ptr->num_nonlinear; i++) {
		*gsl_matrix_ptr (jacobian, i, i) += 1.0 / bezd1[i];
		gsl_vector_set (f, i, voc[i] - gsl_vector_get (vnew, i));
#ifdef LOG_ARRBEZ
		if (op) {
			fprintf (op, "\tjacobian[%d,%d] = %g\n", i, i, gsl_matrix_get (jacobian, i, i));
			fprintf (op, "\tf[%d] = %g\n", i, gsl_vector_get (f, i));
		}
#endif
	}
	for (i = 0; i < ptr->num_nonlinear; i++) {
		for (j = 0; j < ptr->num_nonlinear; j++) {
			*gsl_vector_ptr (f, i) -= gsl_matrix_get (ptr->Rthev, i, j) * bezval[j];
		}
#ifdef LOG_ARRBEZ
		if (op) {
			fprintf (op, "\tf[%d] = %g\n", i, gsl_vector_get (f, i));
		}
#endif
	}
}

/* solve for pole voltages by back substitution */

void solve_pole (struct pole *ptr)
{
	struct arrbez *aptr;
	int i, k, m, count, signum;
	double errf, errx, vl;
	gsl_vector *inew = ptr->inew;
	gsl_vector *vnew = ptr->vnew;
	gsl_vector *f = ptr->f;
	gsl_matrix *jacobian = ptr->jacobian;
	gsl_permutation *jperm = ptr->jperm;
	static double bezval[10];  //REVISIT - use #define for max number of nonlinears
	static double bezd1[10];   // start indexing these at one
	static double voc[10];
	gsl_vector_view rhs, inj;
	
	rhs = gsl_vector_subvector (ptr->voltage, 1, number_of_nodes);
	inj = gsl_vector_subvector (ptr->injection, 1, number_of_nodes);
	
	if (ptr->solve) {
		gsl_vector_memcpy (&rhs.vector, &inj.vector);
		gsl_linalg_LU_svx (ptr->y, ptr->perm, &rhs.vector);
	}
	// now ptr->voltage has the open-circuit voltage
	if (ptr->num_nonlinear > 0) {
#ifdef LOG_ARRBEZ
		if (op) {
			fprintf (op, "Arrbez iteration begins at %g\n", t);
		}
#endif
		count = 0;
		errx = 2.0 * NR_TOLX;
		errf = 2.0 * NR_TOLF;
		for (i = 0; i < ptr->num_nonlinear; i++) {  // initial voltage guess
			aptr = ptr->backptr[i];
			k = aptr->from;
			m = aptr->to;
			gsl_vector_set (inew, i, aptr->amps);
			voc[i] = 0.0;
			if (k > 0) voc[i] += gsl_vector_get (ptr->voltage, k);
			if (m > 0) voc[i] -= gsl_vector_get (ptr->voltage, m);
			if (aptr->rl > 0.0) {
				voc[i] += aptr->h * aptr->rl;
			}
			gsl_vector_set (vnew, i, voc[i]);
#ifdef LOG_ARRBEZ
			if (op) {
				fprintf (op, "\tvoc[%d] = %g\n", i, voc[i]);
			}
#endif
		}
		for (i = 0; i < ptr->num_nonlinear; i++) {
			aptr = ptr->backptr[i];
			*gsl_matrix_ptr (ptr->Rthev, i, i) += aptr->r;
			for (k = 0; k < ptr->num_nonlinear; k++) {
				*gsl_vector_ptr (vnew, i) -= gsl_matrix_get (ptr->Rthev, i, k) * gsl_vector_get (inew, k);
			}
			bezval[i] = bez_eval (aptr->shape, gsl_vector_get (vnew, i));
			bezd1[i] = bez_d1 (aptr->shape, gsl_vector_get (vnew, i));
#ifdef LOG_ARRBEZ
			if (op) {
				fprintf (op, "\tvnew[%d] = %g\n", i, gsl_vector_get (vnew, i));
				fprintf (op, "\tbezval[%d] = %g\n", i, bezval[i]);
				fprintf (op, "\tbezd1[%d] = %g\n", i, bezd1[i]);
			}
#endif
		}
		while (count < MAX_NR_ITER && errx > NR_TOLX && errf > NR_TOLF) {
			++count;
			++nr_iter;
#ifdef LOG_ARRBEZ
			if (op) {
				fprintf (op, "\tcount %d, iter %d\n", count, nr_iter);
			}
#endif
			errx = errf = 0.0;
			form_newton_raphson (ptr, jacobian, vnew, f, bezval, bezd1, voc);
			for (i = 0; i < ptr->num_nonlinear; i++) {
				errf += fabs (gsl_vector_get (f, i));
			}
			gsl_linalg_LU_decomp (jacobian, jperm, &signum);
			gsl_linalg_LU_svx (jacobian, jperm, f);
			for (i = 0; i < ptr->num_nonlinear; i++) {
				errx += fabs (gsl_vector_get (f, i));
				aptr = ptr->backptr[i];
				*gsl_vector_ptr (vnew, i) += gsl_vector_get (f, i) / bezd1[i];
				bezval[i] = bez_eval (aptr->shape, gsl_vector_get (vnew, i));
				bezd1[i] = bez_d1 (aptr->shape, gsl_vector_get (vnew, i));
#ifdef LOG_ARRBEZ
				if (op) {
					fprintf (op, "\tvnew[%d] = %g\n", i, gsl_vector_get (vnew, i));
					fprintf (op, "\tbezval[%d] = %g\n", i, bezval[i]);
					fprintf (op, "\tbezd1[%d] = %g\n", i, bezd1[i]);
				}
#endif
			}
		}
		for (i = 0; i < ptr->num_nonlinear; i++) {  // convert solved v to injected i
			gsl_vector_set (inew, i, bezval[i]);
			aptr = ptr->backptr[i];
			*gsl_matrix_ptr (ptr->Rthev, i, i) -= aptr->r;
		}
		for (i = 0; i < ptr->num_nonlinear; i++) {  // save currents for next time step
			aptr = ptr->backptr[i];
			k = aptr->from;
			m = aptr->to;
			aptr->amps = gsl_vector_get (inew, i);
			aptr->varr = gsl_vector_get (vnew, i);
			if (aptr->rl > 0.0) {
				vl = aptr->rl * (gsl_vector_get (inew, i) - aptr->h);
				aptr->h += vl * aptr->gl;
			}
			if (k > 0) *gsl_vector_ptr (ptr->injection, k) -= gsl_vector_get (inew, i);
			if (m > 0) *gsl_vector_ptr (ptr->injection, m) += gsl_vector_get (inew, i);
		}
		gsl_vector_memcpy (&rhs.vector, &inj.vector);   //  repeat the solution with compensation
		gsl_linalg_LU_svx (ptr->y, ptr->perm, &rhs.vector);
		if (count > nr_max) {
			nr_max = count;
		}
	}
}

void zero_pole_injection (struct pole *ptr)
{
	gsl_vector_set_zero (ptr->injection);
	gsl_vector_set_zero (ptr->imode);
}

/* convert the modal injection currents to phase coordinates - add to
existing injections */
 /* only for non-network systems */
void inject_pole_imode (struct pole *ptr)
{
	gsl_vector_view rhs;

	if (ptr->solve) {
		rhs = gsl_vector_subvector (ptr->injection, 1, number_of_nodes);
		gsl_blas_dgemv (CblasNoTrans, 1.0, span_head->Ti, ptr->imode, 1.0, &rhs.vector);
	}
}

/* convert node voltages to modal coordinates, for use in calculating
past history terms for the lines */
 /* only for non-network systems */

void calc_pole_vmode (struct pole *ptr)
{
	int i;
	gsl_vector_view rhs;
	
	if (ptr->solve) {
		rhs = gsl_vector_subvector (ptr->voltage, 1, number_of_nodes);
		gsl_blas_dgemv (CblasNoTrans, 1.0, span_head->Tvt, &rhs.vector, 0.0, ptr->vmode);
	} else { /* we aren't solving phase voltages at this pole, so no transformation needed */
		for (i = 0; i < number_of_conductors; i++) {  /* travelling waves "pass through" */
			gsl_vector_set (ptr->vmode, i, gsl_vector_get (ptr->imode, i) * gsl_matrix_get (span_head->Zm, i, i) * 0.5);
		}
	}
}

void do_all_poles (void (*verb) (struct pole *))
{
	pole_ptr = pole_head;
	while (((pole_ptr = pole_ptr->next) != NULL)) {
		verb (pole_ptr);
	}
}

/* construct a new pole with initialized parameters */

struct pole *new_pole (int location)
{
	struct pole *ptr;
  
	if (((ptr = (struct pole *) malloc (sizeof *ptr)) != NULL)) {
		pole_ptr->next = ptr;
		pole_ptr = ptr;
		pole_ptr->location = location;
		pole_ptr->solve = FALSE;
		pole_ptr->dirty = TRUE;
		pole_ptr->num_nonlinear = 0;
		pole_ptr->vmode = gsl_vector_calloc (number_of_nodes);
		pole_ptr->imode = gsl_vector_calloc (number_of_nodes);
		pole_ptr->voltage = gsl_vector_calloc (number_of_nodes + 1);   // [0] is ground
		pole_ptr->injection = gsl_vector_calloc (number_of_nodes + 1); // [0] is ground
		pole_ptr->perm = gsl_permutation_alloc (number_of_nodes);
		pole_ptr->Ybus = gsl_matrix_calloc (number_of_nodes, number_of_nodes);
		pole_ptr->y = gsl_matrix_calloc (number_of_nodes, number_of_nodes);
		pole_ptr->rcols = NULL;
		pole_ptr->Rthev = NULL;
		pole_ptr->backptr = NULL;
		pole_ptr->vnew = NULL;
		pole_ptr->inew = NULL;
		pole_ptr->f = NULL;
		pole_ptr->jperm = NULL;
		pole_ptr->jacobian = NULL;
		pole_ptr->next = NULL;
		return (ptr);
	} else {
		if (logfp) fprintf( logfp, "can't build pole at %d\n", location);
		oe_exit (ERR_MALLOC);
		return (NULL);
	}
}

/* add surge impedance termination, with power-frequency voltages (dc),
 to a pole */

void terminate_pole (struct pole *ptr, struct span *defn)
{
	struct source *s_ptr;
	
	gsl_matrix_add (ptr->Ybus, defn->Yp);
	if (((s_ptr = (struct source *) malloc (sizeof *s_ptr)) != NULL)) {
		if (!(s_ptr->val = gsl_vector_calloc (number_of_nodes))) { /* must cover all nodes */
			if (logfp) fprintf( logfp, "can't allocate source currents\n");
			oe_exit (ERR_MALLOC);
		}
		gsl_blas_dgemv (CblasNoTrans, 1.0, defn->Yp, defn->vp_offset, 0.0, s_ptr->val); 
		s_ptr->parent = ptr;
		s_ptr->next = NULL;
		source_ptr->next = s_ptr;
		source_ptr = s_ptr;
	} else {
		if (logfp) fprintf( logfp, "can't allocate new source\n");
		oe_exit (ERR_MALLOC);
	}
}

int init_pole_list (void)
{
	if (((pole_head = (struct pole *) malloc (sizeof *pole_head)) != NULL)) {
		pole_head->next = NULL;
		pole_head->voltage = NULL;
		pole_head->injection = NULL;
		pole_head->vmode = NULL;
		pole_head->imode = NULL;
		pole_head->perm = NULL;
		pole_head->Ybus = NULL;
		pole_head->y = NULL;
		pole_head->Rthev = NULL;
		pole_head->rcols = NULL;
		pole_head->backptr = NULL;
		pole_head->vnew = NULL;
		pole_head->inew = NULL;
		pole_head->f = NULL;
		pole_head->jperm = NULL;
		pole_head->jacobian = NULL;
		pole_ptr = pole_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize pole list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}
