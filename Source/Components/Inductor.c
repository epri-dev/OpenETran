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
#include "Pole.h"
#include "Line.h"
#include "Inductor.h"

char inductor_token[] = "inductor";

struct inductor *inductor_head, *inductor_ptr;

/* inject magnetic energy storage current at the pole */

void inject_inductor_history (struct inductor *ptr)
{
	gsl_vector *c;
	
	c = ptr->parent->injection;
	*gsl_vector_ptr (c, ptr->from) -= ptr->h;
	*gsl_vector_ptr (c, ptr->to) += ptr->h;
}

/* calculate past history currents from branch voltage */

void update_inductor_history (struct inductor *ptr)
{
	gsl_vector *v;
	
	v = ptr->parent->voltage;
	ptr->h = ptr->zi * ptr->h + ptr->yi * (gsl_vector_get (v, ptr->from) - gsl_vector_get (v, ptr->to));
}

/* past history currents start at zero - unless there is a power-frequency
initial condition.  In that case, there had better be some series R */

void init_inductor_history (struct inductor *ptr)
{
	double vdc, denom;
	struct span *defn;
	int i, j;

	vdc = 0.0;
	i = ptr->from;
	j = ptr->to;
	defn = find_pole_defn (ptr->parent);
	if (i) vdc += gsl_vector_get (defn->vp_offset, i-1);
	if (j) vdc -= gsl_vector_get (defn->vp_offset, j-1);
	if (fabs (vdc) >= V_MIN) {
		denom = 1.0 - ptr->zi;
		if (denom != 0.0) {
			ptr->h = vdc * ptr->yi / denom;
			if (logfp) {
				fprintf( logfp, "Warning!\n");
				fprintf( logfp, "Lossy inductor from %d to %d ", i, j);
				fprintf( logfp, "has an initial dc voltage.\n");
				fprintf( logfp, "Results may be invalid.\n");
			}
		} else {
			if (logfp) {
				fprintf( logfp, "Inductor from %d to %d ", i, j);
				fprintf( logfp, "has initial Vdc = %le, ", vdc);
				fprintf( logfp, "but no resistance.\n");
				fprintf( logfp, "Please add R, or change V.\n");
			}
			oe_exit (ERR_LVDC);
		}
	}
}

int init_inductor_list (void)
{
	if (((inductor_head = (struct inductor *) malloc (sizeof *inductor_head)) != NULL)) {
		inductor_head->next = NULL;
		inductor_ptr = inductor_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize inductor list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}

void do_all_inductors (void (*verb) (struct inductor *))
{
	inductor_ptr = inductor_head;
	while (((inductor_ptr = inductor_ptr->next) != NULL)) {
		verb (inductor_ptr);
	}
}

/* read input parameters for a series RL, set up struct */

int read_inductor (void)
{
	int i, j, k;
	double res, ind, f_y, f_zi, f_yi;
	struct inductor *ptr;
	
	(void) next_double (&res);
	(void) next_double (&ind);
/* integration adjustments for RL - see Dommel's papers */
	f_y = 1.0 / (res + 2.0 * ind / dT);
	f_zi = 1.0 - 2.0 * res * f_y;
	f_yi = 2.0 * f_y * (1.0 - res * f_y);
	(void) read_pairs ();
	(void) read_poles ();
	(void) reset_assignments ();
	while (!next_assignment (&i, &j, &k)) {
		if (((ptr = (struct inductor *) malloc (sizeof *ptr)) != NULL)) {
			ptr->res = res;
			ptr->ind = ind;
			ptr->y = f_y;
			ptr->zi = f_zi;
			ptr->yi = f_yi;
			reset_inductor (ptr);
			ptr->parent = find_pole (i);
			if (!ptr->parent) oe_exit (ERR_BAD_POLE);
			ptr->parent->solve = TRUE;
			add_y (ptr->parent, j, k, f_y);
			ptr->from = j;
			ptr->to = k;
			ptr->next = NULL;
			inductor_ptr->next = ptr;
			inductor_ptr = ptr;
		} else {
			if (logfp) fprintf( logfp, "can't allocate new inductor\n");
			oe_exit (ERR_MALLOC);
		}
	}
	return (0);
}

/* zero the past history current */

void reset_inductor (struct inductor *ptr)
{
	ptr->h = 0.0;
}

