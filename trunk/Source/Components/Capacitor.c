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

#include <wtypes.h>
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
#include "Capacitor.h"

char capacitor_token[] = "capacitor";

struct capacitor *capacitor_head, *capacitor_ptr;

void inject_capacitor_history (struct capacitor *ptr)
{
	gsl_vector *c;
	
	c = ptr->parent->injection;
	*gsl_vector_ptr (c, ptr->from) -= ptr->h;
	*gsl_vector_ptr (c, ptr->to) += ptr->h;
}

/* calculate past history currents from branch voltage */

void update_capacitor_history (struct capacitor *ptr)
{
	gsl_vector *v;
	
	v = ptr->parent->voltage;
	ptr->h = ptr->yc * (gsl_vector_get (v, ptr->to) - gsl_vector_get (v, ptr->from)) - ptr->h;
}

/* past history currents start at zero, unless there is trapped charge */

void init_capacitor_history (struct capacitor *ptr)
{
	double vdc;
	int i, j;
	struct span *defn;
	
	vdc = 0.0;
	i = ptr->from;
	j = ptr->to;
	defn = find_pole_defn (ptr->parent);
	if (i) {
		vdc += gsl_vector_get (defn->vp_offset, i-1);
	}
	if (j) {
		vdc -= gsl_vector_get (defn->vp_offset, j-1);
	}
	if (vdc != 0.0) {
		ptr->h = -vdc * ptr->y;
	}
}

int init_capacitor_list (void)
{
	if ((capacitor_head = (struct capacitor *) malloc (sizeof *capacitor_head)) != NULL) {
		capacitor_head->next = NULL;
		capacitor_ptr = capacitor_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize capacitor list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}

void do_all_capacitors (void (*verb) (struct capacitor *))
{
	capacitor_ptr = capacitor_head;
	while ((capacitor_ptr = capacitor_ptr->next) != NULL) {
		verb (capacitor_ptr);
	}
}

/* read capacitor input data from the file or buffer, then set up
the data struct */

int read_capacitor (void)
{
	int i, j, k;
	double cap, f_y, f_yc;
	struct capacitor *ptr;
	
	(void) next_double (&cap);
	f_y = 2 * cap / dT;
	f_yc = f_y + f_y;
	(void) read_pairs ();
	(void) read_poles ();
	(void) reset_assignments ();
	while (!next_assignment (&i, &j, &k)) {
		if ((ptr = (struct capacitor *) malloc (sizeof *ptr)) != NULL) {
			ptr->y = f_y;
			ptr->yc = f_yc;
			reset_capacitor (ptr);
			ptr->parent = find_pole (i);
			if (!ptr->parent) oe_exit (ERR_BAD_POLE);
			ptr->parent->solve = TRUE;
			add_y (ptr->parent, j, k, f_y);
			ptr->from = j;
			ptr->to = k;
			ptr->next = NULL;
			capacitor_ptr->next = ptr;
			capacitor_ptr = ptr;
		} else {
			if (logfp) fprintf( logfp, "can't allocate new capacitor\n");
			oe_exit (ERR_MALLOC);
		}
	}
	return (0);
}

/* zero out the past history current */

void reset_capacitor (struct capacitor *ptr)
{
	ptr->h = 0.0;
}
