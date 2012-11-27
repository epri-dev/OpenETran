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
#include "Surge.h"

char surge_token[] = "surge";

struct surge *surge_head, *surge_ptr;

/* calculate surge current value, and inject it at the pole */

void inject_surge (struct surge *ptr)
{
	double cosarg, x, i;

	x = t - ptr->tstart; /* shift time for surge starting time */
	if (x > 0.0) {
		if (x > ptr->tailadvance) { /* on the tail, use exponential */
			x -= ptr->tailadvance;
			i = ptr->peak * exp (-x / ptr->tau);
		} else { /* on the front, use 1 - cosine */
			cosarg = x * ptr->cfront;
			i = ptr->peak * 0.5 * (1.0 - cos (cosarg));
		}
		*gsl_vector_ptr (ptr->parent->injection, ptr->from) += i;
		*gsl_vector_ptr (ptr->parent->injection, ptr->to) -= i;
	}
}

int init_surge_list (void)
{
	if ((surge_head = (struct surge *) malloc (sizeof *surge_head)) != NULL) {
		surge_head->next = NULL;
		surge_ptr = surge_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize surge list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}

void do_all_surges (void (*verb) (struct surge *))
{
	surge_ptr = surge_head;
	while ((surge_ptr = surge_ptr->next) != NULL) {
		verb (surge_ptr);
	}
}

/* read a current surge from the file or buffer, and set up the surge struct */

int read_surge (void)
{
	int i, j, k;
	double fpeak, ftf, ftt, ftstart;
	struct surge *ptr;
	
	(void) next_double (&fpeak);
	(void) next_double (&ftf);
	(void) next_double (&ftt);
	(void) next_double (&ftstart);
	(void) read_pairs ();
	(void) read_poles ();
	(void) reset_assignments ();
	while (!next_assignment (&i, &j, &k)) {
		if ((ptr = (struct surge *) malloc (sizeof *ptr)) != NULL) {
			move_surge (ptr, i, j, k, fpeak, ftf, ftt, ftstart);
			ptr->next = NULL;
			surge_ptr->next = ptr;
			surge_ptr = ptr;
		} else {
			if (logfp) fprintf( logfp, "can't allocate new surge\n");
			oe_exit (ERR_MALLOC);
		}
	}
	return (0);
}

/* adjust the surge parameters and location, can be called many times
during flashover critical current iterations */

void move_surge (struct surge *ptr, int i, int j, int k, double fpeak,
	double ftf, double ftt, double ftstart)
{
	double fcfront, fctail, ftailadvance, ftau;

/* set struct parameters to input parameters */
	fcfront = TWOPI / (CFKONST * ftf);
	fctail = TWOPI / (CTKONST * ftt);
	ftailadvance = 0.5 * CFKONST * ftf;
	ftau = ETKONST * (ftt - ftailadvance);
	ptr->front = ftf;
	ptr->tail = ftt;
	ptr->cfront = fcfront;
	ptr->ctail = fctail;
	ptr->tailadvance = ftailadvance;
	ptr->tstart = ftstart;
	ptr->tau = ftau;
	ptr->peak = fpeak;
	ptr->parent = find_pole (i);
	if (!ptr->parent) oe_exit (ERR_BAD_POLE);
	ptr->parent->solve = TRUE;
	ptr->from = j;
	ptr->to = k;
}
