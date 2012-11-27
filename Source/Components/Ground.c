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
#include "../WritePlotFile.h"
#include "Pole.h"
#include "Meter.h"
#include "Ground.h"

char ground_token[] = "ground";

struct ground *ground_head, *ground_ptr;

/* see if the ground resistance is reduced by impulse current flow */

void check_ground (struct ground *ptr)
{
	double It, Vt, Vg, Vl, Imag;
	
	Vt = gsl_vector_get (ptr->parent->voltage, ptr->from) - gsl_vector_get (ptr->parent->voltage, ptr->to);
	It = Vt * ptr->y + ptr->i;  /* total ground current */
	ptr->amps = It;
	Imag = fabs (It);
	ptr->Ri = ptr->R60 / sqrt (1.0 + Imag / ptr->Ig); /* desired impulse resistance */
	Vg = It * ptr->Ri; /* ground voltage rise caused by Ri times It */
/* inject this current into R60 to produce a back emf, so total ground voltage is Vg */
	ptr->i_bias = Vg * (1.0 / ptr->Ri - ptr->y60);
/* update past history of the built-in ground inductance */
	Vl = Vt - Vg;
	if (ptr->zl > 0.0) {
		ptr->h = It + Vl / ptr->zl;
	}
	ptr->i = ptr->h * ptr->yzl + ptr->i_bias * ptr->yr;
}

/* add ground bias current plus inductive past history current at the pole */

void inject_ground (struct ground *ptr)
{
	gsl_vector *c;
	double val;
	
	c = ptr->parent->injection;
	val = ptr->i;
	*gsl_vector_ptr (c, ptr->from) -= val;
	*gsl_vector_ptr (c, ptr->to) += val;
}

int init_ground_list (void)
{
	if ((ground_head = (struct ground *) malloc (sizeof *ground_head)) != NULL) {
		ground_head->next = NULL;
		ground_ptr = ground_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize ground list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}

void do_all_grounds (void (*verb) (struct ground *))
{
	ground_ptr = ground_head;
	while ((ground_ptr = ground_ptr->next) != NULL) {
		verb (ground_ptr);
	}
}

struct ground *find_ground (int at, int from, int to)
{
	ground_ptr = ground_head;
	while ((ground_ptr = ground_ptr->next) != NULL) {
		if ((ground_ptr->parent->location == at) &&
			(ground_ptr->from == from) &&(ground_ptr->from == from)) return ground_ptr;
	}
	return NULL;
}

/* read ground data from file or buffer, set up structs */

int read_ground (void)
{
	int i, j, k;
	double R60, Rho, e0, L, length;
	struct ground *ptr;
	int monitor;
	double *target;
	
	(void) next_double (&R60);
	if (R60 < 0.0) { /* input R60 < 0 means we want an ammeter */
		R60 *= -1.0;
		monitor = 1;
	} else {
		monitor = 0;
	}
	(void) next_double (&Rho);
	(void) next_double (&e0);
	(void) next_double (&L);
	(void) next_double (&length);
	L *= length;
	(void) read_pairs ();
	(void) read_poles ();
	(void) reset_assignments ();
	while (!next_assignment (&i, &j, &k)) {
		ptr = add_ground (i, j, k, R60, Rho, e0, L);
		if (monitor) {
			target = &(ptr->amps);
			(void) add_ammeter (i, j, IPG_FLAG, target);
		}
	}
	return (0);
}

/* reset the ground history parameters */

void reset_ground (struct ground *ptr)
{
	ptr->h = 0.0;
	ptr->i = 0.0;
	ptr->i_bias = 0.0;
	ptr->amps = 0.0;
	ptr->Ri = ptr->R60; 
}

/* add a new ground struct to the linked list, called either by read_ground
or by read_customer */
	
struct ground *add_ground (int i, int j, int k, 
double R60, double Rho, double e0, double L)
{
	struct ground *ptr;

	if ((ptr = (struct ground *) malloc (sizeof *ptr)) != NULL) {
		ptr->R60 = R60;
		ptr->y60 = 1.0 / R60;
		ptr->Ig = e0 * Rho / R60 / R60 / 6.283185;
		ptr->parent = find_pole (i);
		if (!ptr->parent) oe_exit (ERR_BAD_POLE);
		ptr->parent->solve = TRUE;
		ptr->zl = 2.0 * L / dT;
		ptr->y = 1.0 / (R60 + ptr->zl);
		ptr->yr = ptr->y * R60;
		ptr->yzl = ptr->y * ptr->zl;
		add_y (ptr->parent, j, k, ptr->y);
		ptr->from = j;
		ptr->to = k;
		ptr->next = NULL;
		reset_ground (ptr);
		ground_ptr->next = ptr;
		ground_ptr = ptr;
		return (ptr);
	} else {
		if (logfp) fprintf( logfp, "can't allocate new ground\n");
		oe_exit (ERR_MALLOC);
	}
	return (NULL);
}
