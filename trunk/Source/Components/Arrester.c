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
#include "Meter.h"
#include "../WritePlotFile.h"
#include "../Parser.h"
#include "../ReadUtils.h"
#include "Pole.h"

#include "Arrester.h"

char arrester_token[] = "arrester";

struct arrester *arrester_head, *arrester_ptr;

void print_arrester_data (struct arrester *ptr)
{
	if (ptr->t_start > 0.0) {
		fprintf (op, "Arrester at pole %d, from %d to %d ",
			ptr->parent->location, ptr->from, ptr->to);
		fprintf (op, "discharged %le Amperes.\n", ptr->i_peak);
		fprintf (op, "\tTime operated: %le", ptr->t_start);
		fprintf (op, "\tTime of peak: %le\n", ptr->t_peak);
		fprintf (op, "\tCharge: %le\n", ptr->charge);
		fprintf (op, "\tEnergy: %le\n", ptr->energy);
	}
}

/* find the highest discharge duty among all arresters */

void arrester_answers_cleanup (struct arrester *ptr)
{
	if (ptr->conducting) {
		add_y (ptr->parent, ptr->from, ptr->to, -ptr->y);
	}
	if (ptr->energy > energy) {
		energy = ptr->energy;
	}
	if (fabs (ptr->i_peak) > fabs (current)) {
		current = ptr->i_peak;
	}
	if (fabs (ptr->charge) > fabs (charge)) {
		charge = ptr->charge;
	}
}

/* add arrester current to the pole, if conducting */

void inject_arrester (struct arrester *ptr)
{
	gsl_vector *c;
	double val;
	
	if (ptr->conducting) {
		c = ptr->parent->injection;
		val = ptr->i_past;
		*gsl_vector_ptr (c, ptr->from) -= val;
		*gsl_vector_ptr (c, ptr->to) += val;
	}
}

/* see if the arrester gap flashed over, or started conduction.  If
the arrester is gapped, we have one time step at the gap sparkover voltage,
then switch to the VI characteristic.  Arrester clears immediately
after current zero, and can operate any number of times. */

void check_arrester (struct arrester *ptr)
{
	struct pole *p;
	int i, j, pos_now;
	double volts, amps, vl, vr;
	
	p = ptr->parent;
	i = ptr->from;
	j = ptr->to;
	volts = gsl_vector_get (p->voltage, i) - gsl_vector_get (p->voltage, j); /* find arrester voltage and polarity */
	if (volts > 0.0) {
		pos_now = TRUE;
	} else {
		pos_now = FALSE;
	}

	if (ptr->conducting) {
		amps = volts * ptr->y + ptr->i_past;   /* i_past = injection */
		ptr->amps = amps;
		if (pos_now) {
			vr = ptr->r_slope * (amps + ptr->i_bias);
		} else {
			vr = ptr->r_slope * (amps - ptr->i_bias);
		}
		ptr->i_bias = ptr->knee_bias;
		vl = volts - vr;
		ptr->energy += dT * amps * vr;
		ptr->charge += dT * amps;
		if (ptr->zl > 0.0) {
			ptr->h = amps + vl / ptr->zl;
		}
		ptr->i = ptr->h * ptr->yzl;	/* i = new injection */
		if (pos_now) {
			ptr->i -= ptr->yr * ptr->i_bias;
		} else {
			ptr->i += ptr->yr * ptr->i_bias;
		}
		if (fabs (amps) > fabs (ptr->i_peak)) {
			ptr->i_peak = amps;
			ptr->t_peak = t;
		}
		if (fabs (vr) < ptr->v_knee) { /* voltage dropped below knee - stop conduction */
			ptr->conducting = FALSE;
			add_y (p, i, j, -ptr->y);
			ptr->h = ptr->i = 0.0;
		}
	} else {  /* arrester not conducting - check for sparkover */
		if (fabs (volts) > ptr->v_gap) {
			ptr->conducting = TRUE;
			add_y (p, i, j, ptr->y); /* add VI segment slope to pole matrix */
			ptr->i_bias = ptr->gap_bias;
			if (pos_now) {
				ptr->i = -ptr->yr * ptr->i_bias;
			} else {
				ptr->i = ptr->yr * ptr->i_bias;
			}
			solution_valid = FALSE; /* force a re-solve for this time step */
			ptr->i_past = ptr->i;  /* update injection for turn-on */
			if (ptr->t_start < dT) {
				ptr->t_start = t;
			}
		}
	}
}

void update_arrester_history (struct arrester *ptr)
{
	ptr->i_past = ptr->i;
}

int init_arrester_list (void)
{
	if (((arrester_head = (struct arrester *) malloc (sizeof *arrester_head)) != NULL)) {
		arrester_head->next = NULL;
		arrester_ptr = arrester_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize arrester list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}

void do_all_arresters (void (*verb) (struct arrester *))
{
	arrester_ptr = arrester_head;
	while (((arrester_ptr = arrester_ptr->next) != NULL)) {
		verb (arrester_ptr);
	}
}

struct arrester *find_arrester (int at, int from, int to)
{
	arrester_ptr = arrester_head;
	while (((arrester_ptr = arrester_ptr->next) != NULL)) {
		if ((arrester_ptr->parent->location == at) &&
			(arrester_ptr->from == from) &&(arrester_ptr->from == from)) return arrester_ptr;
	}
	return NULL;
}

/* read the arrester data from file or buffer, set up structs */

int read_arrester (void)
{
	int i, j, k;
	double f_knee, f_gap, f_r, knee_bias, gap_bias, f_L, f_length;
	struct arrester *ptr;
	double *target;
	int monitor;
	
	(void) next_double (&f_gap);
	(void) next_double (&f_knee);
	(void) next_double (&f_r);
	(void) next_double (&f_L);
	(void) next_double (&f_length);
	if (f_gap < 0.0) {  /* input Vgap < 0 means we want an ammeter */
		f_gap *= -1.0;
		monitor = 1;
	} else {
		monitor = 0;
	}
	f_L *= f_length;
	if (f_knee < 0.0) {
		f_knee *= -1.0;
	}
	if (f_r < 0.0) {
		f_r *= -1.0;
	}
	if (f_gap < f_knee) {
		f_gap = f_knee;
	}
	knee_bias = f_knee / f_r;
	gap_bias = f_gap / f_r;
	(void) read_pairs ();
	(void) read_poles ();
	(void) reset_assignments ();
	while (!next_assignment (&i, &j, &k)) {
		if (((ptr = (struct arrester *) malloc (sizeof *ptr)) != NULL)) {
			ptr->v_knee = f_knee;
			ptr->v_gap = f_gap;
			ptr->r_slope = f_r;
			ptr->knee_bias = knee_bias;
			ptr->gap_bias = gap_bias;
			ptr->zl = 2.0 * f_L / dT;
			ptr->y = 1.0 / (f_r + ptr->zl);
			ptr->yr = ptr->y * ptr->r_slope;
			ptr->yzl = ptr->y * ptr->zl;
			ptr->parent = find_pole (i);
			if (!ptr->parent) oe_exit (ERR_BAD_POLE);
			ptr->parent->solve = TRUE;
			ptr->from = j;
			ptr->to = k;
			reset_arrester (ptr);
			ptr->next = NULL;
			arrester_ptr->next = ptr;
			arrester_ptr = ptr;
			if (monitor) { /* measure discharge current */
				target = &(ptr->amps);
				(void) add_ammeter (i, j, IARR_FLAG, target);
			}
		} else {
			if (logfp) fprintf( logfp, "can't allocate new arrester\n");
			oe_exit (ERR_MALLOC);
		}
	}
	return (0);
}

/* reset the arrester history parameters */

void reset_arrester (struct arrester *ptr)
{
	ptr->i_bias = ptr->gap_bias;
	ptr->t_start = 0.0;
	ptr->t_peak = 0.0;
	ptr->energy = 0.0;
	ptr->charge = 0.0;
	ptr->i_peak = 0.0;
	ptr->h = 0.0;
	ptr->i = 0.0;
	ptr->i_past = 0.0;
	ptr->amps = 0.0;
	ptr->conducting = FALSE;
}
