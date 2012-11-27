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
#include "../ChangeTimeStep.h"
#include "Pole.h"
#include "Monitor.h"
#include "Insulator.h"

char insulator_token[] = "insulator";

struct insulator *insulator_head, *insulator_ptr;

void print_insulator_data (struct insulator *ptr)
{
	double highest_de;
	
	highest_de = ptr->de_pos;
	if (ptr->de_neg > highest_de) {
		highest_de = ptr->de_neg;
	}
	fprintf (op, "Insulator at pole %d, from %d to %d ",
		ptr->parent->location, ptr->from, ptr->to);
	if (ptr->flashed == 1) {
		fprintf (op, "flashed over at %le seconds\n", ptr->t_flash);
	} else {
		fprintf (op, "per-unit SI = %le\n",
			pow (highest_de / ptr->de_max, 1.0 / ptr->beta));
	}
}

/* find the insulator that had the highest severity index */

void insulator_answers_cleanup (struct insulator *ptr)
{
	double highest_de;
	
	highest_de = ptr->de_pos;
	if (ptr->de_neg > highest_de) {
		highest_de = ptr->de_neg;
	}
	if (ptr->flashed == 1) {
		ptr->SI = 1.0;
		add_y (ptr->parent, ptr->from, ptr->to, -Y_SHORT);
	} else {
		ptr->SI = pow (highest_de / ptr->de_max, 1.0 / ptr->beta);
	}
	if (fabs (ptr->SI) > fabs (SI)) {
		SI = ptr->SI;
	}
}

/* see if an insulator flashed over - if so, modify the pole y matrix,
and possibly set flash_halt = TRUE.  This version of the DE model
"remembers" separate positive and negative "leaders" after polarity
changes. */

void check_insulator (struct insulator *ptr)
{
	struct pole *p;
	int i, j;
	double volts, mag, de_inc;
	
	if (!ptr->flashed && !dT_switched) { /* haven't flashed over yet - disable during second_dT */
		p = ptr->parent;
		i = ptr->from;
		j = ptr->to;
		volts = gsl_vector_get (p->voltage, i) - gsl_vector_get (p->voltage, j);
		mag = fabs (volts) - ptr->vb;
		if (mag > 0.0) {
			de_inc = pow (mag, ptr->beta) * dT; /* integrate destructive effect */
			if (volts >= 0.0) { /* add to either the positive or negative de */
				ptr->de_pos += de_inc;
			} else {
				ptr->de_neg += de_inc;
			}
		} else {
			de_inc = 0.0;
		}
		if ((ptr->de_pos >= ptr->de_max) || 
			(ptr->de_neg >= ptr->de_max)) {
			ptr->flashed = TRUE; /* flashover for either positive or negative polarity */
			if (flash_halt_enabled) {
				flash_halt = TRUE;
			}
			ptr->t_flash = t;
			add_y (p, i, j, Y_SHORT); /* change pole matrix - short out insulator */
		}
	}
}

int init_insulator_list (void)
{
	if ((insulator_head = (struct insulator *) malloc (sizeof *insulator_head)) != NULL) {
		insulator_head->next = NULL;
		insulator_ptr = insulator_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize insulator list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}

void do_all_insulators (void (*verb) (struct insulator *))
{
	insulator_ptr = insulator_head;
	while ((insulator_ptr = insulator_ptr->next) != NULL) {
		verb (insulator_ptr);
	}
}

/* read insulator data from file or buffer, set up structs */

int read_insulator (void)
{
	int i, j, k;
	double f_cfo, f_vb, f_de, f_beta;
	struct insulator *ptr;
	
	(void) next_double (&f_cfo);
	(void) next_double (&f_vb);
	(void) next_double (&f_beta);
	(void) next_double (&f_de);
	if (f_cfo < 0.0) {
		f_cfo *= -1.0;
	}
	f_vb *= f_cfo / 100.0e3;  /* input vb is for CFO = 100 kV, so adjust */
	f_de *= pow (f_cfo / 100.0e3, f_beta);  /* max de to cause flashover, adjusted from CFO = 100 */
	(void) read_pairs ();
	(void) read_poles ();
	(void) reset_assignments ();
	while (!next_assignment (&i, &j, &k)) {
		if ((ptr = (struct insulator *) malloc (sizeof *ptr)) != NULL) {
			ptr->cfo = f_cfo;
			ptr->de_max = f_de;
			ptr->vb = f_vb;
			ptr->beta = f_beta;
			reset_insulator (ptr);
			move_insulator (ptr, i);
			ptr->from = j;
			ptr->to = k;
			ptr->next = NULL;
			insulator_ptr->next = ptr;
			insulator_ptr = ptr;
		} else {
			if (logfp) fprintf( logfp, "can't allocate new insulator\n");
			oe_exit (ERR_MALLOC);
		}
	}
	return (0);
}

/* reset the insulator "memory" parameters */

void reset_insulator (struct insulator *ptr)
{
	ptr->de_pos = 0.0;
	ptr->de_neg = 0.0;
	ptr->t_flash = 0.0;
	ptr->SI = 0.0;
	ptr->flashed = FALSE;
}

/* move insulator to a new pole, keep same node connections.  Used when
simulating strokes to different poles under iteration control */ 

void move_insulator (struct insulator *ptr, int i)
{
	ptr->parent = find_pole (i);
	if (!ptr->parent) oe_exit (ERR_BAD_POLE);
	ptr->parent->solve = TRUE;
}
