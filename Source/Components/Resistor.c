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
#include "Source.h"
#include "Resistor.h"

char resistor_token[] = "resistor";

struct resistor *resistor_head, *resistor_ptr;

int read_resistor (void)
{
	int i, j, k;
	double r, y;
	struct resistor *ptr;
	struct source *s_ptr; // if there is a power frequency source
	struct span *defn;
	double vdc, idc;
	
	(void) next_double (&r);
	if (r != 0.0) {
		y = 1.0 / r;
	} else {
		y = Y_SHORT;
	}
	(void) read_pairs ();
	(void) read_poles ();
	(void) reset_assignments ();
	while (!next_assignment (&i, &j, &k)) {
		if ((ptr = (struct resistor *) malloc (sizeof *ptr)) != NULL) {
			ptr->Rphase = r;
			ptr->parent = find_pole (i);
			if (!ptr->parent) oe_exit (ERR_BAD_POLE);
			defn = find_pole_defn (ptr->parent);
			ptr->parent->solve = TRUE;
			add_y (ptr->parent, j, k, y);
			ptr->from = j;
			ptr->to = k;
			ptr->next = NULL;
			resistor_ptr->next = ptr;
			resistor_ptr = ptr;
			vdc = 0.0;  // power frequency bias
			if (j > 0) {
				vdc += gsl_vector_get (defn->vp_offset, j-1);
			}
			if (k > 0) {
				vdc -= gsl_vector_get (defn->vp_offset, k-1);
			}
			if (vdc != 0.0) {
				if ((s_ptr = (struct source *) malloc (sizeof *s_ptr)) != NULL) {
					if (!(s_ptr->val = gsl_vector_calloc (number_of_nodes))) {
						if (logfp) fprintf( logfp, "can't allocate source currents\n");
						oe_exit (ERR_MALLOC);
					}
					idc = vdc * y;
					if (j > 0) {
						gsl_vector_set (s_ptr->val, j-1, idc);
					}
					if (k > 0) {
						gsl_vector_set (s_ptr->val, k-1, -idc);
					}
					s_ptr->parent = ptr->parent;
					s_ptr->next = NULL;
					source_ptr->next = s_ptr;
					source_ptr = s_ptr;
				} else {
					if (logfp) fprintf( logfp, "can't allocate new source\n");
					oe_exit (ERR_MALLOC);
				}
			} 
		} else {
			if (logfp) fprintf( logfp, "can't allocate new resistor\n");
			oe_exit (ERR_MALLOC);
		}
	}
	return (0);
}

int init_resistor_list (void)
{
	if ((resistor_head = (struct resistor *) malloc (sizeof *resistor_head)) != NULL) {
		resistor_head->next = NULL;
		resistor_ptr = resistor_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize resistor list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}

