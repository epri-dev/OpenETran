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
#include "Line.h"
#include "Source.h"

struct source *source_head, *source_ptr;

void print_source_data (struct source *ptr)
{
	int i;
	struct span *defn;
	
	fprintf (op, "Power-frequency source at pole %d.\n", 
		ptr->parent->location);
	fprintf (op, "   #        vp        vm         i\n");
	defn = find_pole_defn (ptr->parent);
	for (i = 0; i < number_of_nodes; i++) {
		fprintf (op, "%4d%10.1lf%10.1lf%10.1lf\n", 
			i, 
			gsl_vector_get (defn->vp_offset, i), 
			gsl_vector_get (defn->vm, i), 
			gsl_vector_get (ptr->val, i));
	}
}

/* add dc current to the pole injection */

void inject_source (struct source *ptr)
{
	gsl_vector_view rhs = gsl_vector_subvector (ptr->parent->injection, 1, number_of_nodes);
	gsl_vector_add (&rhs.vector, ptr->val);
}

int init_source_list (void)
{
	if ((source_head = (struct source *) malloc (sizeof *source_head)) != NULL) {
		source_head->next = NULL;
		source_head->val = NULL;
		source_ptr = source_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize source list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}

void do_all_sources (void (*verb) (struct source *))
{
	source_ptr = source_head;
	while ((source_ptr = source_ptr->next) != NULL) {
		verb (source_ptr);
	}
}

