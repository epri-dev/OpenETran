/*
  Copyright (c) 1992, 1994, 1998, 2002, 2011, 2012, 
  Electric Power Research Institute, Inc.
  All rights reserved.
  
  This file is part of OpenETran.

  OpenETran is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenETran is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenETran.  If not, see <http://www.gnu.org/licenses/>.
*/

/* This module contains functions to parse the input for transient
simulation in oeengine.c */

#include <wtypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "OETypes.h"
#include "Parser.h"
#include "ChangeTimeStep.h"
#include "ReadUtils.h"
#include "Components/Meter.h"

int assign_i;
int assign_j;
int assign_k;

char pair_token[] = "pairs";
char pole_token[] = "poles";

int read_pole_label (void)
{
	int at;
	char *label;

	next_int (&at);
	label = rest_of_line ();
	set_pole_label (at, label);

	return (0);
}

int read_phase_label (void)
{
	int at;
	char *label;

	next_int (&at);
	label = rest_of_line ();
	set_phase_label (at, label);

	return (0);
}

/* read input line of the form "pairs j1 k1 j2 k2 ... " from the input buffer */

int read_pairs (void)
{
	int j, k;
	char *p;

/* initialize for no pairs marked */
	gsl_matrix_int_set_zero (pairs_used);
	p = first_token ();
	if (!strcmp (p, pair_token)) {
/* check for "pairs", then get the first pair of integer node numbers */
		p = next_token ();
		j = atoi (p);
		if (!next_int (&k)) {
			(void) mark_pair (j, k);
		} else {
			if (logfp) fprintf( logfp, "unmatched pair for %d\n", j);
			oe_exit (ERR_UNMATCHED_PAIR);
		}
/* read additional pairs, if they exist */
		while (!next_int (&j)) {
			if (!next_int (&k)) {
				(void) mark_pair (j, k);
			} else {
				if (logfp) fprintf( logfp, "unmatched pair for %d\n", j);
				oe_exit (ERR_UNMATCHED_PAIR);
			}
		}
	}
	return (0);
}

/* indicate that a branch is connected between nodes j and k */

int mark_pair (int j, int k)
{
	if ((j == 0) && (k > 0) && (k <= number_of_nodes)) {
		gsl_matrix_int_set (pairs_used, k-1, k-1, 1);  /* from k to ground */
	} else if ((k == 0) && (j > 0) && (j <= number_of_nodes)) {
		gsl_matrix_int_set (pairs_used, j-1, j-1, 1);  /* from j to ground */
	} else if ((j > 0) && (k > 0) && (j <= number_of_nodes)
		&& (k <= number_of_nodes)) {
		gsl_matrix_int_set (pairs_used, j-1, k-1, 1);  /* between nodes j and k */
	} else {
		if (logfp) fprintf( logfp, "cannot use pair %d, %d\n", j, k);
		oe_exit (ERR_BAD_PAIR);
	}
	return (0);
}

/* used in preparation for connecting component models to poles and
nodes.  We start looking at pole number 1, node 1 */

int reset_assignments ()
{
	assign_i = assign_j = assign_k = 1;
	return (0);
}

/* in the static variables assign_i, assign_j, and assign_k, keep track of
where we are in the branch connection process.  We just processed pole i,
nodes j and k. */

int update_assignments (int i, int j, int k)
{
	if (k >= number_of_nodes) {
		k = 1; /* went through all node pair columns, go to the next row */
		j++;
	} else {
		k++; /* stay in the same row, look at next column */
	}
	if (j > number_of_nodes) { /* been through all the node pair rows, go to next pole */
		j = 1;
		i++;
	}
	assign_i = i;
	assign_j = j;
	assign_k = k;
	return (0);
}

/* return the next pole (next_i) and node pair (next_j, next_k) for
component model connections.  We use the static arrays poles_used and
pairs_used to pick out the pole and pair numbers in order.  Return 0
if this function should be called again for another branch connection,
or 1 if there are no more branch connections. */

int next_assignment (int *next_i, int *next_j, int *next_k)
{
	int i, j, k;
	
	i = assign_i;
	j = assign_j;
	k = assign_k;
	while (i <= number_of_poles) {
		if (gsl_vector_int_get (poles_used, i-1) > 0) { /* loop over poles until we find one to use */
			while (j <= number_of_nodes) { /* loop over node rows */
				while (k <= number_of_nodes) { /* loop over node columns */
					if (gsl_matrix_int_get (pairs_used, j-1, k-1) > 0) { /* found a match */
 /* remember this point for next call to next_assignment */
						(void) update_assignments (i, j, k);
						*next_i = i;
						*next_j = j;
						if (j == k) {  /* to ground */
							*next_k = 0;
						} else {
							*next_k = k; /* node j to node k */
						}
						return (0);
					}
					k++; /* next column */
				}
				j++;  /* go to next row, set column to 1 */
				k = 1;
			}
		}
		i++;  /* next pole, set row to 1 */
		j = 1;
	}
	*next_i = *next_j = *next_k = 0;   /* no more assignments */
	return (1);
}

/* read an input line of the form "poles .... " from the input buffer */

int read_poles (void)
{
	int i;
	char *p;
	
/* initialize for no poles used */
	gsl_vector_int_set_zero (poles_used);
	p = first_token ();
	if (!strcmp (p, pole_token)) {
/* check for "poles", then look at the next token, which may not be integer */
		p = next_token ();
		if (!strcmp (p, "all")) {
			for (i = 0; i < number_of_poles; i++) {
				gsl_vector_int_set (poles_used, i, 1);  /* flag all poles as used, done */
			}
		} else if (!strcmp (p, "even")) {
			for (i = 1; i < number_of_poles; i += 2) {
				gsl_vector_int_set (poles_used, i, 1);  /* use all the even-number poles, done */
			}
		} else if (!strcmp (p, "odd")) {
			for (i = 0; i < number_of_poles; i += 2) {
				gsl_vector_int_set (poles_used, i, 1); /* use all the odd-number poles, done */
			}
		} else { /* selecting individual poles */
			i = atoi (p); /* should be an integer, mark that pole as used */
			if ((i > 0) && (i <= number_of_poles)) {
				gsl_vector_int_set (poles_used, i-1, 1);
			} else {
				if (logfp) fprintf( logfp, "bad pole: %d\n", i);
				oe_exit (ERR_BAD_POLE);
			}
/* read integers from the rest of that line, mark those poles as used */
			while (!next_int (&i)) { 
				if ((i > 0) && (i <= number_of_poles)) {
					gsl_vector_int_set (poles_used, i-1, 1);
				} else {
					if (logfp) fprintf( logfp, "bad pole: %d\n", i);
					oe_exit (ERR_BAD_POLE);
				}
			}
		}
	}
	return (0);
}

void set_pole_label (int pole_number, char *label)
{
	if (pole_number < 0 || pole_number > number_of_poles) return;

	if (pole_labels[pole_number]) {
		free (pole_labels[pole_number]);
	}
	pole_labels[pole_number] = (char *) malloc (strlen (label) + 1);
	strcpy (pole_labels[pole_number], label);
}

void set_phase_label (int phase_number, char *label)
{
	if (phase_number < 0 || phase_number > number_of_nodes) return;

	if (phase_labels[phase_number]) {
		free (phase_labels[phase_number]);
	}
	phase_labels[phase_number] = (char *) malloc (strlen (label) + 1);
	strcpy (phase_labels[phase_number], label);
}
