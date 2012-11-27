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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "OETypes.h"
#include "Parser.h"
#include "ChangeTimeStep.h"
#include "OEEngine.h"
#include "OERead.h"
#include "ReadUtils.h"
#include "AllComponents.h"

#define DEFAULT_LABEL_SIZE  10

char **pole_labels; /* 0..number_of_poles labels for SuperTran graphs */
char **phase_labels; /* 0..number_of_phases labels for SuperTran graphs */
char pole_label_token[] = "labelpole";
char phase_label_token[] = "labelphase";

/*  &&&&   functions to read input file */

/* this function supervises reading from either a file or buffer */

int readfile (void)
{
	char *p;
	int i;

/* init_parser may be found in parser.c   In the DOS version, we read the
whole input file into a char buffer.  In LPDW, this char buffer was set
up and passed in from driver or xfmr.  Either way, we get input data
by parsing a char buffer in memory. */
	init_parser (sn);
	number_of_poles = number_of_conductors = number_of_nodes = 0;
	first_dT = second_dT = dT_switch_time = 0.0;
	using_second_dT = dT_switched = FALSE;
	using_multiple_span_defns = FALSE;
/* read the first "line" of input - simulation control parameters */
	p = first_token ();
	if (!p) {
		if (logfp) fprintf( logfp, "elt error trying to read number of phases\n");
		oe_exit (ERR_PHASE_READ);
	}
	if (!strcmp (p, time_token)) {  /*using network span and line input */
		using_network = TRUE;
		(void) next_int (&number_of_nodes);
		(void) next_double (&dT);
		(void) next_double (&Tmax);
	} else if (!strcmp (p, change_dt_token)) {  /* using a second time step */
		using_network = FALSE;
		using_second_dT = TRUE;
		(void) next_int (&number_of_nodes);
		(void) next_int (&number_of_poles);
		(void) next_double (&span_length);
		(void) next_int (&left_end_z);
		(void) next_int (&right_end_z);
		(void) next_double (&dT);
		(void) next_double (&Tmax);
		if (number_of_poles > 0) {
			poles_used = gsl_vector_int_calloc (number_of_poles);
			for (i = 1; i <= number_of_poles; i++) {
				(void) new_pole (i);
			}
		} else {
			if (logfp) fprintf( logfp, "bad number of poles: %d\n", number_of_poles);
			oe_exit (ERR_NPOLES);
		}
		(void) next_double (&dT_switch_time);
	} else {
		using_network = FALSE;
		number_of_nodes = atoi (p);
		(void) next_int (&number_of_poles);
		(void) next_double (&span_length);
		(void) next_int (&left_end_z);
		(void) next_int (&right_end_z);
		(void) next_double (&dT);
		(void) next_double (&Tmax);
		if (number_of_poles > 0) {
			poles_used = gsl_vector_int_calloc (number_of_poles);
			for (i = 1; i <= number_of_poles; i++) {
				(void) new_pole (i);
			}
		} else {
			if (logfp) fprintf( logfp, "bad number of poles: %d\n", number_of_poles);
			oe_exit (ERR_NPOLES);
		}
	}
	first_dT = dT;  /* so we can reset dT under loop control */
/* set up arrays needed for branch connections */
	if (number_of_nodes > 0) {
		pairs_used = gsl_matrix_int_calloc (number_of_nodes, number_of_nodes);
	} else {
		if (logfp) fprintf( logfp, "bad number of nodes: %d\n", number_of_nodes);
		oe_exit (ERR_NPHASES);
	}
/* the conductor data must follow the first line */
	number_of_conductors = number_of_nodes; /* unless fewer specified on conductor cards */
	if (using_network) {
		read_spans ();
		read_lines (); /* also determines number_of_poles */
		if (number_of_poles > 0) {
			poles_used = gsl_vector_int_calloc (number_of_poles);
		} else {
			if (logfp) fprintf( logfp, "bad number of poles: %d\n", number_of_poles);
			oe_exit (ERR_NPOLES);
		}
	} else {
		read_conductors (span_ptr);
	}
	reset_lines ();

/* set up default pole and phase labels */
	pole_labels = (char **) malloc ((size_t) (number_of_poles + 1) * sizeof (char *));
	for (i = 0; i <= number_of_poles; i++) {
		pole_labels[i] = (char *) malloc (DEFAULT_LABEL_SIZE * sizeof (char));
		sprintf (pole_labels[i], "%d", i);
	}
	phase_labels = (char **) malloc ((size_t) (number_of_nodes + 1) * sizeof (char *));
	for (i = 0; i <= number_of_nodes; i++) {
		phase_labels[i] = (char *) malloc (DEFAULT_LABEL_SIZE * sizeof (char));
		sprintf (phase_labels[i], "%d", i);
	}

/* after reading the conductor data, go through the rest of the input
looking for string tokens that identify various model components.  As
each token is identified, dispatch to the appropriate reading function.
We stop when there are no more tokens in the buffer */
	while ((p = first_token ()) != NULL) {
		if (!strcmp (p, ground_token)) {
			(void) read_ground ();
		} else if (!strcmp (p, arrester_token)) {
			(void) read_arrester ();
		} else if (!strcmp (p, pipegap_token)) {
			(void) read_pipegap ();
		} else if (!strcmp (p, meter_token)) {
			(void) read_meter ();
		} else if (!strcmp (p, pole_label_token)) {
			(void) read_pole_label ();
		} else if (!strcmp (p, phase_label_token)) {
			(void) read_phase_label ();
		} else if (!strcmp (p, surge_token)) {
			(void) read_surge ();
		} else if (!strcmp (p, insulator_token)) {
			(void) read_insulator ();
		} else if (!strcmp (p, resistor_token)) {
			(void) read_resistor ();
		} else if (!strcmp (p, inductor_token)) {
			(void) read_inductor ();
		} else if (!strcmp (p, customer_token)) {
			(void) read_customer ();
		} else if (!strcmp (p, capacitor_token)) {
			(void) read_capacitor ();
		} else if (!strcmp (p, arrbez_token)) {
			(void) read_arrbez ();
		} else if (!strcmp (p, lpm_token)) {
			(void) read_lpm ();
		} else if (!strcmp (p, steepfront_token)) {
			(void) read_steepfront ();
		}
	}
	return (0);
}

void reset_system (void)
{
	restore_time_step ();
	do_all_grounds (reset_ground);
	do_all_arresters (reset_arrester);
	do_all_arrbezs (reset_arrbez);
	do_all_meters (reset_meter);
	do_all_insulators (reset_insulator);
	do_all_lpms (reset_lpm);
	do_all_inductors (reset_inductor);
	do_all_customers (reset_customer);
	do_all_capacitors (reset_capacitor);
	reset_lines ();
	do_all_lines (init_line_history);
	do_all_inductors (init_inductor_history);
	do_all_capacitors (init_capacitor_history);
	do_all_poles (triang_pole);
}

/* free memory */

int cleanup (void)
{
	int i;
	if (pole_labels) {
		for (i = 0; i <= number_of_poles; i++) {
			free (pole_labels[i]);
		}
		free (pole_labels);
	}
	if (phase_labels) {
		for (i = 0; i <= number_of_nodes; i++) {
			free (phase_labels[i]);
		}
		free (phase_labels);
	}
	if (poles_used) {
		gsl_vector_int_free (poles_used);
	}
	if (pairs_used) {
		gsl_matrix_int_free (pairs_used);
	}
	while (span_head) {
		span_ptr = span_head->next;
		if (number_of_conductors > 0) {
			gsl_matrix_free (span_head->Ti);
			gsl_matrix_free (span_head->Tit);
			gsl_matrix_free (span_head->Tv);
			gsl_matrix_free (span_head->Tvt);
			gsl_matrix_free (span_head->Zp);
			gsl_matrix_free (span_head->Zm);
			gsl_matrix_free (span_head->Yp);
			gsl_matrix_free (span_head->Ym);
			gsl_vector_free (span_head->vm);
			gsl_vector_free (span_head->vp_offset);
		}
		free (span_head);
		span_head = span_ptr;
	}
	while (surge_head) {
		surge_ptr = surge_head->next;
		free (surge_head);
		surge_head = surge_ptr;
	}
	while (source_head) {
		source_ptr = source_head->next;
		if (source_head->val) {
			gsl_vector_free (source_head->val);
		}
		free (source_head);
		source_head = source_ptr;
	}
	while (ground_head) {
		ground_ptr = ground_head->next;
		free (ground_head);
		ground_head = ground_ptr;
	}
	while (resistor_head) {
		resistor_ptr = resistor_head->next;
		free (resistor_head);
		resistor_head = resistor_ptr;
	}
	while (inductor_head) {
		inductor_ptr = inductor_head->next;
		free (inductor_head);
		inductor_head = inductor_ptr;
	}
	while (capacitor_head) {
		capacitor_ptr = capacitor_head->next;
		free (capacitor_head);
		capacitor_head = capacitor_ptr;
	}
	while (customer_head) {
		customer_ptr = customer_head->next;
		free (customer_head);
		customer_head = customer_ptr;
	}
	while (insulator_head) {
		insulator_ptr = insulator_head->next;
		free (insulator_head);
		insulator_head = insulator_ptr;
	}
	while (arrester_head) {
		arrester_ptr = arrester_head->next;
		free (arrester_head);
		arrester_head = arrester_ptr;
	}
	while (pipegap_head) {
		pipegap_ptr = pipegap_head->next;
		free (pipegap_head);
		pipegap_head = pipegap_ptr;
	}
	while (meter_head) {
		meter_ptr = meter_head->next;
		free (meter_head);
		meter_head = meter_ptr;
	}
	while (line_head) {
		line_ptr = line_head->next;
		if (line_head->hist_left) {
			gsl_matrix_free (line_head->hist_left);
		}
		if (line_head->hist_right) {
			gsl_matrix_free (line_head->hist_right);
		}
		free (line_head);
		line_head = line_ptr;
	}
	while (pole_head) {
		pole_ptr = pole_head->next;
		if (pole_head->vmode) gsl_vector_free (pole_head->vmode);  
		if (pole_head->imode) gsl_vector_free (pole_head->imode);
		if (pole_head->voltage) gsl_vector_free (pole_head->voltage);
		if (pole_head->injection) gsl_vector_free (pole_head->injection);
		if (pole_head->Ybus) gsl_matrix_free (pole_head->Ybus);
		if (pole_head->rcols)  gsl_matrix_free (pole_head->rcols);
		if (pole_head->Rthev) gsl_matrix_free (pole_head->Rthev);
		if (pole_head->backptr) free (pole_head->backptr);
		if (pole_head->vnew) gsl_vector_free (pole_head->vnew);
		if (pole_head->inew) gsl_vector_free (pole_head->inew);
		if (pole_head->f) gsl_vector_free (pole_head->f);
		if (pole_head->jperm) gsl_permutation_free (pole_head->jperm);
		if (pole_head->jacobian) gsl_matrix_free (pole_head->jacobian);
		if (pole_head->y) gsl_matrix_free (pole_head->y);
		if (pole_head->perm) gsl_permutation_free (pole_head->perm);
		free (pole_head);
		pole_head = pole_ptr;
	}
	while (arrbez_head) {
		arrbez_ptr = arrbez_head->next;
		if (arrbez_head->shape) {
			free_bezier_fit (arrbez_head->shape);
			free (arrbez_head->shape);
		}
		free (arrbez_head);
		arrbez_head = arrbez_ptr;
	}
	while (lpm_head) {
		lpm_ptr = lpm_head->next;
		if (lpm_head->pts) {
			free (lpm_head->pts);
		}
		free (lpm_head);
		lpm_head = lpm_ptr;
	}
	while (steepfront_head) {
		steepfront_ptr = steepfront_head->next;
		if (steepfront_head->shape) {
			free_bezier_fit (steepfront_head->shape);
			free (steepfront_head->shape);
		}
		free (steepfront_head);
		steepfront_head = steepfront_ptr;
	}
	if (sp) {
		free (sp);
	}
	return (0);
}

