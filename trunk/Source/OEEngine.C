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
#include <gsl/gsl_roots.h>

#include "OETypes.h"
#include "OERead.h"
#include "OEEngine.h"
#include "ChangeTimeStep.h"
#include "WritePlotFile.h"
#include "AllComponents.h"

/* Cigre lightning stroke parameters */
#define Q_MEDIAN_FIRST  4.65
#define I_MEDIAN_FIRST 31.10
#define T3090_FIRST     3.83
#define Q_MEDIAN_SUBS   0.938
#define I_MEDIAN_SUBS  12.30
#define T3090_SUBS      0.67

#define MIN_STROKE 3.0e3
#define MAX_STROKE 500.0e3
#define MAX_ITER 200
#define ITER_TOL 1.0

#undef LOG_POLES_AND_LINES
#undef LOG_ARRBEZ

/* definitions for externs in ltengine.h */

char *sp;
char *sn;

int header_stop[3] = {0, 0, 0};

int using_network;  /* used time, span, and line tokens */
int using_multiple_span_defns;
int number_of_nodes;
int number_of_conductors;
int number_of_poles;

double span_length;
double dT;
double Tmax;
double t;
int step;
int solution_valid;
int left_end_z;
int right_end_z;

gsl_vector_int *poles_used;
gsl_matrix_int *pairs_used;

double predischarge;  /* maximum predischarge current in pipegaps */
double SI, energy, current, charge; /* maximum insulator severity index, and
	maximum arrester energy, current, and charge, of all components in the simulation */
int flash_halt, flash_halt_enabled; /* flags to stop simulation if an insulator flashes over */
int want_si_calculation;  /* set 1 for solution by bisection, 0 for an estimate */

int gi_iteration_mode;

/* main simulation function */

int lt (LPLTINSTRUCT lt_input, LPLTOUTSTRUCT answers)
{
	unsigned int bytes = 0;
	int i, j;
	
	gi_iteration_mode = lt_input->iteration_mode;

	if (lt_input->fp) { /* input comes from a file */
		sp = sn = (char *) malloc (BUFFER_LENGTH);
		bytes = fread (sn, 1, BUFFER_LENGTH_LESS_1, lt_input->fp);
		sn [bytes] = '\0';
	} else {
		if (logfp) fprintf( logfp, "No input available for lt simulation\n");
		oe_exit (ERR_BUFFER_MISSING);
	}
	op = lt_input->op; /* text output */
	bp = lt_input->bp; /* plot file */
	if (lt_input->stop_on_flashover) {
		flash_halt_enabled = TRUE;
	} else {
		flash_halt_enabled = FALSE;
	}
	if (gi_iteration_mode == ONE_SHOT) {
		want_si_calculation = TRUE;
	} else {
		want_si_calculation = TRUE; /* need SI when iterating for critical current */
	}
/* set up head pointers for the linked lists */
	(void) init_surge_list ();
	(void) init_source_list ();
	(void) init_meter_list ();
	(void) init_pole_list ();
	(void) init_span_list ();
	(void) init_line_list ();
	(void) init_ground_list ();
	(void) init_resistor_list ();
	(void) init_inductor_list ();
	(void) init_capacitor_list ();
	(void) init_customer_list ();
	(void) init_insulator_list ();
	(void) init_arrester_list ();
	(void) init_pipegap_list ();
	(void) init_lpm_list ();
	(void) init_arrbez_list ();
	(void) init_steepfront_list ();
/* read input from either the file or the memory buffer */
	(void) readfile ();
	if (op && (gi_iteration_mode == ONE_SHOT)) { /* DOS only */
		fprintf (op, "  N   span     dT   Tmax\n");
		fprintf (op, "%3d %6.2f %.4g %.4g\n", 
			number_of_nodes, span_length, dT, Tmax);
		fprintf (op, "Z-phase\n");
		for (i = 0; i < number_of_nodes; i++) {
			for (j = 0; j < number_of_nodes; j++) {
				fprintf (op, " %14.5e", gsl_matrix_get (span_head->Zp, i, j));
			}
			fprintf (op, "\n");
		}
		fprintf (op, "Z-modal\n");
		for (i = 0; i < number_of_nodes; i++) {
			for (j = 0; j < number_of_nodes; j++) {
				fprintf (op, " %14.5e", gsl_matrix_get (span_head->Zm, i, j));
			}
			fprintf (op, "\n");
		}
		fprintf (op, "Modal Transformation (Ti)\n");
		for (i = 0; i < number_of_nodes; i++) {
			for (j = 0; j < number_of_nodes; j++) {
				fprintf (op, " %14.5e", gsl_matrix_get (span_head->Ti, i, j));
			}
			fprintf (op, "\n");
		}
	}
/* complete set-up of the system model, with past history currents for
storage elements */
	if (!using_network) {
		connect_lines ();
	}
	do_all_lines (init_line_history);
	do_all_inductors (init_inductor_history);
	do_all_capacitors (init_capacitor_history);
/* make sure we solve for voltages at each end of the circuit, and add
surge impedance terminations */
	pole_ptr = find_pole (1);
	if (pole_ptr) {
		pole_ptr->solve = TRUE;
		if (left_end_z) {
			terminate_pole (pole_ptr, span_head);
		}
	}
	pole_ptr = find_pole (number_of_poles);
	if (pole_ptr) {
		pole_ptr->solve = TRUE;
		if (right_end_z) {
			terminate_pole (pole_ptr, span_head);
		}
	}
	if (op && (gi_iteration_mode == ONE_SHOT)) {
		do_all_sources (print_source_data);
	}
/* perform initial y matrix factoring at each pole - now ready to start */
	do_all_poles (triang_pole);
	Tmax += 0.5 * dT;
/* there are three running modes for the transient simulation: */
	if (gi_iteration_mode == FIND_CRITICAL_CURRENT) {
/* critical flashover current iteration - as called by driver */
		if (logfp) fprintf( logfp, "lt in location control mode\n");
		loop_control (lt_input, answers);
	} else {
/* single-shot run, as called by the DOS version */
		if (logfp) fprintf( logfp, "lt in stand-alone mode\n");
		time_step_loops (answers);
	}
	if (op) { /* print results, DOS only */
		if (logfp) fprintf( logfp, "\n");
		if (gi_iteration_mode == ONE_SHOT) {
			do_all_meters (print_meter_data);
			do_all_insulators (print_insulator_data);
			do_all_lpms (print_lpm_data);
			do_all_customers (print_customer_data);
			do_all_arresters (print_arrester_data);
			do_all_arrbezs (print_arrbez_data);
			do_all_pipegaps (print_pipegap_data);
		} else {
			fprintf (op, "\nAverage Critical Currents, Poles %d to %d\n", 
				lt_input->first_pole_hit, lt_input->last_pole_hit);
			for (i = 0; i < MAX_WIRES_HIT; ++i) {
				if (lt_input->wire_struck[i] > 0) {
					fprintf (op, "wire %2d: %4e\n", i+1, answers->icritical[i]); 
				}
			}
		}
	}
	if (op) {
		fprintf (op, "nr_iter = %ld, nr_max = %d\n", nr_iter, nr_max);
	} else if (logfp) {
		fprintf (logfp, "nr_iter = %ld, nr_max = %d\n", nr_iter, nr_max);
	}
	(void) cleanup ();
	return (0);
}

void run_loop_case (int pole_number, int wire_number, double i_pk, double ftf, double ftt, 
					LPLTOUTSTRUCT answers)
{
	reset_system ();  /* reset initial conditions for each simulation */
	surge_ptr = surge_head->next;
	steepfront_ptr = steepfront_head->next;
/* re-set the surge component parameters */
	if (surge_ptr) {
		move_surge (surge_ptr, pole_number, wire_number, 0, i_pk, ftf, ftt, 0.0);
	} else if (steepfront_ptr) {
		move_steepfront (steepfront_ptr, pole_number, wire_number, 
			0, i_pk, ftf, ftt, 0.0, steepfront_ptr->pu_si);
	}
/* run the transient simulation */
	time_step_loops (answers);
}

double icrit_function (double i_pk, void *params)
{
	struct icrit_params *p = (struct icrit_params *) params;
	double ftt = Q_MEDIAN_FIRST / I_MEDIAN_FIRST / 1000.0 / ETKONST; 
	double ftf = 1.0e-6 * T3090_FIRST;
	double ret;

	run_loop_case (p->pole_number, p->wire_number, i_pk, ftf, ftt, p->answers);
	ret = p->answers->SI - 1.0;
	if (ret >= 0.0) {
		ret += (Tmax - t) * 1.0e5;
	}

	return ret;
}

/* this function simulates a stroke to each pole and exposed wire,
at each histogram midpoint value */

void loop_control (LPLTINSTRUCT lt_input, LPLTOUTSTRUCT answers)
{
	int wire_idx, pole_number, wire_number;
	int total_cases, case_number, wires_hit;
	int insulators_at_one_pole, first_ins_pole;
	double num_poles;
	int has_arresters;

	struct icrit_params params;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	gsl_function F;
	double i_pk, i_lo, i_hi;
	int status, iter;

/* set up the GSL root finder */
	F.function = &icrit_function;
	F.params = &params;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);

/* zero out the answer arrays */
	wires_hit = 0;
	for (wire_idx = 0; wire_idx < MAX_WIRES_HIT; wire_idx++) {
		answers->icritical[wire_idx] = 0.0;
		if (lt_input->wire_struck[wire_idx] > 0) ++wires_hit;
	}
	
	has_arresters = FALSE;
	if (arrbez_head->next) has_arresters = TRUE;
	if (arrester_head->next) has_arresters = TRUE;
	if (logfp) fprintf (logfp, "has_arresters = %d\n", has_arresters);

	total_cases = (lt_input->last_pole_hit - lt_input->first_pole_hit - 1) * wires_hit;
	num_poles = lt_input->last_pole_hit - lt_input->first_pole_hit + 1.0;

	case_number = 0;
	params.answers = answers;
/* check all of the requested poles */
	for (pole_number = lt_input->first_pole_hit; pole_number <= lt_input->last_pole_hit; pole_number++) {
		params.pole_number = pole_number;
/*  if there are insulators at just one pole, we want to move them with
    the surge.  If insulators at more than one pole, leave them in place.
    First, look through the insulators for presence of different poles: */
		insulators_at_one_pole = TRUE;
		first_ins_pole = 0;
		insulator_ptr = insulator_head;
		while ((insulator_ptr = insulator_ptr->next) != NULL) {
			if (first_ins_pole == 0) {
				first_ins_pole =
					insulator_ptr->parent->location;
			}
			if (first_ins_pole!=insulator_ptr->parent->location) {
				insulators_at_one_pole = FALSE;
			}
		}
		lpm_ptr = lpm_head;
		while ((lpm_ptr = lpm_ptr->next) != NULL) {
			if (first_ins_pole == 0) {
				first_ins_pole = lpm_ptr->parent->location;
			}
			if (first_ins_pole!=lpm_ptr->parent->location) {
				insulators_at_one_pole = FALSE;
			}
		}
/*  Now move the insulators if only one insulator pole was found */
		if (insulators_at_one_pole == TRUE) {
			insulator_ptr = insulator_head;
			while ((insulator_ptr = insulator_ptr->next) != NULL) {
				move_insulator (insulator_ptr, pole_number);
			}
			lpm_ptr = lpm_head;
			while ((lpm_ptr = lpm_ptr->next) != NULL) {
				move_lpm (lpm_ptr, pole_number);
			}
		}
/* check all of the exposed conductors at this pole.  The set of exposed
conductors was originally determined in egm, passed in by driver */
		for (wire_idx = 0; wire_idx < MAX_WIRES_HIT; wire_idx++) {
			if (lt_input->wire_struck[wire_idx] > 0) {
				wire_number = wire_idx + 1;
				params.wire_number = wire_number;
				iter = 0;
				status = GSL_SUCCESS;
				if (icrit_function (MIN_STROKE, &params) >= 0.0) { /* always have a flashover */
					answers->icritical[wire_idx] += (MIN_STROKE / num_poles);
				} else if (icrit_function (MAX_STROKE, &params) <= 0.0) { /* never have a flashover */
					answers->icritical[wire_idx] += (MAX_STROKE / num_poles);
				} else { /* iterate for critical current */
					gsl_root_fsolver_set (s, &F, MIN_STROKE, MAX_STROKE);
					do {
						++iter;
						status = gsl_root_fsolver_iterate (s);
						i_pk = gsl_root_fsolver_root (s);
						i_lo = gsl_root_fsolver_x_lower (s);
						i_hi = gsl_root_fsolver_x_upper (s);
						status = gsl_root_test_interval (i_lo, i_hi, ITER_TOL, 0.0);
						if (status == GSL_SUCCESS) {
							answers->icritical[wire_idx] += (i_pk / num_poles);
						}
					} while (status == GSL_CONTINUE && iter < MAX_ITER);
				}
				++case_number;
				if (logfp) {
					fprintf (logfp, "case %d, pole %d, wire %d, i_pk = %G, ftf = %G, ftt = %G, SI = %G, Energy = %G, iter = %d, status = %d\n",
						case_number, pole_number, wire_number, 0.001 * answers->icritical[wire_idx], T3090_FIRST, 
						1000.0 * Q_MEDIAN_FIRST / I_MEDIAN_FIRST / ETKONST, 
						answers->SI, answers->energy, iter, status);
					fflush (logfp);
				}
			}
		}
	}
	gsl_root_fsolver_free (s);
}

/* run a complete simulation, assuming the initial conditions have been
set properly */

void time_step_loops (LPLTOUTSTRUCT answers)
{
	int last_pct = 0;

	t = 0.0;
	step = 0;
	flash_halt = FALSE;
	if (op) {
//		printf ("%le <- Tmax\n", Tmax);
//		printf ("%le <- Time Now\r", t);
	}
	if (bp) {  /* initialize plot file */
		InitializePlotOutput (meter_head, dT, Tmax);
	}
	do_all_monitors (find_monitor_links);
	do { /* keep going till we hit Tmax, or an insulator flashes over when flash_halt_enabled */
		solution_valid = FALSE;
		while (!solution_valid) { /* get a valid solution for this step - no arrester state changes */
/* solve for voltages at this step */
			do_all_poles (zero_pole_injection);
			do_all_surges (inject_surge);
			do_all_steepfronts (inject_steepfront);
			do_all_sources (inject_source);
			do_all_grounds (inject_ground);
			if (using_multiple_span_defns) {
				do_all_lines (inject_line_iphase);
			} else {
				do_all_lines (inject_line_imode);
				do_all_poles (inject_pole_imode);
			}
			do_all_arresters (inject_arrester);
			do_all_pipegaps (inject_pipegap);
			do_all_inductors (inject_inductor_history);
			do_all_capacitors (inject_capacitor_history);
			do_all_poles (triang_pole);
			do_all_poles (solve_pole);
			solution_valid = TRUE; /* see if an arrester changed state - need to resolve */
			do_all_arresters (check_arrester);
			do_all_pipegaps (check_pipegap);
		}
/* update the non-linear and energy-storage history terms for the next step */
		do_all_grounds (check_ground);
		do_all_insulators (check_insulator);  /* may set flash_halt */
		do_all_lpms (check_lpm);
		do_all_inductors (update_inductor_history);
		do_all_arresters (update_arrester_history);
		do_all_arrbezs (update_arrbez_history);
		do_all_capacitors (update_capacitor_history);
		do_all_customers (update_customer_history);
		if (using_multiple_span_defns) {
			do_all_lines (update_vmode_and_history);
		} else {
			do_all_poles (calc_pole_vmode);
			do_all_lines (update_line_history);
		}
		if (bp) {  
			WritePlotTimeStep (meter_head, t);
		} else {
			do_all_meters (update_meter_peaks);
		}
		do_all_monitors (update_monitor_pts);
#ifdef LOG_POLES_AND_LINES
		do_all_poles (print_pole_data);
		do_all_lines (print_line_history);
#endif
		if (op) {  /* progress report */
			if (step % 10 == 0) {
//				printf ("%le\r", t);
			}
		}
		if (using_second_dT && !dT_switched) {
			if (t >= dT_switch_time) {
				change_time_step();
			}
		}
		t += dT;  /* advance the time step */
		++step;
	} while (t <= Tmax && !flash_halt);
	if (logfp) fprintf( logfp, "\n");
/* set up "quick answers" for DOS version */
	SI = energy = charge = current = predischarge = 0.0;
	do_all_arresters (arrester_answers_cleanup);
	do_all_pipegaps (pipegap_answers_cleanup);
	do_all_insulators (insulator_answers_cleanup);
	do_all_lpms (lpm_answers_cleanup);
	do_all_arrbezs (arrbez_answers_cleanup);
	answers->SI = SI; /* severity index for insulator closest to flashing over */
	answers->energy = energy;  /* highest arrester discharge parameters */
	answers->charge = charge;
	answers->current = current;
	answers->predischarge = predischarge;
	do_all_monitors (update_monitor_summary);
	if (bp) {
		FinalizePlotHeader (t, step);
	}
}

