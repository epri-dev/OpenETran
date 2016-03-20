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
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "../OETypes.h"
#include "../Parser.h"
#include "../ReadUtils.h"
#include "Meter.h"
#include "../WritePlotFile.h"
#include "../ChangeTimeStep.h"
#include "Pole.h"
#include "Line.h"

char span_token[] = "span";
char line_token[] = "line";
char end_token[] = "end";
char conductor_token[] = "conductor";
char node_token[] = "node";
char cable_token[] = "cable";

struct line *line_head, *line_ptr;
struct span *span_head, *span_ptr;

/* supervisory function to add all the line sections between poles */
 /* only for non-network systems */

void connect_lines (void)
{
	int left_pole, right_pole, travel_steps;

/* round the span length to nearest integer number of time steps */
	travel_steps = (int) (0.5 + span_length / span_head->wave_velocity / dT);
	if (using_second_dT) {
		second_dT = dT * travel_steps;
	}
	for (left_pole = 1; left_pole < number_of_poles; ++left_pole) {
		right_pole = left_pole + 1;
		insert_line (left_pole, right_pole, span_head, travel_steps);
	}
}

/* build a line, and connect it to poles at each end.  Usually,
right_pole == left_pole + 1 */

void insert_line (int left_pole, int right_pole, struct span *defn, int travel_steps)
{
	struct line *ptr;

	if (((ptr = (struct line *) malloc (sizeof *ptr)) != NULL)) {
		ptr->left = find_pole (left_pole); /* terminal pole pointers */
		if (!ptr->left) oe_exit (ERR_BAD_POLE);
		ptr->right = find_pole (right_pole);
		if (!ptr->right) oe_exit (ERR_BAD_POLE);
		ptr->steps = ptr->alloc_steps = travel_steps;
		ptr->defn = defn;
		if (!(ptr->hist_left = gsl_matrix_calloc (number_of_conductors, travel_steps))) {
			if (logfp) fprintf( logfp, "can't allocate history space\n");
			oe_exit (ERR_MALLOC);
		}
		if (!(ptr->hist_right = gsl_matrix_calloc (number_of_conductors, travel_steps))) {
			if (logfp) fprintf( logfp, "can't allocate history space\n");
			oe_exit (ERR_MALLOC);
		}
/* add surge impedances to the terminal poles */
		gsl_matrix_add (ptr->left->Ybus, defn->Yp);
		gsl_matrix_add (ptr->right->Ybus, defn->Yp);
		ptr->next = NULL;
		line_ptr->next = ptr;
		line_ptr = ptr;
	} else {
		if (logfp) fprintf( logfp, "can't allocate new line\n");
		oe_exit (ERR_MALLOC);
	}
}

/* past history currents start at zero, unless there is a trapped charge */
void init_line_history (struct line *ptr)
{
	int i, j;
	double idc;
	
	for (i = 0; i < number_of_conductors; i++) {
		 /* dc current to maintain initial voltage in modal coordinates */
		idc = -gsl_matrix_get (ptr->defn->Ym, i, i) * gsl_vector_get (ptr->defn->vm, i);
		for (j = 0; j < ptr->steps; j++) {
			gsl_matrix_set (ptr->hist_left, i, j, idc);
			gsl_matrix_set (ptr->hist_right, i, j, idc);
		}
	}
}

/* only for network systems */
/* converts line history currents to phase coordinates, adds to pole injections */
void inject_line_iphase (struct line *ptr)
{
	gsl_vector *im;
	gsl_matrix *h, *Ti;
	int i, k;
	gsl_vector_view ip;

	k = step % ptr->steps;  /* cycle through the past history array in circular fashion */
	Ti = ptr->defn->Ti;

	im = ptr->left->imode;  /* add to left pole */
	ip = gsl_vector_subvector (ptr->left->injection, 1, number_of_nodes);
	h = ptr->hist_left;
	for (i = 0; i < number_of_conductors; i++) {
		gsl_vector_set (im, i, -gsl_matrix_get (h, i, k));    /* using pole's storage buffer, so don't accumulate imode */
	}
	gsl_blas_dgemv (CblasNoTrans, 1.0, Ti, im, 1.0, &ip.vector);

	im = ptr->right->imode;  /* add to right pole */
	ip = gsl_vector_subvector (ptr->right->injection, 1, number_of_nodes);
	h = ptr->hist_right;
	for (i = 0; i < number_of_conductors; i++) {
		gsl_vector_set (im, i, -gsl_matrix_get (h, i, k));    /* using pole's storage buffer, so don't accumulate imode */
	}
	gsl_blas_dgemv (CblasNoTrans, 1.0, Ti, im, 1.0, &ip.vector);
}

/* only for network systems */
/* convert pole phase voltages to modes for this line's span, update history currents */
void update_vmode_and_history (struct line *ptr)
{
	gsl_vector *vl, *vr;
	gsl_matrix *Tvt, *hl, *hr;
	int i, k;
	double y, irl, ilr;
	gsl_vector_view vp_left, vp_right;
	
	k = step % ptr->steps;
	Tvt = ptr->defn->Tvt;

	vp_left = gsl_vector_subvector (ptr->left->voltage, 1, number_of_conductors);
	vl = ptr->left->vmode;
	hl = ptr->hist_left;

	vp_right = gsl_vector_subvector (ptr->right->voltage, 1, number_of_conductors);
	vr = ptr->right->vmode;
	hr = ptr->hist_right;

/* calculate vmode at each end, using pole's local storage */
	gsl_blas_dgemv (CblasNoTrans, 1.0, Tvt, &vp_left.vector, 0.0, vl);
	gsl_blas_dgemv (CblasNoTrans, 1.0, Tvt, &vp_right.vector, 0.0, vr);

/* update line history at each end */
	for (i = 0; i < number_of_conductors; i++) {
		y = gsl_matrix_get (ptr->defn->Ym, i, i);
		ilr = gsl_vector_get (vl, i) * y + gsl_matrix_get (hl, i, k);
		irl = gsl_vector_get (vr, i) * y + gsl_matrix_get (hr, i, k);
		gsl_matrix_set (hl, i, k, -gsl_vector_get (vr, i) * y - irl);
		gsl_matrix_set (hr, i, k, -gsl_vector_get (vl, i) * y - ilr);
	}
}

/* add the modal current injections to each terminal pole.  These are converted
  to phase coordinates later */
/* only for non-network systems */
void inject_line_imode (struct line *ptr)
{
	gsl_vector *c;
	gsl_matrix *h;
	int i, k;
	
	k = step % ptr->steps;  /* cycle through the past history array in circular fashion */
	c = ptr->left->imode;  /* add to left pole */
	h = ptr->hist_left;
	for (i = 0; i < number_of_conductors; i++) {
		*gsl_vector_ptr (c, i) -= gsl_matrix_get (h, i, k);
	}
	c = ptr->right->imode;  /* add to right pole */
	h = ptr->hist_right;
	for (i = 0; i < number_of_conductors; i++) {
		*gsl_vector_ptr (c, i) -= gsl_matrix_get (h, i, k);
	}
}

/* calculate the line past history currents in modal coordinates, based on
the solved pole voltages in modal coordinates */
 /* only for non-network systems */
void update_line_history (struct line *ptr)
{
	gsl_vector *vl, *vr;
	gsl_matrix *hl, *hr;
	double y, irl, ilr;
	int i, k;
	
	k = step % ptr->steps;
	vl = ptr->left->vmode;
	vr = ptr->right->vmode;
	hl = ptr->hist_left;
	hr = ptr->hist_right;
	for (i = 0; i < number_of_conductors; i++) {
		y = gsl_matrix_get (ptr->defn->Ym, i, i);
		ilr = gsl_vector_get (vl, i) * y + gsl_matrix_get (hl, i, k);
		irl = gsl_vector_get (vr, i) * y + gsl_matrix_get (hr, i, k);
		gsl_matrix_set (hl, i, k, -gsl_vector_get (vr, i) * y - irl);
		gsl_matrix_set (hr, i, k, -gsl_vector_get (vl, i) * y - ilr);
	}
}

void print_line_history (struct line *ptr)
{
	int i, j;
	if (op) {
		fprintf (op, "line %d-%d, step %d, t = %g\n", ptr->left->location, ptr->right->location, step, t);
		fprintf (op, "\tHist Left\n");
		for (i = 0; i < number_of_conductors; i++) {
			fprintf (op, "\t");
			for (j = 0; j < ptr->alloc_steps; j++) {
				fprintf (op, " %14.5e", gsl_matrix_get (ptr->hist_left, i, j));
			}
			fprintf (op, "\n");
		}
		fprintf (op, "\tHist Right\n");
		for (i = 0; i < number_of_conductors; i++) {
			fprintf (op, "\t");
			for (j = 0; j < ptr->alloc_steps; j++) {
				fprintf (op, " %14.5e", gsl_matrix_get (ptr->hist_right, i, j));
			}
			fprintf (op, "\n");
		}
	}
}

int init_span_list (void)
{
    if ((span_head = (struct span *) malloc (sizeof *span_head))) {
        span_head->next = NULL;
        span_head->Zm = NULL;
        span_head->Ym = NULL;
        span_head->Zp = NULL;
        span_head->Yp = NULL;
        span_head->Ti = NULL;
        span_head->Tit = NULL;
        span_head->Tv = NULL;
        span_head->Tvt = NULL;
        span_head->vm = NULL;
        span_head->wave_velocity = LIGHT;
        span_head->span_id = 1;
        span_ptr = span_head;
        return (0);
    }
    if (logfp) fprintf (logfp, "can't initialize span list\n");
    oe_exit (ERR_MALLOC);
    return (1);
}

void do_all_spans (void (*verb) (struct span *))
{
    span_ptr = span_head;  /* unlike other lists, span_head is used */
    while (span_ptr) {
        verb (span_ptr);
        span_ptr = span_ptr->next;
    }
}

/* set initial voltages in modal coordinates, based on phase-coordinate
voltages (which may have had a power-frequency initial condition added) */

void reset_span (struct span *ptr)
{
	gsl_blas_dgemv (CblasNoTrans, 1.0, ptr->Tvt, ptr->vp_offset, 0.0, ptr->vm); 
}

struct span *find_span (int span_id)
{
    span_ptr = span_head;
    while (span_ptr) {
        if (span_ptr->span_id == span_id) {
            return (span_ptr);
        }
        span_ptr = span_ptr->next;
    }
    return (NULL);
}

void read_spans (void)
{
    struct span *ptr;
    char *p;
    int span_id;

    while (((p = first_token()) != NULL) && !strcmp (p, span_token)) {
        next_int (&span_id);
        ptr = find_span (span_id);
        if (!ptr) {
            if (((ptr = (struct span *) malloc (sizeof *ptr)) != NULL)) {
                ptr->Zm = NULL;
                ptr->Ym = NULL;
                ptr->Zp = NULL;
                ptr->Yp = NULL;
                ptr->Ti = NULL;
                ptr->Tit = NULL;
                ptr->Tv = NULL;
                ptr->Tvt = NULL;
                ptr->vm = NULL;
				ptr->vp_offset = NULL;
                ptr->wave_velocity = LIGHT;
                ptr->span_id = span_id;
                ptr->next = NULL;
                span_ptr = span_head;
                while (span_ptr->next) {
                    span_ptr = span_ptr->next;
                }
                span_ptr->next = ptr;
                span_ptr = ptr;
                using_multiple_span_defns = TRUE;
            } else {
                if (logfp) fprintf (logfp, "can't allocate new span\n");
                oe_exit (ERR_MALLOC);
            }
        }
        read_conductors (ptr);
    }
}

void read_lines (void)
{
    int from, to, term_left, term_right, span_id, line_steps;
    double length;
    char *p;
    struct pole *left, *right;
    struct line *ptr;
    struct span *defn;

    while (((p = first_token()) != NULL) && !strcmp (p, line_token)) {
        next_int (&from);
        next_int (&to);
        next_int (&span_id);
        next_double (&length);
        next_int (&term_left);
        next_int (&term_right);
        defn = find_span (span_id);
        line_steps = (int) (0.5 + length / defn->wave_velocity / dT);
        if (from > number_of_poles) {
            number_of_poles = from;
        }
        if (to > number_of_poles) {
            number_of_poles = to;
        }
        if (NULL == (left = find_pole (from))) {
            pole_ptr = pole_head;
            while (pole_ptr->next) {
                pole_ptr = pole_ptr->next;
            }
            left = new_pole (from);
        }
        if (NULL == (right = find_pole (to))) {
            pole_ptr = pole_head;
            while (pole_ptr->next) {
                pole_ptr = pole_ptr->next;
            }
            right = new_pole (to);
        }
        left->solve = TRUE;
        right->solve = TRUE;
        if (((ptr = (struct line *) malloc (sizeof *ptr)) != NULL)) {
            ptr->left = left;
            ptr->right = right;
            ptr->steps = ptr->alloc_steps = line_steps;
            ptr->defn = defn;
            if (!(ptr->hist_left = gsl_matrix_calloc (number_of_conductors, line_steps))) {
                if (logfp) fprintf( logfp, "can't allocate history space\n");
                oe_exit (ERR_MALLOC);
            }
            if (!(ptr->hist_right = gsl_matrix_calloc (number_of_conductors, line_steps))) {
                if (logfp) fprintf( logfp, "can't allocate history space\n");
                oe_exit (ERR_MALLOC);
            }
/* add surge impedances to the terminal poles */
            gsl_matrix_add (left->Ybus, defn->Yp);
            gsl_matrix_add (right->Ybus, defn->Yp);
            if (term_left) {
                terminate_pole (left, defn);
            }
            if (term_right) {
                terminate_pole (right, defn);
            }
            ptr->next = NULL;
            line_ptr->next = ptr;
            line_ptr = ptr;
        } else {
            if (logfp) fprintf (logfp, "can't allocate new line\n");
            oe_exit (ERR_MALLOC);
        }
    }
}

int read_cables (struct span *defn)
{
	int i;
	double z_surge, v_prop, vpf;
	
	(void) next_int (&i);
	if (i > number_of_nodes || i < 0) {
		if (logfp) fprintf (logfp, "bad cable number %d\n", i);
		oe_exit (ERR_CONDUCTOR_N);
	}
	(void) next_double (&z_surge);  /* read Z and V directly */
	(void) next_double (&v_prop);
	(void) next_double (&vpf);
	defn->wave_velocity = v_prop;
/* the matrix definitions are real simple for uncoupled case */
	--i; // zero-based indexing
	gsl_matrix_set (defn->Ti, i, i, 1.0);
	gsl_matrix_set (defn->Tit, i, i, 1.0);
	gsl_matrix_set (defn->Tv, i, i, 1.0);
	gsl_matrix_set (defn->Tvt, i, i, 1.0);
	gsl_matrix_set (defn->Zp, i, i, z_surge);
	gsl_matrix_set (defn->Zm, i, i, z_surge);
	gsl_matrix_set (defn->Yp, i, i, 1.0 / z_surge);
	gsl_matrix_set (defn->Ym, i, i, 1.0 / z_surge);
	gsl_vector_set (defn->vp_offset, i, vpf);
	return (0);
}

/*  &&&&  read conductors and set up modal transformation   */

static void transform_conductors (struct span *defn, gsl_vector *x, gsl_vector *y, gsl_vector *r, gsl_vector *v)
{
	gsl_matrix *zcopy;
	gsl_vector *lambda;
	gsl_eigen_symmv_workspace *w;
	gsl_permutation *p;
	int signum;
	double dx, dy, hs, yi, yj, xi, xj, ri;
	int n = number_of_conductors;
	int i, j;

/* temp matrices for eigenvalue solution */
	zcopy = gsl_matrix_calloc (n, n);
	lambda = gsl_vector_calloc (n);

	gsl_vector_memcpy (defn->vp_offset, v);

/*   simple error checking on the input  */
	for (i = 0; i < n; i++) {
		if (gsl_vector_get (r, i) <= 0.0) {
			if (logfp) fprintf( logfp, "bad radius: %le on conductor %d\n", gsl_vector_get (r, i), i);
			oe_exit (ERR_RADIUS);
		}
		if (gsl_vector_get (y, i) <= 0.0) {
			if (logfp) fprintf( logfp, "bad height: %le on conductor %d\n", gsl_vector_get (y, i), i);
			oe_exit (ERR_HEIGHT);
		}
	}

/*   calculate self and mutual surge impedances in phase coordinates   */

	for(i = 0; i < n; i++) {
		xi = gsl_vector_get (x, i);
		yi = gsl_vector_get (y, i);
		ri = gsl_vector_get (r, i);
		gsl_matrix_set (defn->Zp, i, i, 60.0 * log (2.0 * yi / ri));
		for(j = (i+1); j < n; j++) {
			xj = gsl_vector_get (x, j);
			yj = gsl_vector_get (y, j);
			dx = xi - xj;
			dy = yi - yj;
			if (fabs (dx) < 0.001 && fabs (dy) < 0.001) {
				if (logfp) fprintf( logfp, "wires %d and %d overlap\n",	i, j);
				oe_exit (ERR_OVERLAP);
			}
			hs = yi + yj;
			gsl_matrix_set (defn->Zp, i, j, 60.0*log(sqrt(dx*dx + hs*hs) 
				/ sqrt(dx*dx + dy*dy)));
		}
	}

/*   fill in the lower triangle     */

	if (n>1) {
		for(i = 1; i < n; i++) {
			for (j = 0; j<i; j++) {
				gsl_matrix_set (defn->Zp, i, j, gsl_matrix_get (defn->Zp, j, i));
			}
		}
	}

/*   solve eigenvalue problem on a copy of Zp   */

	gsl_matrix_memcpy (zcopy, defn->Zp);
	w = gsl_eigen_symmv_alloc (n);
	gsl_eigen_symmv (zcopy, lambda, defn->Ti, w);
	gsl_eigen_symmv_sort (lambda, defn->Ti, GSL_EIGEN_SORT_VAL_ASC);
	gsl_eigen_symmv_free (w);
	gsl_matrix_transpose_memcpy (defn->Tit, defn->Ti);

	p = gsl_permutation_alloc (n);
	gsl_matrix_memcpy (zcopy, defn->Ti);
	gsl_linalg_LU_decomp (zcopy, p, &signum);
	gsl_linalg_LU_invert (zcopy, p, defn->Tvt);
	gsl_matrix_transpose_memcpy (defn->Tv, defn->Tvt);

/*   calculate Zm = Tit * Zp * Ti, and then Ym and Yp    */

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, defn->Zp, defn->Ti, 0.0, zcopy);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, defn->Tit, zcopy, 0.0, defn->Zm);
	gsl_matrix_memcpy (zcopy, defn->Zp);
	gsl_linalg_LU_decomp (zcopy, p, &signum);
	gsl_linalg_LU_invert (zcopy, p, defn->Yp);
	gsl_matrix_memcpy (zcopy, defn->Zm);
	gsl_linalg_LU_decomp (zcopy, p, &signum);
	gsl_linalg_LU_invert (zcopy, p, defn->Ym);
	gsl_permutation_free (p);

/*  release temporary matrices and vectors for modal transformation  */
  
	gsl_matrix_free (zcopy);
	gsl_vector_free (lambda);
}

void allocate_definition_memory (struct span *defn, int n)
{
	defn->Ti = gsl_matrix_calloc (n, n);
	defn->Tit = gsl_matrix_calloc (n, n);
	defn->Tv = gsl_matrix_calloc (n, n);
	defn->Tvt = gsl_matrix_calloc (n, n);
	defn->Zm = gsl_matrix_calloc (n, n);
	defn->Ym = gsl_matrix_calloc (n, n);
	defn->Zp = gsl_matrix_calloc (n, n);
	defn->Yp = gsl_matrix_calloc (n, n);
	defn->vm = gsl_vector_calloc (n);
	defn->vp_offset = gsl_vector_calloc (n);

/*   assume waves travel at speed of light on overhead lines    */
	defn->wave_velocity = LIGHT;
}

int read_conductors (struct span *defn)
{
	gsl_vector *height, *radius, *location, *voltage;
	int i, n;
	int using_geometry;
	int count_conductors;
	char *p;

	n = number_of_nodes;
	if (n <= 0) {
		if (logfp) fprintf( logfp, "bad n = %d in read_conductors\n", n);
		oe_exit (ERR_NPHASES);
	}
/* conductor geometry data */
	height = gsl_vector_calloc (n);
	radius = gsl_vector_calloc (n);
	location = gsl_vector_calloc (n); /* horizontal displacement from pole centerline */
	voltage = gsl_vector_calloc (n);

/*  read n conductors from fp   */

	using_geometry = TRUE;
	n = 0;
	count_conductors = 0;

	/* read up to number_of_phases conductor cards, till reaching "end" token */
	while (n < number_of_nodes && (NULL != (p = first_token())) && strcmp (p, end_token)) {
		n++;
		if (!strcmp (p, cable_token)) {
			using_geometry = FALSE;
			if (!defn->Ti) {
				allocate_definition_memory (defn, number_of_nodes);
			}
			read_cables (defn);
		} else if (!strcmp (p, conductor_token)) {
			if (!using_geometry) {
				if (logfp) fprintf (logfp, "can't mix cable and conductor input in the same span\n");
				oe_exit (ERR_MIXED_LINES);
			}
/* read the geometry and initial voltage for each overhead conductor */
			next_int (&i);
			if (i > 0 && i <= number_of_nodes) {
				++count_conductors;
				next_double (gsl_vector_ptr (height, i-1));
				next_double (gsl_vector_ptr (location, i-1));
				next_double (gsl_vector_ptr (radius, i-1));
				next_double (gsl_vector_ptr (voltage, i-1));
			} else {
				if (logfp) fprintf( logfp, "bad conductor number %d\n", i);
				oe_exit (ERR_CONDUCTOR_N);
			}
		} else if (!strcmp (p, node_token)) {
			if (!using_geometry) {
				if (logfp) fprintf (logfp, "can't mix cable and conductor input in the same span\n");
				oe_exit (ERR_MIXED_LINES);
			}
			next_int (&i);
			if (i > 0 && i <= number_of_nodes) {
			} else {
				if (logfp) fprintf( logfp, "bad node number %d\n", i);
				oe_exit (ERR_CONDUCTOR_N);
			}
		}
	}

	if (using_geometry) {
		number_of_conductors = count_conductors;
		if (count_conductors < 1) {
			if (logfp) fprintf( logfp, "input file lt.dat missing conductors\n");
			oe_exit (ERR_MISSING_CONDUCTOR);
		} else if (number_of_conductors > number_of_nodes) {
			if (logfp) fprintf( logfp, "input file lt.dat has too many conductors\n");
			oe_exit (ERR_CABLE_PHASES);
		}

		if (!defn->Ti) {
			allocate_definition_memory (defn, number_of_conductors);
		}

		transform_conductors (defn, location, height, radius, voltage);
	}

	gsl_vector_free (height);
	gsl_vector_free (radius);
	gsl_vector_free (location);
	gsl_vector_free (voltage);

	return (0);
}

/* set initial voltages in modal coordinates, based on phase-coordinate
voltages (which may have had a power-frequency initial condition added) */

void reset_lines (void)
{
	do_all_spans (reset_span);
}

void do_all_lines (void (*verb) (struct line *))
{
	line_ptr = line_head;
	while (((line_ptr = line_ptr->next) != NULL)) {
		verb (line_ptr);
	}
}

int init_line_list (void)
{
	if (((line_head = (struct line *) malloc (sizeof *line_head)) != NULL)) {
		line_head->next = NULL;
		line_head->defn = NULL;
		line_head->hist_left = NULL;
		line_head->hist_right = NULL;
		line_ptr = line_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize line list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}
