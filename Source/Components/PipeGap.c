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

#include <wtypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
#include "PipeGap.h"

char pipegap_token[] = "pipegap";

struct pipegap *pipegap_head, *pipegap_ptr;

int init_pipegap_list (void)
{
    if (pipegap_head = (struct pipegap *) malloc (sizeof *pipegap_head)) {
        pipegap_head->next = NULL;
        pipegap_ptr = pipegap_head;
        return (0);
    }
    if (logfp) fprintf (logfp, "can't initialize pipegap list\n");
    oe_exit (ERR_MALLOC);
    return (1);
}

void do_all_pipegaps (void (*verb) (struct pipegap *))
{
    pipegap_ptr = pipegap_head;
    while (pipegap_ptr = pipegap_ptr->next) {
        verb (pipegap_ptr);
    }
}

struct pipegap *find_pipegap (int at, int from, int to)
{
	pipegap_ptr = pipegap_head;
	while ((pipegap_ptr = pipegap_ptr->next) != NULL) {
		if ((pipegap_ptr->parent->location == at) &&
			(pipegap_ptr->from == from) &&(pipegap_ptr->from == from)) return pipegap_ptr;
	}
	return NULL;
}

void print_pipegap_data (struct pipegap *ptr)
{
    if (ptr->i_peak > 0.0) {
        fprintf (op, "pipegap at pole %d, from %d to %d ",
            ptr->parent->location, ptr->from, ptr->to);
        fprintf (op, "discharged %le Amperes.\n", ptr->i_peak);
    }
}

void pipegap_answers_cleanup (struct pipegap *ptr)
{
    if (ptr->conducting) {
        add_y (ptr->parent, ptr->from, ptr->to, -ptr->y);
    }
    if (fabs (ptr->i_peak) > fabs (predischarge)) {
        predischarge = ptr->i_peak;
    }
}

void inject_pipegap (struct pipegap *ptr)
{
    gsl_vector *v;
    double val;

    if (ptr->conducting) {
        v = ptr->parent->injection;
        val = ptr->i_past;
        gsl_vector_set (v, ptr->from, gsl_vector_get (v, ptr->from) - val);
        gsl_vector_set (v, ptr->to, gsl_vector_get (v, ptr->to) + val);
    }
}

void check_pipegap (struct pipegap *ptr)
{
    struct pole *p;
    int i, j, pos_now;
    double volts;

    p = ptr->parent;
    i = ptr->from;
    j = ptr->to;
    volts = gsl_vector_get (p->voltage, i) - gsl_vector_get (p->voltage, j);
    if (volts > 0.0) {
        pos_now = TRUE;
    } else {
        pos_now = FALSE;
    }

    if (ptr->conducting) {
        ptr->amps = volts * ptr->y + ptr->i_past;
        if (fabs (ptr->amps) > fabs (ptr->i_peak)) {
            ptr->i_peak = ptr->amps;
        }
        if (fabs (volts) < ptr->v_knee) {
            ptr->conducting = FALSE;
            add_y (p, i, j, -ptr->y);
            ptr->i_past = 0.0;
        }
    } else {
        if (fabs (volts) > ptr->v_knee) {
            ptr->conducting = TRUE;
            add_y (p, i, j, ptr->y);
            if (pos_now) {
                ptr->i_past = -ptr->i_bias;
            } else {
                ptr->i_past = ptr->i_bias;
            }
            solution_valid = FALSE;
        }
    }
}

int read_pipegap (void)
{
    int i, j, k;
    double f_knee, f_r;
    struct pipegap *ptr;
    double *target;
    int monitor;

    (void) next_double (&f_knee);
    (void) next_double (&f_r);
    if (f_knee < 0.0) {
        f_knee *= -1.0;
        monitor = 1;
    } else {
        monitor = 0;
    }
    if (f_r < 0.0) {
        f_r *= -1.0;
    }
    (void) read_pairs ();
    (void) read_poles ();
    (void) reset_assignments ();
    while (!next_assignment (&i, &j, &k)) {
        if (ptr = (struct pipegap *) malloc (sizeof *ptr)) {
            ptr->v_knee = f_knee;
            ptr->r_slope = f_r;
            ptr->i_bias = f_knee / f_r;
            ptr->y = 1.0 / f_r;
            ptr->parent = find_pole (i);
			if (!ptr->parent) oe_exit (ERR_BAD_POLE);
            ptr->parent->solve = TRUE;
            ptr->from = j;
            ptr->to = k;
            reset_pipegap (ptr);
            ptr->next = NULL;
            pipegap_ptr->next = ptr;
            pipegap_ptr = ptr;
            if (monitor) {
                target = &(ptr->amps);
                (void) add_ammeter (i, j, IPD_FLAG, target);
            }
        } else {
            if (logfp) fprintf (logfp, "can't allocate new pipegap\n");
            oe_exit (ERR_MALLOC);
        }
    }
    return (0);
}

void reset_pipegap (struct pipegap *ptr)
{
    ptr->i_peak = 0.0;
    ptr->i_past = 0.0;
    ptr->amps = 0.0;
    ptr->conducting = FALSE;
}

