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

/* This file contains functions to increase the time step during
a simulation. */

#include <wtypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "ChangeTimeStep.h"
#include "OETypes.h"
#include "OEEngine.h"
#include "Components/Arrbez.h"
#include "Components/Arrester.h"
#include "Components/Capacitor.h"
#include "Components/Customer.h"
#include "Components/Ground.h"
#include "Components/Inductor.h"
#include "Components/Line.h"
#include "Components/Pole.h"

char time_token[] = "time";
char change_dt_token[] = "2dt";

void change_customer_time_step (struct customer *ptr);
void change_capacitor_time_step (struct capacitor *ptr);
void change_inductor_time_step (struct inductor *ptr);
void change_arrester_time_step (struct arrester *ptr);
void change_arrbez_time_step (struct arrbez *ptr);
void change_line_time_step (struct line *ptr);
void change_ground_time_step (struct ground *ptr);

void restore_customer_time_step (struct customer *ptr);
void restore_capacitor_time_step (struct capacitor *ptr);
void restore_inductor_time_step (struct inductor *ptr);
void restore_arrester_time_step (struct arrester *ptr);
void restore_arrbez_time_step (struct arrbez *ptr);
void restore_line_time_step (struct line *ptr);
void restore_ground_time_step (struct ground *ptr);

double first_dT;
double second_dT;
double dT_switch_time;
int dT_switched;
int using_second_dT;

/* these three functions are called at the system level */

void restore_time_step (void)
{
	if (!dT_switched) return;

	dT = first_dT;

	do_all_arrbezs (restore_arrbez_time_step);
	do_all_arresters (restore_arrester_time_step);
	do_all_capacitors (restore_capacitor_time_step);
	do_all_customers (restore_customer_time_step);
	do_all_grounds (restore_ground_time_step);
	do_all_inductors (restore_inductor_time_step);
	do_all_lines (restore_line_time_step);

	dT_switched = FALSE;
}

void change_time_step (void)
{
	dT = second_dT;

	do_all_arrbezs (change_arrbez_time_step);
	do_all_arresters (change_arrester_time_step);
	do_all_capacitors (change_capacitor_time_step);
	do_all_customers (change_customer_time_step);
	do_all_grounds (change_ground_time_step);
	do_all_inductors (change_inductor_time_step);
	do_all_lines (change_line_time_step);

	do_all_poles (triang_pole);

	dT_switched = TRUE;
}

/* arrbez functions to support second_dT */

void change_arrbez_time_step (struct arrbez *ptr)
{
	double vl = ptr->rl * (ptr->amps - ptr->h);

	ptr->rl *= (first_dT / second_dT);
	ptr->gl *= (second_dT / first_dT);
    ptr->r = ptr->rl + ptr->rgap + 1.0 / ptr->g;

	ptr->h = ptr->amps - 0.5 * ptr->gl * vl;
}

void restore_arrbez_time_step (struct arrbez *ptr)
{
	ptr->rl *= (second_dT / first_dT);
	ptr->gl *= (first_dT / second_dT);
    ptr->r = ptr->rl + ptr->rgap + 1.0 / ptr->g;
}

/* arrester functions to support second_dT */

void change_arrester_time_step (struct arrester *ptr)
{
	double old_y = ptr->y;
	double volts, vl, vr;
	int pos_now;

	ptr->zl *= (first_dT / second_dT);
	ptr->y = 1.0 / (ptr->r_slope + ptr->zl);
	ptr->yr = ptr->y * ptr->r_slope;
	ptr->yzl = ptr->y * ptr->zl;
	if (ptr->conducting) {
		add_y (ptr->parent, ptr->from, ptr->to, ptr->y - old_y);
		volts = gsl_vector_get (ptr->parent->voltage, ptr->from) - gsl_vector_get (ptr->parent->voltage, ptr->to);
		if (volts > 0.0) {
			pos_now = TRUE;
			vr = ptr->r_slope * (ptr->amps + ptr->i_bias);
		} else {
			pos_now = FALSE;
			vr = ptr->r_slope * (ptr->amps - ptr->i_bias);
		}
		vl = volts - vr;
		if (ptr->zl > 0.0) {
			ptr->h = ptr->amps + vl / ptr->zl;
		}
		ptr->i = ptr->h * ptr->yzl;
		if (pos_now) {
			ptr->i -= ptr->yr * ptr->i_bias;
		} else {
			ptr->i += ptr->yr * ptr->i_bias;
		}
		ptr->i_past = ptr->i;
	} 
}

void restore_arrester_time_step (struct arrester *ptr)
{
	double old_y = ptr->y;

	ptr->zl *= (second_dT / first_dT);
	ptr->y = 1.0 / (ptr->r_slope + ptr->zl);
	ptr->yr = ptr->y * ptr->r_slope;
	ptr->yzl = ptr->y * ptr->zl;
	if (ptr->conducting) {
//		add_y (ptr->parent, ptr->from, ptr->to, ptr->y - old_y);
	}
}

/* capacitor functions to support second_dT */

void change_capacitor_time_step (struct capacitor *ptr)
{
	double old_y = ptr->y;

	ptr->y *= (first_dT / second_dT);
	ptr->yc = ptr->y + ptr->y;
	add_y (ptr->parent, ptr->from, ptr->to, ptr->y - old_y);

	ptr->h *= (first_dT / second_dT);
}

void restore_capacitor_time_step (struct capacitor *ptr)
{
	double old_y = ptr->y;

	ptr->y *= (second_dT / first_dT);
	ptr->yc = ptr->y + ptr->y;
	add_y (ptr->parent, ptr->from, ptr->to, ptr->y - old_y);
}

/* customer functions to support second_dT */

void change_customer_time_step (struct customer *ptr)
{
	ptr->Kv *= (second_dT / first_dT);
}

void restore_customer_time_step (struct customer *ptr)
{
	ptr->Kv *= (first_dT / second_dT);
}

/* ground functions to support second_dT */

void change_ground_time_step (struct ground *ptr)
{
	double old_y = ptr->y;
	double Vl, Vg, Vt;

	ptr->zl *= (first_dT / second_dT);
	ptr->y = 1.0 / (ptr->R60 + ptr->zl);
	ptr->yr = ptr->y * ptr->R60;
	ptr->yzl = ptr->y * ptr->zl;
	add_y (ptr->parent, ptr->from, ptr->to, ptr->y - old_y);

	Vt = gsl_vector_get (ptr->parent->voltage, ptr->from) - gsl_vector_get (ptr->parent->voltage, ptr->to);
	Vg = ptr->amps * ptr->Ri;
	Vl = Vt - Vg;
	if (ptr->zl > 0.0) {
		ptr->h = ptr->amps + Vl / ptr->zl;
	} else {
		ptr->h = 0.0;
	}
	ptr->i = ptr->h * ptr->yzl + ptr->i_bias * ptr->yr;
}

void restore_ground_time_step (struct ground *ptr)
{
	double old_y = ptr->y;

	ptr->zl *= (second_dT / first_dT);
	ptr->y = 1.0 / (ptr->R60 + ptr->zl);
	ptr->yr = ptr->y * ptr->R60;
	ptr->yzl = ptr->y * ptr->zl;
	add_y (ptr->parent, ptr->from, ptr->to, ptr->y - old_y);
}

/* inductor functions to support second_dT */

void change_inductor_time_step (struct inductor *ptr)
{
	double old_y = ptr->y;
	double Vt, It;

	Vt = gsl_vector_get (ptr->parent->voltage, ptr->from) - gsl_vector_get (ptr->parent->voltage, ptr->to);
	It = old_y * Vt + ptr->h;

	ptr->y = 1.0 / (ptr->res + 2.0 * ptr->ind / dT);
	ptr->yi = 2.0 * ptr->y * (1.0 - ptr->res * ptr->y);
	ptr->zi = 1.0 - 2.0 * ptr->res * ptr->y;
	add_y (ptr->parent, ptr->from, ptr->to, ptr->y - old_y);

	ptr->h = ptr->y * ((2.0 * ptr->ind / dT - ptr->res) * It + Vt);
}

void restore_inductor_time_step (struct inductor *ptr)
{
	double old_y = ptr->y;

	ptr->y = 1.0 / (ptr->res + 2.0 * ptr->ind / dT);
	ptr->yi = 2.0 * ptr->y * (1.0 - ptr->res * ptr->y);
	ptr->zi = 1.0 - 2.0 * ptr->res * ptr->y;
	add_y (ptr->parent, ptr->from, ptr->to, ptr->y - old_y);
}

/* line functions to support second_dT */

void change_line_time_step (struct line *ptr)
{
	int i, k;

	k = step % ptr->steps;
	for (i = 0; i < number_of_conductors; i++) {
		gsl_matrix_set (ptr->hist_left, i, 0, gsl_matrix_get (ptr->hist_left, i, k));
		gsl_matrix_set (ptr->hist_right, i, 0, gsl_matrix_get (ptr->hist_right, i, k));
	}
	ptr->steps = 1;
}

void restore_line_time_step (struct line *ptr)
{
	ptr->steps = ptr->alloc_steps;
}
