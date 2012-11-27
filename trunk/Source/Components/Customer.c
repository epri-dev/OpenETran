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
#include "Ground.h"
#include "Meter.h"
#include "Customer.h"

char customer_token[] = "customer";

struct customer *customer_head, *customer_ptr;

void print_customer_data (struct customer *ptr)
{
	fprintf (op, "Customer at pole %d, from %d to %d.\n",
		ptr->parent->location, ptr->from, ptr->to);
	fprintf (op, "\tMax Vp  = %le volts\n", ptr->Vp);
	fprintf (op, "\tMax Ihg = %le amps\n", ptr->Ihg);
	fprintf (op, "\tMax Ix2 = %le amps\n", ptr->Ix2_peak);
}

/* perform the primary voltage integration, and keep track of the peak
voltage and current values */

void update_customer_history (struct customer *ptr)
{
	double v, i, i_new;
	
	v = gsl_vector_get (ptr->parent->voltage, ptr->from) - gsl_vector_get (ptr->parent->voltage, ptr->to);
	i = ptr->in->amps;
	ptr->integral += v * ptr->Kv;
	i_new = ptr->Ki * i + ptr->integral; /* Ix2 at this time step */
	if (fabs (i) > fabs (ptr->Ihg)) {
		ptr->Ihg = i;
	}
	if (fabs (v) > fabs (ptr->Vp)) {
		ptr->Vp = v;
	}
	if (fabs (i_new) > fabs (ptr->Ix2_peak)) {
		ptr->Ix2_peak = i_new;
	}
	ptr->Ix2 = i_new;
}

int init_customer_list (void)
{
	if ((customer_head = (struct customer *) malloc (sizeof *customer_head)) != NULL) {
		customer_head->next = NULL;
		customer_ptr = customer_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize customer list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}

void do_all_customers (void (*verb) (struct customer *))
{
	customer_ptr = customer_head;
	while ((customer_ptr = customer_ptr->next) != NULL) {
		verb (customer_ptr);
	}
}

struct customer *find_customer (int at, int from, int to)
{
	customer_ptr = customer_head;
	while ((customer_ptr = customer_ptr->next) != NULL) {
		if ((customer_ptr->parent->location == at) &&
			(customer_ptr->from == from) &&(customer_ptr->from == from)) return customer_ptr;
	}
	return NULL;
}

/* read the extensive input for a pole-top transformer, service drop,
and house ground */

int read_customer (void)
{
	int i, j, k;
	double Rhg, rho, e0, Lhg, Dhg, N, Lp, Ls1, Ls2, ra, rn, Dan, Daa, l;
	double La, Ln, Laa, Lan, Lcm, Lfw, Ki, Kv, Denom;
	struct customer *ptr;
	double *target;
	
	(void) next_double (&Rhg);
	(void) next_double (&rho);
	(void) next_double (&e0);
	(void) next_double (&Lhg);
	(void) next_double (&Dhg);
	Lhg *= Dhg;
	(void) next_double (&N);
	(void) next_double (&Lp);
	(void) next_double (&Ls1);
	(void) next_double (&Ls2);
	(void) next_double (&Lcm);
	(void) next_double (&ra);
	(void) next_double (&rn);
	(void) next_double (&Dan);
	(void) next_double (&Daa);
	(void) next_double (&l);
	Lcm *= l;
/* process the service drop parameters into secondary inductances */
	La = PRIM_L * l * (log (2.0 * l / ra) - 1.0);
	Ln = PRIM_L * l * (log (2.0 * l / rn) - 1.0);
	Laa = PRIM_L * l * (log (2.0 * l / Daa) - 1.0);
	Lan = PRIM_L * l * (log (2.0 * l / Dan) - 1.0);
	Lfw = 4.0 * Lp / N / N + Ls1 + Ls2;
	Denom = 0.5 * (Ls1 + Ls2) + La + 2.0 * Ln + Laa - 4.0 * Lan
		- 0.5 * (Ls1 - Ls2) * (Ls1 - Ls2) / 
			(Lfw + 2.0 * La - 2.0 * Laa);
/* coefficients for calculating Ix2 - see Dave Smith's IEEE papers */
	Ki = (Ln - Lan) / Denom;
	Kv = (Ls2 - Ls1) / N / (Lfw + 2.0 * La - 2.0 * Laa) / Denom;
	(void) read_pairs ();
	(void) read_poles ();
	(void) reset_assignments ();
	while (!next_assignment (&i, &j, &k)) {
		if ((ptr = (struct customer *) malloc (sizeof *ptr)) != NULL) {
			ptr->parent = find_pole (i);
			if (!ptr->parent) oe_exit (ERR_BAD_POLE);
			ptr->parent->solve = TRUE;
			ptr->in = add_ground (i, k, 0, 
				Rhg, rho, e0, Lcm + Lhg);
			ptr->Ki = 2.0 * Ki;
			ptr->Kv = 2.0 * Kv * dT;
			reset_customer (ptr);
			ptr->from = j;
			ptr->to = k;
			ptr->next = NULL;
			customer_ptr->next = ptr;
			customer_ptr = ptr;
/* add meters for the house ground and X2 currents automatically,
mainly useful for the DOS version */
			target = &(ptr->in->amps);
			(void) add_ammeter (i, j, IHG_FLAG, target);
			target = &(ptr->Ix2);
			(void) add_ammeter (i, j, IX2_FLAG, target);
		} else {
			if (logfp) fprintf( logfp, "can't allocate new customer\n");
			oe_exit (ERR_MALLOC);
		}
	}
	return (0);
}

/* reset max values and integral for coordination current iterations */

void reset_customer (struct customer *ptr)
{
	ptr->Ix2 = 0.0;
	ptr->Ihg = 0.0;
	ptr->Vp = 0.0;
	ptr->integral = 0.0;
	ptr->Ix2_peak = 0.0;
}
