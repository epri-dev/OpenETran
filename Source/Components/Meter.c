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
#include "Arrester.h"
#include "Arrbez.h"
#include "Customer.h"
#include "Ground.h"
#include "PipeGap.h"
#include "Meter.h"

char meter_token[] = "meter";

struct meter *meter_head, *meter_ptr;

/* we have 5 different kinds of meters, depending on the "to" node:
	non-negative integer	=> voltmeter
    -1						=> arrester (or arrbez) current
	-2						=> pole ground current
	-3						=> house ground current
	-4						=> transformer X2 current
	-5						=> pipegap current
*/

void print_meter_data (struct meter *ptr)
{
	int flag;
	
	flag = ptr->to;
	fprintf (op, "Meter at pole %d, ", ptr->at);
	if (flag == -1) {
		fprintf (op, "arrester max I      ");
	} else if (flag == -2) {
		fprintf (op, "pole ground max I   ");
	} else if (flag == -3) {
		fprintf (op, "house ground max I  ");
	} else if (flag == -4) {
		fprintf (op, "transformer X2 max I");
	} else {
		fprintf (op, "from %d to %d, max V  ", ptr->from, ptr->to);
	}
	fprintf (op, " = %le\n", ptr->vmax);
}

void update_meter_peaks (struct meter *ptr)
{
	double volts;
	
	volts = *(ptr->v_from) - *(ptr->v_to);
	if (fabs (volts) > fabs (ptr->vmax)) {
		ptr->vmax = volts;
	}
}

int init_meter_list (void)
{
	if ((meter_head = (struct meter *) malloc (sizeof *meter_head)) != NULL) {
		meter_head->next = NULL;
		meter_ptr = meter_head;
		return (0);
	}
	if (logfp) fprintf( logfp, "can't initialize meter list\n");
	oe_exit (ERR_MALLOC);
	return (1);
}

void do_all_meters (void (*verb) (struct meter *))
{
	meter_ptr = meter_head;
	while ((meter_ptr = meter_ptr->next) != NULL) {
		verb (meter_ptr);
	}
}

/* read meters from the file or buffer, and add them as voltmeters.
ammeters are added with ground, arrester, or customer input */

int read_meter (void)
{
	int i, j, k, mtype;
	struct arrester *aptr;
	struct arrbez *bptr;
	struct ground *gptr;
	struct customer *cptr;
	struct pipegap *pptr;
	
    (void) next_int (&mtype);
	mtype *= -1;
	(void) read_pairs ();
	(void) read_poles ();
	(void) reset_assignments ();
	while (!next_assignment (&i, &j, &k)) {
		switch (mtype) {
		case VOLT_FLAG:
			(void) add_voltmeter (i, j, k);
			break;
		case IARR_FLAG:
			if ((aptr = find_arrester (i,j,k)) != NULL) {
				(void) add_ammeter (i, j, IARR_FLAG, &(aptr->amps));
			} else if ((bptr = find_arrbez (i,j,k)) != NULL) {
				(void) add_ammeter (i, j, IARR_FLAG, &(bptr->amps));
			}
			break;
		case IPG_FLAG:
			if ((gptr = find_ground (i,j,k)) != NULL) {
				(void) add_ammeter (i, j, IPG_FLAG, &(gptr->amps));
			}
			break;
		case IHG_FLAG:
			if ((cptr = find_customer (i,j,k)) != NULL) {
				(void) add_ammeter (i, j, IHG_FLAG, &(cptr->in->amps));
			}
			break;
		case IX2_FLAG:
			if ((cptr = find_customer (i,j,k)) != NULL) {
				(void) add_ammeter (i, j, IX2_FLAG, &(cptr->Ix2));
			}
			break;
		case IPD_FLAG:
			if ((pptr = find_pipegap (i,j,k)) != NULL) {
				(void) add_ammeter (i, j, IPD_FLAG, &(pptr->amps));
			}
			break;
		default: break;
		}
	}
	return (0);
}

/* add a meter to measure voltage from node j to k, at pole i */

struct meter *add_voltmeter (int i, int j, int k)
{
	struct meter *ptr;
	struct pole *pptr;
	
	if ((ptr = (struct meter *) malloc (sizeof *ptr)) != NULL) {
		pptr = find_pole (i);
		if (!pptr) oe_exit (ERR_BAD_POLE);
		ptr->v_from = gsl_vector_ptr (pptr->voltage, j);
		ptr->v_to = gsl_vector_ptr (pptr->voltage, k);
		pptr->solve = TRUE;
		ptr->at = i;
		ptr->from = j;
		ptr->to = k;
		reset_meter (ptr);
		ptr->next = NULL;
		meter_ptr->next = ptr;
		meter_ptr = ptr;
		return (ptr);
	} else {
		if (logfp) fprintf( logfp, "can't allocate new voltmeter\n");
		oe_exit (ERR_MALLOC);
	}
	return (NULL);
}

/* add a meter to measure current in a branch.  Since there can be more
than one branch between 2 nodes, we need to pass in a pointer to the
actual monitored value */

static double ground_voltage = 0.0;

struct meter *add_ammeter (int i, int j, int k, double *target)
{
	struct meter *ptr;
	
	if ((ptr = (struct meter *) malloc (sizeof *ptr)) != NULL) {
		ptr->v_from = target;
		ptr->v_to = &ground_voltage;
		ptr->at = i;
		ptr->from = j;
		ptr->to = k;
		reset_meter (ptr);
		ptr->next = NULL;
		meter_ptr->next = ptr;
		meter_ptr = ptr;
		return (ptr);
	} else {
		if (logfp) fprintf( logfp, "can't allocate new ammeter\n");
		oe_exit (ERR_MALLOC);
	}
	return (NULL);
}

/* reset maximum voltage to zero */

void reset_meter (struct meter *ptr)
{
	ptr->vmax = 0.0;
}

