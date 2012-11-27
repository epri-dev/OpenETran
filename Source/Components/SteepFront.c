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
#include "Pole.h"
#include "SteepFront.h"

#define MAX_SF_PTS   25
#define DX_LOW    0.300
#define DX_HIGH   0.005
#define DKNOT     1.005

char steepfront_token[] = "steepfront";

struct steepfront *steepfront_head, *steepfront_ptr;

int init_steepfront_list (void)
{
    if (steepfront_head = (struct steepfront *) malloc (sizeof *steepfront_head)) {
        steepfront_head->next = NULL;
        steepfront_head->shape = NULL;
        steepfront_ptr = steepfront_head;
        return (0);
    }
    if (logfp) fprintf (logfp, "can't initialize steepfront list\n");
    oe_exit (ERR_MALLOC);
    return (1);
}

void do_all_steepfronts (void (*verb) (struct steepfront *))
{
    steepfront_ptr = steepfront_head;
    while (steepfront_ptr = steepfront_ptr->next) {
        verb (steepfront_ptr);
    }
}

int read_steepfront (void)
{
    int i, j, k;
    double fpeak, ftf, ftt, ftstart, fsi;
    struct steepfront *ptr;

    (void) next_double (&fpeak);
    (void) next_double (&ftf);
    (void) next_double (&ftt);
    (void) next_double (&ftstart);
    (void) next_double (&fsi);
    (void) read_pairs ();
    (void) read_poles ();
    (void) reset_assignments ();
    while (!next_assignment (&i, &j, &k)) {
        if (ptr = (struct steepfront *) malloc (sizeof *ptr)) {
            ptr->shape = NULL;
            move_steepfront (ptr, i, j, k, fpeak, ftf, ftt, ftstart, fsi);
            ptr->next = NULL;
            steepfront_ptr->next = ptr;
            steepfront_ptr = ptr;
        } else {
            if (logfp) fprintf (logfp, "can't allocate new steepfront\n");
            oe_exit (ERR_MALLOC);
        }
    }
    return (0);
}

void move_steepfront (struct steepfront *ptr, int i, int j, int k, double fpeak,
    double ftf, double ftt, double ftstart, double pu_si)
{
    double xpts[MAX_SF_PTS], ypts[MAX_SF_PTS];
    double x, dx, t50, t10, t30, t90, tau, si, xstart;
    int npts = 0;

    if (ptr->shape) {
        free_bezier_fit (ptr->shape);
        free (ptr->shape);
    }
    ptr->front = ftf;
    ptr->tail = ftt;
    ptr->tstart = ftstart;
    ptr->peak = fpeak;
    ptr->pu_si = pu_si;
    si = pu_si * fpeak / ftf;
    ptr->si = si;

    t10 = 0.78 * ftf;
    t30 = 1.16 * ftf;
    t90 = 1.76 * ftf;
    xpts[npts] = 0.0;     ypts[npts++] = 0.00;
    xpts[npts] = t10;     ypts[npts++] = 0.10 * fpeak;
    xpts[npts] = t30;     ypts[npts++] = 0.30 * fpeak;
    xpts[npts] = t30 * DKNOT; ypts[npts++] = 0.30 * fpeak * DKNOT;
    dx = DX_LOW * fpeak / si;
    xpts[npts] = t90 - dx;    ypts[npts++] = (0.90 - DX_LOW) * fpeak;
    xpts[npts] = t90;     ypts[npts++] = 0.90 * fpeak;
    dx = DX_HIGH * fpeak / si;
    xpts[npts] = t90 + dx;    ypts[npts++] = (0.90 + DX_HIGH) * fpeak;
    x = t90 + dx * 0.1 / DX_HIGH;
    xpts[npts] = x;       ypts[npts++] = 1.00 * fpeak;   x *= 1.2;
    xpts[npts] = x;       ypts[npts++] = 1.00 * fpeak;
    xstart = x;
    t50 = ftt - xstart;
    tau = ETKONST * t50;
    dx = 0.5 * tau;
    x += dx;
    xpts[npts] = x;       ypts[npts++] = fpeak * exp (-(x - xstart) / tau);    x += dx;
    xpts[npts] = x;       ypts[npts++] = fpeak * exp (-(x - xstart) / tau);    x += dx;
    xpts[npts] = x;      ypts[npts++] = fpeak * exp (-(x - xstart) / tau);   x += dx;
    xpts[npts] = x;      ypts[npts++] = fpeak * exp (-(x - xstart) / tau);   x += dx;
    xpts[npts] = x;      ypts[npts++] = fpeak * exp (-(x - xstart) / tau);   x += dx;
    xpts[npts] = x;      ypts[npts++] = fpeak * exp (-(x - xstart) / tau);
    x *= 10.0;
    xpts[npts] = x;      ypts[npts++] = fpeak * exp (-(x - xstart) / tau);

    ptr->shape = build_bezier (xpts, ypts, npts, FALSE);

    ptr->parent = find_pole (i);
	if (!ptr->parent) oe_exit (ERR_BAD_POLE);
    ptr->parent->solve = TRUE;
    ptr->from = j;
    ptr->to = k;
}

void inject_steepfront (struct steepfront *ptr)
{
    double x, i;
	gsl_vector *v = ptr->parent->injection;

    x = t - ptr->tstart;
    if (x > 0.0) {
        i = bez_eval (ptr->shape, x);
        gsl_vector_set (v, ptr->from, gsl_vector_get (v, ptr->from) + i);
        gsl_vector_set (v, ptr->to, gsl_vector_get (v, ptr->to) - i);
    }
}
