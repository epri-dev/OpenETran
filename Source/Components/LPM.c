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
#include "Meter.h"
#include "../WritePlotFile.h"
#include "../ChangeTimeStep.h"
#include "Pole.h"
#include "Monitor.h"
#include "LPM.h"

#define SI_FOR_FO_STARTED  0.9999
#define SCALE_TOLERANCE  0.0001
#define MAX_SCALE     100.0
#define MIN_SCALE      0.01

static double lpm_si_counter = 0.0;  /* for progress feedback to SDW */

char lpm_token[] = "lpm";

struct lpm *lpm_head, *lpm_ptr;

int init_lpm_list (void)
{
    if ((lpm_head = (struct lpm *) malloc (sizeof *lpm_head))) {
        lpm_head->next = NULL;
        lpm_head->pts = NULL;
        lpm_ptr = lpm_head;
        return (0);
    }
    if (logfp) fprintf (logfp, "can't initialize lpm list\n");
    oe_exit (ERR_MALLOC);
    return (1);
}

void do_all_lpms (void (*verb) (struct lpm *))
{
    lpm_ptr = lpm_head;
    while ((lpm_ptr = lpm_ptr->next)) {
        verb (lpm_ptr);
    }
}

/*  &&&&  lpm functions  */

int read_lpm (void)
{
    int i, j, k;
    double f_cfo, f_e0, f_k;
    int flash_mode = LPM_NOT_FLASHED;
    struct lpm *ptr;

    (void) next_double (&f_cfo);
    (void) next_double (&f_e0);
    (void) next_double (&f_k);
    if (f_cfo < 0.0) {
        f_cfo *= -1.0;
        flash_mode = LPM_DISABLE_FLASH;
    }
    (void) read_pairs ();
    (void) read_poles ();
    (void) reset_assignments ();
    while (!next_assignment (&i, &j, &k)) {
        if ((ptr = (struct lpm *) malloc (sizeof *ptr))) {
            ptr->cfo = f_cfo;
            ptr->e0 = f_e0;
            ptr->k = f_k;
            ptr->pts = NULL;
            ptr->flash_mode = flash_mode;
            reset_lpm (ptr);
            move_lpm (ptr, i);
            ptr->from = j;
            ptr->to = k;
            ptr->next = NULL;
            lpm_ptr->next = ptr;
            lpm_ptr = ptr;
        } else {
            if (logfp) fprintf (logfp, "can't allocate new lpm\n");
            oe_exit (ERR_MALLOC);
        }
    }
    return (0);
}

void reset_lpm (struct lpm *ptr)
{
    int nsteps = (int) (Tmax / dT) + 2;
	int i;

    ptr->d = ptr->cfo / 560.0e3;
    ptr->xpos = ptr->d;
    ptr->xneg = ptr->d;
    ptr->t_flash = 0.0;
    ptr->vpk_pos = 0.0;
    ptr->vpk_neg = 0.0;
	ptr->SI = 0.0;
    if (ptr->flash_mode != LPM_DISABLE_FLASH) {
        ptr->flash_mode = LPM_NOT_FLASHED;
    }
    if (ptr->pts) {
        free (ptr->pts);
		ptr->pts = NULL;
    }
    ptr->pts = (float *) malloc (nsteps * sizeof (float));
	for (i = 0; i < nsteps; i++) {
		ptr->pts[i] = 0.0;
	}
	lpm_si_counter = 0.0;
}

void move_lpm (struct lpm *ptr, int i)
{
    ptr->parent = find_pole (i);
	if (!ptr->parent) oe_exit (ERR_BAD_POLE);
    ptr->parent->solve = TRUE;
}

static int lpm_flashes_over (struct lpm *ptr, double scale, int nsteps)
{
    int i, sign=0;
    double volts, ds, ds2;
    double x=1.0, dx;
    double xpos = ptr->d;
    double xneg = ptr->d;

    for (i = 0; i < nsteps ; i++) {
        volts = scale * ptr->pts[i];
        if (volts > 0.0) {
            sign = 1;
            x = xpos;
        } else if (volts < 0.0) {
            sign = -1;
            x = xneg;
        }
        volts = fabs (volts);
        ds = volts * ptr->k * dT;
        ds2 = ds * volts / x;
        ds *= ptr->e0;
        dx = ds2 - ds;
        if (sign > 0 && dx > 0.0) {
            xpos -= dx;
        } else if (dx > 0.0) {
            xneg -= dx;
        }
        if (xpos <= 0.0 || xneg <= 0.0) {
            return 1;
        }
    }
    return 0;
}

double calculate_lpm_si (struct lpm *ptr)
{
    int nsteps = (int) (Tmax / dT) + 1;
    double scale_low, scale_high, scale_mid;

    if (ptr->flash_mode == LPM_FLASHED) {
        return 1.0;
    }
	if (ptr->vpk_pos <= 0.0 && ptr->vpk_neg <= 0.0) {
		return 0.0;
	}
    // first bracket the root
    scale_low = scale_high = 1.0;
    while (scale_low > MIN_SCALE && lpm_flashes_over (ptr, scale_low, nsteps)) {
        scale_low *= 0.5;
    }
    while (scale_high < MAX_SCALE && !lpm_flashes_over (ptr, scale_high, nsteps)) {
        scale_high *= 2.0;
    }
    // now find the SI using bisection
    while (scale_high - scale_low > SCALE_TOLERANCE) {
        scale_mid = 0.5 * (scale_high + scale_low);
        if (lpm_flashes_over (ptr, scale_mid, nsteps)) {
            scale_high = scale_mid;
        } else {
            scale_low = scale_mid;
        }
    }
    scale_mid = 0.5 * (scale_high + scale_low);
    return 1.0 / scale_mid;
}

double estimate_lpm_si (struct lpm *ptr)
{
    double si_pos = 0.0;
    double si_neg = 0.0;

    if (ptr->flash_mode == LPM_FLASHED) {
        return 1.0;
    }
    if (ptr->xpos < ptr->d) {
        si_pos = SI_FOR_FO_STARTED;
    } else if (ptr->vpk_pos > 0.0) {
        si_pos = ptr->vpk_pos / ptr->cfo;
    }
    if (ptr->xneg < ptr->d) {
        si_neg = SI_FOR_FO_STARTED;
    } else if (ptr->vpk_neg > 0.0) {
        si_neg = ptr->vpk_neg / ptr->cfo;
    }
    return si_pos > si_neg ? si_pos : si_neg;
}

void print_lpm_data (struct lpm *ptr)
{
    fprintf (op, "insulator at pole %d, from %d to %d ",
        ptr->parent->location, ptr->from, ptr->to);
    if (ptr->flash_mode == LPM_FLASHED) {
        fprintf (op, "flashed over at %le seconds\n", ptr->t_flash);
    } else {
        fprintf (op, "per-unit SI = %le\n", calculate_lpm_si (ptr));
    }
}

void lpm_answers_cleanup (struct lpm *ptr)
{
    if (ptr->flash_mode == LPM_FLASHED) {
        ptr->SI = 1.0;
        add_y (ptr->parent, ptr->from, ptr->to, -Y_SHORT);
    } else if (want_si_calculation) {  /* didn't flashover, or was disabled */
        ptr->SI = calculate_lpm_si (ptr);
    } else {
        ptr->SI = estimate_lpm_si (ptr);
    }
    if (ptr->SI > SI) {
        SI = ptr->SI;
    }
}

void check_lpm (struct lpm *ptr)
{
    struct pole *p;
    int i, j, sign;
    double volts, ds, ds2;
    double x, dx;

	if (dT_switched) return;  /* disable flashover, and pts writing, after switching dT */
	/* CAUTION - pts would be overwritten unless it's resized after changing dT */
    if (ptr->flash_mode != LPM_FLASHED) {
        p = ptr->parent;
        i = ptr->from;
        j = ptr->to;
        volts = gsl_vector_get (p->voltage, i) - gsl_vector_get (p->voltage, j);
        ptr->pts[step] = (float) volts;
        if (volts > 0.0) {
            sign = 1;
            x = ptr->xpos;
        } else if (volts < 0.0) {
            sign = -1;
            x = ptr->xneg;
        } else {  // no voltage means no leader propagation
            return;
        }
        volts = fabs (volts);
        ds = volts * ptr->k * dT;
        ds2 = ds * volts / x;
        ds *= ptr->e0;
        dx = ds2 - ds;
        if (sign > 0) {
            if (dx > 0.0) {  // leader moves only when pushed the right direction
                ptr->xpos -= dx;
            }
            if (volts > ptr->vpk_pos) {
                ptr->vpk_pos = volts;
            }
        } else {
            if (dx > 0.0) {
                ptr->xneg -= dx;
            }
            if (volts > ptr->vpk_neg) {
                ptr->vpk_neg = volts;
            }
        }
        if (ptr->flash_mode == LPM_DISABLE_FLASH) { // need to keep waveshape, vpk_pos, and vpk_neg
            return;
        }
        if ((ptr->xpos <= 0.0) ||
            (ptr->xneg <= 0.0)) {
            ptr->flash_mode = LPM_FLASHED;
            if (flash_halt_enabled) {
                flash_halt = TRUE;
            }
            ptr->t_flash = t;
            add_y (p, i, j, Y_SHORT);
        }
    }
}

