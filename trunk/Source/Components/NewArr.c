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
#include "Meter.h"
#include "../WritePlotFile.h"
#include "../Parser.h"
#include "../ReadUtils.h"
#include "Pole.h"

#include "ArrBez.h"

#include "NewArr.h"

char newarr_token[] = "newarr";

#define SHORT_CIRCUIT_G   1.0e6
#define OPEN_CIRCUIT_G    1.0e-7
#define newarr_TREF      80.0
#define newarr_IREF       5.4e3

#define MAX_ARR_PTS  20  /* >= number of rows in double[][3] arrays below */

/*
struct bezier_fit *build_arrester (double v10, enum arr_size_type arr_size,
    enum arr_char_type arr_char, enum arr_minmax_type arr_minmax, int use_linear)
*/
/*  &&&&  newarr functions  */

struct newarr *newarr_head, *newarr_ptr;

int init_newarr_list (void)
{
    if ((newarr_head = (struct newarr *) malloc (sizeof *newarr_head))) {
        newarr_head->next = NULL;
        newarr_head->shape = NULL;
        newarr_ptr = newarr_head;
        return (0);
    }
    if (logfp) fprintf (logfp, "can't initialize newarr list\n");
    oe_exit (ERR_MALLOC);
    return (1);
}

void do_all_newarrs (void (*verb) (struct newarr *))
{
    newarr_ptr = newarr_head;
    while ((newarr_ptr = newarr_ptr->next)) {
        verb (newarr_ptr);
    }
}

struct newarr *find_newarr (int at, int from, int to)
{
	newarr_ptr = newarr_head;
	while (((newarr_ptr = newarr_ptr->next) != NULL)) {
		if ((newarr_ptr->parent->location == at) &&
			(newarr_ptr->from == from) &&(newarr_ptr->from == from)) return newarr_ptr;
	}
	return NULL;
}

int read_newarr (void)
{
    int i, j, k;
    double f_v10, f_vgap, f_L, f_length, f_Uref;
    struct newarr *ptr;
    double *target;
    int monitor;
    int use_linear = FALSE;

    (void) next_double (&f_vgap);
    (void) next_double (&f_v10);
    (void) next_double (&f_Uref);
    (void) next_double (&f_L);
    (void) next_double (&f_length);
    (void) next_int (&monitor);
    if (f_v10 < 0.0) {
        use_linear = TRUE;
        f_v10 = -f_v10;
    }
    f_L *= f_length;
    (void) read_pairs ();
    (void) read_poles ();
    (void) reset_assignments ();
    while (!next_assignment (&i, &j, &k)) {
        if ((ptr = (struct newarr *) malloc (sizeof *ptr))) {
            ptr->vgap = f_vgap;
            ptr->v10 = f_v10;
            ptr->Uref = f_Uref * f_v10;
            ptr->rl = 2.0 * f_L / dT;
            if (ptr->rl > 0.0) {
                ptr->gl = dT / f_L;
            } else {
                ptr->gl = 0.0;
            }
            ptr->parent = find_pole (i);
			if (!ptr->parent) oe_exit (ERR_BAD_POLE);
            ptr->parent->solve = TRUE;
            ptr->parent->num_nonlinear += 1;
            ptr->pole_num_nonlinear = ptr->parent->num_nonlinear;
            ptr->from = j;
            ptr->to = k;
            ptr->shape = build_arrester (f_v10,
                f_v10 > 140.0e3 ? arrsize_54_to_360 : arrsize_2pt7_to_48,
                arr_char_8x20, arr_use_vmax, use_linear);
            reset_newarr (ptr);
            ptr->next = NULL;
            newarr_ptr->next = ptr;
            newarr_ptr = ptr;
            if (monitor) {
                target = &(ptr->amps);
                (void) add_ammeter (i, j, IARR_FLAG, target);
            }
        } else {
            if (logfp) fprintf (logfp, "can't allocate new newarr\n");
            oe_exit (ERR_MALLOC);
        }
    }
    return (0);
}

void reset_newarr (struct newarr *ptr)
{
    ptr->t_start = 0.0;
    ptr->t_peak = 0.0;
    ptr->energy = 0.0;
    ptr->charge = 0.0;
    ptr->i_peak = 0.0;
    ptr->amps = 0.0;
    ptr->varr = 0.0;
    ptr->h = 0.0;
    if (ptr->vgap > 0.0) {
        ptr->rgap = ptr->vgap / 1.0e-3;
    } else {
        ptr->rgap = 0.0;
		ptr->t_start = dT;
    }
    if (ptr->Uref > 0.0) {
        ptr->Gref = 34.0 / (ptr->v10 / 1000.0);
        ptr->g = OPEN_CIRCUIT_G;
    } else {
        ptr->Gref = 0.0;
        ptr->g = SHORT_CIRCUIT_G;
    }
    ptr->r = ptr->rl + ptr->rgap + 1.0 / ptr->g;
}

void print_newarr_data (struct newarr *ptr)
{
    if (ptr->t_start > 0.0) {
        fprintf (op, "Arrester at pole %d, from %d to %d ",
            ptr->parent->location, ptr->from, ptr->to);
        fprintf (op, "discharged %le Amperes.\n", ptr->i_peak);
        fprintf (op, "\tTime operated: %le", ptr->t_start);
        fprintf (op, "\tTime of peak: %le\n", ptr->t_peak);
        fprintf (op, "\tCharge: %le\n", ptr->charge);
        fprintf (op, "\tEnergy: %le\n", ptr->energy);
    }
}

void newarr_answers_cleanup (struct newarr *ptr)
{
    if (ptr->energy > energy) {
        energy = ptr->energy;
    }
    if (fabs (ptr->i_peak) > fabs (current)) {
        current = ptr->i_peak;
    }
    if (fabs (ptr->charge) > fabs (charge)) {
        charge = ptr->charge;
    }
}

void update_newarr_history (struct newarr *ptr)
{
    double dCharge, dG, Ipu, Vpu, Gpu, Vgap;
    struct pole *p;
    int i, j;

    p = ptr->parent;
    i = ptr->from;
    j = ptr->to;
    Vgap = gsl_vector_get (p->voltage, i) - gsl_vector_get (p->voltage, j);

    if (ptr->rgap > 0.0) {  // gap did not sparkover yet
        if (fabs (Vgap) > fabs (ptr->vgap)) {  // start conducting next time step
            ptr->rgap = 0.0;
            ptr->t_start = t;
            ptr->r = ptr->rl + ptr->rgap + 1.0 / ptr->g;
        }
        return;
    }
    if (ptr->Uref > 0.0 && ptr->g < SHORT_CIRCUIT_G) {  // has the Cigre dynamics
        Ipu = ptr->amps / newarr_IREF;
        Vpu = fabs (Vgap) / ptr->Uref;
        Gpu = ptr->g / ptr->Gref;
        dG = (ptr->Gref / newarr_TREF) * (1.0 + Gpu) * (1.0 + Gpu * Ipu * Ipu) * exp (Vpu);
        ptr->g += dG * dT;
        ptr->r = ptr->rl + ptr->rgap + 1.0 / ptr->g;
//          printf ("t, Ur, Ipu, Vpu, Gpu, dG/dT, g, r = %le %le %le %le %le %le %le %le\n",
//              t, ptr->Uref, Ipu, Vpu, Gpu, dG, ptr->g, ptr->r);
    }
    dCharge = dT * ptr->amps;
    ptr->charge += dCharge;
    ptr->energy += dCharge * ptr->varr;
    if (fabs (ptr->amps) > fabs (ptr->i_peak)) {
        ptr->i_peak = ptr->amps;
        ptr->t_peak = t;
    }
}

