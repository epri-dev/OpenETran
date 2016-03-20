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

#include <stdlib.h>
#include <math.h>

#include "BezUtils.h"

int get_bez_size (int size)
{
    return (size-1)*3+1;
}

struct bezier_fit *allocate_bezier (int npts)
{
    struct bezier_fit *p;

    p = (struct bezier_fit *) malloc (sizeof *p);

    p->npts = npts;
    p->bez_size = get_bez_size (npts);
    p->xarr = (double *) malloc (npts * sizeof (double));
    p->ybez_arr = (double *) malloc (p->bez_size * sizeof (double));
    p->next = NULL;

    return p;
}

void fill_bezier (struct bezier_fit *p, double *xpts, double *ypts, int use_linear)
{
    int ii, i, jj, kk;
    int npts = p->npts;

    for (ii = 0; ii < npts; ii++) {
        p->xarr[ii] = xpts[ii];
        p->ybez_arr[3*ii] = ypts[ii];
    }
    p->start_slope = (ypts[1] - ypts[0]) / (xpts[1] - xpts[0]);
    p->end_slope = (ypts[npts-1] - ypts[npts-2])
        / (xpts[npts-1] - xpts[npts-2]);

/* Place control points at equal y spacings between "main" y points */
    for( i = 0; i < npts -1; i++){
        ii = 3*i;
        kk = ii +3;
        jj = ii +1;
        p->ybez_arr[jj]   = (p->ybez_arr[ii]*2.0+p->ybez_arr[kk])/3.0;
        p->ybez_arr[jj+1] = (p->ybez_arr[ii]+p->ybez_arr[kk]*2.0)/3.0;
    }

    if (use_linear) {  /* piecewise linear characteristic */
        return;
    }

  /* readjust each of the "main" y points so that it lies on a line */
  /* between the two nearest control points */
    for (i = 1 ; i < npts-1; i++){
        jj = 3*i;
        p->ybez_arr[jj] =
            p->ybez_arr[jj-1] +
            (p->ybez_arr[jj+1] - p->ybez_arr[jj-1]) /
            (p->xarr[i+1] - p->xarr[i-1]) *
            (p->xarr[i] - p->xarr[i-1]);
    }
}

struct bezier_fit *build_bezier (double *xpts, double *ypts, int npts, int use_linear)
{
    struct bezier_fit *p = allocate_bezier (npts);

    fill_bezier (p, xpts, ypts, use_linear);
    return p;
}

double bez_eval (struct bezier_fit *p, double xx)
{
    double x1, x4, y1, y2, y3, y4, z, c1, c2, c3;
    int    i, j1;
    double y = 0.0;
    double sign = 1.0;

    if (xx < p->xarr[0]) {
        xx = 2.0 * p->xarr[0] - xx;
        sign = -1.0;
    }
    if (xx <= p->xarr[0]) {
        y = p->ybez_arr[0]
            + p->start_slope * (xx - p->xarr[0]);
        return sign * y;
    } else if (xx >= p->xarr[p->npts-1]) {
        y = p->ybez_arr[p->bez_size-1]
            + p->end_slope * (xx - p->xarr[p->npts-1]);
        return sign * y;
    }
    for (i = 0; i < p->npts - 1; i++) {
        x1 = p->xarr[i];
        x4 = p->xarr[i+1];
        if (xx <= x4) {
            j1 = 3*i;
            y1 = p->ybez_arr[j1];
            y2 = p->ybez_arr[++j1];
            y3 = p->ybez_arr[++j1];
            y4 = p->ybez_arr[++j1];
            z = (xx - x1) / (x4 - x1);
            c1 = 3.0*(y2 - y1);
            c2 = 3.0*(y3 - y2);
            c3 = y4 - y1 - c2;
            c2 -= c1;
            y = y1 + z*(c1 + z*(c2 + z*c3));
            return sign * y;
        }
    }
    return 0.0;
}

double bez_d1 (struct bezier_fit *p, double xx)
{
    double x1, x4, y1, y2, y3, y4, z, c1, c2, c3;
    int    i, j1;
    double dz_dx;

    if (xx < p->xarr[0]) {
        xx = 2.0 * p->xarr[0] - xx;
    }
    if (xx <= p->xarr[0]) {
        return p->start_slope;
    } else if (xx >= p->xarr[p->npts-1]) {
        return p->end_slope;
    }
    for (i = 0; i < p->npts - 1; i++) {
        x1 = p->xarr[i];
        x4 = p->xarr[i+1];
        if (xx <= x4) {
            dz_dx = 1.0 / (x4 - x1);
            j1 = 3*i;
            y1 = p->ybez_arr[j1];
            y2 = p->ybez_arr[++j1];
            y3 = p->ybez_arr[++j1];
            y4 = p->ybez_arr[++j1];
            z = (xx - x1) / (x4 - x1);
            c1 = 3.0*(y2 - y1);
            c2 = 3.0*(y3 - y2);
            c3 = y4 - y1 - c2;
            c2 -= c1;
            return dz_dx * (c1 + z*(2.0*c2 + 3.0*c3*z));
        }
    }
    return 0.0;
}

double bez_d2 (struct bezier_fit *p, double xx)
{
    double x1, x4, y1, y2, y3, y4, z, c1, c2, c3;
    int    i, j1;
    double sign = 1.0;
    double dz_dx;

    if (xx < p->xarr[0]) {
        xx = 2.0 * p->xarr[0] - xx;
        sign = -1.0;
    }
    if (xx <= p->xarr[0]) {
        return 0.0;
    } else if (xx >= p->xarr[p->npts-1]) {
        return 0.0;
    }
    for (i = 0; i < p->npts - 1; i++) {
        x1 = p->xarr[i];
        x4 = p->xarr[i+1];
        if (xx <= x4) {
            dz_dx = 1.0 / (x4 - x1);
            j1 = 3*i;
            y1 = p->ybez_arr[j1];
            y2 = p->ybez_arr[++j1];
            y3 = p->ybez_arr[++j1];
            y4 = p->ybez_arr[++j1];
            z = (xx - x1) / (x4 - x1);
            c1 = 3.0*(y2 - y1);
            c2 = 3.0*(y3 - y2);
            c3 = y4 - y1 - c2;
            c2 -= c1;
            return sign * dz_dx * dz_dx *
                (6.0 * c3 * z + 2.0 * c2);
        }
    }
    return 0.0;
}

void free_bezier_fit (struct bezier_fit *b)
{
    if (b->xarr) free (b->xarr);
    if (b->ybez_arr) free (b->ybez_arr);
}
