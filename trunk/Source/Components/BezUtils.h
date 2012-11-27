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

#ifndef bezutils_included
#define bezutils_included

struct bezier_fit {
	int npts;
	int bez_size;
	double start_slope;
	double end_slope;
	double *xarr;
	double *ybez_arr;
	struct bezier_fit *next;
};

int get_bez_size (int size);

struct bezier_fit *allocate_bezier (int npts);

void fill_bezier (struct bezier_fit *p, double *xpts, double *ypts, int use_linear);

double bez_d1 (struct bezier_fit *p, double xx);

double bez_d2 (struct bezier_fit *p, double xx);

double bez_eval (struct bezier_fit *p, double xx);

struct bezier_fit *build_bezier (double *xpts, double *ypts, int npts, int use_linear);

void free_bezier_fit (struct bezier_fit *b);

#endif