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

#ifndef arrbez_included
#define arrbez_included

#include "BezUtils.h"

extern char arrbez_token[];

enum arr_size_type {
	arrsize_2pt7_to_48,
	arrsize_54_to_360
};

enum arr_char_type {
	arr_char_fow,
	arr_char_8x20,
	arr_char_36x90,
	arr_char_long
};

enum arr_minmax_type {
	arr_use_vmin,
	arr_use_vmax
};
	
struct bezier_fit *build_arrester (double v10, enum arr_size_type arr_size,
	enum arr_char_type arr_char, enum arr_minmax_type arr_minmax, int use_linear);

struct arrbez {
	double v10;
	double vgap;
	double Uref; // in per-unit of v10
	struct bezier_fit *shape;
	double charge;
	double i_peak;
	double energy;
	double t_start;
	double t_peak;
	double rgap;
	double Gref; // for Cigre model
	double g;    // Cigre turn-on conductance
	double gl;   // dT / L for history update
	double r;    // 1/g + rl + rgap
	double h;    // inductor history current
	double rl;   // 2L / dT
	double amps; // arrester current
	double varr; // arrester voltage, across the block
	int from;
	int to;
	int pole_num_nonlinear;
	struct pole *parent;
	struct arrbez *next;
};

extern struct arrbez *arrbez_head, *arrbez_ptr;

int init_arrbez_list (void);
int read_arrbez (void);
void do_all_arrbezs (void (*verb) (struct arrbez *));
void print_arrbez_data (struct arrbez *ptr);
void update_arrbez_history (struct arrbez *ptr);
void arrbez_answers_cleanup (struct arrbez *ptr);
void reset_arrbez (struct arrbez *ptr);
struct arrbez *find_arrbez (int at, int from, int to);

#endif