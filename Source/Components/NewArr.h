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

#ifndef newarr_included
#define newarr_included

#include "BezUtils.h"

extern char newarr_token[];

struct bezier_fit *build_arrester (double v10, enum arr_size_type arr_size,
	enum arr_char_type arr_char, enum arr_minmax_type arr_minmax, int use_linear);

struct newarr {
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
	struct newarr *next;
};

extern struct newarr *newarr_head, *newarr_ptr;

int init_newarr_list (void);
int read_newarr (void);
void do_all_newarrs (void (*verb) (struct newarr *));
void print_newarr_data (struct newarr *ptr);
void update_newarr_history (struct newarr *ptr);
void newarr_answers_cleanup (struct newarr *ptr);
void reset_newarr (struct newarr *ptr);
struct newarr *find_newarr (int at, int from, int to);

#endif