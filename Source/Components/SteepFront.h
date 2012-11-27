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

#ifndef steepfront_included
#define steepfront_included

#include "BezUtils.h"

extern char steepfront_token[];

struct steepfront {
	double peak;
	double front;
	double tail;
	double tstart;
	double pu_si; /* max steepness in pu of s30-90 */
	double si;  /* max steepness in amps/sec */
	struct bezier_fit *shape;
	int from;
	int to;
	struct pole *parent;
	struct steepfront *next;
};

extern struct steepfront *steepfront_head, *steepfront_ptr;

int init_steepfront_list (void);
int read_steepfront (void);
void do_all_steepfronts (void (*verb) (struct steepfront *));
void inject_steepfront (struct steepfront *ptr);
void move_steepfront (struct steepfront *ptr, int i, int j, int k, double fpeak,
	double ftf, double ftt, double ftstart, double pu_si);

#endif