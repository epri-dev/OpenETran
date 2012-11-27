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

#ifndef surge_included
#define surge_included

extern char surge_token[];

struct surge { /* current surge with 1-cosine front, with exponential tail */
	double peak; /* peak current */
	double front; /* 30-90 front time */
	double tail;  /* time to half value on the tail */
	double cfront;  /* constant for calculating voltage on the front from cos function */
	double ctail; /* constant for calculating voltage on the tail from exp function */
	double tailadvance; /* offset time for starting the exponential tail */
	double tstart; /* surge start time (usually 0) */
	double tau; /* exponential time constant to simulate the tail */
	int from;
	int to;
	struct pole *parent;
	struct surge *next;
};

extern struct surge *surge_head, *surge_ptr;

int init_surge_list (void);
void do_all_surges (void (*verb) (struct surge *));
void inject_surge (struct surge *ptr);
void move_surge (struct surge *ptr, int i, int j, int k, double fpeak,
	double ftf, double ftt, double ftstart);
int read_surge (void);

#endif
