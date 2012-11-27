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

#ifndef inductor_included
#define inductor_included

extern char inductor_token[];

struct inductor {
	double y;  /* L's admittance for trapezoidal integration */
	double yi; /* adjustments for series R - see Dommel's papers */
	double zi;
	double h;  /* past-history current */
	double ind; /* input inductance - needed to change time step */
	double res; /* input resistance - needed to change time step */
	int from;
	int to;
	struct pole *parent;
	struct inductor *next;
};

extern struct inductor *inductor_head, *inductor_ptr;

int init_inductor_list (void);
void do_all_inductors (void (*verb) (struct inductor *));
void update_inductor_history (struct inductor *ptr);
void inject_inductor_history (struct inductor *ptr);
void init_inductor_history (struct inductor *ptr);
void reset_inductor (struct inductor *ptr);
int read_inductor (void);

#endif
