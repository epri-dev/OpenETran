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

#ifndef capacitor_included
#define capacitor_included

extern char capacitor_token[];

struct capacitor {
	double y;  /* C's admittance for trapezoidal integration */
	double yc; /* adjustment for "recursive" solution - see Dommel's papers */
	double h; /* past history current */
	int from;
	int to;
	struct pole *parent;
	struct capacitor *next;
};

extern struct capacitor *capacitor_head, *capacitor_ptr;

int init_capacitor_list (void);
void do_all_capacitors (void (*verb) (struct capacitor *));
void update_capacitor_history (struct capacitor *ptr);
void inject_capacitor_history (struct capacitor *ptr);
void init_capacitor_history (struct capacitor *ptr);
void reset_capacitor (struct capacitor *ptr);
int read_capacitor (void);

#endif
