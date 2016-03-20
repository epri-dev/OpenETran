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

#ifndef line_included
#define line_included

extern char conductor_token[];
extern char node_token[];
extern char cable_token[];
extern char line_token[];

struct line {
	struct span *defn;   /* has the impedances and transformations */
						 /* hist matrices dimensioned number_of_conductors x steps. */
	gsl_matrix *hist_left;  /* history currents for waves traveling left to right */
	gsl_matrix *hist_right; /* history currents for waves traveling right to left */
	int alloc_steps;     /* number of allocated time steps in the pole span */
	int steps;           /* number of time steps used in the pole span */
	struct pole *left;   /* line sections have a pole at each end */
	struct pole *right;
	struct line *next;
};

extern struct line *line_head, *line_ptr;
extern struct span *span_head, *span_ptr;

int init_line_list (void);
void do_all_lines (void (*verb) (struct line *));
void update_line_history (struct line *ptr); /* only for non-network systems */
void inject_line_imode (struct line *ptr); /* only for non-network systems */
void inject_line_iphase (struct line *ptr); /* for network systems */
void update_vmode_and_history (struct line *ptr); /* for network systems */
void init_line_history (struct line *ptr);
void connect_lines (void); /* only for non-network systems */
void insert_line (int left_pole, int right_pole, struct span *defn, int travel_steps);
void reset_lines (void);
void print_line_history (struct line *ptr);
int readfile (void);
int read_conductors (struct span *defn);
int read_cables (struct span *defn);

int init_span_list (void);
void do_all_spans (void (*verb) (struct span *));
void reset_span (struct span *ptr);
void read_spans (void);
struct span *find_span (int span_id);
void read_lines (void);

#endif
