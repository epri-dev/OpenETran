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

#ifndef pipegap_included
#define pipegap_included

extern char pipegap_token[];

struct pipegap {
	double v_knee;
	double i_bias;
	double r_slope;
	double i_peak;
	double y;
	double i_past;
	double amps;
	int conducting;
	int from;
	int to;
	struct pole *parent;
	struct pipegap *next;
};

extern struct pipegap *pipegap_head, *pipegap_ptr;

int init_pipegap_list (void);
int read_pipegap (void);
void do_all_pipegaps (void (*verb) (struct pipegap *));
void print_pipegap_data (struct pipegap *ptr);
void check_pipegap (struct pipegap *ptr);
void inject_pipegap (struct pipegap *ptr);
void pipegap_answers_cleanup (struct pipegap *ptr);
void reset_pipegap (struct pipegap *ptr);
struct pipegap *find_pipegap (int at, int from, int to);

#endif