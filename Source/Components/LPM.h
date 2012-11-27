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

#ifndef lpm_included
#define lpm_included

extern char lpm_token[];

#define LPM_NOT_FLASHED    0
#define LPM_FLASHED        1
#define LPM_DISABLE_FLASH  2

struct lpm {
	double cfo;
	double e0;
	double k;
	double d;
	double xpos;
	double xneg;
	double t_flash;
	double vpk_neg;
	double vpk_pos;
	double SI;
	float *pts;
	int flash_mode;  /* if set to -1, won't flashover */
	int from;
	int to;
	struct pole *parent;
	struct lpm *next;
};

extern struct lpm *lpm_head, *lpm_ptr;

int init_lpm_list (void);
int read_lpm (void);
void do_all_lpms (void (*verb) (struct lpm *));
void check_lpm (struct lpm *ptr);
void print_lpm_data (struct lpm *ptr);
void lpm_answers_cleanup (struct lpm *ptr);
void reset_lpm (struct lpm *ptr);
void move_lpm (struct lpm *ptr, int i);
double estimate_lpm_si (struct lpm *ptr);
double calculate_lpm_si (struct lpm *ptr);

#endif