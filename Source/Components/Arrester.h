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

#ifndef arrester_included
#define arrester_included

extern char arrester_token[];

struct arrester { /* models a single-segment VI characteristic, with optional gap */
	double v_knee; /* voltage where arrester starts conduction */
	double i_bias; /* parallel current source to create arrester discharge voltage */
	double knee_bias; /* i_bias that will create v_knee when flowing through r_slope */
	double gap_bias; /* i_bias that will create v_gap when flowing through r_slope */
	double v_gap; /* gap sparkover voltage (= v_knee for gapless) */
	double r_slope; /* slope of discharge characteristic */
	double charge; /* coulombs discharged */
	double i_peak; /* peak discharge current */
	double energy; /* arrester discharge energy */
	double t_start; /* time arrester started conduction */
	double t_peak; /* time of peak arrester current */
	double y; /* admittance added to the pole y matrix */
	double h; /* past-history current for built-in lead inductance */
	double i; /* prospective arrester current this time step */
	double i_past;  /* injection current for this arrester, previous time step */
	double yr;  /* admittance for r_slope */
	double zl;  /* admittance parameters for the built-in lead inductance */
	double yzl;
	double amps;  /* most recent arrester current */
	int conducting; /* TRUE if arrester conducting */
	int from;
	int to;
	struct pole *parent;
	struct arrester *next;
};

extern struct arrester *arrester_head, *arrester_ptr;

int init_arrester_list (void);
void do_all_arresters (void (*verb) (struct arrester *));
void print_arrester_data (struct arrester *ptr);
void check_arrester (struct arrester *ptr);
void inject_arrester (struct arrester *ptr);
void update_arrester_history (struct arrester *ptr);
void arrester_answers_cleanup (struct arrester *ptr);
void reset_arrester (struct arrester *ptr); /* ltread.c */
int read_arrester (void);
struct arrester *find_arrester (int at, int from, int to);

#endif
