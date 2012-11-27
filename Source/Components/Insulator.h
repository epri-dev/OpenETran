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

#ifndef insulator_included
#define insulator_included

extern char insulator_token[];

struct insulator {
	double cfo; /* critical flashover voltage, volts */
	double vb; /* minimum breakdown voltage, usually 0 */
	double beta; /* exponent for destructive effect integration */
	double de_pos; /* integrated destructive effect for positive polarity voltage */
	double de_neg; /* integrated de for negative polarity */
	double de_max; /* max of de_pos and de_neg */
	double t_flash; /* time flashover occurred */
	double SI;
	int flashed;  /* TRUE if conducting */
	int from;
	int to;
	struct pole *parent;
	struct insulator *next;
};

extern struct insulator *insulator_head, *insulator_ptr;

int init_insulator_list (void);
void do_all_insulators (void (*verb) (struct insulator *));
void check_insulator (struct insulator *ptr);
void print_insulator_data (struct insulator *ptr);
void insulator_answers_cleanup (struct insulator *ptr);
void reset_insulator (struct insulator *ptr);
void move_insulator (struct insulator *ptr, int i);
int read_insulator (void);

#endif
