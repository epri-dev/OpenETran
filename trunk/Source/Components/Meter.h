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

#ifndef meter_included
#define meter_included

extern char meter_token[];

extern char **pole_labels; /* 0..number_of_poles labels for ELT graphs */
extern char **phase_labels; /* 0..number_of_nodes labels for ELT graphs */

void set_pole_label (int pole_number, char *label);
void set_phase_label (int phase_number, char *label);

struct meter {
	int from;
	int to; /* if this is negative, it's an ammeter */
	int at; /* pole number - no pointer to the pole struct */
	double vmax; /* maximum voltage (or current) recorded */
	double *v_from; /* points to node voltage monitored at pole "at" */
	double *v_to;
	struct meter *next;
};

extern struct meter *meter_head, *meter_ptr;

int init_meter_list (void);
void do_all_meters (void (*verb) (struct meter *));
void print_meter_data (struct meter *ptr);
void update_meter_peaks (struct meter *ptr);
void reset_meter (struct meter *ptr);
int read_meter (void);
struct meter *add_voltmeter (int i, int j, int k);
struct meter *add_ammeter (int i, int j, int k, double *target);

#endif
