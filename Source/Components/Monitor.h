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

#ifndef monitor_included
#define monitor_included

extern struct monitor *monitor_head;

struct monitor {
	int from;
	int to;
	int pole;
	double peak;
	double SI;
	int npts;
	double *pts;
	struct meter *mtr;
	struct insulator *ins_de;
	struct lpm *ins_lpm;
	struct monitor *next;
};

/* public interface for cflash and sdw */

void clear_monitors (void);
int init_monitor_list (void);
struct monitor *add_monitor (int pole, int from, int to, int npts);
double get_monitor_peak (int pole, int from, int to);
double get_monitor_si (int pole, int from, int to);
double *get_monitor_pts (int pole, int from, int to);

/* these are only called within lt */

void do_all_monitors (void (*verb) (struct monitor *));
void find_monitor_links (struct monitor *ptr);
void update_monitor_pts (struct monitor *ptr);
void update_monitor_summary (struct monitor *ptr);
void write_monitor_to_logfile (struct monitor *ptr);

#endif