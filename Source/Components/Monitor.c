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

#include <wtypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "../OETypes.h"
#include "../Parser.h"
#include "../ReadUtils.h"
#include "../WritePlotFile.h"
#include "Pole.h"
#include "Meter.h"
#include "LPM.h"
#include "Insulator.h"
#include "Monitor.h"

struct monitor *monitor_head = NULL;
static struct monitor *monitor_ptr;

static struct monitor *find_monitor (int pole, int from, int to)
{
	struct monitor *ptr = monitor_head;

	if (!ptr) return NULL;

	while ((ptr = ptr->next) != NULL) {
		if (ptr->pole == pole) {
			if (ptr->from == from && ptr->to == to) {
				return (ptr);
			}
			if (ptr->from == to && ptr->to == from) {
				return (ptr);
			}
		}
	}
	return (NULL);
}

/* ammeters have "to" < 0, so passing in "to" >= 0 will match only voltmeters */

static struct meter *find_voltmeter (int pole, int from, int to)
{
	struct meter *ptr = meter_head;
	if (logfp) fprintf (logfp, "find a voltmeter at %d %d-%d\n", pole, from, to);
	while ((ptr = ptr->next) != NULL) {
		if (logfp) fprintf (logfp, "  check meter at %d %d-%d\n", ptr->at, ptr->from, ptr->to);
		if (ptr->at == pole) {
			if (ptr->from == from && ptr->to == to) {
				if (logfp) fprintf (logfp, "   matched\n");
				return (ptr);
			}
			if (ptr->from == to && ptr->to == from) {
				if (logfp) fprintf (logfp, "   matched\n");
				return (ptr);
			}
		}
	}
	return (NULL);
}

static struct insulator *find_insulator (int pole, int from, int to)
{
	struct insulator *ptr = insulator_head;
	while ((ptr = ptr->next) != NULL) {
		if (ptr->parent->location == pole) {
			if (ptr->from == from && ptr->to == to) {
				return (ptr);
			}
			if (ptr->from == to && ptr->to == from) {
				return (ptr);
			}
		}
	}
	return (NULL);
}

static struct lpm *find_lpm (int pole, int from, int to)
{
	struct lpm *ptr = lpm_head;
	while ((ptr = ptr->next) != NULL) {
		if (ptr->parent->location == pole) {
			if (ptr->from == from && ptr->to == to) {
				return (ptr);
			}
			if (ptr->from == to && ptr->to == from) {
				return (ptr);
			}
		}
	}
	return (NULL);
}

/* public interface for cflash and sdw */

void clear_monitors (void)
{
	struct monitor *ptr;

	while (monitor_head) {
		ptr = monitor_head->next;
		if (monitor_head->pts) {
			free (monitor_head->pts);
		}
		free (monitor_head);
		monitor_head = ptr;
	}

	monitor_head = NULL;
}

int init_monitor_list (void)
{
	clear_monitors ();
	
	if ((monitor_head = (struct monitor *) malloc (sizeof *monitor_head)) != NULL) {
		monitor_head->next = NULL;
		monitor_head->pts = NULL;
		monitor_ptr = monitor_head;
		return (0);
	}
	
	if (logfp) fprintf (logfp, "can't initialize monitor list\n");
	oe_exit (ERR_MALLOC);
	
	return (1);
}

struct monitor *add_monitor (int pole, int from, int to, int npts)
{
	struct monitor *ptr;
	
	if ((ptr = (struct monitor *) malloc (sizeof *ptr)) != NULL) {
		ptr->from = from;
		ptr->to = to;
		ptr->pole = pole;
		ptr->peak = ptr->SI = 0.0;;
		ptr->mtr = NULL;
		ptr->ins_de = NULL;
		ptr->ins_lpm = NULL;
		ptr->npts = npts;
		if (npts > 0) {
			ptr->pts = (double *) malloc (npts * sizeof (double));
		} else {
			ptr->pts = NULL;
		}
		ptr->next = NULL;
		monitor_ptr->next = ptr;
		monitor_ptr = ptr;
		return (ptr);
	} else {
		if (logfp) fprintf( logfp, "can't allocate new monitor\n");
		oe_exit (ERR_MALLOC);
	}
	return (NULL);
}

double get_monitor_peak (int pole, int from, int to)
{
	struct monitor *ptr = find_monitor (pole, from, to);
	if (ptr) return ptr->peak;
	return 0.0;
}

double get_monitor_si (int pole, int from, int to)
{
	struct monitor *ptr = find_monitor (pole, from, to);
	if (ptr) return ptr->SI;
	return 0.0;
}

double *get_monitor_pts (int pole, int from, int to)
{
	struct monitor *ptr = find_monitor (pole, from, to);
	if (ptr) return ptr->pts;
	return (NULL);
}

/* these are only called within lt */

void do_all_monitors (void (*verb) (struct monitor *))
{
    struct monitor *ptr = monitor_head;
	
	if (!ptr) return;

    while (ptr = ptr->next) {
        verb (ptr);
    }
}

void find_monitor_links (struct monitor *ptr)
{
	ptr->ins_lpm = find_lpm (ptr->pole, ptr->from, ptr->to);
	ptr->ins_de = find_insulator (ptr->pole, ptr->from, ptr->to);
	ptr->mtr = find_voltmeter (ptr->pole, ptr->from, ptr->to);
}

void update_monitor_pts (struct monitor *ptr)
{
	if (ptr->mtr && ptr->pts && step < ptr->npts) {
		ptr->pts[step] = *(ptr->mtr->v_from) - *(ptr->mtr->v_to);
	}
}

void update_monitor_summary (struct monitor *ptr)
{
	if (ptr->mtr) ptr->peak = ptr->mtr->vmax;
	if (ptr->ins_lpm) {
		ptr->SI = ptr->ins_lpm->SI;
	} else if (ptr->ins_de) {
		ptr->SI = ptr->ins_lpm->SI;
	}
}

void write_monitor_to_logfile (struct monitor *ptr)
{
	if (!logfp) return;
	fprintf (logfp, "monitor %d %d-%d peak kV = %G, SI = %G\n",
		ptr->pole, ptr->from, ptr->to, ptr->peak / 1000.0, ptr->SI);
}
