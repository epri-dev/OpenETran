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

#ifndef readutils_included
#define readutils_included

extern gsl_vector_int *poles_used; /* array sized number_of_poles, flags poles used for a model's branch connection */
extern gsl_matrix_int *pairs_used; /* array sized number_of_nodes x number_of_nodes, for node pairs used in branch connections */

/* functions for adding branch connections based on the "pairs ..." and
"poles ... " input lines - ltaux.c */

int reset_assignments (void);
int update_assignments (int i, int j, int k);
int next_assignment (int *i, int *j, int *k);

struct pole *new_pole (int location);
void terminate_pole (struct pole *ptr, struct span *defn);

int read_pole_label (void);
int read_phase_label (void);

int read_poles (void);
int read_pairs (void);
int mark_pair (int j, int k);

#endif