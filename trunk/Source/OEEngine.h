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

#ifndef oeengine_included
#define oeengine_included

/* time step loop and iteration control functions */

void time_step_loops (LPLTOUTSTRUCT answers);

struct icrit_params {
	int pole_number;
	int wire_number;
	LPLTOUTSTRUCT answers;
};
double icrit_function (double i_pk, void *params);
void run_loop_case (int pole_number, int wire_number, double i_pk, double ftf, double ftt, 
	LPLTOUTSTRUCT answers);
void loop_control (LPLTINSTRUCT lt_input, LPLTOUTSTRUCT answers);

#endif