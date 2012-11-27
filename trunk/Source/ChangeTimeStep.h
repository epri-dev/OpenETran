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

#ifndef changetimestep_included
#define changetimestep_included

extern char time_token[];
extern char change_dt_token[];

extern double first_dT;
extern double second_dT;
extern double dT_switch_time;
extern int using_second_dT;
extern int dT_switched;

void restore_time_step (void);  /* back to the first_dT */
void change_time_step (void);   /* to the second_dT */

#endif