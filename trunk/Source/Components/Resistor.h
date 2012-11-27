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

#ifndef resistor_included
#define resistor_included

extern char resistor_token[];

struct resistor {
	double Rphase;  /* 60-Hz resistance */
	int from;
	int to;
	struct pole *parent;
	struct resistor *next;
};


extern struct resistor *resistor_head, *resistor_ptr;

int init_resistor_list (void);
int read_resistor (void);

#endif
