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

#ifndef parser_included
#define parser_included

extern char *sp; /* input buffer */
extern char *sn; /* buffer for a single line of input */

/*  functions to parse input character strings - parser.c */
	     
void init_parser (char *sn);
char *first_token (void);  /* lower case */
char *next_token (void);  /* lower case */
int next_int (int *value);
int next_double (double *value);
int first_int (char *sn, int *value);
int first_double (char *sn, double *value);
char *rest_of_line (void); /* case sensitive */

#endif