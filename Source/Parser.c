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

/* this module is used to parse transient simulation input from a char
buffer in memory */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Parser.h"

static char seps[] = " \t\n\v\f\r\a\b";  /* for the library function strtok */
static char eat_line_seps[] = "\n\r";  /* for the library function strtok */

static char *t;  /* points to a retrieved token */
static char *ps;  /* points to current char location in the buffer sn */
static char ns [64];  /* maximum size of a single token */

/*  &&&&  input file parsing functions  */

void init_parser (char *sn)
{
	ps = sn;
	while (*ps != '\0' && isspace (*ps)) {
		++ps;
	}
}

/* form a "line" of input starting at the current buffer position, and
stopping at a CR/LF.  Then return the first token from this "line". */
	
char *first_token ()
{
	int i, imax;
	char *ts;

	t = NULL;
	ts = strtok (ps, "\n\r");  /* now ts is a NULL-terminated input line */
	if (ts) {
		while (*ts != '\0' && isspace (*ts)) {
			++ts;
		}
		imax = strlen (ts);
		ps += (imax + 1);      /* get ready for the next input line */
		while (*ps != '\0' && isspace (*ps)) {
			++ps;
		}
		t = strtok (ts, seps); /* initialize library function strtok */
		if (t) {
			if (t[0] == '*') { /* comment line - try again */
				return first_token ();
			}
			imax = strlen (t);
			for (i = 0; i < imax; i++) {
				t[i] = (char) tolower (t[i]);  /* convert token to lower case */
			}
		} else {
			return first_token ();
		}
	}
	return (t);
}

/* pull the next char token out of the buffer, converted to lower case */

char *next_token (void)
{
	int i, imax;
	
	t = strtok (NULL, seps);
	if (t) {
		imax = strlen (t);
		for (i = 0; i < imax; i++) {
			t[i] = (char) tolower (t[i]);
		}
	}
	return (t);
}

/* read the rest of the line, case sensitive - usually a text label */

char *rest_of_line (void)
{
	t = strtok (NULL, eat_line_seps);
	return (t);
}

/* read the next int from buffer - this assumes that the next token
is supposed to be an int - and if it isn't, there is an input error.
Note program will eventually crash if next_token exists but is not an int.
This doesn't happen in LPDWCALC, if input is properly prepared in driver.c
*/

int next_int (int *value)
{
	t = next_token ();
	if (t) {
		strcpy (ns, t);
		*value = atoi (ns);  /* convert copy of token to int */
		return (0);
	} else {
		*value = 0;  /* no more data */
		return (1);
	}
}

/* get the next double out of the buffer - again assuming that the next
token is supposed to be a double
Note program will eventually crash if next_token exists but is not a double.
This doesn't happen in LPDWCALC, if input is properly prepared in driver.c
*/

int next_double (double *value)
{
	t = next_token ();
	if (t) {
		strcpy (ns, t);
		*value = atof (ns);  /* convert copy of token to double */
		return (0);
	} else {
		*value = 0.0;  /* no more data */
		return (1);
	}
}

/* get the first integer out of the buffer - as called from ltread.c, this
also initializes the parser.  Otherwise, identical to next_int. */

int first_int (char *sn, int *value)
{
	init_parser (sn);
	t = first_token ();
	if (t) {
		strcpy (ns, t);
		*value = atoi (ns);
		return (0);
	} else {
		*value = 0;
		return (1);
	}
}

/* get the first double out of the buffer - as called from ltread.c, this
also initializes the parser.  Otherwise, identical to next_double. */

int first_double (char *sn, double *value)
{
	init_parser (sn);
	t = first_token ();
	if (t) {
		strcpy (ns, t);
		*value = atof (ns);
		return (0);
	} else {
		*value = 0.0;
		return (1);
	}
}
