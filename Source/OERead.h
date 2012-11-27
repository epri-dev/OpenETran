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

#ifndef oeread_included
#define oeread_included

int readfile (void);

void reset_system (void); /* reset past history terms for critical current iterations - ltaux.c */
int cleanup (void); /* free memory in the linked lists - ltaux.c */

#endif