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

#ifndef writeplotfile_included
#define writeplotfile_included

#define VOLT_FLAG  0
#define IARR_FLAG -1
#define IPG_FLAG  -2
#define IHG_FLAG  -3
#define IX2_FLAG  -4
#define IPD_FLAG  -5

/* also re-orders meters by voltage, then current */
void InitializePlotOutput (struct meter *head, double dT, double Tmax);

/* also needs to maintain absolute vmax on all meters */
void WritePlotTimeStep (struct meter *head, double t);

/* because the simulation may stop early on a flashover */
void FinalizePlotHeader (double t, int step);

/* must be called after InitializePlotOutput */
void FinalizePlotTitles (char *line1, char *line2, char *line3, char *line4, char *line5);

#endif