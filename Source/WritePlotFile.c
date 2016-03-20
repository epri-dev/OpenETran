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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "OETypes.h"
#include "Components/Meter.h"
#include "WritePlotFile.h"

static char delim = ',';

typedef unsigned short USHORT;

#define STO_VNAME_SIZE       9
#define STO_INAME_SIZE       9
#define STO_TITLE_SIZE       80
#define STO_SIGNATURE_SIZE   16

#pragma pack(2)

struct OutputFileHeader
    {
    USHORT      size;               //  size of this header
    char        szSignature [STO_SIGNATURE_SIZE];
    USHORT      idVersionMajor;
    USHORT      idVersionMinor;
    double      dFBase; // Frequency base (Hertz)
    double      dVBase; // Default voltage base (Volts)

//    time_t      tStart;
//    time_t      tFinish;
// sizeof time_t increased in VS2005
//    _int32      tStart;
//    _int32      tFinish;
    int         tStart;
    int         tFinish;

    double      dTStart;
    double      dTFinish;
    double      dDeltaT;

    USHORT      nStep;

    USHORT      nVoltage;
    USHORT      nCurrent;

    USHORT      sizeVoltageName;    //  size of voltage name
    USHORT      sizeCurrentName;    //  size of current name

    long        idxVoltageNames;    //  beginning of voltage name list
    long        idxCurrentNames;    //  beginning of current name list
    long        idxBaseData;        //  beginning of base voltage info (0 if none)
    long        idxData;            //  beginning of output data

    char        szTitle1 [STO_TITLE_SIZE];
    char        szTitle2 [STO_TITLE_SIZE];
    char        szTitle3 [STO_TITLE_SIZE];
    char        szTitle4 [STO_TITLE_SIZE];
    char        szTitle5 [STO_TITLE_SIZE];
    };

static struct OutputFileHeader ofh;

static struct meter *CopyMeter (struct meter *ptr,
	struct meter **first_new, struct meter **last_mtr)
{
	struct meter *next_mtr;

	next_mtr = (struct meter *) malloc (sizeof *next_mtr);
	if (!next_mtr) {
		printf ("can't allocate new voltmeter\n");
		exit (EXIT_FAILURE);
	}
	next_mtr->at = ptr->at;
	next_mtr->from = ptr->from;
	next_mtr->to = ptr->to;
	next_mtr->vmax = ptr->vmax;
	next_mtr->v_from = ptr->v_from;
	next_mtr->v_to = ptr->v_to;
	next_mtr->next = NULL;
	
	if (!(*first_new)) {
		*first_new = next_mtr;
	}
	if (*last_mtr) {
		(*last_mtr)->next = next_mtr;
	}
	*last_mtr = next_mtr;
	return next_mtr;
}

static void SortMetersAndFigureIndices (struct meter *head)
{
	struct meter *ptr = head;
	struct meter *first_new = NULL;
	struct meter *last_mtr = NULL;
	struct meter *next_mtr;

/* build a new list of meter in two passes, first for voltage, then current */
	while ((ptr = ptr->next)) {
		if (ptr->to >= 0) {
			next_mtr = CopyMeter (ptr, &first_new, &last_mtr);
		}
	}
	ptr = head;
	while ((ptr = ptr->next)) {
		if (ptr->to < 0) {
			next_mtr = CopyMeter (ptr, &first_new, &last_mtr);
		}
	}

/* clean up and patch in the re-ordered list of meters */
	ptr = head->next;
	while (ptr) {
		next_mtr = ptr->next;
		free (ptr);
		ptr = next_mtr;
	}
	head->next = first_new;
}

void InitializeSTOOutput (struct meter *head, double dT, double Tmax)
{
	struct meter *ptr = head;

	ofh.size = sizeof ofh;
    strncpy (ofh.szSignature, "OpenETran 1.00", STO_SIGNATURE_SIZE);
    ofh.idVersionMajor = 2;
    ofh.idVersionMinor = 0;
    ofh.dFBase = 376.999;
    ofh.dVBase = 1.0;

    ofh.tStart  = 0; // time (NULL);
    ofh.tFinish = 0; // time (NULL);

    ofh.dTStart  = 0.0;
    ofh.dTFinish = Tmax;
    ofh.dDeltaT  = dT;

    ofh.nStep = 0;

    ofh.nVoltage = 0;
    ofh.nCurrent = 0;

    ofh.sizeVoltageName = STO_VNAME_SIZE;
    ofh.sizeCurrentName = STO_INAME_SIZE;

    ofh.idxVoltageNames = 0;
    ofh.idxCurrentNames = 0;
    ofh.idxBaseData = 0;
    ofh.idxData = 0;

    ofh.szTitle1[0] = 0;
    ofh.szTitle2[0] = 0;
    ofh.szTitle3[0] = 0;
    ofh.szTitle4[0] = 0;
    ofh.szTitle5[0] = 0;

    strcpy (ofh.szTitle1, "EPRI OpenETran Transient Simulation");

/* figure out the file positions */

	while ((ptr = ptr->next)) {
		if (ptr->to >= 0) {
			++ofh.nVoltage;
		} else {
			++ofh.nCurrent;
		}
	}

	ofh.idxVoltageNames = ofh.size;
	ofh.idxCurrentNames =
		ofh.idxVoltageNames + ofh.sizeVoltageName * ofh.nVoltage;
	ofh.idxData         =
		ofh.idxCurrentNames + ofh.sizeCurrentName * ofh.nCurrent;
}

static void trim_buf (char *buf, unsigned int max_len)
{
	char *pc;

	if (strlen (buf) > max_len) {
		pc = strchr (buf, '_');
		while (pc && *pc) {
			*pc = ' ';
			++pc;
		}
	}
	buf [max_len] = '\0';
}

void WriteSTOHeader (struct meter *head)
{
	struct meter *ptr = head;
	USHORT count = 0;
	char buf [STO_INAME_SIZE + STO_VNAME_SIZE];
	
    fwrite (&ofh, sizeof ofh, 1, bp);
	while ((ptr = ptr->next)) {
		memset (buf, ' ', sizeof buf);
		if (count < ofh.nVoltage) {
			sprintf (buf, "V %s_%s%s",
				pole_labels [ptr->at], phase_labels [ptr->from], phase_labels [ptr->to]);
			trim_buf (buf, STO_VNAME_SIZE - 1);
			fwrite (buf, STO_VNAME_SIZE, 1, bp);
		} else {  /* a current meter */
			if (ptr->to == IARR_FLAG) {
				sprintf (buf, "Ia %s_%s",
					pole_labels [ptr->at], phase_labels [ptr->from]);
			} else if (ptr->to == IPG_FLAG) {
				sprintf (buf, "PG %s_%s",
					pole_labels [ptr->at], phase_labels [ptr->from]);
			} else if (ptr->to == IHG_FLAG) {
				sprintf (buf, "HG %s",
					pole_labels [ptr->at]);
			} else if (ptr->to == IX2_FLAG) {
				sprintf (buf, "X2 %s",
					pole_labels [ptr->at]);
			} else if (ptr->to == IPD_FLAG) {
				sprintf (buf, "PD %s_%s",
					pole_labels [ptr->at], phase_labels [ptr->from]);
			} else {
				sprintf (buf, "Ib %s_%s",
					pole_labels [ptr->at], phase_labels [ptr->from]);
			}
			trim_buf (buf, STO_INAME_SIZE - 1);
			fwrite (buf, STO_INAME_SIZE, 1, bp);
		}
		++count;
	}
}

void FinalizeSTOHeader (double t, int step)
{
	ofh.tFinish = 0; // time (NULL);
	ofh.dTFinish = t;
	if (step > USHRT_MAX) {
		ofh.nStep = USHRT_MAX;
	} else {
		ofh.nStep = step;
	}
	rewind (bp);
	fwrite (&ofh, sizeof ofh, 1, bp);
}

void FinalizeSTOTitles (char *line1, char *line2, char *line3, char *line4, char *line5)
{
    ofh.szTitle1[0] = 0;
    ofh.szTitle2[0] = 0;
    ofh.szTitle3[0] = 0;
    ofh.szTitle4[0] = 0;
    ofh.szTitle5[0] = 0;
    if (line1) strncpy (ofh.szTitle1, line1, STO_TITLE_SIZE - 1);
    if (line2) strncpy (ofh.szTitle2, line2, STO_TITLE_SIZE - 1);
    if (line3) strncpy (ofh.szTitle3, line3, STO_TITLE_SIZE - 1);
    if (line4) strncpy (ofh.szTitle4, line4, STO_TITLE_SIZE - 1);
    if (line5) strncpy (ofh.szTitle5, line5, STO_TITLE_SIZE - 1);
	rewind (bp);
	fwrite (&ofh, sizeof ofh, 1, bp);
}

void WriteSTOTimeStep (struct meter *head, double t)
{
	double volts;
	struct meter *ptr = head;

	fwrite (&t, sizeof t, 1, bp);
	while ((ptr = ptr->next)) {
		volts = *(ptr->v_from) - *(ptr->v_to);
		if (fabs (volts) > fabs (ptr->vmax)) {
			ptr->vmax = volts;
		}
		fwrite (&volts, sizeof volts, 1, bp);
	}
}

// tab or csv delimited plot functions

void WriteTextHeader (struct meter *head)
{
	struct meter *ptr = head;

	fprintf (bp, "Time%c", delim);
	while ((ptr = ptr->next)) {
		if (ptr->to >= 0) {
			fprintf (bp, "P%d:%d-%d", ptr->at, ptr->from, ptr->to);
		} else if (ptr->to == IARR_FLAG) {
			fprintf (bp, "P%d:%d-IARR", ptr->at, ptr->from);
		} else if (ptr->to == IPG_FLAG) {
			fprintf (bp, "P%d:%d-IPG", ptr->at, ptr->from);
		} else if (ptr->to == IHG_FLAG) {
			fprintf (bp, "P%d:%d-IHG", ptr->at, ptr->from);
		} else if (ptr->to == IX2_FLAG) {
			fprintf (bp, "P%d:%d-IX2", ptr->at, ptr->from);
		} else if (ptr->to == IPD_FLAG) {
			fprintf (bp, "P%d:%d-IPIPE", ptr->at, ptr->from);
		}
		if (ptr->next) {
			fputc (delim, bp);
		} else {
			fputc ('\n', bp);
		}
	}
}

void WriteTextTimeStep (struct meter *head, double t)
{
	double volts;
	struct meter *ptr = head;

	fprintf (bp, "%e%c", t, delim);
	while ((ptr = ptr->next)) {
		volts = *(ptr->v_from) - *(ptr->v_to);
		if (fabs (volts) > fabs (ptr->vmax)) {
			ptr->vmax = volts;
		}
		fprintf (bp, "%e", volts);
		if (ptr->next) {
			fputc (delim, bp);
		} else {
			fputc ('\n', bp);
		}
	}
}

// public plot functions

void InitializePlotOutput (struct meter *head, double dT, double Tmax)
{
    SortMetersAndFigureIndices (head);
	if (bp) {
		if (plot_type == PLT_ELT) {
			InitializeSTOOutput (head, dT, Tmax);
			WriteSTOHeader (head);
		} else {
			if (plot_type == PLT_TAB) delim = '\t';
			WriteTextHeader (head);
		}
	}
}

void FinalizePlotHeader (double t, int step)
{
	if (bp) {
		if (plot_type == PLT_ELT) {
			FinalizeSTOHeader (t, step);
		}
	}
}

void FinalizePlotTitles (char *line1, char *line2, char *line3, char *line4, char *line5)
{
	if (bp) {
		if (plot_type == PLT_ELT) {
			FinalizeSTOTitles (line1, line2, line3, line4, line5);
		}
	}
}

void WritePlotTimeStep (struct meter *head, double t)
{
	if (bp) {
		if (plot_type == PLT_ELT) {
			WriteSTOTimeStep (head, t);
		} else {
			WriteTextTimeStep (head, t);
		}
	}
}
