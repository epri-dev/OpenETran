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

#ifndef oetypes_included
#define oetypes_included

/*    defined error codes for oe_exit.  */

#define	ERR_OVERLAP		     1
#define	ERR_MALLOC 		     2
#define	ERR_BUFFER_MISSING   3
#define	ERR_LVDC		     4
#define	ERR_PHASE_READ  	 5
#define	ERR_NPHASES		     6
#define	ERR_NPOLES 		     7
#define	ERR_CABLE_PHASES	 8
#define	ERR_CONDUCTOR_N 	 9
#define	ERR_MISSING_CONDUCTOR   10
#define	ERR_RADIUS 		    11
#define	ERR_HEIGHT 		    12
#define	ERR_UNMATCHED_PAIR  13
#define	ERR_BAD_PAIR    	14
#define	ERR_BAD_POLE    	15
#define	ERR_LT_STOPPED  	16
#define	ERR_ARRTAU_STOPPED  17
#define	ERR_MATH_ALLOC  	18
#define	ERR_MATH_CALC   	19
#define	ERR_OUT_OF_RANGE   	20
#define ERR_BAD_ARR_VI      21
#define ERR_MIXED_LINES     22

void oe_exit (int i);

extern FILE *logfp;

#define ONE_SHOT		0  /* ltengine iteration modes */
#define FIND_CRITICAL_CURRENT	1

#define INVALID_PHASE -1
#define MAX_CIRCUITS 3
#define MAX_PHASES 3   /* Phases per circuit */
#define MAX_GNDS 2     /* Grounds per circuit */
#define MAX_GND_WIRES 6
#define MAX_PHS_WIRES 9
#define MAX_WIRES_HIT 15  /* MAX_PHS_WIRES+MAX_GND_WIRES */
#define MAX_POLE_NODES 16 /* MAX_WIRES_HIT + 1 */
#define MAX_CFOS 45       /* MAX_PHS_WIRES*(MAX_PHSWIRES+1)/2 */

#define MAX_REPORT        3 /*array index values defined below */
#define SINGLE_CONDUCTOR  0 /*flashovers on just one conductor in this circuit */
#define MULTI_CONDUCTOR   1 /*flashovers on two or three conductors in this circuit */
#define MULTI_CIRCUIT     2 /*flashover on at least one conductor in this circuit and at least one other circuit on the same pole */

#define TWOPI         6.2831853
#define EPSILON       0.00001     /* tolerance */
#define V_MIN         1.0e-3
#define Y_SHORT       1.0e3       /* admittance for a "short circuit" */
#define LIGHT         3.0e8       /* speed of light in m/sec */
#define I_GROUND_MIN  100.0       /* minimum current for impulse ground modeling */
#define PRIM_L        2.0e-7  /* primitive inductance coefficient, H/m */
#define DEFAULT_LEAD_R	0.00635  /* ground/arrester lead radius, m */
#define CFKONST       2.815863   /* constants for 1 - cosine wave front */
#define CTKONST       4.0
#define ETKONST       1.442695  /* time constant = ETKONST * tail time */
#define TOKEN_LENGTH	256
#define BUFFER_LENGTH	10000  /* buffer size for transient simulation input */
#define BUFFER_LENGTH_LESS_1	9999
#define LOW_START	0    /* histogram bin numbers for starting critical current iterations */
#define BACKFLASH_START	10
#define ARRESTER_START	19

#define THREE 3

enum plot_type {
	PLT_NONE,
	PLT_CSV,
	PLT_TAB,
	PLT_ELT,
	PLT_MAT }; // MatLab not implemented yet

typedef struct tagLTINSTRUCT {
	double ic;  /* critical current to cause flashover */
	FILE *fp;   /* input file */
	FILE *op;   /* text output file, for DOS only */
	FILE *bp;   /* plot output file */
	int plot_type;
	int stop_on_flashover;
	int iteration_mode;  /* none, critical current, transformer case */
	int first_pole_hit;  /* indicates which poles and wires to find critical current for */
	int last_pole_hit;
	int wire_struck [MAX_WIRES_HIT];  /* >0 if wire is exposed to direct stroke */
} LTINSTRUCT;

typedef LTINSTRUCT *LPLTINSTRUCT;

typedef struct tagLTOUTSTRUCT {
	double SI;      /* highest insulator severity index */
	double energy;  /* highest arrester discharge energy */
	double current;  /* highest arrester discharge current */
	double charge;   /* highest arrester charge */
	double predischarge;  /* highest predischarge current */
	double icritical [MAX_WIRES_HIT]; /* critical current for direct stroke to each wire */
} LTOUTSTRUCT;

typedef LTOUTSTRUCT *LPLTOUTSTRUCT;

extern long nr_iter;
extern int nr_max;

extern FILE *op; /* text output file */
extern FILE *bp; /* plot file */

extern int plot_type;

/* global simulation parameters */

extern int using_network;  /* used time, span, and line tokens */
extern int using_multiple_span_defns;  /* can't just work with span_head */
extern int number_of_nodes;  /* number of nodes at each pole */
extern int number_of_conductors;  /* number of conductors at each pole (<= number_of_nodes) */
extern int number_of_poles; /* number of poles in circuit */
extern double span_length; /* pole span in meters */
extern double dT;  /* simulation time step, seconds */
extern double t;  /* current simulation time, advances by dT */
extern double Tmax;  /* simulation stop time, seconds */
extern int step;  /* current simulation step number */
extern int solution_valid; /* flag for solving arresters at each time step */
extern int left_end_z; /* TRUE if left end of circuit has surge impedance terminations */
extern int right_end_z;

extern double predischarge;  /* maximum predischarge current in pipegaps */
extern double SI, energy, current, charge; /* maximum insulator severity index, and
	maximum arrester energy, current, and charge, of all components in the simulation */
extern int flash_halt, flash_halt_enabled; /* flags to stop simulation if an insulator flashes over */
extern int want_si_calculation;  /* set 1 for solution by bisection, 0 for an estimate */
extern int gi_iteration_mode;

/* transient simulation module */
int lt (LPLTINSTRUCT lt_input, LPLTOUTSTRUCT answers);

#endif