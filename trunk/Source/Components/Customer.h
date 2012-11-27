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

#ifndef customer_included
#define customer_included

extern char customer_token[];

struct customer {
	int from;
	int to;
	double Vp; /* voltage on transformer primary */
	double Ihg;  /* current in house ground, this time step */
	double Ix2; /* current in transformer X2 terminal, this time step */
	double Ix2_peak; /* maximum X2 current */
	double integral; /* integral of Vp */
	double Ki; /* Ix2 = Ki times Ihg, plus */
	double Kv; /* ... Kv times integral */
	struct pole *parent; /* pole ground should be connected to this parent */
	struct ground *in;  /* points to the house ground */
	struct customer *next;
};

extern struct customer *customer_head, *customer_ptr;

int init_customer_list (void);
void do_all_customers (void (*verb) (struct customer *));
void update_customer_history (struct customer *ptr);
void print_customer_data (struct customer *ptr);
void reset_customer (struct customer *ptr);
int read_customer (void);
struct customer *find_customer (int at, int from, int to);

#endif
