Electric Power System Transient Simulator, OpenETran

Copyright (c) 1992, 1994, 1998, 2002, 2011, 2012,
Electric Power Research Institute, Inc.
All rights reserved.

Version 1.0.0.2

Introduction
============

During the period from 1990 through 2002, EPRI funded the development of a 
Lightning Protection Design Workstation (LPDW), which was used by many 
utilities to assess the lightning performance of distribution lines.  
Since about 2002, this program has not been available.  EPRI decided to 
release the simulation kernel of LPDW under an open-source license (GPL 
v3), so it may be incorporated into IEEE Flash and other projects.  

OpenETran can presently simulate multi-conductor power lines, insulators, 
surge arresters, non-linear grounds, and lightning strokes.  It 
efficiently calculates energy and charge duty on surge arresters, and 
iterates to find the critical lightning current causing flashover on one 
or more phases.  It is also suitable for use in substation insulation 
coordination.  Capacitor switching, TRV, and other applicitions may be 
added.  

New Features in Version 1.0
===========================

1 - EPRI originally had permission to use code from the Numerical Recipes 
book in LPDW.  These routines have been removed in favor of GSL.  

2 - Added options for comma-delimited and tab-delimted plot files.  

3 - Added options to define current meters apart from the component 
inputs, for arresters, grounds, customers, and pipegaps.  

4 - Provided an optional Excel spreadsheet interface using Visual Basic 
for Applications (VBA) to manage inputs, outputs, and plot data.  

5 - Added critical current iteration mode from the command line and 
spreadsheet interfaces.  

New Features in Version 1.0.0.2
===============================

1 - FIX: arrbez components with no gap now create text outputs.

Third-party Components
======================

This program uses the GNU Scientific Library (GSL), 
available from http://www.gnu.org/software/gsl/

Installation
============

Unzip the contents of this file into a directory of your choice, such as 
c:\openetran.

Source Code
===========

OpenETran source code is available from the following SVN repository:

https://openetran.svn.sourceforge.net/svnroot/openetran 

License
=======

Use of this software is subject to the GPL version 3 open source license.  
The terms are in a file called "license.txt" distributed with the 
software.  
