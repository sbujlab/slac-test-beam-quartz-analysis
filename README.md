# slac-test-beam-quartz-analysis

Gives code and explanation for analyzing Root files containing quartz cerenkov detector data.

This directory is for storing data recorded at SLAC for the benchmarking of Showermax detector among others.

---------------------------------------------------------------------------

First, before using the analyzer, you must open load_files.h, get to 
line 43, and modify the file path pointing to your rootfiles directory.

Or, alternately, you can make this change from the command line with:
 -  sed -i "43s@rootfiles@/home/user/slacdir/rootfiles@" load_files.h

---------------------------------------------------------------------------

The quartz analysis macro is called pmt_analyzer.c 
and can be called from command line:
 - root -l 'pmt_analyzer.c(177)'
 - root -l 'pmt_analyzer.c(int runNum, float sigSize = 100.0)'

The analyzer can handle up to 10 input run numbers:
 - root -l 'pmt_analyzer.c(177, 155.5, 178, 179, 180, 181, ...)'

pmt_analyzer OUTPUT: 
 - Histogram of ADC distribution with 9-parameter fit including 
   the signal size above pedestal, signal width, electron rate, etc.
 - Detector resolution is printed.
 - # of PE's per electron is printed and returned.
 - Relative peak proportions are printed (calculated AND measured)
   for checking the Poissonian shape of the distribution.


the_real_deal_yx.h:  
 - This file contains the fitting functions for the analyzer.
 - The main fitting function calls two other functions to calculate
   the background and signal.
 - The last two parameters (9 and 10) are fixed numbers that tell the 
   algorithm which peaks to consider in the Poissonian distribution 
   (min and max).

read_csv.h: 
 - This file gives the analyzer access to the run parameters
   through function definintions for reading csv files.

values_by_run.csv:
 - This file holds run number, high voltage, pmt,
   detector, beam energy, pedestal value, signal size,
   and signal rms for all runs.

signal_by_hv_and_pmt.csv:
 - This file holds the signal size for one PE
   for a given PMT and HV.
 - The signal size for one PE with units channels per PE 
   allows for conversion from channels to PEs.
