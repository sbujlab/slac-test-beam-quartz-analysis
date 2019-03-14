# slac-test-beam-quartz-analysis

This repo provides code and explanation for analyzing Root files containing quartz cerenkov detector data. The data was acquired at SLAC (Stanford Linear Accelerator Center) for the benchmarking of Showermax detector among others. The data is not included in this repo.

---------------------------------------------------------------------------

## First, before using the analyzer, 
you must open load_files.h, get to 
line 43, and modify the file path pointing to your rootfiles directory.

Or, alternately, you can make this change from the command line with:
 -  sed -i "43s@rootfiles@/home/user/slacdir/rootfiles@" load_files.h

---------------------------------------------------------------------------

### The quartz analysis macro pmt_analyzer.c can be called from command line:
 - root -l 'pmt_analyzer.c(177)'
 - root -l 'pmt_analyzer.c(int runNum, float sigSize = -1.0)'
 - Using -1.0 as the initial signal size input causes the analysis
   to read in the initial signal value from values_by_run.csv.

### The analyzer can handle up to 10 input run numbers:
 - root -l 'pmt_analyzer.c(177, 155.5, 178, 179, 180, 181, ...)'
 - root -l 'pmt_analyzer.c(177, -1.0, 178, 179, 180, 181, ...)'
 - All the data is loaded into a single histogram for the analysis.
 - This functionality introduces a bug somewhere that makes the number
   of events incorrect in the fit display.

### Analyzer output: 
 - Histogram of ADC distribution with 9-parameter fit including 
   the signal size above pedestal, signal width, electron rate, etc.
 - Detector resolution is printed.
 - # of PE's per electron is printed and returned.
 - Relative peak proportions are printed (calculated AND measured)
   for checking the Poissonian shape of the distribution.

--------------------------------------------------------------------------

### A modified version of the analyzer is required for the tandem detectors:
 - root -l 'pmt_analyzer_tandem.c(177)'
 - root -l 'pmt_analyzer_tandem.c(int runNum, float upstreamSig = -1.0, float downstreamSig = -1.0, int run2 = 0, int run3 = 0, ...)'
 - This analyzer works identically to pmt_analyzer.c while supporting
   two PMTs recorded on their own ADC channel. 
 - This analyzer also takes in two signal inputs after the first run number.

### A modified version of the analyzer is required for very large signal runs
 - root -l 'pmt_analyzer_stack.c(318, -1.0, 1250, 2300)
 - root -l 'pmt_analyzer_stack.c(int runNum, float initialSig = -1.0, int fit_min = 0, int fit_max = 0)
 - This analyzer performs a simple gaussian fit to the pedestal and a simple gaussian fit 
   to the signal. 
 - The "fit_min" and "fit_max" parameters should encompass the pedestal and first PE peak.

--------------------------------------------------------------------------

### A shell script has been created for analyzing data in a loop (automatically chooses macro)
 - Execute the analyze_quartz.sh script with a list of runs to look at individually:
    * ./analyze_quartz.sh 137 138 139 140 220 330 440 
 - Execute the analyze_quartz.sh script with NO inputs to loop through all runs:
    * ./analyze_quartz.sh
 - While looping, the user must quit Root in order for the next fit to appear.
 - There is a 0.5 second delay after quitting root to allow the user to kill the script:
    * Simply press (Ctrl-c) after entering the (.q) command into Root.

--------------------------------------------------------------------------

## Other files:

### the_real_deal_yx.h:  
 - This file contains the fitting functions for the analyzer.
 - The main fitting function calls two other functions to calculate
   the background and signal.
 - The last two parameters (9 and 10) are fixed numbers that tell the 
   algorithm which peaks to consider in the Poissonian distribution 
   (min and max).

### read_csv.h: 
 - This file gives the analyzer access to the run parameters
   through function definintions for reading csv files.

### values_by_run.csv:
 - This file holds run number, high voltage, pmt,
   detector, beam energy, pedestal value, signal size,
   and signal rms for all runs.
 - This file makes it possible to get a good fit to all the data runs
   while only needing to enter the run number as input (user-friendly).
 - For tandem runs, the last two values are downstream signal mean and rms.
 - For stack runs with tungsten, the last two values are fit low and high bounds.

### signal_by_hv_and_pmt.csv:
 - This file holds the signal size for one PE
   for a given PMT and HV.
 - The signal size for one PE with units channels per PE 
   allows for conversion from channels to PEs.
