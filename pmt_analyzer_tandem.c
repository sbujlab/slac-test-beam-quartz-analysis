// Brady Lowe // lowebra2@isu.edu

#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMath.h>
#include <TChain.h>
#include <TH1F.h>
#include <TF1.h>

#include "Math/SpecFunc.h"
#include "read_csv.h"
#include "the_real_deal_yx.h"
#include "load_files.h"

#include <TMinuit.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLine.h>


double poisson_peak_calculator(int n, double mu){
	double numerator = TMath::Exp(-mu) * TMath::Power(mu, n);
	double denominator = TMath::Factorial(n);
	return numerator / denominator;
}

// This function reads the PMT root file from a hard-coded directory and fits the distribution
// to a Poisson convoluted with some gaussians and an exponential.
float pmt_analyzer_tandem(int runNum, float initialSigUpstream = -1.0, float initialSigDownstream = -1.0, int run2 = 0, int run3 = 0, int run4 = 0, int run5 = 0, int run6 = 0, int run7 = 0, int run8 = 0, int run9 = 0, int run10 = 0 ){

	// Define histogram numbers
	Int_t binWidth = 1;
	const int MAX_BIN = 4096;
	int adc_range = 1;

	// Define ADC channels used and which PMTs
	int chanUpstream = 2;
	int chanDownstream = 0;
	Int_t pmt_upstream, pmt_downstream;
	if (runNum < 425) {
		pmt_upstream = 2;
		pmt_downstream = 1;
	} else {
		pmt_upstream = 1; 
		pmt_downstream = 2;
	}
        
	// Grab initial values from csv files or use defaults
	Float_t initialMu = 1.0;
        Float_t initialPed = (float)(GetPedestalFromRun(runNum));
        if (initialSigUpstream < 0.0) initialSigUpstream = (float)(GetSignalFromRun(runNum));
        if (initialSigDownstream < 0.0) initialSigDownstream = (float)(GetDownstreamSignalFromRun(runNum));
        Float_t initialSigRmsUpstream = (float)(GetSignalRmsFromRun(runNum));
        Float_t initialSigRmsDownstream = (float)(GetDownstreamSignalRmsFromRun(runNum));
        if (initialPed < 0.0) {
                if (runNum < 177) initialPed = 877.0;
                else initialPed = 1298.0;
        }
        if (initialSigUpstream < 0.0) {
                initialSigUpstream = 40.0;
        }
        if (initialSigDownstream < 0.0) {
                initialSigDownstream = 40.0;
        }
        if (initialSigRmsUpstream < 0.0) {
                initialSigRmsUpstream = 8.0;
        }
        if (initialSigRmsDownstream < 0.0) {
                initialSigRmsDownstream = 8.0;
        }

        // Create histogram
        Int_t bins = MAX_BIN / binWidth + 1;
	TH1F* qdc_upstream = load_files(adc_range, binWidth, chanUpstream, runNum, run2, run3, run4, run5, run6, run7, run8, run9, run10); 
	TH1F* qdc_downstream = load_files(adc_range, binWidth, chanDownstream, runNum, run2, run3, run4, run5, run6, run7, run8, run9, run10); 
	Int_t nentries = qdc_upstream->GetEntries();
	if (nentries < 3) return -1.0;
	else printf("Total entries: %d\n", nentries);

	// Define all variables to be used
	Int_t low, high, sum, ndf, nfitpoints;
	TFitResultPtr res;
	Double_t back[11], peak_ratio[5];
	TF1 *fis_from_fit_pe[10];
	TF1 *fis_from_fit_bg;
	TF1 *fit_func;
	Double_t chi, integral_width, this_pe_mean, wout, pedout, pedrmsout, alphaout, muout, sigout, sigrmsout, injout, realout, wouterr, pedouterr, pedrmsouterr, alphaouterr, muouterr, sigouterr, sigrmsouterr, injouterr, realouterr;
	Int_t hv = GetHvFromRun(runNum);
	// Grab gain conversion factors for pmts and hv's used
	float onePEsig_upstream = GetSignalFromPmtAndHV(pmt_upstream, hv);
	float onePEsig_downstream = GetSignalFromPmtAndHV(pmt_downstream, hv);
	float nPE_upstream, nPE_downstream;
	Double_t total_integral = 0.0;
	string printString = "";

	// Create canvas 
	TCanvas *can = new TCanvas("can","can", 1600, 500);
	can->Divide(2);
	gStyle->SetOptFit(1);
	TGaxis::SetMaxDigits(3);

	//////////////////////////////////////////////////////////////
	//////// Perform these functions for upstream data 
	//////////////////////////////////////////////////////

	can->cd(1);
	// Grab fit bounds from user-defined thresholds
	low = qdc_upstream->FindFirstBinAbove(150) * binWidth - 20;
	high = qdc_upstream->FindLastBinAbove(2) * binWidth + 20;
	printf("range: %d, %d\n", low, high);
	
	// If we are overflowing, just don't even run the fit
	if (high >= MAX_BIN) return -2.0; 

        // Normalize integral of histogram to 1 (same scaling for error bar)
        sum = qdc_upstream->GetSum();
        for (Int_t curBin = 0; curBin < bins; curBin++) {
		Float_t curVal = qdc_upstream->GetBinContent(curBin); 
                qdc_upstream->SetBinContent(curBin, curVal / (sum * (float)(binWidth)));
                qdc_upstream->SetBinError(curBin, sqrt(curVal) / (sum * (float)(binWidth)));
        }

	// Define fitting function
	fit_func = new TF1("fit_func", the_real_deal_yx, 0, MAX_BIN, 11); 
	fit_func->SetLineColor(4); // Dark blue
	fit_func->SetNpx(2000);
	fit_func->SetLineWidth(2);
	// 11 parameters
	fit_func->SetParName(0, "W");
	fit_func->SetParName(1, "Q0");
	fit_func->SetParName(2, "S0");
	fit_func->SetParName(3, "alpha");
	fit_func->SetParName(4, "mu");
	fit_func->SetParName(5, "Q1");
	fit_func->SetParName(6, "S1");
	fit_func->SetParName(7, "inj");
	fit_func->SetParName(8, "real");
	// These two parameters are fixed, just used as memory for fit
	fit_func->SetParName(9, "MIN_PE");
	fit_func->SetParName(10, "MAX_PE");
	// Min PE must be at least one
	fit_func->FixParameter(9, 1);
	fit_func->FixParameter(10, 10);

	// Initialize the pedestal + rms, and signal + rms
	fit_func->SetParameter(0, 0.1);		// w
	fit_func->SetParameter(1, initialPed);	// pedestal mean
	fit_func->SetParameter(2, 2.0);		// pedestal rms
	fit_func->SetParameter(3, 0.1);		// alpha
	fit_func->SetParameter(4, initialMu);	// mu
	fit_func->SetParameter(5, initialSigUpstream);	// signal mean
	fit_func->SetParameter(6, initialSigRmsUpstream); // signal rms
	fit_func->SetParameter(7, 0.0);		// injected proportion
	fit_func->SetParameter(8, 1.0);		// real porportion
	
	// Set bounds on parameters where necessary
	fit_func->SetParLimits(0, 0.0, 1.0); // w
	fit_func->SetParLimits(7, 0.0, 1.0); // inj
	fit_func->SetParLimits(8, 0.0, 1.0); // real
	fit_func->SetParLimits(1, initialPed - 5.0, initialPed + 5.0);
	fit_func->SetParLimits(5, initialSigUpstream * 0.8, initialSigUpstream * 1.2);
	fit_func->SetParLimits(6, initialSigRmsUpstream * 0.8, initialSigRmsUpstream * 1.2);

	// Setup histogram for printing
        qdc_upstream->GetXaxis()->SetTitle("ADC channels");
        qdc_upstream->GetYaxis()->SetTitle("Counts (integral normalized)");
        qdc_upstream->SetTitle("Upstream quartz spectrum");
	qdc_upstream->SetLineColor(1);
	qdc_upstream->SetMarkerSize(0.7);
	qdc_upstream->SetMarkerStyle(20);
	qdc_upstream->GetXaxis()->SetTitle("QDC channel");
	qdc_upstream->GetYaxis()->SetTitle("Normalized yield");
	qdc_upstream->GetXaxis()->SetRangeUser(low, high);
	qdc_upstream->Draw("same");
	res = qdc_upstream->Fit(fit_func, "LRS", "", low, high);
	fit_func->GetParameters(back);

	// Plot deconvoluted signal distribution functions
	for ( int bb = 1; bb < 10; bb++) {
		// Set current PE peak in consideration
		fit_func->GetParameters(back);
		back[9] = (double)(bb);
        	fis_from_fit_pe[bb] = new TF1(Form("pe_fit_%d", bb), the_real_deal_yx_pe, 0, MAX_BIN, 11);
        	fis_from_fit_pe[bb]->SetParameters(back);
        	fis_from_fit_pe[bb]->SetLineStyle(2);
		fis_from_fit_pe[bb]->SetLineColor(2);
		fis_from_fit_pe[bb]->SetNpx(2000);
		fis_from_fit_pe[bb]->Draw("same");
	}

	// Make background distribution function for printing for user
	fis_from_fit_bg = new TF1("fis_from_fit_bg", the_real_deal_yx_bg, 0, MAX_BIN, 11);
	fis_from_fit_bg->SetParameters(back);
        fis_from_fit_bg->SetLineStyle(2);
        fis_from_fit_bg->SetLineColor(7); // Cyan
        fis_from_fit_bg->SetNpx(2000);
	fis_from_fit_bg->Draw("same");
	can->Update();

	// Grab some stats info from the fit
	chi = fit_func->GetChisquare();
	ndf = fit_func->GetNDF();
	nfitpoints = fit_func->GetNumberFitPoints();

	// Get output values
	wout = back[0];
	pedout = back[1];
	pedrmsout = back[2];
	alphaout = back[3];
	muout = back[4];
	sigout = back[5];
	sigrmsout = back[6];
	injout = back[7];
	realout = back[8];
	wouterr = fit_func->GetParError(0);
	pedouterr = fit_func->GetParError(1);
	pedrmsouterr = fit_func->GetParError(2);
	alphaouterr = fit_func->GetParError(3);
	muouterr = fit_func->GetParError(4);
	sigouterr = fit_func->GetParError(5);
	sigrmsouterr = fit_func->GetParError(6);
	injouterr = fit_func->GetParError(7);
	realouterr = fit_func->GetParError(8);
	nPE_upstream = sigout / onePEsig_upstream;

	// Print out results and a copy of all inputs to check them
	printString = "\n\n=============================\nUPSTREAM ANALYSIS\n---------------------------\n";
	printString += "HV:  " + std::to_string(hv) + "\n";
	printString += "PMT: " + std::to_string(pmt_upstream) + "\n";
	printString += "Ped: " + std::to_string(pedout) + " +/- " + std::to_string(pedrmsout) + "\n";
	printString += "Sig: " + std::to_string(sigout) + " +/- " + std::to_string(sigrmsout) + "\n";
	printString += "#PEs: " + std::to_string(nPE_upstream) + " +/- " + std::to_string(sigrmsout / onePEsig_upstream) + "\n";
	printString += "Resolution: " + std::to_string(sigrmsout / sigout) + "\n";
	
	integral_width = sigout / 2.0;
	for (int i = 0; i < 5; i++) {
		this_pe_mean = pedout + i * sigout;
		peak_ratio[i] = qdc_upstream->Integral(this_pe_mean - integral_width, this_pe_mean + integral_width);
	}
	printString += "\nPE peak contributions:  calculated, measured\n";
	printString += "------------------------------------------------------\n";
	printString += "0-PE: \t\t\t" + std::to_string(poisson_peak_calculator(0, muout)) + ", " + std::to_string(peak_ratio[0]) + "\n";
	printString += "1-PE: \t\t\t" + std::to_string(poisson_peak_calculator(1, muout)) + ", " + std::to_string(peak_ratio[1]) + "\n";
	printString += "2-PE: \t\t\t" + std::to_string(poisson_peak_calculator(2, muout)) + ", " + std::to_string(peak_ratio[2]) + "\n";
	printString += "3-PE: \t\t\t" + std::to_string(poisson_peak_calculator(3, muout)) + ", " + std::to_string(peak_ratio[3]) + "\n";
	printString += "4-PE: \t\t\t" + std::to_string(poisson_peak_calculator(4, muout)) + ", " + std::to_string(peak_ratio[4]) + "\n";
	printString += "------------------------------------------------------\n";

	//////////////////////////////////////////////////////////////
	//////// Perform these functions for downstream data 
	////////////////////////////////////////////////////////

	can->cd(2);
	// Grab fit bounds from user-defined thresholds
	low = qdc_downstream->FindFirstBinAbove(150) * binWidth - 20;
	high = qdc_downstream->FindLastBinAbove(2) * binWidth + 20;
	printf("range: %d, %d\n", low, high);
	
	// If we are overflowing, just don't even run the fit
	if (high >= MAX_BIN) return -2.0; 

        // Normalize integral of histogram to 1 (same scaling for error bar)
        sum = qdc_downstream->GetSum();
        for (Int_t curBin = 0; curBin < bins; curBin++) {
		Float_t curVal = qdc_downstream->GetBinContent(curBin); 
                qdc_downstream->SetBinContent(curBin, curVal / (sum * (float)(binWidth)));
                qdc_downstream->SetBinError(curBin, sqrt(curVal) / (sum * (float)(binWidth)));
        }

	initialPed = 1268;

	// Define fitting function
	fit_func=new TF1("fit_func", the_real_deal_yx, 0, MAX_BIN, 11); 
	fit_func->SetLineColor(4); // Dark blue
	fit_func->SetNpx(2000);
	fit_func->SetLineWidth(2);
	// 11 parameters
	fit_func->SetParName(0, "W");
	fit_func->SetParName(1, "Q0");
	fit_func->SetParName(2, "S0");
	fit_func->SetParName(3, "alpha");
	fit_func->SetParName(4, "mu");
	fit_func->SetParName(5, "Q1");
	fit_func->SetParName(6, "S1");
	fit_func->SetParName(7, "inj");
	fit_func->SetParName(8, "real");
	// These two parameters are fixed, just used as memory for fit
	fit_func->SetParName(9, "MIN_PE");
	fit_func->SetParName(10, "MAX_PE");
	// Min PE must be at least one
	fit_func->FixParameter(9, 1);
	fit_func->FixParameter(10, 10);

	// Initialize the pedestal + rms, and signal + rms
	fit_func->SetParameter(0, 0.1);		// w
	fit_func->SetParameter(1, initialPed);	// pedestal mean
	fit_func->SetParameter(2, 2.0);		// pedestal rms
	fit_func->SetParameter(3, 0.01);	// alpha
	fit_func->SetParameter(4, initialMu);	// mu
	fit_func->SetParameter(5, initialSigDownstream);	// signal mean
	fit_func->SetParameter(6, initialSigRmsDownstream); // signal rms
	fit_func->SetParameter(7, 0.0);		// injected proportion
	fit_func->SetParameter(8, 1.0);		// real porportion
	
	// Set bounds on parameters where necessary
	fit_func->SetParLimits(0, 0.0, 1.0); // w
	fit_func->SetParLimits(7, 0.0, 1.0); // inj
	fit_func->SetParLimits(8, 0.0, 1.0); // real
	fit_func->SetParLimits(1, initialPed - 5.0, initialPed + 5.0);
	fit_func->SetParLimits(5, initialSigDownstream * 0.8, initialSigDownstream * 1.2);
	fit_func->SetParLimits(6, initialSigRmsDownstream * 0.8, initialSigRmsDownstream * 1.2);

	// Setup histogram for printing
        qdc_downstream->GetXaxis()->SetTitle("ADC channels");
        qdc_downstream->GetYaxis()->SetTitle("Counts (integral normalized)");
        qdc_downstream->SetTitle("Downstream quartz spectrum");
	qdc_downstream->SetLineColor(1);
	qdc_downstream->SetMarkerSize(0.7);
	qdc_downstream->SetMarkerStyle(20);
	qdc_downstream->GetXaxis()->SetTitle("QDC channel");
	qdc_downstream->GetYaxis()->SetTitle("Normalized yield");
	qdc_downstream->GetXaxis()->SetRangeUser(low, high);
	qdc_downstream->Draw("same");
	res = qdc_downstream->Fit(fit_func, "LRS", "", low, high);
	fit_func->GetParameters(back);

	// Plot deconvoluted signal distribution functions
	for ( int bb = 1; bb < 10; bb++) {
		// Set current PE peak in consideration
		fit_func->GetParameters(back);
		back[9] = (double)(bb);
        	fis_from_fit_pe[bb] = new TF1(Form("pe_fit_%d", bb), the_real_deal_yx_pe, 0, MAX_BIN, 11);
        	fis_from_fit_pe[bb]->SetParameters(back);
        	fis_from_fit_pe[bb]->SetLineStyle(2);
		fis_from_fit_pe[bb]->SetLineColor(2);
		fis_from_fit_pe[bb]->SetNpx(2000);
		fis_from_fit_pe[bb]->Draw("same");
	}

	// Make background distribution function for printing for user
	fis_from_fit_bg = new TF1("fis_from_fit_bg", the_real_deal_yx_bg, 0, MAX_BIN, 11);
	fis_from_fit_bg->SetParameters(back);
        fis_from_fit_bg->SetLineStyle(2);
        fis_from_fit_bg->SetLineColor(7); // Cyan
        fis_from_fit_bg->SetNpx(2000);
	fis_from_fit_bg->Draw("same");
	can->Update();

	// Grab some stats info from the fit
	chi = fit_func->GetChisquare();
	ndf = fit_func->GetNDF();
	nfitpoints = fit_func->GetNumberFitPoints();

	// Get output values
	wout = back[0];
	pedout = back[1];
	pedrmsout = back[2];
	alphaout = back[3];
	muout = back[4];
	sigout = back[5];
	sigrmsout = back[6];
	injout = back[7];
	realout = back[8];
	wouterr = fit_func->GetParError(0);
	pedouterr = fit_func->GetParError(1);
	pedrmsouterr = fit_func->GetParError(2);
	alphaouterr = fit_func->GetParError(3);
	muouterr = fit_func->GetParError(4);
	sigouterr = fit_func->GetParError(5);
	sigrmsouterr = fit_func->GetParError(6);
	injouterr = fit_func->GetParError(7);
	realouterr = fit_func->GetParError(8);
	nPE_downstream = sigout / onePEsig_downstream;

	// Print out results and a copy of all inputs to check them
	printString += "=============================\nDOWNSTREAM ANALYSIS\n---------------------------\n";
	printString += "HV:  " + std::to_string(hv) + "\n";
	printString += "PMT: " + std::to_string(pmt_downstream) + "\n";
	printString += "Ped: " + std::to_string(pedout) + " +/- " + std::to_string(pedrmsout) + "\n";
	printString += "Sig: " + std::to_string(sigout) + " +/- " + std::to_string(sigrmsout) + "\n";
	printString += "#PEs: " + std::to_string(nPE_downstream) + " +/- " + std::to_string(sigrmsout / onePEsig_downstream) + "\n";
	printString += "Resolution: " + std::to_string(sigrmsout / sigout) + "\n";
	
	integral_width = sigout / 2.0;
	for (int i = 0; i < 5; i++) {
		this_pe_mean = pedout + i * sigout;
		peak_ratio[i] = qdc_downstream->Integral(this_pe_mean - integral_width, this_pe_mean + integral_width);
	}
	printString += "\nPE peak contributions:  calculated, measured\n";
	printString += "------------------------------------------------------\n";
	printString += "0-PE: \t\t\t" + std::to_string(poisson_peak_calculator(0, muout)) + ", " + std::to_string(peak_ratio[0]) + "\n";
	printString += "1-PE: \t\t\t" + std::to_string(poisson_peak_calculator(1, muout)) + ", " + std::to_string(peak_ratio[1]) + "\n";
	printString += "2-PE: \t\t\t" + std::to_string(poisson_peak_calculator(2, muout)) + ", " + std::to_string(peak_ratio[2]) + "\n";
	printString += "3-PE: \t\t\t" + std::to_string(poisson_peak_calculator(3, muout)) + ", " + std::to_string(peak_ratio[3]) + "\n";
	printString += "4-PE: \t\t\t" + std::to_string(poisson_peak_calculator(4, muout)) + ", " + std::to_string(peak_ratio[4]) + "\n";
	printString += "=====================================\n\n";

	printf("%s", printString.c_str());

	return nPE_upstream;
}

