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
float pmt_analyzer(int runNum, float initialSig = -1.0, int run2 = 0, int run3 = 0, int run4 = 0, int run5 = 0,
						int run6 = 0, int run7 = 0, int run8 = 0, int run9 = 0, int run10 = 0){

	// Define histogram numbers
	Int_t binWidth = 1;
	const int MAX_BIN = 4096;
	int adc_range = 1;

	// Grab initial values from csv files or use defaults
	Float_t initialPed = (float)(GetPedestalFromRun(runNum));
	if (initialSig < 0.0) initialSig = (float)(GetSignalFromRun(runNum));
	Float_t initialSigRms = (float)(GetSignalRmsFromRun(runNum));
	if (initialPed < 0.0) {
		if (runNum < 177) initialPed = 877.0;
		else initialPed = 1298.0;
	}
	if (initialSig < 0.0) {
		initialSig = 155.0;
	}
	if (initialSigRms < 0.0) {
		initialSigRms = sqrt(initialSig);
	}

	// Update bin width depending on signal size and rms
	if (initialSig > 200.0 || initialSigRms > 65.0) binWidth = 2;
	if (initialSig > 300.0 || initialSigRms > 100.0) binWidth = 4;
	if (initialSig > 500.0 || initialSigRms > 150.0) binWidth = 8;

	// Define ADC channels used
	int chanUpstream = 2;
	int chanDownstream = 0;

        // Create histogram
        Int_t bins = MAX_BIN / binWidth + 1;
	TH1F* h_QDC = load_files(adc_range, binWidth, chanUpstream, runNum, run2, run3, run4, run5, run6, run7, run8, run9, run10); 
	Int_t nentries = h_QDC->GetEntries();
	if (nentries < 3) return -1.0;
	else printf("Total entries: %d\n", nentries);

	// Grab fit bounds from user-defined thresholds
	Int_t low = h_QDC->FindFirstBinAbove(2) * binWidth - 20;
	Int_t high = h_QDC->FindLastBinAbove(2) * binWidth + 20;
	printf("range: %d, %d\n", low, high);
	
	// If we are overflowing, just don't even run the fit
	if (high >= MAX_BIN) return -2.0; 

        // Normalize integral of histogram to 1 (same scaling for error bar)
        Int_t sum = h_QDC->GetSum();
        for (Int_t curBin = 0; curBin < bins; curBin++) {
		Float_t curVal = h_QDC->GetBinContent(curBin); 
                h_QDC->SetBinContent(curBin, curVal / (sum * (float)(binWidth)));
                h_QDC->SetBinError(curBin, sqrt(curVal) / (sum * (float)(binWidth)));
        }
	//TF1 *fit_gaus_ped = new TF1("fit_gaus_ped", "gaus", initialPed - 2.0, initialPed + 2.0);
	//h_QDC->Fit(fit_gaus_ped, "RN", "");

	// Define fitting function
	TF1 *fit_func=new TF1("fit_func", the_real_deal_yx, 0, MAX_BIN, 11); 
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
	fit_func->SetParameter(0, 0.1);
	fit_func->SetParameter(1, initialPed);
	fit_func->SetParameter(2, 2.0);
	fit_func->SetParameter(3, 0.1);
	fit_func->SetParameter(4, 1.0);
	fit_func->SetParameter(5, initialSig);
	fit_func->SetParameter(6, initialSigRms);
	fit_func->SetParameter(7, 0.0);
	fit_func->SetParameter(8, 1.0);
	
	// Set bounds on parameters where necessary
	fit_func->SetParLimits(0, 0.0, 1.0);
	fit_func->SetParLimits(7, 0.0, 1.0);
	fit_func->SetParLimits(8, 0.0, 1.0);
	fit_func->SetParLimits(1, initialPed - 5.0, initialPed + 5.0);
	fit_func->SetParLimits(5, initialSig * 0.8, initialSig * 1.2);
	fit_func->SetParLimits(6, initialSigRms * 0.8, initialSigRms * 1.2);

	// Setup histogram for printing
        h_QDC->GetXaxis()->SetTitle("ADC channels");
        h_QDC->GetYaxis()->SetTitle("Counts (integral normalized)");
        h_QDC->SetTitle("Low-light PE fit");
	h_QDC->SetLineColor(1);
	h_QDC->SetMarkerSize(0.7);
	h_QDC->SetMarkerStyle(20);
	h_QDC->GetXaxis()->SetTitle("QDC channel");
	h_QDC->GetYaxis()->SetTitle("Normalized yield");
	h_QDC->GetXaxis()->SetRangeUser(low, high);
	
	// Create canvas and draw histogram
	TCanvas *can = new TCanvas("can","can");
	can->cd();
	gStyle->SetOptFit(1);
	TGaxis::SetMaxDigits(3);
	h_QDC->Draw("same");

	// Define results pointer 
	TFitResultPtr res = h_QDC->Fit(fit_func, "LRS", "", low, high);

	// Create vector and grab return parameters
	Double_t back[11];
	fit_func->GetParameters(back);

	// Plot deconvoluted signal distribution functions
	TF1 *fis_from_fit_pe[10];
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
	TF1 *fis_from_fit_bg = new TF1("fis_from_fit_bg", the_real_deal_yx_bg, 0, MAX_BIN, 11);
	fis_from_fit_bg->SetParameters(back);
        fis_from_fit_bg->SetLineStyle(2);
        fis_from_fit_bg->SetLineColor(7); // Cyan
        fis_from_fit_bg->SetNpx(2000);
	fis_from_fit_bg->Draw("same");
	can->Update();

	// Grab some stats info from the fit
	Double_t chi = fit_func->GetChisquare();
	Int_t ndf = fit_func->GetNDF();
	Int_t nfitpoints = fit_func->GetNumberFitPoints();

	// Get output values
	Double_t wout = back[0];
	Double_t pedout = back[1];
	Double_t pedrmsout = back[2];
	Double_t alphaout = back[3];
	Double_t muout = back[4];
	Double_t sigout = back[5];
	Double_t sigrmsout = back[6];
	Double_t injout = back[7];
	Double_t realout = back[8];
	Double_t wouterr = fit_func->GetParError(0);
	Double_t pedouterr = fit_func->GetParError(1);
	Double_t pedrmsouterr = fit_func->GetParError(2);
	Double_t alphaouterr = fit_func->GetParError(3);
	Double_t muouterr = fit_func->GetParError(4);
	Double_t sigouterr = fit_func->GetParError(5);
	Double_t sigrmsouterr = fit_func->GetParError(6);
	Double_t injouterr = fit_func->GetParError(7);
	Double_t realouterr = fit_func->GetParError(8);

	// Grab the parameters for this run
	int pmt = GetPmtFromRun(runNum);
	int hv = GetHvFromRun(runNum);
	float onePEsig = GetSignalFromPmtAndHV(pmt, hv);
	if (adc_range == 0) onePEsig = onePEsig * 8.00;
	float nPE = sigout / onePEsig;
	// Print out results and a copy of all inputs to check them
	//printf("Run: %s\n", rootFile.c_str());
	printf("HV:  %d\n", hv);
	printf("PMT: %d\n", pmt);
	//printf("Detector:  %s\n", detector.c_str());
	//printf("Energy:  %.1f\n", energy);
	printf("Ped: %.5f +/- %.5f\n", pedout, pedrmsout);
	printf("Sig: %.5f +/- %.5f\n", sigout, sigrmsout);
	printf("#PEs: %.5f +/- %.5f\n", nPE, sigrmsout / onePEsig);
	printf("Resolution: %.5f\n", sigrmsout / sigout);
	
	// Calculate and measure the relative contributions to dist. from PE peaks
	double peak_ratio[5];
	double total_integral = 0.0;
	double integral_width = sigout / 2.0;
	double this_pe_mean;
	for (int i = 0; i < 5; i++) {
		this_pe_mean = pedout + i * sigout;
		peak_ratio[i] = h_QDC->Integral(this_pe_mean - integral_width, this_pe_mean + integral_width);
	}
	//for (int i = 0; i < 5; i++) { total_integral += peak_ratio[i]; }
	//for (int i = 0; i < 5; i++) { peak_ratio[i] = peak_ratio[i] / total_integral; }
	printf("\nPE peak contributions:  calculated, measured\n");
	printf("------------------------------------------------------\n");
	printf("0-PE: \t\t\t%.4f, %.4f\n", poisson_peak_calculator(0, muout), peak_ratio[0]);
	printf("1-PE: \t\t\t%.4f, %.4f\n", poisson_peak_calculator(1, muout), peak_ratio[1]);
	printf("2-PE: \t\t\t%.4f, %.4f\n", poisson_peak_calculator(2, muout), peak_ratio[2]);
	printf("3-PE: \t\t\t%.4f, %.4f\n", poisson_peak_calculator(3, muout), peak_ratio[3]);
	printf("4-PE: \t\t\t%.4f, %.4f\n", poisson_peak_calculator(4, muout), peak_ratio[4]);
	printf("------------------------------------------------------\n\n");

	printf("### for runs_with_signal.csv\n");
	printf("%d,%.2f,%.2f,%.2f,%.6f,%.6f\n\n", runNum, pedout, sigout, nPE, sigrmsout / sigout, muout);

	return nPE;
}

