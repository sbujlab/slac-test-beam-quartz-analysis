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
float pmt_analyzer_stack(int runNum, float initialSig = -1.0, int low = 0, int high = 0, int run2 = 0, int run3 = 0, int run4 = 0, int run5 = 0,
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
	if (initialSig < 0.0) initialSig = 155.0;
	if (initialSigRms < 0.0) initialSigRms = sqrt(initialSig);
	if (low == 0) low = GetLowFromRun(runNum);
	if (high == 0) high = GetHighFromRun(runNum);

	// Update bin width depending on signal size and rms
	if (initialSig > 100.0 || initialSigRms > 35.0) binWidth = 2;
	if (initialSig > 200.0 || initialSigRms > 65.0) binWidth = 4;
	if (initialSig > 300.0 || initialSigRms > 100.0) binWidth = 8;
	if (initialSig > 500.0 || initialSigRms > 150.0) binWidth = 40;
	if (initialSig > 800.0 || initialSigRms > 200.0) binWidth = 20;

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
	if (low == 0) low = h_QDC->FindFirstBinAbove(2) * binWidth - 20;
	if (high == 0) high = h_QDC->FindLastBinAbove(2) * binWidth + 20;
	printf("range: %d, %d\n", low, high);
	
	// If we are overflowing, just don't even run the fit
	if (high >= MAX_BIN) return -2.0; 

        // Normalize integral of histogram to 1 (same scaling for error bar)
//        Int_t sum = h_QDC->GetSum();
//        for (Int_t curBin = 0; curBin < bins; curBin++) {
//		Float_t curVal = h_QDC->GetBinContent(curBin); 
//                h_QDC->SetBinContent(curBin, curVal / (sum * (float)(binWidth)));
//                h_QDC->SetBinError(curBin, sqrt(curVal) / (sum * (float)(binWidth)));
//        }
	//TF1 *fit_gaus_ped = new TF1("fit_gaus_ped", "gaus", initialPed - 2.0, initialPed + 2.0);
	//h_QDC->Fit(fit_gaus_ped, "RN", "");

	// Define fitting function
	TF1 *ped_func=new TF1("ped_func", "gaus", low, initialPed + 120.0); 
	ped_func->SetLineColor(4); // Dark blue
	ped_func->SetNpx(2000);
	ped_func->SetLineWidth(2);
	ped_func->SetParName(0, "Norm0");
	ped_func->SetParName(1, "Q0");
	ped_func->SetParName(2, "S0");

	TF1 *sig_func=new TF1("sig_func", "gaus", initialPed + 120.0, high); 
	sig_func->SetLineColor(3); 
	sig_func->SetNpx(2000);
	sig_func->SetLineWidth(2);
	sig_func->SetParName(0, "Norm1");
	sig_func->SetParName(1, "Q1");
	sig_func->SetParName(2, "S1");

	// Initialize the pedestal + rms, and signal + rms
	ped_func->SetParameter(0, 2000.0);
	ped_func->SetParameter(1, initialPed);
	ped_func->SetParameter(2, 2.0);
	sig_func->SetParameter(0, 200.0);
	sig_func->SetParameter(1, initialPed + initialSig);
	sig_func->SetParameter(2, initialSigRms);
	
	// Set bounds on parameters where necessary
	ped_func->SetParLimits(1, initialPed - 5.0, initialPed + 5.0);
	sig_func->SetParLimits(1, (initialPed + initialSig) * 0.8, (initialPed + initialSig) * 1.2);
	sig_func->SetParLimits(2, initialSigRms * 0.8, initialSigRms * 1.2);

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
	can->SetLogy();
can->Print(Form("images/stack_%d.png", runNum));

	// Define results pointer 
	TFitResultPtr ped_res = h_QDC->Fit(ped_func, "RS", "", low, initialPed + 40.0);
	TFitResultPtr sig_res = h_QDC->Fit(sig_func, "RS", "", initialPed + 40.0, high);

	// Create vector and grab return parameters
	Double_t ped_back[3];
	Double_t sig_back[3];
	ped_func->GetParameters(ped_back);
	sig_func->GetParameters(sig_back);

	// Grab some stats info from the fit
	Double_t chi = sig_func->GetChisquare();
	Int_t ndf = sig_func->GetNDF();
	Int_t nfitpoints = sig_func->GetNumberFitPoints();

	// Get output values
	Double_t norm1out = ped_back[0];
	Double_t pedout = ped_back[1];
	Double_t pedrmsout = ped_back[2];
	Double_t norm2out = sig_back[0];
	Double_t sigout = sig_back[1] - pedout;
	Double_t sigrmsout = sig_back[2];
	Double_t norm1outerr = ped_func->GetParError(0);
	Double_t pedouterr = ped_func->GetParError(1);
	Double_t pedrmsouterr = ped_func->GetParError(2);
	Double_t norm2outerr = sig_func->GetParError(0);
	Double_t sigouterr = sig_func->GetParError(1);
	Double_t sigrmsouterr = sig_func->GetParError(2);

	// Grab the parameters for this run
	int pmt = GetPmtFromRun(runNum);
	int hv = GetHvFromRun(runNum);
	float onePEsig = GetSignalFromPmtAndHV(pmt, hv);
	if (adc_range == 0) onePEsig = onePEsig * 8.00;
	float nPE = sigout / onePEsig;
	// Print out results and a copy of all inputs to check them
	//printf("Run: %s\n", rootFile.c_str());
	printf("Run:  %d\n", runNum);
	printf("HV:  %d\n", hv);
	printf("PMT: %d\n", pmt);
	//printf("Detector:  %s\n", detector.c_str());
	//printf("Energy:  %.1f\n", energy);
	printf("Ped: %.5f +/- %.5f\n", pedout, pedrmsout);
	printf("Sig: %.5f +/- %.5f\n", sigout, sigrmsout);
	printf("#PEs: %.5f +/- %.5f\n", nPE, sigrmsout / onePEsig);
	printf("Resolution: %.5f\n", sigrmsout / sigout);
	
	return nPE;
}

