
float twopi = 2 * 3.14159265359;

// This function takes in the same parameters from the above fit (+1 extra), and returns
// the signal value from ONLY A SINGLE pe contribution (par[9] = n).
// This will allow us to draw a deconvoluted picture of the fit.
Double_t the_real_deal_yx_pe(Double_t *x, Double_t *par){

	// Grab current x value, current pe, and prepare variables
	Double_t xx = x[0];
	Int_t n = (int)(par[9]);
	Double_t qn, sigma_n, term_1, term_11, term_2, term_3, igne, gn, igne_is;

	// Initialize
	qn = 0.;
	sigma_n = 0.;
	term_1 = 0.;
	term_11 = 0.;
	term_2 = 0.;
	term_3 = 0.;
	igne = 0.;
	gn = 0.;

	// Calculate values to be used for this PE 
	qn = par[1] + n * par[5];				// mean of this PE dist
	sigma_n = sqrt(pow(par[2],2) + n * pow(par[6],2));	// sigma of this PE dist
	term_1 = xx - qn - par[3] * pow(sigma_n,2);		// expo and erf argument
	term_11 = xx - qn - par[3] * pow(sigma_n,2) / 2.0;	// Error or correction
	term_2 = par[1] - qn - par[3] * pow(sigma_n,2);		// Erf argument
	term_3 = xx - qn;					// xx shifted by PE mean

	// Calculate igne
	// Depending on which side of the PE distribution we are on, add or subtract in parentheses
	if (term_1 >= 0.) {
		igne = 	par[3] / 2.0 * exp(-par[3] * term_11) * 
			(
				TMath::Erf(fabs(term_2) / sqrt(2.0) / sigma_n) + 
				TMath::Erf(fabs(term_1) / sqrt(2.0) / sigma_n)
			);
	} else {
		igne = 	par[3] / 2.0 * exp(-par[3] * term_11) * 
		(
			TMath::Erf(fabs(term_2) / sqrt(2.0) / sigma_n) - 
			TMath::Erf(fabs(term_1) / sqrt(2.0) / sigma_n)
		);
	}

	// Calculate gn
	gn = exp(-pow(term_3, 2) / 2.0 / pow(sigma_n, 2)) / (sqrt(twopi) * sigma_n);

	// Put it all together and return
	return TMath::PoissonI(n, par[4]) * par[8] * ((1 - par[0]) * gn + par[0] * igne); 
  
}

// This function returns the background (pedestal) distribution
// which includes real 0-pe events, injected 0-pe events, as
// well as exponential decay distribution of discrete background events.
Double_t the_real_deal_yx_bg(Double_t *x, Double_t *par){

	// Initialize variables, grab current x value
	Double_t xx = x[0];
	Double_t qn, sigma_n, term_1, term_2, term_3, igne, igne_is;
	Double_t poisson_is = exp(-par[4]);
	Double_t gaus_is = exp(-pow(xx - par[1], 2) / 2.0 / pow(par[2], 2)) / par[2] / sqrt(twopi);
  
	// If we are to the right of the pedestal, include the exponential.
	if(xx >= par[1]){
	// Use this line to use noExpoPed fit model
	//if(false){
		igne_is = par[3] * exp(-par[3] * (xx - par[1]));
	} else {
		igne_is = 0.;
	}
	
	// Calculate background portion
	Double_t s_bg = poisson_is * par[8] * ((1 - par[0]) * gaus_is + par[0] * igne_is); 
  
	// Add in clock contribution
	Double_t s_clock = par[7] * ((1-par[0]) * gaus_is + par[0] * igne_is);

	// Sum and return
	return s_bg + s_clock;

}

// This function was written by YuXiang Zhao.
// I had to modify it some amount to get it working on our setup.
// The model used here is described in a NIM publication in the "docs" directory.
Double_t the_real_deal_yx(Double_t *x, Double_t *par){
	Double_t s_real_sum = 0.;
	Double_t initial_par9 = par[9];
	for (int i = (int)(par[9]); i < (int)(par[10]); i++) {
		par[9] = (double)(i);
		s_real_sum += the_real_deal_yx_pe(x, par);
	}
	par[9] = initial_par9;
	Double_t s_bg = the_real_deal_yx_bg(x, par);
	return s_real_sum + s_bg;
}
