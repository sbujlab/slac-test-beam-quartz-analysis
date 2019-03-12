
float twopi = 2 * 3.14159265359;

Double_t poisson_model_pe(Double_t *x, Double_t *par){
	Double_t xx = x[0];
	Int_t n = (int)(par[5]);		// Current poisson peak
	Double_t qn = par[0] + n * par[2]; 	// mean of this PE dist
	Double_t sigma_n = sqrt(pow(par[1], 2) + n * pow(par[3], 2));
	Double_t norm1 = 1.0 / (2.0 * pow(sigma_n, 2));
	Double_t norm2 = 1.0 / (sqrt(twopi) * sigma_n);
	Double_t gaus = norm2 * exp(-pow(xx - qn, 2) * norm1);
	return TMath::PoissonI(n, par[4]) * gaus * (1.0 - par[5]); 
}

Double_t poisson_model_bg(Double_t *x, Double_t *par){
	Double_t norm1 = 1.0 / (2.0 * pow(par[1], 2));
	Double_t norm2 = 1.0 / (sqrt(twopi) * par[1]);
	Double_t gaus = norm2 * exp(-pow(x[0] - par[0], 2) * norm1 );
	Double_t real = (1.0 - par[5]) * gaus * TMath::PoissonI(0, par[4]);
	Double_t injected = par[5] * gaus;
	return real + injected;
}

Double_t poisson_model(Double_t *x, Double_t *par){
	Double_t sum = 0.0;
	Double_t initial_par6 = par[6];
	for (int i = (int)(par[6]); i < (int)(par[7]); i++) {
		par[6] = (double)(i);
		sum += poisson_model_pe(x, par);
	}
	par[6] = initial_par6;
	return sum + poisson_model_bg(x, par);
}
