// Brady Lowe // lowebra2@isu.edu


// This function takes in at least one run number (up to 10)
// and opens all corresponding files and reads the data
// from the given ADC channel into a single histogram, then
// the histogram is returned
TH1F* load_files(int adc_range = 1, int binWidth = 1, int chan = 2, int run0 = 0, 
	int run1 = 0, int run2 = 0, int run3 = 0, 
	int run4 = 0, int run5 = 0, int run6 = 0,
	int run7 = 0, int run8 = 0, int run9 = 0)
{

	// Parse input
	int count = 0;
	int runs[10] = {run0, run1, run2, run3, run4, run5, 
			run6, run7, run8, run9};
	for (int i = 0; i < 10; i++)
		if (runs[i] > 0) count++;

	// Define histogram numbers
	const int MAX_BIN = 4096;

        // Initialize all variables for creating and filling histograms
	Int_t nentries;
        Int_t bins = MAX_BIN / binWidth + 1;
        Float_t minR = -0.5 * (Float_t) binWidth;
        Float_t maxR = MAX_BIN + 0.5 * (Float_t) binWidth;
	TFile* file;
	string filename;
	TTree* tree;
	TLeaf* leaf;
	TBranch* branch;
        char chanStr[32];
	if (adc_range < 1) sprintf(chanStr, "sbs.sbuscint.ladc%d", chan);
        else sprintf(chanStr, "sbs.sbuscint.hadc%d", chan); 
	TH1F* h_QDC = new TH1F("h_QDC", "QDC spectrum", bins, minR, maxR);

	// Loop over all files to fill the histogram
	for (int fileIndex = 0; fileIndex < count; fileIndex++) {

        	// Open .root data file
		filename = Form("../gems/rootfiles/run_%d.root", runs[fileIndex]);
        	printf("Opening file %s\n", filename.c_str());
        	file = new TFile(filename.c_str());
        	if ( file->IsZombie() ) return h_QDC;

        	// Extract data from file, read into memory
        	tree = (TTree*) file->Get("T");
        	if (tree == NULL) {
        	        printf("ERROR ==> Couldn't find \"T\"\n");
        	        return h_QDC;
        	}

        	// Get number of events
        	nentries = (Int_t) tree->GetEntries();
        	printf("Number of entries:  %d\n", nentries);

		// If this is tandem detector data, get data from 2 channels
		// string check = GetDetectorFromRun(runs[i]);
		// if(check.compare("tandem")) {  }

        	// Create leaf and branch from tree
		leaf = tree->GetLeaf(chanStr);
        	branch = leaf->GetBranch();


        	// Fill histogram
        	for (Int_t entry = 0; entry < branch->GetEntries(); entry++) {
        	        branch->GetEntry(entry);
        	        h_QDC->Fill(leaf->GetValue());
        	}
		file->Close();
	}

	return h_QDC;
}

