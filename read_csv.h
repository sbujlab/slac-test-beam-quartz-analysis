
// Reads from a csv file and checks for correct pmt and hv.
// Returns number of channels for a single photo-electron for the given pmt and hv.
float GetSignalFromPmtAndHV(int pmt, int hv) {
	// Reading a float from file
	string pmtString, hvString, valueString;
	int integerValue;
	float floatValue = -1.0;

	// Open file and loop
	ifstream file;
	file.open("signal_by_hv_and_pmt.csv");
	while (file.good()) {
		// Check to see if this is the right line
		getline(file, pmtString, ',');
		getline(file, hvString, ',');
		getline(file, valueString);
		if (hvString.compare("") && std::stoi(hvString) == hv && std::stoi(pmtString) == pmt) {
			return std::stof(valueString);
		}
	}
	// If we don't find a value to return, return 0
	return 0.0;
}

// This function reads the (index)^th value from the line matching (run)
// Index = 1: high voltage
// Index = 2: PMT
// Index = 3: beam energy
// Index = 4: beam energy
// Index = 5: pedestal mean
int GetIntegerFromRun(int run, int index = 1) {
	// Reading run, hv, and pmt integers from file
	string runString, valueString;
	ifstream file;
	file.open("values_by_run.csv");
	while (file.good()) {
		// Find the right line
		getline(file, runString, ',');
		if (runString.compare("") && std::stoi(runString) == run) {
			// Find the right value in this line
			for (int i = 0; i < index; i++) {
				getline(file, valueString, ',');
			}
			file.close();
			return std::stoi(valueString);
		} else getline(file, valueString);
	}
	// If we don't find a value, return -2
	file.close();
	return -2;
}

// This function returns the high voltage corresponding to the run number (run)
int GetHvFromRun(int run) {
	return GetIntegerFromRun(run, 1);
}

// This function returns the PMT index corresponding to the run number (run)
// PMT 1, 2, 3, 4 are Hamamatsu 2" for tandem mount
// PMT 5, 6 are 3" PMTs for showermax (mainly used 5)
// PMT 7 is ET-Tube for Ring 5
// PMT 8 is Hamamatsu-Tube for Ring 5
// PMT 9 is for the SAM
int GetPmtFromRun(int run) {
	return GetIntegerFromRun(run, 2);
}

int GetBeamEnergyFromRun(int run) {
	return GetIntegerFromRun(run, 4);
}

int GetPedestalFromRun(int run) {
	return GetIntegerFromRun(run, 5);
}

int GetSignalFromRun(int run) {
	return GetIntegerFromRun(run, 6);
}

int GetSignalRmsFromRun(int run) {
	return GetIntegerFromRun(run, 7);
}

int GetDownstreamSignalFromRun(int run) {
	return GetIntegerFromRun(run, 8);
}

int GetDownstreamSignalRmsFromRun(int run) {
	return GetIntegerFromRun(run, 9);
}

int GetLowFromRun(int run) {
	return GetIntegerFromRun(run, 8);
}

int GetHighFromRun(int run) {
	return GetIntegerFromRun(run, 9);
}
