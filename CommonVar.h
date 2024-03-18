const Int_t numPart = 7;
const Int_t numChoice = 5; // mean, sigma, purity, yield, efficiency for MC
const Int_t numPtBins = 4;
const Int_t numCent = 3;
Float_t PtBins[numPtBins + 1] = {0, 1.2, 2, 3, 4};
Int_t CentFT0C[numCent + 1] = {10, 30, 50, 90};
Float_t ParticleMassPDG[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};

//---------------------------------------------------------
Bool_t ChosenParticleXi = kTRUE; //Xi, put false for Omega
TString SinputFileName = "16March_New";