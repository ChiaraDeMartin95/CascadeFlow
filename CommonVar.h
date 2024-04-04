const Int_t numPart = 2;
const Int_t numChoice = 5; // mean, sigma, purity, yield, efficiency for MC
const Int_t numPtBins = 8;
const Int_t numCent = 4;
Float_t PtBins[numPtBins + 1] = {0.6, 1.2, 1.6, 2, 2.5, 3, 3.5, 4, 5};
Int_t CentFT0C[numCent + 1] = {10, 20, 30, 40, 50};
Float_t ParticleMassPDG[numPart] = {1.32171, 1.67245};
TString ParticleName[numPart] = {"Xi", "Omega"};
TString IsOneOrTwoGauss[2] = {"_OneGaussFit", ""};
TString SIsBkgParab[3] = {"_BkgRetta", "_BkgParab", "_BkgPol3"};
/// FT0C resolutions 0.514595, 0.7228, 0.760156, 0.733402, 0.659964, 0.540407, 0.383689, 0.218501
float ftcReso[7] = {0.7228, 0.760156, 0.733402, 0.659964, 0.540407, 0.383689, 0.218501}; /// skipping the first as we start from 10%

//---------------------------------------------------------
Bool_t ChosenParticleXi = 1; //Xi, put false for Omega
TString SinputFileName = "LHC23_PbPb_pass2_Train190305_New";
Bool_t ExtrBkgType = 1; //0: pol1, 1:pol2, 2:pol3
Bool_t ExtrUseTwoGauss = 1;
Int_t ExtrParticle = !ChosenParticleXi;