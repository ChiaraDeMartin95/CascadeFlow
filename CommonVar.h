const Int_t numPart = 2;
const Int_t numPtBins = 13; //8
const Int_t numCent = 8;
const Int_t numChoice = 5; // mean, sigma, purity, yield, v2

Int_t ColorMult[] = {634, 628, 807, kOrange - 4, 797, 815, 418, 429, 867, 856, 601, kViolet, kPink + 9, kPink + 1, 1};
Float_t SizeMult[] = {2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8};
Float_t SizeMultRatio[] = {1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8};
Int_t MarkerMult[] = {20, 21, 33, 34, 29, 20, 21, 33, 34, 29, 20, 21, 33, 34, 29};
Float_t ScaleFactor[] = {256, 128, 64, 32, 16, 8, 4, 2, 1};

//Float_t PtBins[numPtBins + 1] = {0.6, 1.2, 1.6, 2, 2.5, 3, 3.5, 4, 5};
Float_t PtBins[numPtBins + 1] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5};
Int_t CentFT0C[numCent + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
float ftcReso[numCent] = {0.514595, 0.7228, 0.760156, 0.733402, 0.659964, 0.540407, 0.383689, 0.218501}; 

Float_t ParticleMassPDG[numPart] = {1.32171, 1.67245};
TString ParticleName[numPart] = {"Xi", "Omega"};
TString IsOneOrTwoGauss[2] = {"_OneGaussFit", ""};
TString SIsBkgParab[4] = {"_BkgRetta", "_BkgParab", "_BkgPol3", "_BkgExpo"};

Float_t MinPt[numPart] = {0.8, 0.8};
Float_t MaxPt[numPart] = {5, 5};
TString TypeHisto[numChoice] = {"Mean", "Sigma", "Purity", "Yield", "V2"};
TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "v2"};
TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";

//---------------------------------------------------------
Bool_t ChosenParticleXi = 1; //Xi, put false for Omega
TString SinputFileName = "LHC23_PbPb_pass2_Train190305";
Bool_t ExtrBkgType = 1; //0: pol1, 1:pol2, 2:pol3
Bool_t ExtrUseTwoGauss = 1;
Int_t ExtrParticle = !ChosenParticleXi;