const Int_t numPart = 2;
const Int_t numCharges = 3; // 0: all, 1: positive, -1: negative
bool isRun2Binning = 0;
const Int_t numPtBins = 15;
// const Int_t numPtBins = 6; //Run2 binning
const Int_t numCent = 8;
const Int_t numChoice = 5; // mean, sigma, purity, yield, v2

Int_t ColorPart[numPart] = {kPink + 9, kAzure + 7};
Int_t MarkerPart[numPart] = {20, 33};
Float_t MarkerPartSize[numPart] = {1.5, 2.};
Int_t ColorMult[] = {634, 628, 807, kOrange - 4, 797, 815, 418, 429, 867, 856, 601, kViolet, kPink + 9, kPink + 1, 1};
Float_t SizeMult[] = {2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8};
Float_t SizeMultRatio[] = {1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8};
Int_t MarkerMult[] = {20, 21, 33, 34, 29, 20, 21, 33, 34, 29, 20, 21, 33, 34, 29};
Float_t ScaleFactor[] = {256, 128, 64, 32, 16, 8, 4, 2, 1};

// Double_t PtBins[numPtBins + 1] = {0.8, 1.4, 2, 2.5, 3, 4, 6}; //Run 2 binning
Double_t PtBins[numPtBins + 1] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8};
Int_t CentFT0C[numCent + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
Float_t fCentFT0C[numCent + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
float ftcReso[numCent] = {0.514595, 0.7228, 0.760156, 0.733402, 0.659964, 0.540407, 0.383689, 0.218501}; // values from Sourav
float ftcResoSP[numCent] = {0.805015, 1.06586, 1.10463, 1.07285, 0.984069, 0.82759, 0.602314, 0.34722};

Float_t ParticleMassPDG[numPart] = {1.32171, 1.67245};
TString ParticleName[numPart] = {"Xi", "Omega"};
TString ChargeName[numCharges] = {"Neg", "", "Pos"};
TString ParticleNameLegend[numPart] = {"#Xi^{#pm}", "#Omega^{#pm}"};
TString IsOneOrTwoGauss[2] = {"_OneGaussFit", ""};
TString SIsBkgParab[4] = {"_BkgRetta", "_BkgParab", "_BkgPol3", "_BkgExpo"};

Float_t MinPt[numPart] = {0.8, 1.};
Float_t MaxPt[numPart] = {8, 8};
TString TypeHisto[numChoice] = {"Mean", "Sigma", "Purity", "Yield", "V2Mixed"};
TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "v2"};
TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";

//---------------------------------------------------------
Bool_t ChosenParticleXi = 1;                                 // 1 for Xi, 0 for Omega
Int_t ExtrParticle = !ChosenParticleXi;
TString SinputFileNameSyst = "LHC23_PbPb_pass3_Train218607"; /*"LHC23_PbPb_pass3_Train207098";//"LHC23_PbPb_pass2_Train192773"; // 190305 --> ok for Xi, not ok for Omegas */
TString SinputFileName = "LHC23_PbPb_pass3_Train218607";
//TString SinputFileName = "LHC23zzh_pass3_Train224930";
Int_t ExtrCharge = 0;   // 0: all, 1: positive, -1: negative
Bool_t ExtrBkgType = 1; // 0: pol1, 1:pol2, 2:pol3, 3:expo
Bool_t ExtrUseTwoGauss = 1;
Bool_t isApplyWeights = 1; // weights to flatten the phi distribution of cascades
Int_t v2type = 2;          // 0: v2 - old task version before train 224930, 1: v2 SP, 2: v2 EP
const bool useCommonBDTValue = 1; // common BDT cut for all centralities, set to DefaultBDTscoreCut
const float bdtCut[numCent] = {0.98, 0.98, 0.98, 0.97, 0.97, 0.95, 0.95, 0.95};
const float LimitForV2woFit = 0.97; //purity limit to decide whether to extract v2 from fit or not
//---------------------------------------------------------

// systematic studies
bool ExtrisSysMultTrial = 0; // 1 for systematic studies, 0 for default analysis
const int trialsBDT = 10; // number of trials for the systematic studies related to BDTscore
const int nsigmaBarlow = 2;
const float DefaultBDTscoreCut = 0.98;
const float UpperlimitBDTscoreCut = 1;
const float LowerlimitBDTscoreCut = 0.9;

TString SEtaSysChoice[3] = {"", "_Etagt0", "_Etasm0"}; // all eta, eta > 0, eta < 0
Int_t ExtrEtaSysChoice = 0;                            // 0: all eta, 1: eta > 0, 2: eta < 0

TString SIRChoice[6] = {"", "_544013", "_544392", "_544098", "_544032", "_544184"};
TString SIRValue[6] = {"", "6 kHz", "12 kHz", "18 kHz", "23 kHz", "33 kHz"};
TString inputFileNameIR = "Train207098";

//weights files
TString weightFileName = "PhiWeights/FixedWeights_LHC23_PbPb_pass3_Train218607_Xi.root";

//resolution files
TString ComputeResoFileName = "TreeForAnalysis/AnalysisResults_LHC23zzh_pass3_Train224930.root";