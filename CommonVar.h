Bool_t isV2 = 0;              // 0 for polarization, 1 for v2
Int_t ChosenParticle = 0;     // 0: Xi, 1: Omega, 2: Xi-, 3: Xi+, 4: Omega-, 5: Omega+
Bool_t ExtrisRapiditySel = 0; // 0: |eta| < 0.8, 1: |y| < 0.5 (for Pzs2)

TString STHN[2] = {"", "_FromTHN"};
Bool_t ExtrisFromTHN = 1; // 0: process the tree, 1: process the THnSparse

const Int_t numPart = 6; // Xi+-, Omega+-, Xi-, Xi+, Omega-, Omega+
bool isRun2Binning = 0;
const Int_t numPtBins = 15;
// const Int_t numPtBins = 6; // Run2 binning
const Int_t numPtBinsEff = 15; // for efficiency
const Int_t numPsiBins = 6;    // bins into which Pz (longitudinal polarization) is computed
const Int_t numCent = 8;
const Int_t numChoice = 10; // mean, sigma, purity, yield, v2, Pzs2, Pzs2 from lambda, Cos2Theta, Cos2Theta from lambda, V2MixedCorr
TString NameAnalysis[2] = {"V2", "Pzs2"};
TString RapidityCoverage[2] = {"Eta08", "Y05"};

Int_t ColorPart[numPart] = {kPink + 9, kAzure + 7, kPink + 1, kPink - 9, kAzure + 3, kAzure - 3};
Int_t MarkerPart[numPart] = {20, 33, 20, 20, 33, 33};
Float_t MarkerPartSize[numPart] = {1.5, 2., 1.5, 1.5, 2., 2.};
Int_t ColorMult[] = {634, 628, 807, kOrange - 4, 797, 815, 418, 429, 867, 856, 601, kViolet, kPink + 9, kPink + 1, 1};
Float_t SizeMult[] = {2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8};
Float_t SizeMultRatio[] = {1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8};
Int_t MarkerMult[] = {20, 21, 33, 34, 29, 20, 21, 33, 34, 29, 20, 21, 33, 34, 29};
Float_t ScaleFactor[] = {256, 128, 64, 32, 16, 8, 4, 2, 1};

// Double_t PtBins[numPtBins + 1] = {0.8, 1.4, 2, 2.5, 3, 4, 6}; // Run 2 binning
Double_t PtBins[numPtBins + 1] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8};
Double_t PtBinsEff[numPtBinsEff + 1] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8};
Int_t CentFT0C[numCent + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
Double_t fCentFT0C[numCent + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
Double_t dNdEtaAbhi[numCent] = {(2080. + 1697.) / 2, 1274, 862, 566, 355, 208, 112, 54}; // values from Abhi
Double_t dNdEtaAbhiErr[numCent] = {63, 40, 27, 19, 13, 8, 5, 3};
Double_t v2PubRun2[numCent] = {(0.02839 + 0.04566) / 2, 0.06551, 0.08707, 0.0991, 0.10414, 0.10286, 0.09746, 0.08881}; // values from Run2 https://arxiv.org/pdf/1602.01119

float ftcResoSourav[numCent] = {0.514595, 0.7228, 0.760156, 0.733402, 0.659964, 0.540407, 0.383689, 0.218501}; // values from Sourav
// float ftcResoSP[numCent] = {0.805015, 1.06586, 1.10463, 1.07285, 0.984069, 0.82759, 0.602314, 0.34722};

Float_t ParticleMassPDG[numPart] = {1.32171, 1.67245, 1.32171, 1.32171, 1.67245, 1.67245}; // Xi+-, Omega+-, Xi-, Xi+, Omega-, Omega+
TString ParticleName[numPart] = {"Xi", "Omega", "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};
// from PDG 2024 //for Xi+ and Omega+-, set to 1 as it has no meaning:
// Float_t AlphaH[numPart] = {1, 1, -0.390, 0.371, 0.0154, -0.0181};
Float_t AlphaH[numPart] = {1, 1, -0.390, 0.371, 0.0154, -0.018};
// Float_t AlphaHErrors[numPart] = {1, 1, 0.006, sqrt(pow(0.007,2) + pow(0.002,2)), 0.0020, sqrt(pow(0.0028,2) + pow(0.0026,2))};
Float_t AlphaHErrors[numPart] = {1, 1, 0.007, 0.007, 0.0020, 0.004};
Float_t CXiToLambda = 0.925;
Float_t AlphaLambda[numPart] = {1, 1, 0.747, -0.757, 0.747, -0.757}; // decay parameter for Lambda -> p pi
TString ParticleNameLegend[numPart] = {"#Xi^{#pm}", "#Omega^{#pm}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};
TString IsOneOrTwoGauss[2] = {"_OneGaussFit", ""};
TString SIsBkgParab[4] = {"_BkgRetta", "_BkgParab", "_BkgPol3", "_BkgExpo"};

Float_t MinPt[numPart] = {0.8, 1., 0.8, 0.8, 1., 1.};
Float_t MaxPt[numPart] = {8, 8, 8, 8, 8, 8};
// Float_t MaxPt[numPart] = {10, 10, 10, 10, 10, 10};
//  Float_t MaxPt[numPart] = {6, 6, 6, 6, 6, 6}; // Run 2 binning
TString TypeHisto[numChoice] = {"Mean", "SigmaWeighted", "Purity", "Yield", "V2Mixed", "Pzs2Mixed", "Pzs2LambdaFromCMixed", "Cos2ThetaNoFit", "Cos2ThetaLambdaFromC", "V2MixedCorr"};
TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "v2", "Pz,s2", "Pz,s2", "#LTcos^{2}(#theta)#GT", "#LTcos^{2}(#theta)#GT", "v2, corr"};
TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";
TString TitleXCent = "Centrality (%)";
TString TitleYPzs = "P_{z,s2}";

//---------------------------------------------------------
// TString SinputFileNameSyst = "LHC23_PbPb_pass3_Train218607"; OLD
// TString SinputFileName = "LHC23_PbPb_pass4_Train300924";
// TString SinputFileName = "LHC23_PbPb_pass4_Train321006"; //All PbPb statistics

// TString SinputFileName = "LHC23_PbPb_pass4_Train332687"; //Partial 2023 PbPb statistics, pass4 training, BDT > 0.8
// TString SinputFileName = "LHC23_PbPb_pass4_Train331632"; //Partial 2023 PbPb statistics, pass4 training
// TString SinputFileName = "LHC23_PbPb_pass4_Train333596"; //All PbPb stata, pass4 training, bug in Pz fixed
// TString SinputFileName = "TestTHN";
// TString SinputFileName = "LHC23_PbPb_pass4_medium_Train346163"; //THN only
// TString SinputFileName = "LHC23_PbPb_pass4_Train347929"; // THN only, BDT score > 0.4, Pzs2 from Lambda
TString SinputFileName = "LHC23_PbPb_pass4_Train354079"; // THN only, BDT score > 0.2, Pzs2 from Lambda
//TString SinputFileName = "LHC23_PbPb_pass4_Train361757"; // THN only, for Pz,s2 of Xi (direct measurement), acceptance applied on the fly

// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train333596";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_medium_Train346163";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train347929";
TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train354079";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train361757";

TString SinputFileNameEff = "LHC24g3_pass4_Train331315";
TString SinputFileNameEffSyst = "LHC24g3_pass4_Train331315";

Bool_t ExtrBkgType = 1; // 0: pol1, 1:pol2, 2:pol3, 3:expo
Bool_t ExtrUseTwoGauss = 1;
Bool_t isApplyWeights = 0;        // weights to flatten the phi distribution of cascades
Bool_t ExtrisApplyEffWeights = 0; // weights to take into account efficiency dependence on multiplciity (for v2 only)
Int_t v2type = 2;                 // 0: v2 - old task version before train 224930, 1: v2 SP, 2: v2 EP

// BDT selections---------------------------------------------------------
const float DefaultBDTscoreCut = 0.96;
const bool useCommonBDTValue = 1; // common BDT cut for all centralities, set to DefaultBDTscoreCut
const float bdtCut[numCent + 1] = {0.95, 0.95, 0.95, 0.85, 0.85, 0.85, 0.85, 0.85, 0.95};

// Variabls used in FitV2OrPol.C macro----------------------
const bool isApplyAcceptanceCorrection = 0; // for recent files, acceptance correction is applied on the fly
const bool useMixedBDTValueInFitMacro = 1;  // variable used in FitV2OrPol.C macro
// if = 1: pt and multiplicity dependent value defined in:
//   - the function DefineMixedBDTValue (for the pt differential measurement) or
//   - BDTscoreCutPtInt (for the integrated pt measurement)
//   - BDTscoreCutPtIntLoosest (for the integrated pt measurement) if isTightMassCut=1
// if = 0: common BDT value defined by DefaultBDTscoreCut
const double BDTscoreCutPtInt[numCent + 1] = {0.8, 0.8, 0.6, 0.52, 0.44, 0.2, 0.2, 0.2, 0.56}; // BDT cut for integrated pt measurement
const float LimitForV2woFit = 0.97;                                                           // purity limit to decide whether to extract v2 from fit or not
bool isTightMassCut = 1;                                                                      // 1 for tight mass cut, 0 for loose mass cut
float Extrsigmacentral[2] = {4.2, 2.1};                                                       // 2.1
const double BDTscoreCutPtIntLoosest[numCent + 1] = {0.96, 0.92, 0.88, 0.76, 0.52, 0.4, 0.24, 0.2, 0.92};
// BDT cut for integrated pt measurement, loosest cut that give a purity > 0.95 within Extrsigmacentral[1];

//---------------------------------------------------------
// BDT scores applied to produce the acceptance correction (chosen in order to have large purity)
bool isProducedAcceptancePlots = 0; // 1 for acceptance production, 0 for default analysis
const float BDTscoreCutAcceptance[numCent + 1] = {0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96};

// systematic studies
bool ExtrisSysMultTrial = 0; // 1 for systematic studies, 0 for default analysis
const int trialsBDT = 20;    // number of trials for the systematic studies related to BDTscore
const float nsigmaBarlow = 1;
const float UpperlimitBDTscoreCut = 1;
const float LowerlimitBDTscoreCut = 0.2;
// const float MinBDTscorePtInt[numCent + 1] = {0.4, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.4}; // minimum BDT value for syst. evaluation
//const float MaxBDTscorePtInt[numCent + 1] = {0.96, 0.96, 0.8, 0.8, 0.8, 0.6, 0.48, 0.48, 0.96}; // maximum BDT value for syst. evaluation
const double MinBDTscorePtInt[numCent + 1] = {0.96, 0.92, 0.88, 0.76, 0.52, 0.4, 0.24, 0.2, 0.92};
const double MaxBDTscorePtInt[numCent + 1] = {0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.8, 0.76, 0.96}; 

//systematic studies on mass cut
const int trialsMassCut = 7; 
bool ExtrisSysMassCut = 0; // 1 for systematic studies, 0 for default analysis
float ExtrLowLimitSysXi[trialsMassCut] = {1.316, 1.317, 1.318, 1.316, 1.316, 1.317, 1.318};
float ExtrUpLimitSysXi[trialsMassCut] = {1.327, 1.326, 1.325, 1.326, 1.325, 1.327, 1.327};

TString SEtaSysChoice[3] = {"", "_Etagt0", "_Etasm0"}; // all eta, eta > 0, eta < 0
Int_t ExtrEtaSysChoice = 0;                            // 0: all eta, 1: eta > 0, 2: eta < 0

TString SIRChoice[6] = {"", "_544013", "_544392", "_544098", "_544032", "_544184"};
TString SIRValue[6] = {"", "6 kHz", "12 kHz", "18 kHz", "23 kHz", "33 kHz"};
TString inputFileNameIR = "Train207098";
//---------------------------------------------------------

// Files to compute resolution
// TString ComputeResoFileNameLF = "TreeForAnalysis/AnalysisResults_LHC23zzh_pass3_Train224930.root";      //<-- LF framework
// TString ComputeResoFileNameLF = "TreeForAnalysis/AnalysisResults_LHC23zzh_pass4_test5_Train235645.root";      //<-- LF framework
TString inputFileResoCFW = "LHC23zzh_pass3_Train226234_CFW";
TString inputFileResoLF = SinputFileName;

// Files with stored resolution
TString ResoFileName_EPLF = "Resolution/Resolution_EP_LF_" + inputFileResoLF;
TString ResoFileName_EPCFW = "Resolution/Resolution_EP_CFW_" + inputFileResoCFW;
TString ResoFileName_SPLF = "Resolution/Resolution_SP_LF_" + inputFileResoLF;
TString ResoFileName_SPCFW = "Resolution/Resolution_SP_CFW_" + inputFileResoCFW;