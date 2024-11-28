Bool_t isV2 = 1;          // 0 for polarization, 1 for v2
Int_t ChosenParticle = 0; // 0: Xi, 1: Omega, 2: Xi-, 3: Xi+, 4: Omega-, 5: Omega+

const Int_t numPart = 6; // Xi+-, Omega+-, Xi-, Xi+, Omega-, Omega+
bool isRun2Binning = 0;
const Int_t numPtBins = 15;
// const Int_t numPtBins = 6; // Run2 binning
const Int_t numPtBinsEff = 16; // for efficiency
const Int_t numPsiBins = 6;    // bins into which Pz (longitudinal polarization) is computed
const Int_t numCent = 8;
const Int_t numChoice = 9; // mean, sigma, purity, yield, v2, Pzs2, Pzs2 from lambda, Cos2Theta, Cos2Theta from lambda
TString NameAnalysis[2] = {"V2", "Pzs2"};

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
Double_t PtBinsEff[numPtBinsEff + 1] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8};
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
Float_t AlphaH[numPart] = {1, 1, -0.390, 0.371, 0.0154, -0.0181};
// Float_t AlphaHErrors[numPart] = {1, 1, 0.006, sqrt(pow(0.007,2) + pow(0.002,2)), 0.0020, sqrt(pow(0.0028,2) + pow(0.0026,2))};
Float_t AlphaHErrors[numPart] = {1, 1, 0.006, 0.0073, 0.0020, 0.0038};
Float_t CXiToLambda = 0.944;
Float_t AlphaLambda[numPart] = {1, 1, 0.747, -0.757, 0.747, -0.757}; // decay parameter for Lambda -> p pi
TString ParticleNameLegend[numPart] = {"#Xi^{#pm}", "#Omega^{#pm}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};
TString IsOneOrTwoGauss[2] = {"_OneGaussFit", ""};
TString SIsBkgParab[4] = {"_BkgRetta", "_BkgParab", "_BkgPol3", "_BkgExpo"};

Float_t MinPt[numPart] = {0.8, 1., 0.8, 0.8, 1., 1.};
Float_t MaxPt[numPart] = {8, 8, 8, 8, 8, 8};
//Float_t MaxPt[numPart] = {10, 10, 10, 10, 10, 10};
// Float_t MaxPt[numPart] = {6, 6, 6, 6, 6, 6}; // Run 2 binning
TString TypeHisto[numChoice] = {"Mean", "Sigma", "Purity", "Yield", "V2Mixed", "Pzs2Mixed", "Pzs2LambdaFromCMixed", "Cos2ThetaMixed", "Cos2ThetaLambdaFromCMixed"};
TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "v2", "Pz,s2", "Pz,s2", "#LTcos^{2}(#theta)#GT", "#LTcos^{2}(#theta)#GT"};
TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";
TString TitleXCent = "Centrality (%)";
TString TitleYPzs = "P_{z,s2}";

//---------------------------------------------------------
TString SinputFileNameSyst = "LHC23_PbPb_pass3_Train218607";
// TString SinputFileName = "LHC23_PbPb_pass4_Train268801";
TString SinputFileName = "LHC23_PbPb_pass4_Train300924" /*"LHC23_pass4_QC1_Train268802"*/;
TString SinputFileNameEff = "LHC24g3_pass4_Train292305"; /*"LHC24g3_pass4_Train290168";*/
// TString SinputFileName = "LHC23zzh_pass4_test5_Train235645"; //<-- test on pass4 (test5)
// TString SinputFileName = "LHC23zzh_pass4_test3_Train232412"; //<-- test on pass4 (test3)
//  TString SinputFileName = "LHC23zzh_pass4_test3_Train232126"; //<--test on pass4 with occupancy cut
//  TString SinputFileName = "LHC23_PbPb_pass3_Train231308"; //<--the largest dataset so far, new occupancy selection applied
//  TString SinputFileName = "LHC23_PbPb_pass3_Train218607"; //<--large dataset, no info about SP
//    TString SinputFileName = "LHC23zzh_pass3_Train226234_CFW"; //<-- event plane from central FW
//  TString SinputFileName = "LHC23zzh_pass3_Train224930"; //<-- contains info about EP and SP from LF
//   TString SinputFileName = "LHC23zzh_pass3_Train225737_BDT0.7"; //<-- BDT selection > 0.7

Int_t ExtrCharge = 0;   // 0: all, 1: positive, -1: negative
Bool_t ExtrBkgType = 1; // 0: pol1, 1:pol2, 2:pol3, 3:expo
Bool_t ExtrUseTwoGauss = 1;
Bool_t isApplyWeights = 0;        // weights to flatten the phi distribution of cascades
Bool_t ExtrisApplyEffWeights = 0; // weights to take into account efficiency dependence on multiplciity (for v2 only)
Int_t v2type = 2;                 // 0: v2 - old task version before train 224930, 1: v2 SP, 2: v2 EP
const bool useCommonBDTValue = 0; // common BDT cut for all centralities, set to DefaultBDTscoreCut
const float bdtCut[numCent] = {0.98, 0.98, 0.98, 0.97, 0.97, 0.96, 0.96, 0.96};
const float LimitForV2woFit = 0.97; // purity limit to decide whether to extract v2 from fit or not
//---------------------------------------------------------

// systematic studies
bool ExtrisSysMultTrial = 0; // 1 for systematic studies, 0 for default analysis
const int trialsBDT = 10;    // number of trials for the systematic studies related to BDTscore
const int nsigmaBarlow = 2;
const float DefaultBDTscoreCut = 0.98; // 0.98
const float UpperlimitBDTscoreCut = 1;
const float LowerlimitBDTscoreCut = 0.9;

TString SEtaSysChoice[3] = {"", "_Etagt0", "_Etasm0"}; // all eta, eta > 0, eta < 0
Int_t ExtrEtaSysChoice = 0;                            // 0: all eta, 1: eta > 0, 2: eta < 0

TString SIRChoice[6] = {"", "_544013", "_544392", "_544098", "_544032", "_544184"};
TString SIRValue[6] = {"", "6 kHz", "12 kHz", "18 kHz", "23 kHz", "33 kHz"};
TString inputFileNameIR = "Train207098";

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