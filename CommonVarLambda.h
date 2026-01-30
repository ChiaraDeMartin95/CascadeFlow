Bool_t isV2 = 0;              // 0 for polarization, 1 for v2
Int_t ChosenParticle = 6;     // 0: Xi, 1: Omega, 2: Xi-, 3: Xi+, 4: Omega-, 5: Omega+, 6: Lambda + ALambda
Bool_t ExtrisRapiditySel = 0; // 0: |eta| < 0.8, 1: |y| < 0.5 (for Pzs2),
Bool_t ExtrisPartialEta = 0;  // 1: select only 0 < eta < 0.8 (opposite to FT0C)
Int_t ExtrBkgType = 1;        // 0: pol1, 1:pol2, 2:pol3, 3:expo
Int_t ExtrBkgTypeSyst = 1;    // for syst. uncertainty: 0: pol1, 1:pol2, 2:pol3, 3:expo
Bool_t ExtrUseTwoGauss = 1;
Bool_t isApplyWeights = 0;          // weights to flatten the phi distribution of cascades
Bool_t isApplyCentWeight = 1;       // 0 for acceptance
Bool_t ExtrisApplyEffWeights = 0;   // weights to take into account efficiency dependence on multiplciity (for v2 only)
Bool_t ExtrisApplyResoOnTheFly = 1; // 0 for acceptance
Int_t v2type = 2;                   // 0: v2 - old task version before train 224930, 1: v2 SP, 2: v2 EP
Bool_t ExtrisFromTHN = 0;           // 0: process the tree, 1: process the THnSparse
Bool_t isReducedPtBins = 1;
Bool_t isOOCentrality = 1;

const Int_t numPart = 9; // Xi+-, Omega+-, Xi-, Xi+, Omega-, Omega+, Lambda + ALambda, Lambda, AntiLambda
bool isRun2Binning = 0;
// const Int_t numPtBins = 15;
//  const Int_t numPtBins = 6; // Run2 binning
const Int_t numPtBinsEff = 15; // for efficiency
const Int_t numPsiBins = 6;    // bins into which Pz (longitudinal polarization) is computed
const Int_t numCent = 8;
const Int_t numCentLambdaOO = 10;
const Int_t commonNumCent = 10; // the maximum of the two above (?) (numCent for Xi, numCentLambdaOO for Lambda in OO)
// const Int_t numCent_PtDiff = 3; // for pt differential measurement
const Int_t numChoice = 12; // mean, sigma, purity, yield, v2, Pzs2, Pzs2 from lambda, Cos2Theta, Cos2Theta from lambda, V2MixedCorr, Cos2ThetaFromLambdaVsPtLambda

TString sPolFromLambda[2] = {"", "LambdaFromC"};
TString STHN[2] = {"", "_FromTHN"};
TString V2FromFit[2] = {"NoFit", ""};
TString NameAnalysis[2] = {"V2", "Pzs2"};
TString RapidityCoverage[2] = {"Eta08", "Y05"};
TString IsOneOrTwoGauss[2] = {"_OneGaussFit", ""};
TString SIsBkgParab[4] = {"_BkgRetta", "_BkgParab", "_BkgPol3", "_BkgExpo"};
Float_t ParticleMassPDG[numPart] = {1.32171, 1.67245, 1.32171, 1.32171, 1.67245, 1.67245, 1.115683, 1.115683, 1.115683}; // Xi+-, Omega+-, Xi-, Xi+, Omega-, Omega+, Lambda + ALambda
TString ParticleName[numPart] = {"Xi", "Omega", "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus", "Lambda", "LambdaPart", "AntiLambda"};
TString ParticleNameLegend[numPart] = {"#Xi^{#pm}", "#Omega^{#pm}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}", "#Lambda + #overline{#Lambda}", "#Lambda", "#overline{#Lambda}"};
TString TypeHisto[numChoice] = {"Mean", "SigmaWeighted", "Purity", "Yield", "V2Mixed", "Pzs2Mixed", "Pzs2LambdaFromCMixed", "Cos2ThetaNoFit", "Cos2ThetaLambdaFromC", "V2MixedCorr", "Cos2ThetaLambdaFromCVsPt", "Cos2ThetaLambdaFromCVsEta"};
TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "v2", "Pz,s2", "Pz,s2", "#LTcos^{2}(#theta*_{#Lambda})#GT", "#LTcos^{2}(#theta*_{p})#GT", "v2, corr", "#LTcos^{2}(#theta*_{p})#GT", "#LTcos^{2}(#theta*_{p})#GT"};
TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";
TString TitleXCent = "Centrality (%)";
TString TitleYPzs = "#it{P}_{z,s2}";

// Centrality
Int_t CentFT0C[numCent + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80}; //{0, 30, 50, 80}; // for pt differential measurement
Double_t fCentFT0C[numCent + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
Int_t CentFT0CLambdaOO[numCentLambdaOO + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
Double_t fCentFT0CLambdaOO[numCentLambdaOO + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
Double_t dNdEtaAbhi[numCent] = {(2047. + 1668.) / 2, 1253, 848, 559, 351, 205, 110, 53}; // values from Abhi
Double_t dNdEtaAbhiErr[numCent] = {(54. + 42.)/2, 33, 25, 19, 14, 11, 8, 5};
Double_t dNdEtaOOPrel[numCentLambdaOO] = {(126.6660 + 106.8340) / 2, 87.2877, 67.1562, 51.1201, 37.8919, 26.9060, 19.03, 13.22, 8.50, 0}; // approved up to 60%
Double_t dNdEtaOOErrPrel[numCentLambdaOO] = {(4.23 + 3.44) / 2, 2.76, 2.08, 1.55, 1.13, 0.81, 0.57, 0.39, 0.23, 0};         // approved up to 60%
Double_t dNdEtaOOErrPrelSyst[numCentLambdaOO] = {(4.1009 + 3.3921) / 2, 2.8101, 2.2507, 1.8847, 1.6370, 1.4909, 0.57, 0.39, 0.23, 0};         // approved up to 60%
Double_t dNdEtaOO[numCentLambdaOO] = {(126.6660 + 106.8340) / 2, 87.2877, 67.1562, 51.1201, 37.8919, 26.9060, 19.91, 14.87, 11.11, 0};    // from 60% to 90%, extrapolated with MultVsCent.C macro
Double_t dNdEtaOOErr[numCentLambdaOO] = {(0.0335 + 0.0250) / 2, 0.0143, 0.0117, 0.0099, 0.0084, 0.0069, 0.58, 0.54, 0.48, 0};             // from 60% to 90%, extrapolated with MultVsCent.C macro
Double_t dNdEtaOOErrSyst[numCentLambdaOO] = {(4.1009 + 3.3921) / 2, 2.8101, 2.2507, 1.8847, 1.6370, 1.4909, 0.58, 0.54, 0.48, 0};             // from 60% to 90%, extrapolated with MultVsCent.C macro
Double_t dNdEtaNeNe[2] = {105.59, 20.63};                                                                                // for Junlee results. Averages computed from analysis note (0-40%, 40-90% even if polarization uses 40-100%; multiplicity available only up to 90%)
Double_t dNdEtaNeNeErr[2] = {3.52, 0.69};                                                                                // random reasonable errors assigned

Double_t v2PubRun2[numCent] = {(0.02839 + 0.04566) / 2, 0.06551, 0.08707, 0.0991, 0.10414, 0.10286, 0.09746, 0.08881}; // values from Run2 https://arxiv.org/pdf/1602.01119

// Pt bins
const Int_t numPtBins = 7;
const Int_t numPtBinsReduced = 7;
Double_t PtBinsEff[numPtBinsEff + 1] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8};
Double_t PtBins[numPtBins + 1] = {0.5, 1.0, 1.5, 2, 2.5, 3, 4, 8};
Float_t MinPt[numPart] = {0.8, 1., 0.8, 0.8, 1., 1., 0.5, 0.5, 0.5};
Float_t MaxPt[numPart] = {8, 8, 8, 8, 8, 8, 8, 8, 8};

// Acceptance correction
const Int_t numEtaBins = 16; // was 8
const Int_t numPtBinsLambda = 9;
//Double_t EtaBins[numEtaBins + 1] = {-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8};
Double_t EtaBins[numEtaBins + 1] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
Double_t PtBinsLambda[numPtBinsLambda + 1] = {0.4, 0.8, 1.2, 1.6, 2, 2.5, 3, 4, 6, 10};

// Colors and markers
Int_t ColorPart[numPart] = {kPink + 9, kAzure + 7, kPink + 1, kPink - 9, kAzure + 3, kAzure - 3, kOrange, kOrange - 2, kOrange - 4};
Int_t MarkerPart[numPart] = {20, 33, 20, 20, 33, 33, 33, 33, 33};
Float_t MarkerPartSize[numPart] = {1.5, 2., 1.5, 1.5, 2., 2., 2., 2., 2};
Int_t ColorMult[] = {634, 628, 807, kOrange - 4, 797, 815, 418, 429, 867, 856, 601, kViolet, kPink + 9, kPink + 1, 1, kRed + 1, kGray + 1};
Float_t SizeMult[] = {2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8};
Float_t SizeMultRatio[] = {1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8};
Int_t MarkerMult[] = {20, 21, 33, 34, 29, 20, 21, 33, 34, 29, 20, 21, 33, 34, 29};
Float_t ScaleFactor[] = {256, 128, 64, 32, 16, 8, 4, 2, 1};

// Decay parameters
Float_t AlphaH[numPart] = {1, 1, -0.390, 0.371, 0.0154, -0.018, 1, 1, 1}; // from PDG 2024, for Xi+ and Omega+-, set to 1 as it has no meaning
// Float_t AlphaHErrors[numPart] = {1, 1, 0.006, sqrt(pow(0.007,2) + pow(0.002,2)), 0.0020, sqrt(pow(0.0028,2) + pow(0.0026,2))};
Float_t AlphaHErrors[numPart] = {1, 1, 0.007, 0.007, 0.0020, 0.004, 1, 1, 1};
Float_t CXiToLambda = 0.925;
Float_t AlphaLambda[numPart] = {1, 1, 0.746, -0.758, 0.746, -0.758, 1, 1, 1};     // decay parameter for Lambda -> p pi
Float_t AlphaLambdaErrors[numPart] = {1, 1, 0.008, 0.005, 0.008, 0.005, 1, 1, 1}; // decay parameter for Lambda -> p pi

// File names
// TString SinputFileNameSyst = "LHC23_PbPb_pass3_Train218607"; OLD
// TString SinputFileName = "LHC23_PbPb_pass4_Train300924";
// TString SinputFileName = "LHC23_PbPb_pass4_Train321006"; //All PbPb statistics
// TString SinputFileName = "LHC23_PbPb_pass4_Train332687"; //Partial 2023 PbPb statistics, pass4 training, BDT > 0.8
// TString SinputFileName = "LHC23_PbPb_pass4_Train331632"; //Partial 2023 PbPb statistics, pass4 training
// TString SinputFileName = "LHC23_PbPb_pass4_Train333596"; //All PbPb stata, pass4 training, bug in Pz fixed
// TString SinputFileName = "TestTHN";
// TString SinputFileName = "LHC23_PbPb_pass4_medium_Train346163"; //THN only
// TString SinputFileName = "LHC23_PbPb_pass4_Train347929"; // THN only, BDT score > 0.4, Pzs2 from Lambda
// TString SinputFileName = "LHC23_PbPb_pass4_Train354079"; // THN only, BDT score > 0.2, Pzs2 from Lambda
// TString SinputFileName = "LHC23_PbPb_pass4_Train361757"; // THN only, for Pz,s2 of Xi (direct measurement), acceptance applied on the fly
// TString SinputFileName = "LHC23_PbPb_pass4_Train365784"; // THN only, for Pz,s2 of Xi (direct measurement), acceptance not applied
// TString SinputFileName = "LHC23_PbPb_pass4_Train366446"; // proton acceptance calculation (vs pt Lambda)
// TString SinputFileName = "LHC23_PbPb_pass4_Train368064_ProtonAcc"; // proton acceptance calculation (vs pt and eta of Lambda)
// TString SinputFileName = "LHC23_PbPb_pass4_Train369742"; // Pzs2 from Lambda, proton acceptance vs pt on the fly
// TString SinputFileName = "LHC23_PbPb_pass4_Train370610_ProtonAcc"; // PRELIMINARIES: Pzs2 from Lambda, proton acceptance vs pt on the fly, proton acceptance vs pt and eta of Lambda
// TString SinputFileName = "LHC23_PbPb_pass5_Train456578_ProtonAcc"; // Xi polarization, proton acceptance vs pt and eta of Lambda
// TString SinputFileName = "LHC23_PbPb_pass5_Train456579_ProtAccFromPass4"; // Pzs2 of Xi from Lambda, proton acceptance vs pt and eta of Lambda from PASS4
// TString SinputFileName = "LHC23_PbPb_pass5_Train463978_PrimaryProtonAcceptance"; //proton acceptance for primary lambdas
// TString SinputFileName = "LHC23_PbPb_pass5_Train463979_ProtAcceptanceFromSecondayLambdas"; //Lambda polarization, proton acceptance for secondary lambdas
// TString SinputFileName = "464640_NewEP";
// TString SinputFileName = "TestLFEP";
// TString SinputFileName = "LHC23_PbPb_pass5_small_testEP"; // test with LF EP
// TString SinputFileName = "LHC23_PbPb_pass5_Train481586"; // test with Lambdas
// TString SinputFileName = "LHC25_OO_Train487953"; // test with Lambdas
// TString SinputFileName = "LHC25_OO_LambdaPol_Train491711"; // Pzs2 of Lambda
// TString SinputFileName = "LHC25_OO_pass2_Train503805"; // Pzs2 of Lambda

// TString SinputFileName = "LHC25_OO_pass2_Train510678_CorrectReso"; //Pzs2 of Lambda up to 100%
// TString SinputFileName = "LHC25_OO_pass2_Train510678"; //Pzs2 of Lambda up to 100%
// TString SinputFileName = "LHC25_OO_pass2_Train562132_wTHN";
// TString SinputFileName = "LHC25_OO_pass2_Train562850"; //latest
// TString SinputFileName = "LHC25_OO_pass2_Train562132_wTHN";
// TString SinputFileName = "LHC25_OO_pass2_SecondaryProtonAcc_Train508938"; //used for acceptance
// TString SinputFileName = "LHC25_OO_pass2_Train576495"; //tree with |eta| < 0.8 -- the DEFAULT ONE FOR THE MOMENT
// TString SinputFileName = "LHC25_OO_pass2_Train575744"; //THN
// TString SinputFileName = "LHC25_OO_pass2_Train576496"; //tree with |eta| < 0.8 and |z| < 8 cm
// TString SinputFileName = "LHC25_OO_pass2_Train589559"; //THN larger range pzs2
//TString SinputFileName = "LHC25_OO_pass2_Train589711"; // tree with new acceptance
//TString SinputFileName = "LHC25_OO_pass2_Train589711Bis"; // tree with new acceptance
//TString SinputFileName = "LHC25_OO_pass2_Train597528_NewAcc"; // THN with new acceptance (wrt previou: |etaDau| < 0.8)
//TString SinputFileName = "LHC25_OO_pass2_Train597527_NewAcc_EtaPos"; // THN with new acceptance (wrt previou: |etaDau| < 0.8)
//TString SinputFileName = "LHC25_OO_pass2_Train597526_NewAcc_EtaNeg"; // THN with new acceptance (wrt previou: |etaDau| < 0.8)
TString SinputFileName = "LHC25_OO_pass2_Train598890"; 
//TString SinputFileName = "LHC25_OO_pass2_Train598891_EtaPos"; 
//TString SinputFileName = "LHC25_OO_pass2_Train598892_EtaNeg"; 

// TString SinputFileNameAR = "LHC25_OO_pass2_Train510678";
// TString SinputFileNameAR = "LHC25_OO_pass2_Train562132_wTHN";
// TString SinputFileNameAR = "LHC25_OO_pass2_Train562850";
//TString SinputFileNameAR = "LHC25_OO_pass2_Train589711";
TString SinputFileNameAR = "LHC25_OO_pass2_Train598890";
//TString SinputFileNameAR = "LHC25_OO_pass2_Train598891_EtaPos";
//TString SinputFileNameAR = "LHC25_OO_pass2_Train598892_EtaNeg";
// TString SinputFileNameAR = "LHC25_OO_pass2_Train567017"; //latest resolution

// TString SinputFileNameCentWeight = "LHC25_OO_pass2_Train503805";
TString SinputFileNameCentWeight = "LHC25_OO_pass2_Train562132_wTHN";

// TString SinputFileNameResoWeight = "Resolution_SP_CFW_LHC25_OO_pass2_Train510916.root";
// TString SinputFileNameResoWeight = "Resolution_EP_CFW_LHC25_OO_pass2_Train557787_T0CShiftCorr_TPCCorr_WithT0A.root";
TString SinputFileNameResoWeight = "Resolution_EP_CFW_LHC25_OO_pass2_Train567017.root"; // the most recent ones, compatible with Train557787_T0CShiftCorr_TPCCorr_WithT0A
// TString SinputFileNameResoWeight = "Resolution_EP_CFW_LHC25_OO_pass2_Train562132_wTHN.root";

// File names for systematics
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train333596";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_medium_Train346163";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train347929";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train354079";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train361757";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train365784";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train366446";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train368064_ProtonAcc";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train369742";
// TString SinputFileNameSyst = "LHC23_PbPb_pass4_Train370610_ProtonAcc";
// TString SinputFileNameSyst = "LHC23_PbPb_pass5_Train456578_ProtonAcc";
// TString SinputFileNameSyst = "LHC23_PbPb_pass5_Train456579_ProtAccFromPass4";
// TString SinputFileNameSyst = "LHC23_PbPb_pass5_Train481586"; // test with Lambdas
// TString SinputFileNameSyst = "LHC25_OO_Train487953";
// TString SinputFileNameSyst = "LHC23_PbPb_pass5_Train456579_ProtAccFromPass4";
// TString SinputFileNameSyst = "LHC25_OO_LambdaPol_Train491711";
// TString SinputFileNameSyst = "LHC25_OO_pass2_Train503805";
// TString SinputFileNameSyst = "LHC25_OO_pass2_Train510678_CorrectReso";
// TString SinputFileNameSyst = "LHC23_PbPb_pass5_Train463979_ProtAcceptanceFromSecondayLambdas";
//TString SinputFileNameSyst = "LHC25_OO_pass2_Train562850";
//TString SinputFileNameSyst = "LHC25_OO_pass2_Train589711";
TString SinputFileNameSyst = "LHC25_OO_pass2_Train589711";

// File name for efficiency correction (if ExtrisApplyEffWeights == 1)
TString SinputFileNameEff = "LHC24g3_pass4_Train331315";
TString SinputFileNameEffSyst = "LHC24g3_pass4_Train331315";

// MC file for Lambda feed-down fraction
TString SinputFileNameFDFraction = "LHC25h3b_pass2_Train591313";

// BDT selections---------------------------------------------------------
const float DefaultBDTscoreCut = 0.96;
const bool useCommonBDTValue = 1; // common BDT cut for all centralities, set to DefaultBDTscoreCut
const float bdtCut[numCent + 1] = {0.95, 0.95, 0.95, 0.85, 0.85, 0.85, 0.85, 0.85, 0.95};

//---------------------------------------------------------
// BDT scores applied to produce the acceptance correction (chosen in order to have large purity)
bool isProducedAcceptancePlots = 0; // 1 for acceptance production, 0 for default analysis
const float BDTscoreCutAcceptance[numCent + 1] = {0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96};
//---------------------------------------------------------

// Variabls used in FitV2OrPol.C macro ----------------------
const bool isApplyAcceptanceCorrection = 0;                        // for recent files, acceptance correction is applied on the fly
const bool isAcceptanceFromExternalFile = 0;                       // 1 for acceptance from external file, 0 for acceptance from the same file
TString SAcceptanceFile = "../AcceptancePlots/Acceptance_Xi.root"; // file where acceptance is taken from if isAcceptanceFromExternalFile == 1
const bool useMixedBDTValueInFitMacro = 0;                         // variable used in FitV2OrPol.C macro
// if = 1: pt and multiplicity dependent value defined in:
//   - the function DefineMixedBDTValue (for the pt differential measurement) or
//   - BDTscoreCutPtInt (for the integrated pt measurement)
//   - BDTscoreCutPtIntLoosest (for the integrated pt measurement) if isTightMassCut=1
// if = 0: common BDT value defined by DefaultBDTscoreCut
const double BDTscoreCutPtInt[numCent + 1] = {0.8, 0.8, 0.6, 0.52, 0.44, 0.2, 0.2, 0.2, 0.56}; // BDT cut for integrated pt measurement
const float LimitForV2woFit = 0.97;                                                            // purity limit to decide whether to extract v2 from fit or not
bool isTightMassCut = 1;                                                                       // 1 for tight mass cut, 0 for loose mass cut
float Extrsigmacentral[2] = {4.2, 2.1};                                                        // 2.1
const double BDTscoreCutPtIntLoosest[numCent + 1] = {0.96, 0.92, 0.88, 0.76, 0.52, 0.4, 0.24, 0.2, 0.92};
// BDT cut for integrated pt measurement, loosest cut that give a purity > 0.95 within Extrsigmacentral[1];

// --------------------------- SYST ------------------------------
const int trialsLambdaTopo = 396; // number of trials for the systematic studies related to Lambda topology
// systematic studies on BDT score variation ----------------------
bool ExtrisSysMultTrial = 0;       // 1 for systematic studies, 0 for default analysis
bool ExtrisSysLambdaMultTrial = 0; // 1 for systematic studies, 0 for default analysis
const int trialsBDT = 201;         // number of trials for the systematic studies related to BDTscore
const float nsigmaBarlow = 0;
const float UpperlimitBDTscoreCut = 1;
const float LowerlimitBDTscoreCut = 0.2;
// const float MinBDTscorePtInt[numCent + 1] = {0.4, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.4}; // minimum BDT value for syst. evaluation
// const float MaxBDTscorePtInt[numCent + 1] = {0.96, 0.96, 0.8, 0.8, 0.8, 0.6, 0.48, 0.48, 0.96}; // maximum BDT value for syst. evaluation
const double MinBDTscorePtInt[numCent + 1] = {0.959, 0.92, 0.879, 0.76, 0.52, 0.4, 0.24, 0.2, 0.92};
const double MaxBDTscorePtInt[numCent + 1] = {0.98, 0.96, 0.96, 0.96, 0.96, 0.96, 0.8, 0.76, 0.96};
const bool isLoosest = 0;
const bool isTightest = 0;

// systematics for Lambda
const float DefaultV0RadiusCut = 1.0;
const float UpperlimitV0RadiusCut = 1.4;
const float LowerlimitV0RadiusCut = 0.9; // derived data limit
const float DefaultDcaV0DauCut = 1.2;
const float UpperlimitDcaV0DauCut = 1.5; // derived data limit
const float LowerlimitDcaV0DauCut = 0.26;
const double DefaultV0CosPA = 0.995;
const double UpperlimitV0CosPA = 0.998;
const double LowerlimitV0CosPA = 0.99; // 0.97 is the derived data limit
const float DefaultDcaNegToPV = 0.06;
const float UpperlimitDcaNegToPV = 0.09;
const float LowerlimitDcaNegToPV = 0.05; // derived data limit
const float DefaultDcaPosToPV = 0.06;
const float UpperlimitDcaPosToPV = 0.09;
const float LowerlimitDcaPosToPV = 0.05; // derived data limit

// Systematic studies on mass cut
const int trialsMassCut = 7;
bool ExtrisSysMassCut = 0; // 1 for systematic studies, 0 for default analysis
float ExtrLowLimitSysXi[trialsMassCut] = {1.316, 1.317, 1.318, 1.316, 1.316, 1.317, 1.318};
float ExtrUpLimitSysXi[trialsMassCut] = {1.327, 1.326, 1.325, 1.326, 1.325, 1.327, 1.327};

// Systematic studies on eta selection
TString SEtaSysChoice[3] = {"", "_Etagt0", "_Etasm0"}; // all eta, eta > 0, eta < 0
Int_t ExtrEtaSysChoice = 0;                            // 0: all eta, 1: eta > 0, 2: eta < 0

// Systematic studies on IR
TString SIRChoice[6] = {"", "_544013", "_544392", "_544098", "_544032", "_544184"};
TString SIRValue[6] = {"", "6 kHz", "12 kHz", "18 kHz", "23 kHz", "33 kHz"};
TString inputFileNameIR = "Train207098";

// -------- Event plane resolution ------------------------------
// float ftcResoSourav[numCent] = {0.514595, 0.7228, 0.760156, 0.733402, 0.659964, 0.540407, 0.383689, 0.218501}; // values from Sourav
// float ftcResoSP[numCent] = {0.805015, 1.06586, 1.10463, 1.07285, 0.984069, 0.82759, 0.602314, 0.34722};

// Files to compute resolution
// TString ComputeResoFileNameLF = "TreeForAnalysis/AnalysisResults_LHC23zzh_pass3_Train224930.root";      //<-- LF framework
// TString ComputeResoFileNameLF = "TreeForAnalysis/AnalysisResults_LHC23zzh_pass4_test5_Train235645.root";      //<-- LF framework
TString inputFileResoCFW = SinputFileNameAR; // OLD: "LHC23zzh_pass3_Train226234_CFW";
TString inputFileResoLF = SinputFileNameAR;

// Files with stored resolution
TString ResoFileName_EPLF = "Resolution/Resolution_EP_LF_" + inputFileResoLF;
TString ResoFileName_EPCFW = "Resolution/Resolution_EP_CFW_" + inputFileResoCFW;
TString ResoFileName_SPLF = "Resolution/Resolution_SP_LF_" + inputFileResoLF;
TString ResoFileName_SPCFW = "Resolution/Resolution_SP_CFW_" + inputFileResoCFW;

// theory predictions
// A.Palermo, Pzs2 of Lambda vs centrality with bulk viscosity
Double_t CentPalermo[9] = {2.5e+00, 7.5e+00, 1.5e+01, 2.5e+01, 3.5e+01, 4.5e+01, 5.5e+01, 6.5e+01, 7.5e+01};
Double_t Pzs2Palermo[9] = {-3.050173509930762550e-05, -7.858818176538021343e-05, -8.719012066015866002e-05, 7.493273322344973971e-06,
                           2.665266717715395481e-04, 7.819382729092881927e-04, 1.527638672979295988e-03, 2.227462747363428888e-03, 2.583982606087120888e-03};

// Published V2 of charged particles in OO collisions (arxiv.org/pdf/2509.06428)
const Int_t numV2OOPubCent = 16;
Double_t V2OOPubCent[numV2OOPubCent + 1] = {0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
Double_t V2OOPubCentMid[numV2OOPubCent] = {130.785, 127.023, 123.37, 119.822, 116.376, 106.621, 92.1451, 79.6348,
                                           68.8231, 59.4792, 51.4039, 44.425, 38.3936, 33.181, 28.6761, 24.7829};
Double_t V2OOPubCentMidErr[numV2OOPubCent] = {3.47982, 3.29583, 3.12024, 2.95276, 2.79311, 2.35845, 1.76754, 1.32212,
                                              1.0029, 0.793745, 0.675824, 0.622972, 0.606406, 0.60396, 0.603009, 0.597924};
Double_t V2OOPub[numV2OOPubCent] = {0.058265, 0.059672, 0.060794, 0.062044, 0.062179, 0.063781, 0.065495, 0.066987,
                                    0.067843, 0.068207, 0.068403, 0.068001, 0.067598, 0.066635, 0.066261, 0.064535};
Double_t V2OOPubErrStat[numV2OOPubCent] = {0.000070, 0.000162, 0.000165, 0.000134, 0.000126, 0.000072, 0.000070, 0.000069,
                                           0.000090, 0.000095, 0.000123, 0.000084, 0.000165, 0.000186, 0.000198, 0.000196};
Double_t V2OOPubErrSyst[numV2OOPubCent] = {0.000226, 0.000232, 0.000236, 0.000241, 0.000241, 0.000248, 0.000254, 0.000260,
                                          0.000263, 0.000265, 0.000266, 0.000264, 0.000263, 0.000259, 0.000257, 0.000251};

//Published V2 of charged particles in PbPb collisions (arxiv.org/pdf/1602.01119)
const Int_t numV2PbPbPubCent = 8;
Double_t V2PbPbPubCent[numV2PbPbPubCent + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
Double_t V2PbPbPub[numV2PbPbPubCent] = {(0.02839 + 0.04566)/2, 0.06551, 0.08707, 0.0991, 0.10414, 0.10286, 0.09746, 0.08881};
Double_t V2PbPbPubErrStat[numV2PbPbPubCent] = {(0.00057 + 0.00064)/2, 0.00037, 0.00044, 0.00055, 0.00073, 0.00107, 0.00186, 0.00438};
Double_t V2PbPbPubErrSys[numV2PbPbPubCent] = {(0.00043 + 0.00069)/2, 0.00098, 0.00131, 0.00149, 0.00156, 0.00154, 0.00146, 0.00133};