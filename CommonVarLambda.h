//To be changed according to the following instructions to produce acceptance
//isApplyCentWeight = 0
//ExtrisApplyResoOnTheFly = 0
//ExtrisFromTHN = 1
//isProducedAcceptancePlots = 1
//SinputFileName --> take the proper input for acceptance calculation

//To be changed according to the following instructions to produce systematic variations in input of MultiTrial.C
//ExtrisSysLambdaMultTrial = 1
//trialsLambdaTopo --> actual number of variations
//SinputFileName --> take the proper input for systematic variations

Bool_t isV2 = 0;              // 0 for polarization, 1 for v2
Int_t ChosenParticle = 6;     // 0: Xi, 1: Omega, 2: Xi-, 3: Xi+, 4: Omega-, 5: Omega+, 6: Lambda + ALambda
Bool_t ExtrisRapiditySel = 0; // 0: |eta| < 0.8, 1: |y| < 0.5 (for Pzs2),
Bool_t ExtrisPartialEta = 0;  // 1: select only 0 < eta < 0.8 (opposite to FT0C)
Int_t ExtrBkgType = 1;        // 0: pol1, 1:pol2, 2:pol3, 3:expo
Int_t ExtrBkgTypeSyst = 1;    // for syst. uncertainty: 0: pol1, 1:pol2, 2:pol3, 3:expo
Bool_t ExtrUseTwoGauss = 1;
Bool_t isApplyWeights = 0;          // weights to flatten the phi distribution of cascades
Bool_t isApplyCentWeight = 1;       // 0 for acceptance
Bool_t ExtrisApplyEffWeights = 1;   
Bool_t ExtrisApplyResoOnTheFly = 1; // 0 for acceptance
Int_t v2type = 2;                   // 0: v2 - old task version before train 224930, 1: v2 SP, 2: v2 EP
Bool_t ExtrisFromTHN = 0;           // 1 for acceptance; 0: process the tree, 1: process the THnSparse
Bool_t isReducedPtBins = 1;         //0 for acceptance
Bool_t isOOCentrality = 1;
Bool_t isRun2Binning = 0;

const Int_t commonNumCent = 10; //= numCentLambdaOO for Lambda in OO

// Pt bins
const Int_t numPtBinsEff = 17; // for efficiency
const Int_t numPsiBins = 6;    // bins into which Pz (longitudinal polarization) is computed
const Int_t numPtBins = 7;
const Int_t numPtBinsReduced = 7;
Double_t PtBinsEff[numPtBinsEff + 1] = {0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8};
Double_t PtBins[numPtBins + 1] = {0.5, 1.0, 1.5, 2, 2.5, 3, 4, 8};
Float_t MinPt[numPart] = {0.8, 1., 0.8, 0.8, 1., 1., 0.5, 0.5, 0.5};
Float_t MaxPt[numPart] = {8, 8, 8, 8, 8, 8, 8, 8, 8};

// Acceptance correction
const Int_t numEtaBins = 16; // was 8
const Int_t numPtBinsLambda = 9;
Double_t EtaBins[numEtaBins + 1] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
Double_t PtBinsLambda[numPtBinsLambda + 1] = {0.4, 0.8, 1.2, 1.6, 2, 2.5, 3, 4, 6, 10};

// File names
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
// TString SinputFileName = "LHC25_OO_pass2_Train576495"; //tree with |eta| < 0.8 
// TString SinputFileName = "LHC25_OO_pass2_Train575744"; //THN
// TString SinputFileName = "LHC25_OO_pass2_Train576496"; //tree with |eta| < 0.8 and |z| < 8 cm
// TString SinputFileName = "LHC25_OO_pass2_Train589559"; //THN larger range pzs2
// TString SinputFileName = "LHC25_OO_pass2_Train589711"; // tree with new acceptance
// TString SinputFileName = "LHC25_OO_pass2_Train589711Bis"; // tree with new acceptance
// TString SinputFileName = "LHC25_OO_pass2_Train597528_NewAcc"; // THN with new acceptance (wrt previou: |etaDau| < 0.8)
// TString SinputFileName = "LHC25_OO_pass2_Train597527_NewAcc_EtaPos"; // THN with new acceptance (wrt previou: |etaDau| < 0.8)
// TString SinputFileName = "LHC25_OO_pass2_Train597526_NewAcc_EtaNeg"; // THN with new acceptance (wrt previou: |etaDau| < 0.8)
//TString SinputFileName = "LHC25_OO_pass2_Train598890"; 
TString SinputFileName = "LHC25_OO_pass2_Train598890_MyEff"; //for PRELIMINARIES 2026
// TString SinputFileName = "LHC25_OO_pass2_Train598890_PositiveEta"; //no sel on daughter tracks eta apart from |etaDau| < 0.8 
// TString SinputFileName = "LHC25_OO_pass2_Train598890_NegativeEta"; //no sel on daughter tracks eta apart from |etaDau| < 0.8
// TString SinputFileName = "LHC25_OO_pass2_Train598891_EtaPos"; //Also 0 < etaDau < 0.8
// TString SinputFileName = "LHC25_OO_pass2_Train598892_EtaNeg"; //Also -0.8 < etaDau < 0

//To get number of analyzed events
TString SinputFileNameAR = "LHC25_OO_pass2_Train598890";

// File name for centrality weights
// TString SinputFileNameCentWeight = "LHC25_OO_pass2_Train503805";
TString SinputFileNameCentWeight = "LHC25_OO_pass2_Train562132_wTHN";

//File name for efficiency weights
//TString SinputFileNameEfficiency = "CorrectedSpectra_Lambda_withEvtLoss_withFeeddown.root"; //Romain Schotter input
TString SinputFileNameEfficiency = "LHC25h3c_pass2_Train621696"; //--> in input of ComputeEff.C macro, to be used for the efficiency correction in FitV2OrPol.C macro
//TString SinputFileNameEfficiencyWeight = "../EfficiencyWeight.root"; //Romain efficiency
TString SinputFileNameEfficiencyWeightLambda = "EfficiencyWeight_LHC25h3c_pass2_Train621696_Lambda_Eta08.root";
TString SinputFileNameEfficiencyWeightAntiLambda = "EfficiencyWeight_LHC25h3c_pass2_Train621696_AntiLambda_Eta08.root";

// File name for resolution weights
// TString SinputFileNameResoWeight = "Resolution_SP_CFW_LHC25_OO_pass2_Train510916.root";
// TString SinputFileNameResoWeight = "Resolution_EP_CFW_LHC25_OO_pass2_Train557787_T0CShiftCorr_TPCCorr_WithT0A.root";
TString SinputFileNameResoWeight = "Resolution_EP_CFW_LHC25_OO_pass2_Train567017.root"; // the most recent ones, compatible with Train557787_T0CShiftCorr_TPCCorr_WithT0A
// TString SinputFileNameResoWeight = "Resolution_EP_CFW_LHC25_OO_pass2_Train562132_wTHN.root";

// File names for systematics (taken in input of MultiTrial.C, SystematicErrorVsCent.C, and PzsVsCentrality.C to plot final results)
// TString SinputFileNameSyst = "LHC25_OO_pass2_Train562850";
// TString SinputFileNameSyst = "LHC25_OO_pass2_Train589711";
TString SinputFileNameSyst = "LHC25_OO_pass2_Train589711";

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

// -------- Event plane resolution (not for Pzs of Lambda) ------------------------------
TString inputFileResoCFW = SinputFileNameAR; 
TString inputFileResoLF = SinputFileNameAR;

// Files with stored resolution
TString ResoFileName_EPLF = "Resolution/Resolution_EP_LF_" + inputFileResoLF;
TString ResoFileName_EPCFW = "Resolution/Resolution_EP_CFW_" + inputFileResoCFW;
TString ResoFileName_SPLF = "Resolution/Resolution_SP_LF_" + inputFileResoLF;
TString ResoFileName_SPCFW = "Resolution/Resolution_SP_CFW_" + inputFileResoCFW;