Bool_t isV2 = 0;              // 0 for polarization, 1 for v2
Int_t ChosenParticle = 0;     // 0: Xi, 1: Omega, 2: Xi-, 3: Xi+, 4: Omega-, 5: Omega+, 6: Lambda + ALambda
Bool_t ExtrisRapiditySel = 0; // 0: |eta| < 0.8, 1: |y| < 0.5 (for Pzs2)
Int_t ExtrBkgType = 1;       // 0: pol1, 1:pol2, 2:pol3, 3:expo
Int_t ExtrBkgTypeSyst = 1; // for syst. uncertainty: 0: pol1, 1:pol2, 2:pol3, 3:expo
Bool_t ExtrUseTwoGauss = 1;
Bool_t isApplyWeights = 0;          // weights to flatten the phi distribution of cascades
Bool_t isApplyCentWeight = 0;       // 1 for OO
Bool_t ExtrisApplyEffWeights = 0;   // weights to take into account efficiency dependence on multiplciity (for v2 only)
Bool_t ExtrisApplyResoOnTheFly = 0; // 1 for OO
Int_t v2type = 2;                   // 0: v2 - old task version before train 224930, 1: v2 SP, 2: v2 EP
Bool_t ExtrisFromTHN = 1;           // 0: process the tree, 1: process the THnSparse
Bool_t isReducedPtBins = 0;         // 1 for Lambda in OO
Bool_t isOOCentrality = 0;          // 1 for Lambda in OO
Bool_t isRun2Binning = 0;

const Int_t commonNumCent = 8; // = numCent for Xi

// Pt bins
// const Int_t numPtBins = 15;
const Int_t numPtBins = 7;
const Int_t numPtBinsReduced = 7;
//  const Int_t numPtBins = 6; // Run2 binning
const Int_t numPtBinsEff = 15; // for efficiency
const Int_t numPsiBins = 6;    // bins into which Pz (longitudinal polarization) is computed
// Double_t PtBins[numPtBins + 1] = {0.8, 1.4, 2, 2.5, 3, 4, 6}; // Run 2 binning for v2
// Double_t PtBins[numPtBins + 1] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8};
Double_t PtBinsEff[numPtBinsEff + 1] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8};
Double_t PtBins[numPtBins + 1] = {0.8, 1.0, 1.5, 2, 2.5, 3, 4, 8};
Float_t MinPt[numPart] = {0.8, 1., 0.8, 0.8, 1., 1., 0.5, 0.5, 0.5};
Float_t MaxPt[numPart] = {8, 8, 8, 8, 8, 8, 8, 8, 8};

// Acceptance correction
const Int_t numEtaBins = 8;
const Int_t numPtBinsLambda = 9;
Double_t EtaBins[numEtaBins + 1] = {-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8};
Double_t PtBinsLambda[numPtBinsLambda + 1] = {0.4, 0.8, 1.2, 1.6, 2, 2.5, 3, 4, 6, 10};

// File names
// TString SinputFileName = "LHC23_PbPb_pass4_Train370610_ProtonAcc"; // PRELIMINARIES: Pzs2 from Lambda, proton acceptance vs pt on the fly, proton acceptance vs pt and eta of Lambda
// TString SinputFileName = "LHC23_PbPb_pass5_Train456578_ProtonAcc"; // proton acceptance vs pt and eta of Lambda for Xi polarization
// TString SinputFileName = "LHC23_PbPb_pass5_Train456579_ProtAccFromPass4"; // Pzs2 of Xi from Lambda, proton acceptance vs pt and eta of Lambda from PASS4
// TString SinputFileName = "LHC23_PbPb_pass5_Train534683"; // Pzs2 of Xi from Lambda, proton acceptance vs pt and eta of Lambda from PASS5
TString SinputFileName = "LHC23_PbPb_pass5_Train540301"; // PAPER PROPOSAL: Pzs2 of Xi from Lambda, proton acceptance vs pt and eta of Lambda from PASS5, event plane FLAT in phi (shift corrected)
// TString SinputFileName = "LHC23_PbPb_pass5_Train541065"; // Pzs2 of Xi from Lambda, proton acceptance vs pt and eta of Lambda from PASS5, event plane FLAT in phi (shift corrected), |z|< 8 cm
// TString SinputFileName = "LHC23_PbPb_pass5_Train563856"; // same settings as Train540301 but histograms for EP resolutions stored (to be used for resolution!)
// TString SinputFileName = "LHC23_PbPb_pass5_Train566502"; // acceptance for run by run studies
//TString SinputFileName = "LHC23_PbPb_pass5_Train567157_OccupancyCut"; // systematic associated with occupancy; FT0COccupancy < 3000 (cuts a lot of events)
//TString SinputFileName = "LHC23_PbPb_pass5_Train568468_OccupancySel20000";
//TString SinputFileName = "LHC23_PbPb_pass5_Train568467_OccupancySel30000";

// Reso tests in OO
// TString SinputFileName = "LHC25_OO_pass2_Train510678"; // Pzs2 of Lambda
// TString SinputFileName = "LHC25_OO_pass2_Train543247_T0CEventPlaneTest"; // Pzs2 of Lambda
// TString SinputFileName = "LHC25_OO_pass2_Train544235_V0AShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train545130_V0ANShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train546065_V0ANoShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train547804_V0AShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train548457_T0AShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train548636_T0CShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train547803_T0ANoShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train548643_T0MNoShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train549121_T0CNoShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train549122_T0MNoShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train549964_T0CShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train549965_T0MShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train556005_T0CShiftCorr_TPCCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train557895_T0AShiftCorr"; //T0A shift corrected (validated)
// TString SinputFileName = "LHC25_OO_pass2_Train557896_V0AShiftCorr"; //V0A shift corrected (validated)
// TString SinputFileName = "LHC25_OO_pass2_Train557787_T0CShiftCorr_TPCCorr_WithT0A"; //T0C shift corrected with corrected TPC event plane and T0A info for resolution
// TString SinputFileName = "LHC25_OO_pass2_Train558333_T0MShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train562133_T0CShiftCorr";
// TString SinputFileName = "LHC25_OO_pass2_Train518384_V0AResolution"; // V0AResolution
// TString SinputFileName = "LHC25_OO_pass2_Train518383_T0MResolution"; // V0AResolution
// TString SinputFileName = "LHC25_OO_pass2_Train523874";
// TString SinputFileName = "LHC25_OO_pass2_Train515730_V0AResolution";
// TString SinputFileName = "LHC25_OO_pass2_Train515731_T0MResolution";
// TString SinputFileName = "LHC25_OO_pass2_Train510916"; //reso in 1% centrality bins
// TString SinputFileName = "LHC25_OO_pass2_SecondaryProtonAcc_Train508938"; // secondary proton acceptance for Lambda pol in OO (possibly used to evaluate systematics)
// TString SinputFileName = "LHC25_OO_pass2_Train510677"; // resolution in 1% centrality bins

TString SinputFileNameReso = "LHC23_PbPb_pass5_Train563856"; //used for PAPER PROPOSAL
//TString SinputFileNameReso = "LHC23_PbPb_pass5_Train567157_OccupancyCut";
TString SinputFileNameAR = SinputFileName;
TString SinputFileNameResoWeight = ""; // empty, not needed for Xi in Pb-Pb

// File names for systematics
TString SinputFileNameSyst = "LHC23_PbPb_pass5_Train534683"; //these systematics are DONE and to be USED for PAPER
//TString SinputFileNameSyst = "LHC23_PbPb_pass5_Train568467_OccupancySel30000";
//TString SinputFileNameSyst = "LHC23_PbPb_pass5_Train567157_OccupancyCut";
// TString SinputFileNameSyst = "LHC23_PbPb_pass5_Train540301"; //these were not run YET
// TString SinputFileNameSyst = "LHC23_PbPb_pass5_Train541065";

// File name for efficiency correction (if ExtrisApplyEffWeights == 1)
TString SinputFileNameEff = "LHC24g3_pass4_Train331315";
TString SinputFileNameEffSyst = "LHC24g3_pass4_Train331315";

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
const bool isApplyAcceptanceCorrection = 0;                     // for recent files, acceptance correction is applied on the fly
const bool isAcceptanceFromExternalFile = 0;                    // 1 for acceptance from external file, 0 for acceptance from the same file
TString SAcceptanceFile = "AcceptancePlots/Acceptance_Xi.root"; // file where acceptance is taken from if isAcceptanceFromExternalFile == 1
const bool useMixedBDTValueInFitMacro = 1;                      // variable used in FitV2OrPol.C macro
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
const int trialsLambdaTopo = 2;
// systematic studies on BDT score variation ----------------------
bool ExtrisSysMultTrial = 0;       // 1 for systematic studies, 0 for default analysis
bool ExtrisSysLambdaMultTrial = 0; // 1 for systematic studies, 0 for default analysis
const int trialsBDT = 20;          // number of trials for the systematic studies related to BDTscore
const float nsigmaBarlow = 0;
const float UpperlimitBDTscoreCut = 1;
const float LowerlimitBDTscoreCut = 0.2;
// const double MinBDTscorePtInt[numCent + 1] = {0.959, 0.92, 0.879, 0.76, 0.52, 0.4, 0.24, 0.2, 0.92};
// const double MaxBDTscorePtInt[numCent + 1] = {0.98, 0.96, 0.96, 0.96, 0.96, 0.96, 0.8, 0.76, 0.96};
const double MinBDTscorePtInt[numCent + 1] = {0.84, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2, 0.2, 0.6};
const double MaxBDTscorePtInt[numCent + 1] = {0.98, 0.96, 0.96, 0.96, 0.96, 0.96, 0.8, 0.76, 0.96};
// const double MaxBDTscorePtInt[numCent + 1] = {0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72};
const bool isLoosest = 0;
const bool isTightest = 0;

// systematics for Lambda
const float DefaultV0RadiusCut = 1.0;
const float UpperlimitV0RadiusCut = 1.2;
const float LowerlimitV0RadiusCut = 0.9; // derived data limit
const float DefaultDcaV0DauCut = 1.2;
const float UpperlimitDcaV0DauCut = 1.5; // derived data limit
const float LowerlimitDcaV0DauCut = 1.0;
const double DefaultV0CosPA = 0.995;
const double UpperlimitV0CosPA = 0.999;
const double LowerlimitV0CosPA = 0.99; // 0.97 is the derived data limit
const float DefaultDcaNegToPV = 0.06;
const float UpperlimitDcaNegToPV = 0.1;
const float LowerlimitDcaNegToPV = 0.05; // derived data limit
const float DefaultDcaPosToPV = 0.06;
const float UpperlimitDcaPosToPV = 0.1;
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
TString inputFileResoCFW = SinputFileNameReso; // OLD: "LHC23zzh_pass3_Train226234_CFW";
TString inputFileResoLF = SinputFileNameReso;

// Files with stored resolution
TString ResoFileName_EPLF = "../Resolution/Resolution_EP_LF_" + inputFileResoLF;
TString ResoFileName_EPCFW = "../Resolution/Resolution_EP_CFW_" + inputFileResoCFW;
TString ResoFileName_SPLF = "../Resolution/Resolution_SP_LF_" + inputFileResoLF;
TString ResoFileName_SPCFW = "../Resolution/Resolution_SP_CFW_" + inputFileResoCFW;