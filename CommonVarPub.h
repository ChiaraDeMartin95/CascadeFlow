const Int_t numPart = 9;    // Xi+-, Omega+-, Xi-, Xi+, Omega-, Omega+, Lambda + ALambda, Lambda, AntiLambda
const Int_t numChoice = 12; // mean, sigma, purity, yield, v2, Pzs2, Pzs2 from lambda, Cos2Theta, Cos2Theta from lambda, V2MixedCorr, Cos2ThetaFromLambdaVsPtLambda

// Titles
TString sPolFromLambda[2] = {"", "LambdaFromC"};
TString STHN[2] = {"", "_FromTHN"};
TString V2FromFit[2] = {"NoFit", ""};
TString NameAnalysis[2] = {"V2", "Pzs2"};
TString RapidityCoverage[2] = {"Eta08", "Y05"};
TString IsOneOrTwoGauss[2] = {"_OneGaussFit", ""};
TString SIsBkgParab[5] = {"_BkgRetta", "_BkgParab", "_BkgPol3", "_BkgExpo", "_BkgCheb"};
Float_t ParticleMassPDG[numPart] = {1.32171, 1.67245, 1.32171, 1.32171, 1.67245, 1.67245, 1.115683, 1.115683, 1.115683}; // Xi+-, Omega+-, Xi-, Xi+, Omega-, Omega+, Lambda + ALambda
TString ParticleName[numPart] = {"Xi", "Omega", "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus", "Lambda", "LambdaPart", "AntiLambda"};
TString ParticleNameLegend[numPart] = {"#Xi^{#pm}", "#Omega^{#pm}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}", "#Lambda + #bar{#Lambda}", "#Lambda", "#overline{#Lambda}"};
TString TypeHisto[numChoice] = {"Mean", "SigmaWeighted", "Purity", "Yield", "V2Mixed", "Pzs2Mixed", "Pzs2LambdaFromCMixed", "Cos2ThetaNoFit", "Cos2ThetaLambdaFromC", "V2MixedCorr", "Cos2ThetaLambdaFromCVsPt", "Cos2ThetaLambdaFromCVsEta"};
TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "v2", "Pz,s2", "Pz,s2", "#LTcos^{2}(#theta*_{#Lambda})#GT", "#LTcos^{2}(#theta*_{p})#GT", "v2, corr", "#LTcos^{2}(#theta*_{p})#GT", "#LTcos^{2}(#theta*_{p})#GT"};
TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";
TString TitleXCent = "Centrality (%)";
TString TitleYPzs = "#it{P}_{z,s2}";
TString TitleYPzsLambda = "#it{P}_{z,s2} (#Lambda + #bar{#Lambda})";

// Mass ranges were puity is calculatedin case of DSCB fit
Float_t ExtrLowLimitDSCB[numPart] = {1.312, 1.672, 1.312, 1.312, 1.672, 1.672, 1.1145, 1.1145, 1.1145};
Float_t ExtrUpLimitDSCB[numPart] = {1.33, 1.672, 1.33, 1.33, 1.672, 1.672, 1.1145, 1.1145, 1.1145};
// 2sigma range approximately
// Float_t ExtrLowLimitDSCB[numPart] = {1.316, 1.672, 1.316, 1.316, 1.672, 1.672, 1.1145, 1.1145, 1.1145};
// Float_t ExtrUpLimitDSCB[numPart] = {1.326, 1.672, 1.326, 1.326, 1.672, 1.672, 1.1145, 1.1145, 1.1145};

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
Float_t AlphaHErrors[numPart] = {1, 1, 0.007, 0.007, 0.0020, 0.004, 1, 1, 1};
Float_t CXiToLambda = 0.925;
Float_t COmegaToLambda = 1;
// Float_t COmegaToLambda = -0.6;
Float_t AlphaLambda[numPart] = {1, 1, 0.746, -0.758, 0.746, -0.758, 1, 1, 1};     // decay parameter for Lambda -> p pi
Float_t AlphaLambdaErrors[numPart] = {1, 1, 0.008, 0.005, 0.008, 0.005, 1, 1, 1}; // decay parameter for Lambda -> p pi

// Centrality Pb-Pb (Lambda - Xi)
Bool_t ExtrisCentXiRed = 0; // 1 for reduced centrality bins for Xi in Pb-Pb
const Int_t numCent = 8;
Double_t CentFT0CMaxPbPb = 80;
Int_t CentFT0C[numCent + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
Double_t fCentFT0C[numCent + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
Double_t dNdEtaAbhi[numCent] = {(2047. + 1668.) / 2, 1253, 848, 559, 351, 205, 110, 53}; // values from https://arxiv.org/pdf/2504.02505; statistical uncertainites negligible
Double_t dNdEtaAbhiErr[numCent] = {(54. + 42.) / 2, 33, 25, 19, 14, 11, 8, 5};           // systematic uncertainties (https://arxiv.org/pdf/2504.02505)
Double_t EccPbPb[numCent] = {(0.074 + 0.110) / 2, 0.172, 0.246, 0.314, 0.377, 0.442, 0.518, 0.610};
Double_t v2PubRun2[numCent] = {(0.02839 + 0.04566) / 2, 0.06551, 0.08707, 0.0991, 0.10414, 0.10286, 0.09746, 0.08881}; // values from Run2 https://arxiv.org/pdf/1602.01119
const Int_t numCentXiRed = 3;
Double_t CentFT0CMaxXiRed = 70;
Int_t CentFT0CXiRed[numCentXiRed + 1] = {0, 30, 50, 70};
Double_t fCentFT0CXiRed[numCentXiRed + 1] = {0, 30, 50, 70};

// Centrality Pb-Pb (Omega)
Bool_t ExtrisCentOmegaRed = 0; // 1 for reduced centrality bins for Omega in Pb-Pb
const Int_t numCentOmega = 8;
Double_t CentFT0CMaxOmega = 80;
Int_t CentFT0COmega[numCentOmega + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
Double_t fCentFT0COmega[numCentOmega + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
const Int_t numCentOmegaRed = 3;
Double_t CentFT0CMaxOmegaRed = 70;
Int_t CentFT0COmegaRed[numCentOmegaRed + 1] = {0, 30, 50, 70};
Double_t fCentFT0COmegaRed[numCentOmegaRed + 1] = {0, 30, 50, 70};

// min and max pt
Float_t MinPt[numPart] = {0.8, 1.4, 0.8, 0.8, 1.4, 1.4, 0.5, 0.5, 0.5};
Float_t MaxPt[numPart] = {8, 8, 8, 8, 8, 8, 8, 8, 8};

// Centrality OO (Lambda)
const Int_t numCentLambdaOO = 10;
Int_t CentFT0CLambdaOO[numCentLambdaOO + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
Double_t fCentFT0CLambdaOO[numCentLambdaOO + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
Double_t CentFT0CMaxLambdaOO = 100;
Double_t dNdEtaOOPrel[numCentLambdaOO] = {(129.6660 + 106.8340) / 2, 87.2877, 67.1562, 51.1201, 37.8919, 26.9060, 0, 0, 0, 0};               // approved up to 60% (version 3 of Marta AN)
Double_t dNdEtaOOErrPrel[numCentLambdaOO] = {(0.0335 + 0.0250) / 2, 0.0143, 0.0117, 0.0099, 0.0084, 0.0069, 0, 0, 0, 0};                     // approved up to 60% (version 3 of Marta AN)
Double_t dNdEtaOOErrPrelSyst[numCentLambdaOO] = {(4.1009 + 3.3921) / 2, 2.8101, 2.2507, 1.8847, 1.6370, 1.4909, 0, 0, 0, 0};                 // approved up to 60% (version 3 of Marta AN)
Double_t dNdEtaOO[numCentLambdaOO] = {(129.6660 + 106.8340) / 2, 87.2877, 67.1562, 51.1201, 37.8919, 26.9060, 21.0754, 15.8103, 11.8605, 0}; // from 60% to 90%, extrapolated with MultVsCent.C macro
Double_t dNdEtaOOErr[numCentLambdaOO] = {(0.0335 + 0.0250) / 2, 0.0143, 0.0117, 0.0099, 0.0084, 0.0069, 0, 0, 0, 0};                         // from 60% to 90%, extrapolated with MultVsCent.C macro
Double_t dNdEtaOOErrSyst[numCentLambdaOO] = {(4.1009 + 3.3921) / 2, 2.8101, 2.2507, 1.8847, 1.6370, 1.4909, 0.868578, 0.793424, 0.70358, 0}; // from 60% to 90%, extrapolated with MultVsCent.C macro
Double_t EccOO[numCentLambdaOO] = {(0.302 + 0.323) / 2, 0.364, 0.421, 0.479, 0.532, 0.577, 0.614, 0.644, 0.672, 0.699};
Double_t dNdEtaOOPrel6090 = 11.6993;
Double_t dNdEtaOOErrPrel6090 = 0.0024;
Double_t dNdEtaOOErrPrelSyst6090 = 1.2348;

// NeNe dNdeta
Double_t dNdEtaNeNe[2] = {105.59, 20.63}; // for Junlee results. Averages computed from analysis note (0-40%, 40-90% even if polarization uses 40-100%; multiplicity available only up to 90%)
Double_t dNdEtaNeNeErr[2] = {3.52, 0.69}; // random reasonable errors assigned

// theory predictions
// A.Palermo, Pzs2 of Lambda vs centrality with bulk viscosity
Double_t CentPalermo[9] = {2.5e+00, 7.5e+00, 1.5e+01, 2.5e+01, 3.5e+01, 4.5e+01, 5.5e+01, 6.5e+01, 7.5e+01};
Double_t Pzs2Palermo[9] = {-3.050173509930762550e-05, -7.858818176538021343e-05, -8.719012066015866002e-05, 7.493273322344973971e-06,
                           2.665266717715395481e-04, 7.819382729092881927e-04, 1.527638672979295988e-03, 2.227462747363428888e-03, 2.583982606087120888e-03};

// Published V2 of charged particles in OO collisions (arxiv.org/pdf/2509.06428)
const Int_t numV2OOPubCent = 16;
Double_t V2OOPubCent[numV2OOPubCent + 1] = {0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
Double_t V2OOPubCentMid[numV2OOPubCent] = {134.579, 130.765, 127.06, 123.46, 119.961, 110.05, 95.3172, 82.5568, 71.5046, 61.9321, 53.641, 46.4599, 40.2402, 34.8531, 30.1872, 26.1459};
// these values are extracted from a fit to dNdeta (preliminary) vs centrality in the MultVsCent.C macro, and then used to compute the weighted average of v2 in each centrality bin
Double_t V2OOPubCentMidErr[numV2OOPubCent] = {3.68493, 3.4772, 3.27966, 3.09199, 2.91387, 2.43402, 1.80162, 1.35779, 1.08181, 0.944225, 0.899436, 0.899721, 0.912353, 0.920963, 0.919386, 0.906505};
// these values are extracted from a fit to dNdeta (preliminary) vs centrality in the MultVsCent.C macro, and then used to compute the weighted average of v2 in each centrality bin
Double_t V2OOPub[numV2OOPubCent] = {0.058265, 0.059672, 0.060794, 0.062044, 0.062179, 0.063781, 0.065495, 0.066987,
                                    0.067843, 0.068207, 0.068403, 0.068001, 0.067598, 0.066635, 0.066261, 0.064535};
Double_t V2OOPubErrStat[numV2OOPubCent] = {0.000070, 0.000162, 0.000165, 0.000134, 0.000126, 0.000072, 0.000070, 0.000069,
                                           0.000090, 0.000095, 0.000123, 0.000084, 0.000165, 0.000186, 0.000198, 0.000196};
Double_t V2OOPubErrSyst[numV2OOPubCent] = {0.000226, 0.000232, 0.000236, 0.000241, 0.000241, 0.000248, 0.000254, 0.000260,
                                           0.000263, 0.000265, 0.000266, 0.000264, 0.000263, 0.000259, 0.000257, 0.000251};

// Published V2 of charged particles in PbPb collisions (arxiv.org/pdf/1602.01119)
const Int_t numV2PbPbPubCent = 8;
Double_t V2PbPbPubCent[numV2PbPbPubCent + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
// Double_t V2PbPbPub[numV2PbPbPubCent] = {(0.02839 + 0.04566) / 2, 0.06551, 0.08707, 0.0991, 0.10414, 0.10286, 0.09746, 0.08881};
Double_t V2PbPbPub[numV2PbPbPubCent] = {(2047. * 0.02839 + 1668. * 0.04566) / (2047. + 1668.), 0.06551, 0.08707, 0.0991, 0.10414, 0.10286, 0.09746, 0.08881};
// Double_t V2PbPbPubErrStat[numV2PbPbPubCent] = {(0.00057 + 0.00064) / 2, 0.00037, 0.00044, 0.00055, 0.00073, 0.00107, 0.00186, 0.00438};
Double_t V2PbPbPubErrStat[numV2PbPbPubCent] = {(2047. * 0.00057 + 1668. * 0.00064) / (2047. + 1668.), 0.00037, 0.00044, 0.00055, 0.00073, 0.00107, 0.00186, 0.00438};
// Double_t V2PbPbPubErrSys[numV2PbPbPubCent] = {(0.00043 + 0.00069) / 2, 0.00098, 0.00131, 0.00149, 0.00156, 0.00154, 0.00146, 0.00133};
Double_t V2PbPbPubErrSys[numV2PbPbPubCent] = {(2047. * 0.00043 + 1668. * 0.00069) / (2047. + 1668.), 0.00098, 0.00131, 0.00149, 0.00156, 0.00154, 0.00146, 0.00133};