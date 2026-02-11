#include <Riostream.h>
#include <string>
#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TLine.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TCutG.h>
#include "TFitResult.h"
#include "TLegend.h"
// #include "CommonVar.h"
#include "CommonVarLambda.h"

void StyleCanvas(TCanvas *canvas, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(LMargin);
  canvas->SetRightMargin(RMargin);
  canvas->SetTopMargin(TMargin);
  canvas->SetBottomMargin(BMargin);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  // gStyle->SetPalette(55, 0);
}

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX,
                TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)
    histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->SetTitle(title);
}
void SetFont(TH1F *histo)
{
  histo->GetXaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetYaxis()->SetLabelFont(43);
}
void SetTickLength(TH1F *histo, Float_t TickLengthX, Float_t TickLengthY)
{
  histo->GetXaxis()->SetTickLength(TickLengthX);
  histo->GetYaxis()->SetTickLength(TickLengthY);
}

void SetHistoTextSize(TH1F *histo, Float_t XSize, Float_t XLabelSize, Float_t XOffset, Float_t XLabelOffset, Float_t YSize, Float_t YLabelSize, Float_t YOffset, Float_t YLabelOffset)
{
  histo->GetXaxis()->SetTitleSize(XSize);
  histo->GetXaxis()->SetLabelSize(XLabelSize);
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetXaxis()->SetLabelOffset(XLabelOffset);
  histo->GetYaxis()->SetTitleSize(YSize);
  histo->GetYaxis()->SetLabelSize(YLabelSize);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetYaxis()->SetLabelOffset(YLabelOffset);
}
void StylePad(TPad *pad, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  pad->SetFillColor(0);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetLeftMargin(LMargin);
  pad->SetRightMargin(RMargin);
  pad->SetTopMargin(TMargin);
  pad->SetBottomMargin(BMargin);
}
void StyleHistoYield(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(TitleX);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(TitleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset); // 1.2
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->SetTitle(title);
}

const Float_t UpperLimitLSBOmega = 1.655;  // upper limit of fit of left sidebands for omega
const Float_t LowerLimitRSBOmega = 1.689;  // lower limit of fit of right sidebands for omega
const Float_t UpperLimitLSBXi = 1.302;     // upper limit of fit of left sidebands for Xi
const Float_t LowerLimitRSBXi = 1.34;      // lower limit of fit of right sidebands for Xi
const Float_t UpperLimitLSBLambda = 1.107; // upper limit of fit of left sidebands for Lambda
const Float_t LowerLimitRSBLambda = 1.124; // lower limit of fit of right sidebands for Lambda
const Int_t numBinsEta = 8;

struct v2fit
{
  double operator()(double *x, double *par)
  {
    int bin = bkgfraction.FindBin(x[0]);
    float BkgFraction = bkgfraction.GetBinContent(bin);
    float SigFraction = 1.f - BkgFraction;
    return par[0] * SigFraction + (par[1] + par[2] * x[0]) * BkgFraction;
  }
  void setBkgFraction(TF1 *bkg, TF1 *total, float min, float max)
  {
    bkgfraction = TH1D(Form("fraction%s_%s", bkg->GetName(), total->GetName()), "", 2000, min, max);
    bkgfraction.Add(bkg);
    bkgfraction.Divide(total);
  }
  TH1D bkgfraction;
};

Double_t v2bkgfit(Double_t *x, Double_t *par)
{
  return par[0] + par[1] * x[0];
}

Bool_t reject = 1;
Double_t fparab(Double_t *x, Double_t *par)
{
  Float_t LimInf = 0;
  Float_t LimSup = 0;
  if (par[3] == 0)
  {
    LimInf = UpperLimitLSBXi;
    LimSup = LowerLimitRSBXi;
  }
  else if (par[3] == 1)
  {
    LimInf = UpperLimitLSBOmega;
    LimSup = LowerLimitRSBOmega;
  }
  else if (par[3] == 2)
  {
    LimInf = UpperLimitLSBLambda;
    LimSup = LowerLimitRSBLambda;
  }
  if (reject && x[0] > LimInf && x[0] < LimSup)
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}

Double_t fpol3(Double_t *x, Double_t *par)
{
  Float_t LimInf = 0;
  Float_t LimSup = 0;
  if (par[4] == 0)
  {
    LimInf = UpperLimitLSBXi;
    LimSup = LowerLimitRSBXi;
  }
  else if (par[4] == 1)
  {
    LimInf = UpperLimitLSBOmega;
    LimSup = LowerLimitRSBOmega;
  }
  else if (par[4] == 2)
  {
    LimInf = UpperLimitLSBLambda;
    LimSup = LowerLimitRSBLambda;
  }
  if (reject && x[0] > LimInf && x[0] < LimSup)
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
}

Double_t fexpo(Double_t *x, Double_t *par)
{
  Float_t LimInf = 0;
  Float_t LimSup = 0;
  if (par[2] == 0)
  {
    LimInf = UpperLimitLSBXi;
    LimSup = LowerLimitRSBXi;
  }
  else if (par[2] == 1)
  {
    LimInf = UpperLimitLSBOmega;
    LimSup = LowerLimitRSBOmega;
  }
  else if (par[2] == 2)
  {
    LimInf = UpperLimitLSBLambda;
    LimSup = LowerLimitRSBLambda;
  }
  if (reject && x[0] > LimInf && x[0] < LimSup)
  {
    TF1::RejectPoint();
    return 0;
  }
  // return par[0] + par[1] * exp(par[2] * x[0]);
  return exp(par[0] + par[1] * x[0]);
}

Double_t fretta(Double_t *x, Double_t *par)
{
  Float_t LimInf = 0;
  Float_t LimSup = 0;
  if (par[2] == 0)
  {
    LimInf = UpperLimitLSBXi;
    LimSup = LowerLimitRSBXi;
  }
  else if (par[2] == 1)
  {
    LimInf = UpperLimitLSBOmega;
    LimSup = LowerLimitRSBOmega;
  }
  else if (par[2] == 2)
  {
    LimInf = UpperLimitLSBLambda;
    LimSup = LowerLimitRSBLambda;
  }
  if (reject && x[0] > LimInf && x[0] < LimSup)
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0];
}

Float_t DefineMixedBDTValue(Int_t mul = 0, Int_t pt = 0)
{
  if (mul == numCent)
  {
    if (PtBins[pt] < 1.5)
    {
      return 0.8;
    }
    else if (PtBins[pt] < 2)
    {
      return 0.8;
    }
    else
      return 0.6;
  }
  if (CentFT0C[mul] >= 60)
  {
    if (PtBins[pt] < 1.5)
    {
      return 0.4;
    }
    else if (PtBins[pt] < 2)
    {
      return 0.28;
    }
    else
      return 0.2;
  }
  else if (CentFT0C[mul] >= 50)
  {
    if (PtBins[pt] < 2)
    {
      return 0.48;
    }
    else
      return 0.2;
  }
  else if (CentFT0C[mul] >= 40)
  {
    if (PtBins[pt] < 1.5)
    {
      return 0.56;
    }
    else if (PtBins[pt] < 2)
    {
      return 0.52;
    }
    else
      return 0.28;
  }
  else if (CentFT0C[mul] >= 20)
  {
    if (PtBins[pt] < 1.5)
    {
      return 0.8;
    }
    else if (PtBins[pt] < 2)
    {
      return 0.72;
    }
    else
      return 0.4;
  }
  else if (CentFT0C[mul] >= 0)
  {
    if (PtBins[pt] < 1.5)
    {
      return 0.96;
    }
    else if (PtBins[pt] < 2)
    {
      return 0.88;
    }
    else
      return 0.8;
  }
  else
  {
    cout << "Error: multiplicity out of range" << endl;
    return -1; // Error case
  }
}

TString titleYield = "1/N_{evt} dN/dp_{T}";
TString titleYieldNN = "dN/dp_{T}";
TString titlePt = "p_{T} (GeV/c)";
TString TitleInvMass[numPart] = {"(#Lambda, #pi)", "(#Lambda, K)", "(#Lambda, #pi^{-})", "(#overline{#Lambda}, #pi^{+})", "(#Lambda, K^{-})", "(#overline{#Lambda}, K^{+})", "(p, #pi)", "(p, #pi^{-})", "(#overline{p}, #pi^{+})"};
TString SInvMass = "invariant mass (GeV/c^{2})";

// fit ranges
Float_t min_range_signal[numPart] = {1.3, 1.65, 1.3, 1.3, 1.65, 1.65, 1.11, 1.11, 1.11}; // gauss fit range
Float_t max_range_signal[numPart] = {1.335, 1.69, 1.335, 1.335, 1.69, 1.69, 1.12, 1.12, 1.12};
Float_t liminf[numPart] = {1.29, 1.63, 1.29, 1.29, 1.63, 1.63, 1.1, 1.1, 1.1}; // bkg and total fit range
Float_t limsup[numPart] = {1.352, 1.71, 1.352, 1.352, 1.71, 1.71, 1.13, 1.13, 1.13};
Float_t liminfBkg[numPart] = {1.29, 1.63, 1.29, 1.29, 1.63, 1.63, 1.1, 1.1, 1.1}; // bkg and total fit range
Float_t limsupBkg[numPart] = {1.352, 1.71, 1.352, 1.352, 1.71, 1.71, 1.13, 1.13, 1.13};
Float_t liminfV2[numPart] = {1.308, 1.63, 1.29, 1.29, 1.63, 1.63, 1.1, 1.1, 1.1}; // v2 fit range
Float_t limsupV2[numPart] = {1.335, 1.71, 1.352, 1.352, 1.71, 1.71, 1.3, 1.3, 1.3};
Float_t XRangeMin[numPart] = {1.301, 1.656, 1.3, 1.3, 1.656, 1.656, 1.1, 1.1, 1.1};
Float_t XRangeMax[numPart] = {1.344, 1.688, 1.343, 1.343, 1.688, 1.688, 1.13, 1.13, 1.13};

// visualisation ranges
Float_t LowMassRange[numPart] = {1.31, 1.655, 1.31, 1.31, 1.655, 1.6551, 1.1, 1.1, 1.1}; // range to compute approximate yield (signal + bkg)
Float_t UpMassRange[numPart] = {1.33, 1.685, 1.33, 1.33, 1.685, 1.685, 1.13, 1.13, 1.13};
Float_t gaussDisplayRangeLow[numPart] = {1.29, 1.63, 1.29, 1.29, 1.63, 1.63, 1.1, 1.1, 1.1}; // display range of gauss functions (from total fit)
Float_t gaussDisplayRangeUp[numPart] = {1.35, 1.71, 1.35, 1.35, 1.71, 1.71, 1.13, 1.13, 1.13};
Float_t bkgDisplayRangeLow[numPart] = {1.29, 1.626, 1.29, 1.29, 1.626, 1.626, 1.1, 1.1, 1.1}; // display range of bkg function (from total fit)
Float_t bkgDisplayRangeUp[numPart] = {1.35, 1.72, 1.35, 1.35, 1.72, 1.72, 1.13, 1.13, 1.13};
Float_t histoMassRangeLow[numPart] = {1.301, 1.626, 1.301, 1.301, 1.626, 1.626, 1.1, 1.1, 1.1}; // display range of mass histograms
Float_t histoMassRangeUp[numPart] = {1.344, 1.72, 1.344, 1.344, 1.72, 1.72, 1.13, 1.13, 1.13};

// Event plane resolution
Float_t ftcReso[commonNumCent + 1] = {0};

void FitV2orPol(
    Bool_t isPtAnalysis = 1,    // 1 for V2 vs pt and Pzs2 vs pt, 0 for Pz vs 2(phi-Psi)
    Bool_t isPolFromLambda = 0, // 0: polarization of cascades computed directly, 1: polarization of cascades computed from polarization of lambdas
    Bool_t isBkgPol = 1,
    Int_t indexMultTrial = 0,
    Int_t mul = 0,
    Bool_t isMassCutForAcceptance = 1, // default is 1, with 0 we have no mass cut so to be able to do the fit of acceptance vs mass; only for LAMBDA
    Bool_t isTighterPzFitRange = 0,
    Int_t isTightMassForAcceptancePurity = 0,
    Int_t indexMassCut = 0,
    Bool_t isRapiditySel = ExtrisRapiditySel,
    Int_t ChosenPart = ChosenParticle,
    TString inputFileName = SinputFileName,
    Int_t EtaSysChoice = ExtrEtaSysChoice,
    Int_t BkgType = ExtrBkgType,
    Bool_t isLogy = 1,
    Bool_t isYAxisMassZoomed = 0,
    Bool_t UseTwoGauss = ExtrUseTwoGauss,
    Bool_t isMeanFixedPDG = 0,
    Bool_t isSysMultTrial = ExtrisSysMultTrial)
{
  if (isOOCentrality && commonNumCent != numCentLambdaOO)
  {
    cout << "O-O centrality is selected, but commonNumCent is not set to numCentLambdaOO. Please check the settings." << endl;
    return;
  }
  if (!isOOCentrality && commonNumCent != numCent)
  {
    cout << "Pb-Pb centrality is selected, but commonNumCent is not set to numCent. Please check the settings." << endl;
    return;
  }
  Int_t numCentMax = 0;
  if (isOOCentrality)
  {
    numCentMax = numCentLambdaOO;
  }
  else
  {
    numCentMax = numCent;
  }
  if (ChosenPart >= 6 && isPolFromLambda)
  {
    cout << "Polarization of lambda cannot be computed from polarization of lambdas :) Set isPolFromLambda = 0" << endl;
    return;
  }
  if (isReducedPtBins && numPtBins != numPtBinsReduced)
  {
    cout << "Reduced pt bins are selected, but numPtBins is not set to numPtBinsReduced. Please check the settings." << endl;
    return;
  }
  Float_t sigmacentral = Extrsigmacentral[isTightMassCut];
  Float_t ExtrLowLimit = 0;
  Float_t ExtrUpLimit = 0;
  if (ExtrisSysMassCut)
  {
    ExtrLowLimit = ExtrLowLimitSysXi[indexMassCut];
    ExtrUpLimit = ExtrUpLimitSysXi[indexMassCut];
  }

  if (useMixedBDTValueInFitMacro && !useCommonBDTValue)
  {
    cout << "If you set useMixedBDTValueInFitMacro = 1, input files are correct only if you set UseCommonBDTValue = 1" << endl;
    return;
  }
  if (useMixedBDTValueInFitMacro && ExtrisSysMultTrial)
  {
    cout << "UseMixedBDTValueInFitMacro must be set to zero to compute systematics" << endl;
    return;
  }
  if (isProducedAcceptancePlots && ExtrisSysMultTrial)
  {
    cout << "Either you produce acceptance plots or you study systematics" << endl;
    return;
  }
  if (isProducedAcceptancePlots && useMixedBDTValueInFitMacro)
  {
    cout << "Either you produce acceptance plots or you use mixed BDT values" << endl;
    return;
  }

  Int_t ParticleType = 0; // 1 for Xi, 0 for Omega, 2 for Lambda
  if (ChosenPart == 0 || ChosenPart == 2 || ChosenPart == 3)
    ParticleType = 1;
  else if (ChosenPart >= 6)
    ParticleType = 2;
  Int_t part = 0;        // Xi
  if (ParticleType == 1) // Omega
    part = 1;
  else if (ParticleType == 2) // Lambda
    part = 2;

  if (isV2 && !isPtAnalysis)
  {
    cout << "V2 analysis not available as a function of 2(phi-Psi)" << endl;
    return;
  }
  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  if (mul == numCent)
  { // 0-80%
    CentFT0CMin = 0;
    CentFT0CMax = 80;
  }
  else
  {
    CentFT0CMin = CentFT0C[mul];
    CentFT0CMax = CentFT0C[mul + 1];
  }
  if (isOOCentrality)
  {
    if (mul == numCentLambdaOO)
    {
      CentFT0CMin = 0;
      CentFT0CMax = CentFT0CMaxLambdaOO;
    }
    else
    {
      CentFT0CMin = CentFT0CLambdaOO[mul];
      CentFT0CMax = CentFT0CLambdaOO[mul + 1];
    }
  }

  cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << endl;
  TString histoNameMassvsV2 = "";
  TString histoNameMassvsCos2Theta = "";
  TString ProfileV2 = "";
  TString AcceptanceHisto = "";

  if (isSysMultTrial)
    inputFileName = SinputFileNameSyst;
  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trialsBDT)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;
  TString SBDT = "";

  Int_t NEvents = 0;
  TString PathInEvents = "../TreeForAnalysis/AnalysisResults_" + inputFileName + ".root";
  if (ChosenPart >= 6)
    PathInEvents = "../TreeForAnalysis/AnalysisResults_" + SinputFileNameAR + ".root";
  TFile *fileEvt = new TFile(PathInEvents, "");
  if (!fileEvt)
  {
    cout << "File Evt does not exist" << endl;
    return;
  }
  TDirectoryFile *dirEvt = (TDirectoryFile *)fileEvt->Get("lf-cascade-flow/histos");
  TH1F *hEvents = (TH1F *)dirEvt->Get("hEventCentrality");
  if (!hEvents)
  {
    cout << "hEvents not available " << endl;
    return;
  }
  for (Int_t b = 1; b <= hEvents->GetNbinsX(); b++)
  {
    if (hEvents->GetBinCenter(b) >= CentFT0CMin && hEvents->GetBinCenter(b) < CentFT0CMax)
    {
      NEvents += hEvents->GetBinContent(b);
    }
  }

  // Get resolution
  TString fileResoName = "";
  if (v2type == 1)
  {
    if (SinputFileName.Index("CFW") != -1)
      fileResoName = "../" + ResoFileName_SPCFW;
    else
      fileResoName = "../" + ResoFileName_SPLF;
  }
  else
  {
    if (SinputFileName.Index("CFW") != -1)
      fileResoName = "../" + ResoFileName_EPCFW;
    else
      fileResoName = "../" + ResoFileName_EPLF;
  }
  fileResoName = ResoFileName_EPCFW;
  fileResoName += ".root";
  if (ChosenPart >= 6)
    fileResoName = "../Resolution/" + SinputFileNameResoWeight;

  TFile *fileResoEP = new TFile(fileResoName, "");
  if (!fileResoEP)
  {
    cout << "File Reso does not exist" << endl;
    return;
  }
  TH1F *hReso = (TH1F *)fileResoEP->Get("hReso");
  TH1F *hReso080 = (TH1F *)fileResoEP->Get("hReso080");
  if (ChosenPart >= 6)
  {
    hReso = (TH1F *)fileResoEP->Get("hResoV0ATPCA");
    hReso080 = (TH1F *)fileResoEP->Get("hResoV0ATPCA080");
  }
  cout << "Reso name: " << fileResoName << endl;
  if (isOOCentrality)
  {
    if (mul == numCentLambdaOO)
      ftcReso[mul] = hReso080->GetBinContent(1);
    else
      ftcReso[mul] = hReso->GetBinContent(hReso->FindBin(CentFT0CLambdaOO[mul] + 0.001));
  }
  else
  {
    if (mul == numCent)
      ftcReso[mul] = hReso080->GetBinContent(1);
    else
      ftcReso[mul] = hReso->GetBinContent(hReso->FindBin(CentFT0C[mul] + 0.001));
  }
  cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << endl;
  cout << "Resolution: " << ftcReso[mul] << endl;

  Float_t UpperLimitLSB = 0;
  Float_t LowerLimitRSB = 0;
  if (ParticleType == 1)
  {
    UpperLimitLSB = UpperLimitLSBXi;
    LowerLimitRSB = LowerLimitRSBXi;
  }
  else if (ParticleType == 0)
  {
    UpperLimitLSB = UpperLimitLSBOmega;
    LowerLimitRSB = LowerLimitRSBOmega;
  }
  else if (ParticleType == 2)
  {
    UpperLimitLSB = UpperLimitLSBLambda;
    LowerLimitRSB = LowerLimitRSBLambda;
  }

  if (mul > numCentMax + 1)
  {
    cout << "Multiplicity out of range" << endl;
    return;
  }

  Int_t numPtBinsVar = numPtBins;
  if (!isPtAnalysis)
    numPtBinsVar = numPsiBins;
  cout << "Number of bins: " << numPtBinsVar << endl;
  if (numPtBinsVar > numPtBins)
  {
    cout << "Number of bins too large" << endl;
    return;
  }

  if (ChosenParticle >= 6)
    PtBins[0] = 0.5; // for Lambda

  if (ChosenPart >= 6)
  {
    if (isTighterPzFitRange)
    {
      // liminfV2[ChosenPart] = 1.103;
      liminfV2[ChosenPart] = 1.102;
      // limsupV2[ChosenPart] = 1.127;
      limsupV2[ChosenPart] = 1.128;
    }
    else
    {
      liminfV2[ChosenPart] = 1.1;
      limsupV2[ChosenPart] = 1.13;
    }
  }

  TString SPt[numPtBins + 1] = {""};
  TH1F *hInvMass[numPtBins + 1];
  TH1F *hInvMassDraw[numPtBins + 1];
  TH1F *hV2[numPtBins + 1];
  TH2F *hmassVsV2C[numPtBins + 1];
  TH1F *hV2MassIntegrated[numPtBins + 1];
  TF1 *fitV2SP[numPtBins + 1];

  TH1F *hCos2Theta[numPtBins + 1];
  TH2F *hmassVsCos2Theta[numPtBins + 1];
  TH1F *hCos2ThetaMassIntegrated[numPtBins + 1];

  const Int_t numCanvas = 4;
  TCanvas *canvas[numCanvas];
  TCanvas *canvasCos2Theta[numCanvas];
  for (Int_t c = 0; c < numCanvas; c++)
  {
    canvas[c] = new TCanvas(Form("canvas_%i", c), Form("canvas%i", c), 1800, 1400);
    canvas[c]->Divide(4, 2);
    StyleCanvas(canvas[c], 0.15, 0.05, 0.05, 0.15);
    canvasCos2Theta[c] = new TCanvas(Form("canvasCos2Theta_%i", c), Form("canvasCos2Theta%i", c), 1800, 1400);
    canvasCos2Theta[c]->Divide(4, 2);
    StyleCanvas(canvasCos2Theta[c], 0.15, 0.05, 0.05, 0.15);
  }

  Double_t BinsVar[numPtBins + 1] = {0};
  for (Int_t i = 0; i <= numPtBinsVar; i++)
  {
    if (isPtAnalysis)
      BinsVar[i] = PtBins[i];
    else
      BinsVar[i] = i * 2 * TMath::Pi() / numPsiBins;
  }

  TH1F *histoMean = new TH1F("histoMean", "histoMean", numPtBinsVar, BinsVar);
  TH1F *histoMeanPtInt = new TH1F("histoMeanPtInt", "histoMeanPtInt", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoMean1 = new TH1F("histoMean1", "histoMean1", numPtBinsVar, BinsVar);
  TH1F *histoMean1PtInt = new TH1F("histoMean1PtInt", "histoMean1PtInt", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoMean2 = new TH1F("histoMean2", "histoMean2", numPtBinsVar, BinsVar);
  TH1F *histoMean2PtInt = new TH1F("histoMean2PtInt", "histoMean2PtInt", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoSigma = new TH1F("histoSigma", "histoSigma", numPtBinsVar, BinsVar);
  TH1F *histoSigmaPtInt = new TH1F("histoSigmaPtInt", "histoSigmaPtInt", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoSigma1 = new TH1F("histoSigma1", "histoSigma1", numPtBinsVar, BinsVar);
  TH1F *histoSigma1PtInt = new TH1F("histoSigma1PtInt", "histoSigma1PtInt", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoSigma2 = new TH1F("histoSigma2", "histoSigma2", numPtBinsVar, BinsVar);
  TH1F *histoSigma2PtInt = new TH1F("histoSigma2PtInt", "histoSigma2PtInt", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoSigmaWeighted = new TH1F("histoSigmaWeighted", "histoSigmaWeighted", numPtBinsVar, BinsVar);
  TH1F *histoSigmaPtIntWeighted = new TH1F("histoSigmaPtIntWeighted", "histoSigmaPtIntWeighted", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoPurity = new TH1F("histoPurity", "histoPurity", numPtBinsVar, BinsVar);
  TH1F *histoPurityPtInt = new TH1F("histoPurityPtInt", "histoPurityPtInt", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoSignificance = new TH1F("histoSignificance", "histoSignificance", numPtBinsVar, BinsVar);
  TH1F *histoSignificancePtInt = new TH1F("histoSignificancePtInt", "histoSignificancePtInt", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoYield = new TH1F("histoYield", "histoYield", numPtBinsVar, BinsVar);
  TH1F *histoYieldPtInt = new TH1F("histoYieldPtInt", "histoYieldPtInt", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoYieldNN = new TH1F("histoYieldNotNormByEvts", "histoYieldNotNormByEvts", numPtBinsVar, BinsVar);
  TH1F *histoYieldFraction = new TH1F("histoYieldFraction", "histoYieldFraction", numPtBinsVar, BinsVar);
  TH1F *histoYieldFractionPtInt = new TH1F("histoYieldFractionPtInt", "histoYieldFractionPtInt", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoRelErrYield = new TH1F("histoRelErrYield", "histoRelErrYield", numPtBinsVar, BinsVar);
  TH1F *histoTot = new TH1F("histoTot", "histoTot", numPtBinsVar, BinsVar);
  TH1F *histoTotPtInt = new TH1F("histoTotPtInt", "histoTotPtInt", 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoB = new TH1F("histoB", "histoB", numPtBinsVar, BinsVar);
  TH1F *histoBPtInt = new TH1F("histoBPtInt", "histoBPtInt", 1, 0, BinsVar[numPtBinsVar]);

  TH1F *histoAppliedBDT = new TH1F("histoAppliedBDT", "histoAppliedBDT", numPtBinsVar, BinsVar);
  TH1F *histoAppliedBDTPtInt = new TH1F("histoAppliedBDTPtInt", "histoAppliedBDTPtInt", 1, 0, BinsVar[numPtBinsVar]);

  histoMeanPtInt->SetLineColor(kMagenta);
  histoMeanPtInt->SetMarkerColor(kMagenta);
  histoSigmaPtInt->SetLineColor(kMagenta);
  histoSigmaPtInt->SetMarkerColor(kMagenta);
  histoSigmaPtIntWeighted->SetLineColor(kMagenta);
  histoSigmaPtIntWeighted->SetMarkerColor(kMagenta);
  histoPurityPtInt->SetLineColor(kMagenta);
  histoPurityPtInt->SetMarkerColor(kMagenta);
  histoSignificancePtInt->SetLineColor(kMagenta);
  histoSignificancePtInt->SetMarkerColor(kMagenta);
  histoYieldPtInt->SetLineColor(kMagenta);
  histoYieldPtInt->SetMarkerColor(kMagenta);
  histoYieldFractionPtInt->SetLineColor(kMagenta);
  histoYieldFractionPtInt->SetMarkerColor(kMagenta);
  histoAppliedBDTPtInt->SetLineColor(kMagenta);
  histoAppliedBDTPtInt->SetMarkerColor(kMagenta);

  TString ShistoV2 = "histoV2";
  if (!isV2)
  { // polarization
    ShistoV2 = "histoPzs2";
    if (isPolFromLambda)
      ShistoV2 = "histoPzs2LambdaFromC";
    if (!isPtAnalysis)
    {
      ShistoV2 = "histoPz";
      if (isPolFromLambda)
        ShistoV2 = "histoPzLambdaFromC";
    }
  }
  TString ShistoCos2Theta = "histoCos2Theta";
  if (isPolFromLambda)
    ShistoCos2Theta = "histoCos2ThetaLambdaFromC";
  TString titlex = "#it{p}_{T} (GeV/#it{c})";
  if (!isPtAnalysis)
    titlex = "2(#varphi-#Psi_{EP})";
  TString TitleCos2Theta = "#LT cos^{2}(#theta_{#Lambda}*) #GT";
  if (isPolFromLambda)
    TitleCos2Theta = "#LT cos^{2}(#theta_{p}*) #GT";
  TString TitleCos2Theta_sig = "#LT cos^{2}(#theta_{#Lambda}*) #GT_{sig}";
  if (isPolFromLambda)
    TitleCos2Theta_sig = "#LT cos^{2}(#theta_{p}*) #GT_{sig}";

  TH1F *histoV2 = new TH1F(ShistoV2, Form(";%s ;#it{v}_{2}", titlex.Data()), numPtBinsVar, BinsVar);
  TH1F *histoV2NoFit = new TH1F(ShistoV2 + "NoFit", Form(";%s ;#it{v}_{2}", titlex.Data()), numPtBinsVar, BinsVar);
  TH1F *histoV2Mixed = new TH1F(ShistoV2 + "Mixed", Form(";%s ;#it{v}_{2}", titlex.Data()), numPtBinsVar, BinsVar);
  TH1F *histoV2MixedCorr = new TH1F(ShistoV2 + "MixedCorr", Form(";%s ;#it{v}_{2}", titlex.Data()), numPtBinsVar, BinsVar);
  TH1F *histoV2Bkg = new TH1F(ShistoV2 + "Bkg", Form(";%s ;#it{v}_{2}", titlex.Data()), numPtBinsVar, BinsVar);

  TH1F *histoV2Err = new TH1F(ShistoV2 + "Err", Form(";%s ;#it{v}_{2}", titlex.Data()), numPtBinsVar, BinsVar);
  TH1F *histoV2NoFitErr = new TH1F(ShistoV2 + "NoFitErr", Form(";%s ;#it{v}_{2}", titlex.Data()), numPtBinsVar, BinsVar);
  TH1F *histoV2MixedErr = new TH1F(ShistoV2 + "MixedErr", Form(";%s ;#it{v}_{2}", titlex.Data()), numPtBinsVar, BinsVar);

  TH1F *histoV2PtInt = new TH1F(ShistoV2 + "PtInt", Form(";%s ;#it{v}_{2}", titlex.Data()), 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoV2PtIntNoFit = new TH1F(ShistoV2 + "PtIntNoFit", Form(";%s ;#it{v}_{2}", titlex.Data()), 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoV2PtIntMixed = new TH1F(ShistoV2 + "PtIntMixed", Form(";%s ;#it{v}_{2}", titlex.Data()), 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoV2BkgPtInt = new TH1F(ShistoV2 + "BkgPtInt", Form(";%s ;#it{v}_{2}", titlex.Data()), 1, 0, BinsVar[numPtBinsVar]);

  TH1F *histoV2PtIntErr = new TH1F(ShistoV2 + "PtIntErr", Form(";%s ;#it{v}_{2}", titlex.Data()), 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoV2PtIntNoFitErr = new TH1F(ShistoV2 + "PtIntNoFitErr", Form(";%s ;#it{v}_{2}", titlex.Data()), 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoV2PtIntMixedErr = new TH1F(ShistoV2 + "PtIntMixedErr", Form(";%s ;#it{v}_{2}", titlex.Data()), 1, 0, BinsVar[numPtBinsVar]);

  TH1F *histoCos2Theta = new TH1F(ShistoCos2Theta, Form(";%s ;%s", titlex.Data(), TitleCos2Theta.Data()), numPtBinsVar, BinsVar);
  TH1F *histoCos2ThetaNoFit = new TH1F(ShistoCos2Theta + "NoFit", Form(";%s ;%s", titlex.Data(), TitleCos2Theta.Data()), numPtBinsVar, BinsVar);
  TH1F *histoCos2ThetaMixed = new TH1F(ShistoCos2Theta + "Mixed", Form(";%s ;%s", titlex.Data(), TitleCos2Theta.Data()), numPtBinsVar, BinsVar);
  TH1F *histoCos2ThetaPeakPos = new TH1F(ShistoCos2Theta + "PeakPos", Form(";%s ;%s", titlex.Data(), TitleCos2Theta.Data()), numPtBinsVar, BinsVar);
  TH1F *histoCos2ThetaPtInt = new TH1F(ShistoCos2Theta + "PtInt", Form(";%s ;%s", titlex.Data(), TitleCos2Theta.Data()), 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoCos2ThetaPtIntNoFit = new TH1F(ShistoCos2Theta + "PtIntNoFit", Form(";%s ;%s", titlex.Data(), TitleCos2Theta.Data()), 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoCos2ThetaPtIntMixed = new TH1F(ShistoCos2Theta + "PtIntMixed", Form(";%s ;%s", titlex.Data(), TitleCos2Theta.Data()), 1, 0, BinsVar[numPtBinsVar]);
  TH1F *histoCos2ThetaPtIntPeakPos = new TH1F(ShistoCos2Theta + "PtIntPeakPos", Form(";%s ;%s", titlex.Data(), TitleCos2Theta.Data()), 1, 0, BinsVar[numPtBinsVar]);

  TH2F *histoCos2ThetaNoFit2D = new TH2F(ShistoCos2Theta + "NoFit2D", Form(";%s ;%s", titlex.Data(), TitleCos2Theta.Data()), numPtBinsVar, BinsVar, 8, -0.8, 0.8);

  histoV2PtInt->SetLineColor(kMagenta);
  histoV2PtInt->SetMarkerColor(kMagenta);
  histoV2PtIntNoFit->SetLineColor(kMagenta);
  histoV2PtIntNoFit->SetMarkerColor(kMagenta);
  histoV2PtIntMixed->SetLineColor(kMagenta);
  histoV2PtIntMixed->SetMarkerColor(kMagenta);
  histoCos2ThetaPtInt->SetLineColor(kMagenta);
  histoCos2ThetaPtInt->SetMarkerColor(kMagenta);
  histoCos2ThetaPtIntNoFit->SetLineColor(kMagenta);
  histoCos2ThetaPtIntNoFit->SetMarkerColor(kMagenta);
  histoCos2ThetaPtIntPeakPos->SetLineColor(kMagenta);
  histoCos2ThetaPtIntPeakPos->SetMarkerColor(kMagenta);

  Double_t PhiBins[numPsiBins + 1];

  TString SPathIn = "";
  TString SPathInPtInt = "";
  TFile *filein;
  Float_t BDTscoreCutPtInt_checkValue = 0;
  for (Int_t pt = 0; pt < numPtBinsVar + 1; pt++)
  {
    // Define mult and pt-dependent BDt selection
    if (useMixedBDTValueInFitMacro)
    {
      BDTscoreCut = DefineMixedBDTValue(mul, pt);
      if (pt == numPtBinsVar)
      {
        BDTscoreCut = BDTscoreCutPtInt[mul];
        if (isTightMassCut)
          BDTscoreCut = BDTscoreCutPtIntLoosest[mul];
      }
    }
    else if (isProducedAcceptancePlots)
    {
      BDTscoreCut = BDTscoreCutAcceptance[mul];
    }
    if (BDTscoreCut != DefaultBDTscoreCut || isSysMultTrial)
      SBDT = Form("_BDT%.3f", BDTscoreCut);
    else
      SBDT = "";

    if (pt == numPtBinsVar)
      BDTscoreCutPtInt_checkValue = BDTscoreCut;

    histoAppliedBDT->SetBinContent(pt + 1, BDTscoreCut);
    if (pt == numPtBinsVar)
      histoAppliedBDTPtInt->SetBinContent(1, BDTscoreCut);

    SPathIn = "../OutputAnalysis/V2_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice];
    if (ChosenPart < 6)
      SPathIn += SBDT;
    if (isApplyWeights)
      SPathIn += "_Weighted";
    if (isApplyCentWeight && !isProducedAcceptancePlots)
      SPathIn += "_CentWeighted";
    if (v2type == 1)
      SPathIn += "_SP";
    if (!useCommonBDTValue)
      SPathIn += "_BDTCentDep";
    if (isRun2Binning)
      SPathIn += "_Run2Binning";
    if (ExtrisApplyEffWeights)
    {
      SPathIn += "_EffW";
    }
    SPathIn += "_WithAlpha";
    if (!isRapiditySel)
      SPathIn += "_Eta08";
    if (isReducedPtBins)
      SPathIn += "_ReducedPtBins";

    if (ChosenPart >= 6 && ExtrisSysLambdaMultTrial)
    {
      if (isLoosest)
        SPathIn += "_isLoosest";
      else if (isTightest)
        SPathIn += "_isTightest";
      else
        SPathIn += Form("_SysMultTrial_%i", indexMultTrial);
      SPathIn += "_isSysLambdaMultTrial";
    }
    SPathIn += STHN[ExtrisFromTHN];
    if (ExtrisApplyResoOnTheFly && !isProducedAcceptancePlots)
      SPathIn += "_ResoOnTheFly";
    // if (ChosenPart >= 6)
    // SPathIn += "_CorrectReso_TestLeassPtBins";
    // SPathIn += "_SystReso";
    if (ChosenPart >= 6 && !isMassCutForAcceptance && isProducedAcceptancePlots)
      SPathIn += "_NoMassCutForAcceptance";
    // SPathIn += "_TestMoreBins";
    SPathIn += ".root";

    if (pt == numPtBinsVar)
      SPathInPtInt = SPathIn;

    filein = new TFile(SPathIn, "");
    if (!filein)
    {
      cout << "FileIn not available" << endl;
      return;
    }

    if (!isPtAnalysis)
    {
      if (pt == numPtBinsVar)
        continue; // skip the integrated
    }
    PhiBins[pt] = pt * 2 * TMath::Pi() / numPsiBins;
    if (ParticleType == 0 && pt == 0)
      continue;
    SPt[pt] = Form("%.2f < p_{T} < %.2f", PtBins[pt], PtBins[pt + 1]);
    if (pt == numPtBinsVar)
    { // integrated
      SPt[pt] = Form("%.2f < p_{T} < %.2f", PtBins[0], PtBins[numPtBins]);
    }
    if (!isPtAnalysis) // psi bins
      SPt[pt] = Form("%.2f < #psi < %.2f", PhiBins[pt], PhiBins[pt] + 2 * TMath::Pi() / numPsiBins - 0.0001);

    if (ChosenPart < 6)
      cout << "\nFor the centrality: " << CentFT0CMin << "-" << CentFT0CMax << " % and the pt: " << SPt[pt] << " the BDT cut is: " << BDTscoreCut << endl;
    cout << "FileIn: " << SPathIn << endl;

    if (isPtAnalysis)
      cout << "Analysed pt interval: " << SPt[pt] << endl;
    else
      cout << "Analysed psi interval: " << SPt[pt] << endl;

    if (isPtAnalysis)
      hInvMass[pt] = (TH1F *)filein->Get(Form("mass_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt));
    else
      hInvMass[pt] = (TH1F *)filein->Get(Form("mass_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, pt));
    if (!hInvMass[pt])
    {
      cout << "Histogram inv. mass not available" << endl;
      return;
    }

    StyleHisto(hInvMass[pt], 0, 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), 1, 20, TitleInvMass[ChosenPart] + " " + SInvMass, "Counts", SPt[pt] + " GeV/#it{c}", 1, histoMassRangeLow[ChosenPart], histoMassRangeUp[ChosenPart], 1.4, 1.6, 0.7);

    if (isV2)
    {
      histoNameMassvsV2 = Form("MassvsV2C_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      ProfileV2 = Form("V2C_cent%i-%i_pt%i_Profile", CentFT0CMin, CentFT0CMax, pt);
      histoNameMassvsCos2Theta = Form("MassvsCos2Theta_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      AcceptanceHisto = Form("Cos2Theta_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
    }
    else
    { // polarization
      histoNameMassvsV2 = Form("MassvsPzs2_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      ProfileV2 = Form("Pzs2_cent%i-%i_pt%i_Profile", CentFT0CMin, CentFT0CMax, pt);
      histoNameMassvsCos2Theta = Form("MassvsCos2Theta_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      AcceptanceHisto = Form("Cos2Theta_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      if (isPolFromLambda)
      { // polarization from lambdas
        histoNameMassvsV2 = Form("MassvsPzs2LambdaFromC_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
        ProfileV2 = Form("Pzs2LambdaFromC_cent%i-%i_pt%i_Profile", CentFT0CMin, CentFT0CMax, pt);
        histoNameMassvsCos2Theta = Form("MassvsCos2ThetaLambdaFromC_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
        AcceptanceHisto = Form("Cos2ThetaLambdaFromC_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      }
      if (!isPtAnalysis)
      {
        histoNameMassvsV2 = Form("MassvsPz_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, pt);
        ProfileV2 = Form("Pz_cent%i-%i_psi%i_Profile", CentFT0CMin, CentFT0CMax, pt);
        histoNameMassvsCos2Theta = Form("MassvsCos2Theta_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, pt);
        AcceptanceHisto = Form("Cos2Theta_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, pt);
        if (isPolFromLambda)
        {
          histoNameMassvsV2 = Form("MassvsPzLambdaFromC_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, pt);
          ProfileV2 = Form("PzLambdaFromC_cent%i-%i_psi%i_Profile", CentFT0CMin, CentFT0CMax, pt);
          histoNameMassvsCos2Theta = Form("MassvsCos2ThetaLambdaFromC_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, pt);
          AcceptanceHisto = Form("Cos2ThetaLambdaFromC_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, pt);
        }
      }
    }

    hmassVsV2C[pt] = (TH2F *)filein->Get(histoNameMassvsV2);
    if (!hmassVsV2C[pt])
    {
      cout << "Histogram " << histoNameMassvsV2 << " not available" << endl;
      return;
    }

    hmassVsCos2Theta[pt] = (TH2F *)filein->Get(histoNameMassvsCos2Theta);
    if (!hmassVsCos2Theta[pt])
    {
      cout << "Histogram " << histoNameMassvsCos2Theta << " not available" << endl;
      return;
    }

    // hV2[pt] = (TH1F *)filein->Get(Form("V2C_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt));
    hV2[pt] = (TH1F *)filein->Get(ProfileV2);
    if (!hV2[pt])
    {
      cout << "Histogram hV2 not available" << endl;
      return;
    }

    Float_t MaxV2 = 0.02;
    Float_t MinV2 = -0.02;
    if (mul > 4)
    {
      MaxV2 = 0.05;
      MinV2 = -0.05;
    }
    StyleHisto(hV2[pt], -MinV2, MaxV2, 1, 20, titlePt, "v_{2}", TitleInvMass[ChosenPart] + " " + SInvMass, 1, 0, 100, 1.4, 1.6, 0.7);
    hV2[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[ChosenPart], histoMassRangeUp[ChosenPart]);

    hCos2Theta[pt] = (TH1F *)filein->Get(AcceptanceHisto);
    if (!hCos2Theta[pt])
    {
      cout << "Histogram hCos2Theta not available" << endl;
      return;
    }
    StyleHisto(hCos2Theta[pt], -1, 1, 1, 20, titlePt, TitleCos2Theta, TitleInvMass[ChosenPart] + " " + SInvMass, 1, 0, 100, 1.4, 1.6, 0.7);
    hCos2Theta[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[ChosenPart], histoMassRangeUp[ChosenPart]);

    if (isYAxisMassZoomed)
    {
      if (ParticleType == 1)
        hInvMass[pt]->GetYaxis()->SetRangeUser(0, 2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(1.65)));
      else if (ParticleType == 0)
        hInvMass[pt]->GetYaxis()->SetRangeUser(0, 2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(1.29)));
      else if (ParticleType == 2)
        hInvMass[pt]->GetYaxis()->SetRangeUser(0, 2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(1.55)));
    }
    if (isLogy)
    {
      if (hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(histoMassRangeUp[ChosenPart])) > hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(histoMassRangeLow[ChosenPart])))
      {
        hInvMass[pt]->SetMinimum(0.8 * hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(histoMassRangeLow[ChosenPart])));
      }
      else
      {
        hInvMass[pt]->SetMinimum(0.8 * hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(histoMassRangeUp[ChosenPart])));
      }
      hInvMass[pt]->SetMaximum(1.2 * hInvMass[pt]->GetMaximum());
    }
    if (pt < 4)
      canvas[0]->cd(pt + 1);
    else if (pt < 8)
      canvas[1]->cd(pt + 1 - 4);
    else if (pt < 12)
      canvas[2]->cd(pt + 1 - 8);
    else if (pt < 16)
      canvas[3]->cd(pt + 1 - 12);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.2);
    if (isLogy)
      gPad->SetLogy();
    hInvMass[pt]->Draw("e same");
  }

  // fits

  TF1 **functionsFirst = new TF1 *[numPtBins + 1];
  TF1 **functionsSecond = new TF1 *[numPtBins + 1];
  TF1 **functions1 = new TF1 *[numPtBins + 1];
  TF1 **functions2 = new TF1 *[numPtBins + 1];
  TF1 **bkg1 = new TF1 *[numPtBins + 1];
  TF1 **bkg2 = new TF1 *[numPtBins + 1];
  TF1 **bkg3 = new TF1 *[numPtBins + 1];
  TF1 **bkg4 = new TF1 *[numPtBins + 1];
  TF1 **bkgretta = new TF1 *[numPtBins + 1]; // initial bkg fit
  TF1 **bkgparab = new TF1 *[numPtBins + 1]; // initial bkg fit
  TF1 **bkgpol3 = new TF1 *[numPtBins + 1];  // initial bkg fit
  TF1 **bkgexpo = new TF1 *[numPtBins + 1];  // initial bkg fit
  TF1 **total = new TF1 *[numPtBins + 1];
  TF1 **totalbis = new TF1 *[numPtBins + 1];
  TF1 **totalSignal = new TF1 *[numPtBins + 1];

  v2fit v2fitarray[numPtBins + 1];
  TF1 **v2FitFunction = new TF1 *[numPtBins + 1];
  TF1 **v2BkgFunction = new TF1 *[numPtBins + 1];
  TF1 **Cos2ThetaFitFunction = new TF1 *[numPtBins + 1];
  TF1 **Cos2ThetaBkgFunction = new TF1 *[numPtBins + 1];

  Double_t parTwoGaussRetta[numPtBins + 1][8];
  Double_t parTwoGaussParab[numPtBins + 1][9];
  Double_t parTwoGaussPol3[numPtBins + 1][10];
  Double_t parTwoGaussExpo[numPtBins + 1][8];
  Double_t parOneGaussRetta[numPtBins + 1][5];
  Double_t parOneGaussParab[numPtBins + 1][6];
  Double_t parOneGaussPol3[numPtBins + 1][7];
  Double_t parOneGaussExpo[numPtBins + 1][5];

  TFitResultPtr fFitResultPtr0[numPtBins + 1];
  TFitResultPtr fFitResultPtr1[numPtBins + 1];
  TFitResultPtr fFitV2Bkg[numPtBins + 1];

  Double_t mean[numPtBins + 1] = {0};
  Double_t errmean[numPtBins + 1] = {0};
  Double_t sigma[numPtBins + 1] = {0};
  Double_t errsigma[numPtBins + 1] = {0};
  Double_t sigmaw[numPtBins + 1] = {0};
  Double_t errsigmaw[numPtBins + 1] = {0};
  Double_t w1[numPtBins + 1] = {0};
  Double_t w2[numPtBins + 1] = {0};
  Double_t I12[numPtBins + 1] = {0};
  Double_t Yield[numPtBins + 1] = {0};
  Double_t ErrYield[numPtBins + 1] = {0};
  Double_t LowLimit[numPtBins + 1] = {0};
  Double_t UpLimit[numPtBins + 1] = {0};
  Double_t b[numPtBins + 1] = {0};
  Double_t errb[numPtBins + 1] = {0};
  Double_t SSB[numPtBins + 1] = {0};
  Double_t errSSB[numPtBins + 1] = {0};
  Double_t Signif[numPtBins + 1] = {0};
  Double_t errSignif[numPtBins + 1] = {0};
  Double_t FitIntegral[numPtBins + 1] = {0};
  Double_t YieldFromFit[numPtBins + 1] = {0};
  Double_t entries_range[numPtBins + 1] = {0};
  Double_t TotYield = 0;
  Double_t TotSigBkg = 0;

  TLine *lineP3Sigma[numPtBins + 1];
  TLine *lineM3Sigma[numPtBins + 1];
  TLine *lineP3SigmaNorm[numPtBins + 1];
  TLine *lineM3SigmaNorm[numPtBins + 1];

  Float_t bTest[numPtBins + 1] = {0};
  Float_t errbTest[numPtBins + 1] = {0};
  Float_t SignalTest[numPtBins + 1] = {0};
  Float_t errSignalTest[numPtBins + 1] = {0};
  Float_t YieldTest[numPtBins + 1] = {0};
  Float_t ErrYieldTest[numPtBins + 1] = {0};
  TH1F *hYieldTest[numPtBins + 1];
  TH1F *hYieldRelErrorTest[numPtBins + 1];
  TH1F *hYieldRelErrorTestRelative[numPtBins + 1];
  Float_t LowLimitTest[numPtBins + 1] = {0};
  Float_t UpLimitTest[numPtBins + 1] = {0};
  Float_t LowBin0[numPtBins + 1] = {0};
  Float_t UpBin0[numPtBins + 1] = {0};
  Float_t LowBin[numPtBins + 1] = {0};
  Float_t UpBin[numPtBins + 1] = {0};
  Int_t numMassInt = 10;
  TF1 *LineAt1 = new TF1("LineAt1", "[0]+[1]*x", 0, numMassInt);
  LineAt1->SetParameter(0, 1);
  TF1 *LineAt995 = new TF1("LineAt995", "[0]+[1]*x", 0, numMassInt);
  LineAt995->SetParameter(0, 0.995);

  TLine *lineBkgLimitA[numPtBins + 1];
  TLine *lineBkgLimitB[numPtBins + 1];
  TLine *lineBkgLimitC[numPtBins + 1];
  TLine *lineBkgLimitD[numPtBins + 1];

  Bool_t isV2FromFit[numPtBins + 1] = {0};

  for (Int_t pt = 0; pt < numPtBinsVar + 1; pt++)
  {
    if (!isPtAnalysis)
    {
      if (pt == numPtBinsVar)
        continue; // skip the integrated
    }

    if (ParticleType == 0 && pt == 0)
      continue;
    if (pt < 4)
      canvas[0]->cd(pt + 1);
    else if (pt < 8)
      canvas[1]->cd(pt + 1 - 4);
    else if (pt < 12)
      canvas[2]->cd(pt + 1 - 8);
    else if (pt < 16)
      canvas[3]->cd(pt + 1 - 12);

    functionsFirst[pt] = new TF1(Form("1f_%i", pt), "gaus", min_range_signal[ChosenPart], max_range_signal[ChosenPart]);
    functionsFirst[pt]->SetLineColor(881);
    functionsFirst[pt]->SetParameter(1, ParticleMassPDG[ChosenPart]);
    functionsFirst[pt]->SetParName(0, "norm");
    functionsFirst[pt]->SetParName(1, "mean");
    functionsFirst[pt]->SetParName(2, "sigma");
    functionsFirst[pt]->SetParLimits(1, min_range_signal[ChosenPart], max_range_signal[ChosenPart]);
    functionsFirst[pt]->SetParLimits(2, 0.001, 0.1);
    functionsFirst[pt]->SetParLimits(0, 0, 1.1 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));

    functionsSecond[pt] = new TF1(Form("2f_%i", pt), "gaus", min_range_signal[ChosenPart], max_range_signal[ChosenPart]);
    functionsSecond[pt]->SetLineColor(867);
    functionsSecond[pt]->SetParameter(1, ParticleMassPDG[ChosenPart]);
    functionsSecond[pt]->SetParName(0, "norm");
    functionsSecond[pt]->SetParName(1, "mean");
    functionsSecond[pt]->SetParName(2, "sigma");
    functionsSecond[pt]->SetParLimits(1, min_range_signal[ChosenPart], max_range_signal[ChosenPart]);
    functionsSecond[pt]->SetParLimits(2, 0.001, 0.15);
    functionsSecond[pt]->SetParLimits(0, 0, 1.1 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));

    functions1[pt] = new TF1(Form("1f_%i_final", pt), "gaus", gaussDisplayRangeLow[ChosenPart], gaussDisplayRangeUp[ChosenPart]);
    functions1[pt]->SetLineColor(kRed); // 867
    functions1[pt]->SetParName(0, "norm");
    functions1[pt]->SetParName(1, "mean");
    functions1[pt]->SetParName(2, "sigma");

    functions2[pt] = new TF1(Form("2f_%i_final", pt), "gaus", gaussDisplayRangeLow[ChosenPart], gaussDisplayRangeUp[ChosenPart]);
    functions2[pt]->SetLineColor(kMagenta); // 891
    functions2[pt]->SetParName(0, "norm");
    functions2[pt]->SetParName(1, "mean");
    functions2[pt]->SetParName(2, "sigma");

    bkg1[pt] = new TF1(Form("bkg1%i", pt), "pol1", bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
    bkg1[pt]->SetLineColor(kGreen + 7);
    bkg1[pt]->SetLineStyle(2);

    bkg2[pt] = new TF1(Form("bkg2%i", pt), "pol2", bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
    bkg2[pt]->SetLineColor(1);
    bkg2[pt]->SetLineStyle(2);

    bkg3[pt] = new TF1(Form("bkg3%i", pt), "pol3", bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
    bkg3[pt]->SetLineColor(kOrange + 7);
    bkg3[pt]->SetLineStyle(2);

    bkg4[pt] = new TF1(Form("bkg4%i", pt), "expo", bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
    bkg4[pt]->SetLineColor(kOrange + 7);
    bkg4[pt]->SetLineStyle(2);

    bkgretta[pt] = new TF1(Form("retta%i", pt), fretta, liminf[ChosenPart], limsup[ChosenPart], 3);
    bkgretta[pt]->SetLineColor(kGreen + 3);
    bkgretta[pt]->FixParameter(2, part);

    bkgparab[pt] = new TF1(Form("parab%i", pt), fparab, liminf[ChosenPart], limsup[ChosenPart], 4);
    bkgparab[pt]->SetLineColor(kAzure + 7);
    bkgparab[pt]->FixParameter(3, part);

    bkgpol3[pt] = new TF1(Form("pol3%i", pt), fpol3, liminf[ChosenPart], limsup[ChosenPart], 5);
    bkgpol3[pt]->SetLineColor(kRed + 7);
    bkgpol3[pt]->FixParameter(4, part);

    bkgexpo[pt] = new TF1(Form("expo%i", pt), fexpo, liminf[ChosenPart], limsup[ChosenPart], 3);
    bkgexpo[pt]->SetLineColor(kGreen + 2);
    bkgexpo[pt]->FixParameter(2, part);
    Bool_t UseTwoGaussUpdated = 1;
    TF1 *totalFunction = nullptr, *bkgFunction = nullptr;
    if (UseTwoGauss)
    {
      cout << "\n\e[35mFit with two gauss \e[39m" << endl;
      if (isPtAnalysis)
        cout << "Analysed pt interval: " << SPt[pt] << endl;
      else
        cout << " Psi: " << PhiBins[pt] << "-" << PhiBins[pt] + 2 * TMath::Pi() / numPsiBins << endl;

      if (BkgType == 0)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+gaus(3)+pol1(6)", liminf[ChosenPart], limsup[ChosenPart]);
      else if (BkgType == 1)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+gaus(3)+pol2(6)", liminf[ChosenPart], limsup[ChosenPart]);
      else if (BkgType == 2)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+gaus(3)+pol3(6)", liminf[ChosenPart], limsup[ChosenPart]);
      else
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+gaus(3)+expo(6)", liminf[ChosenPart], limsup[ChosenPart]);
      total[pt]->SetLineColor(597);
      total[pt]->SetParName(0, "norm");
      total[pt]->SetParName(1, "mean");
      total[pt]->SetParName(2, "sigma");
      total[pt]->SetParName(3, "norm2");
      total[pt]->SetParName(4, "mean2");
      total[pt]->SetParName(5, "sigma2");

      cout << "\n\n fit gauss1 " << endl;
      hInvMass[pt]->Fit(functionsFirst[pt], "R");
      cout << "\n\n fit gauss2 " << endl;
      hInvMass[pt]->Fit(functionsSecond[pt], "RB");

      bkg1[pt]->SetRange(bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
      bkg2[pt]->SetRange(bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
      bkg3[pt]->SetRange(bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
      bkg4[pt]->SetRange(bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
      bkgparab[pt]->SetRange(liminfBkg[ChosenPart], limsupBkg[ChosenPart]);
      bkgretta[pt]->SetRange(liminf[ChosenPart], limsup[ChosenPart]);
      bkgexpo[pt]->SetRange(liminf[ChosenPart], limsup[ChosenPart]);
      bkgpol3[pt]->SetRange(liminf[ChosenPart], limsup[ChosenPart]);
      total[pt]->SetRange(liminf[ChosenPart], limsup[ChosenPart]);

      cout << "\n\n fit bkg " << endl;
      if (BkgType == 0)
        hInvMass[pt]->Fit(bkgretta[pt], "RB0");
      else if (BkgType == 1)
        hInvMass[pt]->Fit(bkgparab[pt], "RB0");
      else if (BkgType == 2)
        hInvMass[pt]->Fit(bkgpol3[pt], "RB0");
      else
        hInvMass[pt]->Fit(bkgexpo[pt], "RB");

      if (BkgType == 0)
      {
        functionsFirst[pt]->GetParameters(&parTwoGaussRetta[pt][0]);
        functionsSecond[pt]->GetParameters(&parTwoGaussRetta[pt][3]);
        bkgretta[pt]->GetParameters(&parTwoGaussRetta[pt][6]);
        total[pt]->SetParameters(parTwoGaussRetta[pt]);
      }
      else if (BkgType == 1)
      {
        functionsFirst[pt]->GetParameters(&parTwoGaussParab[pt][0]);
        functionsSecond[pt]->GetParameters(&parTwoGaussParab[pt][3]);
        bkgparab[pt]->GetParameters(&parTwoGaussParab[pt][6]);
        total[pt]->SetParameters(parTwoGaussParab[pt]);
      }
      else if (BkgType == 2)
      {
        functionsFirst[pt]->GetParameters(&parTwoGaussPol3[pt][0]);
        functionsSecond[pt]->GetParameters(&parTwoGaussPol3[pt][3]);
        bkgpol3[pt]->GetParameters(&parTwoGaussPol3[pt][6]);
        total[pt]->SetParameters(parTwoGaussPol3[pt]);
      }
      else
      {
        functionsFirst[pt]->GetParameters(&parTwoGaussExpo[pt][0]);
        functionsSecond[pt]->GetParameters(&parTwoGaussExpo[pt][3]);
        bkgexpo[pt]->GetParameters(&parTwoGaussExpo[pt][6]);
        total[pt]->SetParameters(parTwoGaussExpo[pt]);
      }

      cout << "\n\n fit total " << endl;
      if (ParticleType == 1) // Xi
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.318, 1.324);
        if (mul == 6)
          total[pt]->SetParLimits(2, 0.0015, 0.010);
        else
          total[pt]->SetParLimits(2, 0.0012, 0.010);
        total[pt]->SetParLimits(3, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin())); // maximum was wothout 0.3
        total[pt]->SetParLimits(4, 1.318, 1.324);
        total[pt]->SetParLimits(5, 0.001, 0.01);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, ParticleMassPDG[ChosenPart]);
          total[pt]->FixParameter(4, ParticleMassPDG[ChosenPart]);
        }
      }
      else if (ParticleType == 0) // Omega
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.66, 1.68);
        total[pt]->SetParLimits(2, 0.002, 0.01);
        total[pt]->SetParLimits(3, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin())); // maximum was wothout 0.3
        total[pt]->SetParLimits(4, 1.66, 1.68);
        total[pt]->SetParLimits(5, 0.001, 0.01);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, ParticleMassPDG[ChosenPart]);
          total[pt]->FixParameter(4, ParticleMassPDG[ChosenPart]);
        }
      }
      else if (ParticleType == 2) // Lambda
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), 0.3 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.1, 1.13);
        total[pt]->SetParLimits(2, 0.002, 0.01);
        total[pt]->SetParLimits(3, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin())); // maximum was wothout 0.3
        total[pt]->SetParLimits(4, 1.1, 1.13);
        total[pt]->SetParLimits(5, 0.001, 0.01);
        // if (mul >= 7)
        // total[pt]->SetParLimits(2, 0.001, 0.010); //0.0015

        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, ParticleMassPDG[ChosenPart]);
          total[pt]->FixParameter(4, ParticleMassPDG[ChosenPart]);
        }
        if (BkgType == 0)
        {
          total[pt]->SetParLimits(1, 1.113, 1.117);
          total[pt]->SetParLimits(4, 1.113, 1.117);
        }
      }

      fFitResultPtr0[pt] = hInvMass[pt]->Fit(total[pt], "SRB+"); // per errore gaussiana, S indica che il risultato del fit e' accessibile da fFitResultPtr0
      // la gaussiana più larga deve esserte quella più bassa
      if (total[pt]->GetParameter(2) > total[pt]->GetParameter(5))
      {
        if (total[pt]->GetParameter(0) > total[pt]->GetParameter(3))
          UseTwoGaussUpdated = kFALSE;
      }
      else
      {
        if (total[pt]->GetParameter(0) < total[pt]->GetParameter(3))
          UseTwoGaussUpdated = kFALSE;
      }

      cout << "UseTwoGauss = " << UseTwoGaussUpdated << endl;

      totalbis[pt] = (TF1 *)total[pt]->Clone(Form("totalbis_pt%i", pt));
      fFitResultPtr1[pt] = fFitResultPtr0[pt];
      totalSignal[pt] = (TF1 *)total[pt]->Clone(Form("totalSignal_pt%i", pt));

      functions1[pt]->FixParameter(0, total[pt]->GetParameter(0));
      functions1[pt]->FixParameter(1, total[pt]->GetParameter(1));
      functions1[pt]->FixParameter(2, total[pt]->GetParameter(2));
      functions2[pt]->FixParameter(0, total[pt]->GetParameter(3));
      functions2[pt]->FixParameter(1, total[pt]->GetParameter(4));
      functions2[pt]->FixParameter(2, total[pt]->GetParameter(5));

      totalbis[pt]->FixParameter(0, 0);
      totalbis[pt]->FixParameter(1, 0);
      totalbis[pt]->FixParameter(2, 0);
      totalbis[pt]->FixParameter(3, 0);
      totalbis[pt]->FixParameter(4, 0);
      totalbis[pt]->FixParameter(5, 0);

      totalSignal[pt]->FixParameter(6, 0);
      totalSignal[pt]->FixParameter(7, 0);
      if (BkgType == 1)
      {
        totalSignal[pt]->FixParameter(8, 0);
      }
      else if (BkgType == 2)
      {
        totalSignal[pt]->FixParameter(8, 0);
        totalSignal[pt]->FixParameter(9, 0);
      }
      totalSignal[pt]->SetLineColor(kOrange + 2);

      totalFunction = (TF1 *)total[pt]->Clone(Form("totalFunction_%i", pt));
      if (BkgType == 0)
      {
        bkg1[pt]->FixParameter(0, total[pt]->GetParameter(6));
        bkg1[pt]->FixParameter(1, total[pt]->GetParameter(7));
        bkgFunction = (TF1 *)bkg1[pt]->Clone(Form("bkgFunction1_%i", pt));
      }
      else if (BkgType == 1)
      {
        bkg2[pt]->FixParameter(0, total[pt]->GetParameter(6));
        bkg2[pt]->FixParameter(1, total[pt]->GetParameter(7));
        bkg2[pt]->FixParameter(2, total[pt]->GetParameter(8));
        bkgFunction = (TF1 *)bkg2[pt]->Clone(Form("bkgFunction2_%i", pt));
      }
      else if (BkgType == 2)
      {
        bkg3[pt]->FixParameter(0, total[pt]->GetParameter(6));
        bkg3[pt]->FixParameter(1, total[pt]->GetParameter(7));
        bkg3[pt]->FixParameter(2, total[pt]->GetParameter(8));
        bkg3[pt]->FixParameter(3, total[pt]->GetParameter(9));
        bkgFunction = (TF1 *)bkg3[pt]->Clone(Form("bkgFunction3_%i", pt));
      }
      else
      {
        bkg4[pt]->FixParameter(0, total[pt]->GetParameter(6));
        bkg4[pt]->FixParameter(1, total[pt]->GetParameter(7));
        bkgFunction = (TF1 *)bkg4[pt]->Clone(Form("bkgFunction4_%i", pt));
      }

      if (UseTwoGaussUpdated)
      {

        if (pt < 4)
          canvas[0]->cd(pt + 1);
        else if (pt < 8)
          canvas[1]->cd(pt + 1 - 4);
        else if (pt < 12)
          canvas[2]->cd(pt + 1 - 8);
        else if (pt < 16)
          canvas[3]->cd(pt + 1 - 12);

        hInvMass[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[ChosenPart], histoMassRangeUp[ChosenPart]);
        hInvMass[pt]->Draw("same e");
        functions1[pt]->Draw("same");
        functions2[pt]->Draw("same");
        if (BkgType == 0)
          bkg1[pt]->Draw("same");
        else if (BkgType == 1)
          bkg2[pt]->Draw("same");
        else if (BkgType == 2)
          bkg3[pt]->Draw("same");
        else
          bkg4[pt]->Draw("same");

        TMatrixDSym cov = fFitResultPtr0[pt]->GetCovarianceMatrix();
        Double_t cov_mean = cov[1][4];
        Double_t cov_sigma = cov[2][5];
        I12[pt] = functions1[pt]->Integral(0, 2) + functions2[pt]->Integral(0, 2);
        cout << "SignalIntegral " << I12[pt] << endl;
        w1[pt] = functions1[pt]->Integral(0, 2) / I12[pt];
        w2[pt] = functions2[pt]->Integral(0, 2) / I12[pt];
        mean[pt] = (functions1[pt]->GetParameter(1) + functions2[pt]->GetParameter(1)) / 2;
        errmean[pt] = (total[pt]->GetParError(1) + total[pt]->GetParError(4)) / 2;
        sigma[pt] = (functions1[pt]->GetParameter(2) + functions2[pt]->GetParameter(2)) / 2;
        errsigma[pt] = sqrt(pow(total[pt]->GetParError(2), 2) + pow(total[pt]->GetParError(5), 2) + 2 * cov_sigma) / 2;
        sigmaw[pt] = (functions1[pt]->GetParameter(2) * functions1[pt]->Integral(0, 2) + functions2[pt]->GetParameter(2) * functions2[pt]->Integral(0, 2)) / I12[pt];
        errsigmaw[pt] = sqrt(pow(w1[pt] / I12[pt] * functions1[pt]->GetParError(2), 2) + pow(w2[pt] / I12[pt] * functions2[pt]->GetParError(2), 2) + 2 * functions2[pt]->Integral(0, 2) * functions1[pt]->Integral(0, 2) * cov_sigma) / I12[pt];
      }
    }
    if (!UseTwoGaussUpdated || !UseTwoGauss)
    {
      cout << "\n\e[36mFit with one gauss only: \e[39m"
           << " Pt: " << PtBins[pt] << "-" << PtBins[pt + 1] << endl;

      if (BkgType == 0)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+pol1(3)", liminf[ChosenPart], limsup[ChosenPart]);
      else if (BkgType == 1)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+pol2(3)", liminf[ChosenPart], limsup[ChosenPart]);
      else if (BkgType == 2)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+pol3(3)", liminf[ChosenPart], limsup[ChosenPart]);
      else
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+expo(3)", liminf[ChosenPart], limsup[ChosenPart]);

      total[pt]->SetLineColor(7);
      total[pt]->SetParName(0, "norm");
      total[pt]->SetParName(1, "mean");
      total[pt]->SetParName(2, "sigma");

      cout << "\n\n fit gauss " << endl;
      hInvMass[pt]->Fit(functionsFirst[pt], "RB");

      bkg1[pt]->SetRange(bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
      bkg2[pt]->SetRange(bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
      bkg3[pt]->SetRange(bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
      bkg4[pt]->SetRange(bkgDisplayRangeLow[ChosenPart], bkgDisplayRangeUp[ChosenPart]);
      bkgparab[pt]->SetRange(liminf[ChosenPart], limsup[ChosenPart]);
      bkgretta[pt]->SetRange(liminf[ChosenPart], limsup[ChosenPart]);
      bkgpol3[pt]->SetRange(liminf[ChosenPart], limsup[ChosenPart]);
      bkgexpo[pt]->SetRange(liminf[ChosenPart], limsup[ChosenPart]);
      total[pt]->SetRange(liminf[ChosenPart], limsup[ChosenPart]);

      cout << "\n\n fit bkg " << endl;
      if (BkgType == 0)
        hInvMass[pt]->Fit(bkgretta[pt], "RB0");
      else if (BkgType == 1)
        hInvMass[pt]->Fit(bkgparab[pt], "RB0");
      else if (BkgType == 2)
        hInvMass[pt]->Fit(bkgpol3[pt], "RB0");
      else
        hInvMass[pt]->Fit(bkgexpo[pt], "RB0");

      if (BkgType == 0)
      {
        functionsFirst[pt]->GetParameters(&parOneGaussRetta[pt][0]);
        bkgretta[pt]->GetParameters(&parOneGaussRetta[pt][3]);
        total[pt]->SetParameters(parOneGaussRetta[pt]);
      }
      else if (BkgType == 1)
      {
        functionsFirst[pt]->GetParameters(&parOneGaussParab[pt][0]);
        bkgparab[pt]->GetParameters(&parOneGaussParab[pt][3]);
        total[pt]->SetParameters(parOneGaussParab[pt]);
      }
      else if (BkgType == 2)
      {
        functionsFirst[pt]->GetParameters(&parOneGaussPol3[pt][0]);
        bkgpol3[pt]->GetParameters(&parOneGaussPol3[pt][3]);
        total[pt]->SetParameters(parOneGaussPol3[pt]);
      }
      else
      {
        functionsFirst[pt]->GetParameters(&parOneGaussExpo[pt][0]);
        bkgexpo[pt]->GetParameters(&parOneGaussExpo[pt][3]);
        total[pt]->SetParameters(parOneGaussExpo[pt]);
      }

      cout << "\n\n fit total " << endl;
      if (ParticleType == 1) // Xi
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.31, 1.335);
        total[pt]->SetParLimits(2, 0.0012, 0.010);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, ParticleMassPDG[ChosenPart]);
        }
      }
      else if (ParticleType == 0) // Omega
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.66, 1.68);
        total[pt]->SetParLimits(2, 0.001, 0.02);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, ParticleMassPDG[ChosenPart]);
        }
      }
      else if (ParticleType == 2) // Lambda
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.1, 1.13);
        total[pt]->SetParLimits(2, 0.001, 0.02);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, ParticleMassPDG[ChosenPart]);
        }
        if (BkgType == 0)
        {
          total[pt]->SetParLimits(1, 1.113, 1.117);
        }
      }

      cout << "max value " << hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()) << endl;
      fFitResultPtr0[pt] = hInvMass[pt]->Fit(total[pt], "SRB0"); // per errore gaussiana, S indica che il risultato del fit e' accessibile da fFitResultPtr0

      totalbis[pt] = (TF1 *)total[pt]->Clone();
      fFitResultPtr1[pt] = fFitResultPtr0[pt];

      functions1[pt]->FixParameter(0, total[pt]->GetParameter(0));
      functions1[pt]->FixParameter(1, total[pt]->GetParameter(1));
      functions1[pt]->FixParameter(2, total[pt]->GetParameter(2));

      totalbis[pt]->FixParameter(0, 0);
      totalbis[pt]->FixParameter(1, 0);
      totalbis[pt]->FixParameter(2, 0);

      totalSignal[pt] = (TF1 *)total[pt]->Clone(Form("totalSignal1Gaus_pt%i", pt));
      totalSignal[pt]->FixParameter(3, 0);
      totalSignal[pt]->FixParameter(4, 0);
      if (BkgType == 1)
      {
        totalSignal[pt]->FixParameter(5, 0);
      }
      else if (BkgType == 2)
      {
        totalSignal[pt]->FixParameter(5, 0);
        totalSignal[pt]->FixParameter(6, 0);
      }
      totalSignal[pt]->SetLineColor(kOrange + 2);

      totalFunction = (TF1 *)total[pt]->Clone(Form("totalFunction_%i", pt));
      if (BkgType == 0)
      {
        bkg1[pt]->FixParameter(0, total[pt]->GetParameter(3));
        bkg1[pt]->FixParameter(1, total[pt]->GetParameter(4));
        bkgFunction = (TF1 *)bkg1[pt]->Clone(Form("bkgFunction1_%i", pt));
      }
      else if (BkgType == 1)
      {
        bkg2[pt]->FixParameter(0, total[pt]->GetParameter(3));
        bkg2[pt]->FixParameter(1, total[pt]->GetParameter(4));
        bkg2[pt]->FixParameter(2, total[pt]->GetParameter(5));
        bkgFunction = (TF1 *)bkg2[pt]->Clone(Form("bkgFunction2_%i", pt));
      }
      else if (BkgType == 2)
      {
        bkg3[pt]->FixParameter(0, total[pt]->GetParameter(3));
        bkg3[pt]->FixParameter(1, total[pt]->GetParameter(4));
        bkg3[pt]->FixParameter(2, total[pt]->GetParameter(5));
        bkg3[pt]->FixParameter(3, total[pt]->GetParameter(6));
        bkgFunction = (TF1 *)bkg3[pt]->Clone(Form("bkgFunction3_%i", pt));
      }
      else
      {
        bkg4[pt]->FixParameter(0, total[pt]->GetParameter(3));
        bkg4[pt]->FixParameter(1, total[pt]->GetParameter(4));
        // bkg4[pt]->FixParameter(2, total[pt]->GetParameter(5));
        bkgFunction = (TF1 *)bkg4[pt]->Clone(Form("bkgFunction4_%i", pt));
      }
      if (pt < 4)
        canvas[0]->cd(pt + 1);
      else if (pt < 8)
        canvas[1]->cd(pt + 1 - 4);
      else if (pt < 12)
        canvas[2]->cd(pt + 1 - 8);
      else if (pt < 16)
        canvas[3]->cd(pt + 1 - 12);

      hInvMass[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[ChosenPart], histoMassRangeUp[ChosenPart]);
      hInvMass[pt]->Draw("same e");
      mean[pt] = total[pt]->GetParameter(1);
      errmean[pt] = total[pt]->GetParError(1);
      sigma[pt] = total[pt]->GetParameter(2);
      errsigma[pt] = total[pt]->GetParError(2);
      sigmaw[pt] = sigma[pt];
      errsigmaw[pt] = errsigma[pt];
    }

    cout << "\nMean: " << mean[pt] << " +/- " << errmean[pt] << endl;
    cout << "Sigma: " << sigma[pt] << " +/- " << errsigma[pt] << endl;

    TLine *linebkgFitLL = new TLine(liminf[ChosenPart], 0, liminf[ChosenPart], hInvMass[pt]->GetMaximum()); // low limit of left SB
    TLine *linebkgFitRR = new TLine(limsup[ChosenPart], 0, limsup[ChosenPart], hInvMass[pt]->GetMaximum()); // upper limit of right SB
    linebkgFitLL->SetLineColor(kBlue);
    linebkgFitRR->SetLineColor(kBlue);
    TLine *linebkgFitLR = new TLine(UpperLimitLSB, 0, UpperLimitLSB, hInvMass[pt]->GetMaximum()); // upper limit of left SB
    TLine *linebkgFitRL = new TLine(LowerLimitRSB, 0, LowerLimitRSB, hInvMass[pt]->GetMaximum()); // lower limit of right SB
    linebkgFitLR->SetLineColor(kBlue);
    linebkgFitRL->SetLineColor(kBlue);
    lineBkgLimitA[pt] = new TLine(liminf[ChosenPart], 0, liminf[ChosenPart], hInvMass[pt]->GetMaximum());
    lineBkgLimitB[pt] = new TLine(UpperLimitLSB, 0, UpperLimitLSB, hInvMass[pt]->GetMaximum());
    lineBkgLimitC[pt] = new TLine(limsup[ChosenPart], 0, limsup[ChosenPart], hInvMass[pt]->GetMaximum());
    lineBkgLimitD[pt] = new TLine(LowerLimitRSB, 0, LowerLimitRSB, hInvMass[pt]->GetMaximum());
    lineBkgLimitA[pt]->SetLineColor(kViolet + 1);
    lineBkgLimitB[pt]->SetLineColor(kViolet + 1);
    lineBkgLimitC[pt]->SetLineColor(kViolet + 1);
    lineBkgLimitD[pt]->SetLineColor(kViolet + 1);
    // lineBkgLimitA[pt]->Draw("same");
    // lineBkgLimitB[pt]->Draw("same");
    // lineBkgLimitC[pt]->Draw("same");
    // lineBkgLimitD[pt]->Draw("same");

    // linebkgFitLL->Draw("same");
    // linebkgFitRR->Draw("same");
    // linebkgFitLR->Draw("same");
    // linebkgFitRL->Draw("same");

    /*
        if (BkgType == 1)
          bkgparab[pt]->Draw("same");
        else
          bkgretta[pt]->Draw("same");
    */
    if (pt < 4)
      canvas[0]->cd(pt + 1);
    else if (pt < 8)
      canvas[1]->cd(pt + 1 - 4);
    else if (pt < 12)
      canvas[2]->cd(pt + 1 - 8);
    else if (pt < 16)
      canvas[3]->cd(pt + 1 - 12);

    LowLimit[pt] = hInvMass[pt]->GetXaxis()->GetBinLowEdge(hInvMass[pt]->GetXaxis()->FindBin(mean[pt] - sigmacentral * sigmaw[pt]));
    UpLimit[pt] = hInvMass[pt]->GetXaxis()->GetBinUpEdge(hInvMass[pt]->GetXaxis()->FindBin(mean[pt] + sigmacentral * sigmaw[pt]));
    if (ExtrisSysMassCut)
    {
      LowLimit[pt] = ExtrLowLimit;
      UpLimit[pt] = ExtrUpLimit;
    }
    if (isTightMassForAcceptancePurity && ChosenPart >= 6)
    {
      // LowLimit[pt] = 1.112;
      // UpLimit[pt] = 1.119;
      //  LowLimit[pt] = 1.1145;
      //  UpLimit[pt] = 1.1158;
      // LowLimit[pt] = 1.1145;
      // UpLimit[pt] = 1.1158;
      LowLimit[pt] = hInvMass[pt]->GetXaxis()->GetBinLowEdge(hInvMass[pt]->GetXaxis()->FindBin(1.1145));
      UpLimit[pt] = hInvMass[pt]->GetXaxis()->GetBinUpEdge(hInvMass[pt]->GetXaxis()->FindBin(1.1158));
    }
    lineP3Sigma[pt] = new TLine(UpLimit[pt], 0, UpLimit[pt], hInvMass[pt]->GetMaximum());
    lineM3Sigma[pt] = new TLine(LowLimit[pt], 0, LowLimit[pt], hInvMass[pt]->GetMaximum());
    lineP3Sigma[pt]->Draw("same");
    lineM3Sigma[pt]->Draw("same");

    b[pt] = 0;
    errb[pt] = 0;
    if (BkgType == 0)
    {
      b[pt] = bkg1[pt]->Integral(LowLimit[pt], UpLimit[pt]);
      errb[pt] = totalbis[pt]->IntegralError(LowLimit[pt], UpLimit[pt], fFitResultPtr1[pt]->GetParams(),
                                             (fFitResultPtr1[pt]->GetCovarianceMatrix()).GetMatrixArray());
    }
    else if (BkgType == 1)
    {
      b[pt] = bkg2[pt]->Integral(LowLimit[pt], UpLimit[pt]);
      errb[pt] = totalbis[pt]->IntegralError(LowLimit[pt], UpLimit[pt], fFitResultPtr1[pt]->GetParams(),
                                             (fFitResultPtr1[pt]->GetCovarianceMatrix()).GetMatrixArray());
    }
    else if (BkgType == 2)
    {
      b[pt] = bkg3[pt]->Integral(LowLimit[pt], UpLimit[pt]);
      errb[pt] = totalbis[pt]->IntegralError(LowLimit[pt], UpLimit[pt], fFitResultPtr1[pt]->GetParams(),
                                             (fFitResultPtr1[pt]->GetCovarianceMatrix()).GetMatrixArray());
    }
    else
    {
      b[pt] = bkg4[pt]->Integral(LowLimit[pt], UpLimit[pt]);
      errb[pt] = totalbis[pt]->IntegralError(LowLimit[pt], UpLimit[pt], fFitResultPtr1[pt]->GetParams(),
                                             (fFitResultPtr1[pt]->GetCovarianceMatrix()).GetMatrixArray());
    }

    b[pt] = b[pt] / hInvMass[pt]->GetBinWidth(1);
    errb[pt] = errb[pt] / hInvMass[pt]->GetBinWidth(1);
    entries_range[pt] = 0;
    for (Int_t l = hInvMass[pt]->GetXaxis()->FindBin(LowLimit[pt] + 0.00001); l <= hInvMass[pt]->GetXaxis()->FindBin(UpLimit[pt] - 0.00001); l++)
    { // I inlcude bins where the limits lie
      entries_range[pt] += hInvMass[pt]->GetBinContent(l);
    }

    Yield[pt] = entries_range[pt] - b[pt];
    ErrYield[pt] = sqrt(entries_range[pt] + pow(errb[pt], 2));
    TotYield += Yield[pt];
    TotSigBkg += entries_range[pt];

    SSB[pt] = (entries_range[pt] - b[pt]) / entries_range[pt];
    errSSB[pt] = SSB[pt] * sqrt(1. / entries_range[pt] + pow(errb[pt] / b[pt], 2));

    Signif[pt] = Yield[pt] / sqrt(entries_range[pt]);
    errSignif[pt] = sqrt(pow(ErrYield[pt] / sqrt(entries_range[pt]), 2) + pow(Yield[pt] / (2 * sqrt(entries_range[pt]) * entries_range[pt]), 2));

    YieldFromFit[pt] = totalSignal[pt]->Integral(LowLimit[pt], UpLimit[pt]);
    FitIntegral[pt] = totalSignal[pt]->Integral(histoMassRangeLow[ChosenPart], histoMassRangeUp[ChosenPart]);
    cout << "FitIntegral " << FitIntegral[pt] << endl;
    cout << "YieldFromFit " << YieldFromFit[pt] << endl;

    //*********************************************
    if (pt < numPtBinsVar)
    {
      histoYield->SetBinContent(pt + 1, Yield[pt] / NEvents / histoYield->GetBinWidth(pt + 1));
      histoYield->SetBinError(pt + 1, ErrYield[pt] / NEvents / histoYield->GetBinWidth(pt + 1));

      histoYieldNN->SetBinContent(pt + 1, Yield[pt] / histoYield->GetBinWidth(pt + 1));
      histoYieldNN->SetBinError(pt + 1, ErrYield[pt] / histoYield->GetBinWidth(pt + 1));

      histoRelErrYield->SetBinContent(pt + 1, ErrYield[pt] / Yield[pt]);
      histoRelErrYield->SetBinError(pt + 1, 0);

      histoYieldFraction->SetBinContent(pt + 1, YieldFromFit[pt] / FitIntegral[pt]);
      histoYieldFraction->SetBinError(pt + 1, 0);

      histoTot->SetBinContent(pt + 1, entries_range[pt] / NEvents / histoYield->GetBinWidth(pt + 1));
      histoTot->SetBinError(pt + 1, sqrt(entries_range[pt]) / NEvents / histoYield->GetBinWidth(pt + 1));

      histoB->SetBinContent(pt + 1, b[pt] / NEvents / histoYield->GetBinWidth(pt + 1));
      histoB->SetBinError(pt + 1, errb[pt] / NEvents / histoYield->GetBinWidth(pt + 1));

      histoMean->SetBinContent(pt + 1, mean[pt]);
      histoMean->SetBinError(pt + 1, errmean[pt]);

      histoMean1->SetBinContent(pt + 1, functions1[pt]->GetParameter(1));
      histoMean1->SetBinError(pt + 1, total[pt]->GetParError(1));

      histoMean2->SetBinContent(pt + 1, functions2[pt]->GetParameter(1));
      histoMean2->SetBinError(pt + 1, total[pt]->GetParError(4));

      histoSigma->SetBinContent(pt + 1, sigma[pt]);
      histoSigma->SetBinError(pt + 1, errsigma[pt]);

      histoSigmaWeighted->SetBinContent(pt + 1, sigmaw[pt]);
      histoSigmaWeighted->SetBinError(pt + 1, errsigmaw[pt]);

      histoSigma1->SetBinContent(pt + 1, functions1[pt]->GetParameter(2));
      histoSigma1->SetBinError(pt + 1, total[pt]->GetParError(2));

      histoSigma2->SetBinContent(pt + 1, functions2[pt]->GetParameter(2));
      histoSigma2->SetBinError(pt + 1, total[pt]->GetParError(5));

      histoPurity->SetBinContent(pt + 1, SSB[pt]);
      histoPurity->SetBinError(pt + 1, errSSB[pt]);

      histoSignificance->SetBinContent(pt + 1, Yield[pt] / ErrYield[pt]);
      histoSignificance->SetBinError(pt + 1, 0);
    }

    if (!bkgFunction || !totalFunction)
      continue;
    v2fitarray[pt].setBkgFraction(bkgFunction, totalFunction, liminf[ChosenPart], limsup[ChosenPart]);
    if (ChosenPart >= 6 && !isTighterPzFitRange)
    {
      if (mul == 8)
      {
        liminfV2[ChosenPart] = 1.101;
        limsupV2[ChosenPart] = 1.129;
      }
    }
    v2FitFunction[pt] = new TF1(Form("v2function%i", pt), v2fitarray[pt], liminfV2[ChosenPart], limsupV2[ChosenPart], 3);
    v2FitFunction[pt]->SetLineColor(kRed + 1);
    if (pt < 4)
      canvas[0]->cd(pt + 4 + 1);
    else if (pt < 8)
      canvas[1]->cd(pt + 4 + 1 - 4);
    else if (pt < 12)
      canvas[2]->cd(pt + 4 + 1 - 8);
    else if (pt < 16)
      canvas[3]->cd(pt + 4 + 1 - 12);
    cout << "Pt " << SPt[pt] << " GeV/c" << endl;
    cout << "Fitting the V2 / polarization... " << endl;
    if (isBkgPol == 0)
    {
      v2FitFunction[pt]->FixParameter(1, 0); // bkg v2 constant
      v2FitFunction[pt]->FixParameter(2, 0); // bkg v2 constant
    }
    // hV2[pt]->Rebin(2);
    // hV2[pt]->Fit(v2FitFunction[pt], "R0");
    fFitV2Bkg[pt] = hV2[pt]->Fit(v2FitFunction[pt], "SRB+");

    v2BkgFunction[pt] = new TF1(Form("v2bkgfunction%i", pt), v2bkgfit, liminf[ChosenPart], limsup[ChosenPart], 2);
    v2BkgFunction[pt]->FixParameter(0, v2FitFunction[pt]->GetParameter(1)); // quota
    v2BkgFunction[pt]->FixParameter(1, v2FitFunction[pt]->GetParameter(2)); // pendenza
    v2BkgFunction[pt]->SetLineColor(kBlack);
    v2BkgFunction[pt]->SetLineStyle(8);
    TMatrixDSym covV2 = fFitV2Bkg[pt]->GetCovarianceMatrix();
    Double_t covV2Bkg = covV2[1][2];
    Double_t PzBkgError = sqrt(pow(mean[pt] * v2FitFunction[pt]->GetParError(2), 2) + pow(v2FitFunction[pt]->GetParError(1), 2) + 2 * covV2Bkg * mean[pt]);
    cout << "Eval Pz,bkg at peak position: " << v2BkgFunction[pt]->Eval(mean[pt]) << endl;
    // cout << "sigmaq : " << v2FitFunction[pt]->GetParError(1) << "sigmam : " << v2FitFunction[pt]->GetParError(2) << " covar: " << covV2Bkg << endl;
    cout << "Error: " << PzBkgError << endl;

    hV2[pt]->GetXaxis()->SetTitle(TitleInvMass[ChosenPart] + " " + SInvMass);
    hV2[pt]->SetTitle(SPt[pt] + " GeV/#it{c}");
    cout << "Pt " << SPt[pt] << " GeV/c" << endl;
    cout << "v2: " << v2FitFunction[pt]->GetParameter(0) << " +- " << v2FitFunction[pt]->GetParError(0) << endl;
    cout << "Rel. error: " << v2FitFunction[pt]->GetParError(0) / v2FitFunction[pt]->GetParameter(0) << endl;
    if (pt < numPtBinsVar)
    {
      histoV2->SetBinContent(pt + 1, v2FitFunction[pt]->GetParameter(0));
      histoV2->SetBinError(pt + 1, v2FitFunction[pt]->GetParError(0));
      histoV2Bkg->SetBinContent(pt + 1, v2BkgFunction[pt]->Eval(mean[pt]));
      histoV2Bkg->SetBinError(pt + 1, PzBkgError);
    }
    else
    {
      histoV2PtInt->SetBinContent(1, v2FitFunction[pt]->GetParameter(0));
      histoV2PtInt->SetBinError(1, v2FitFunction[pt]->GetParError(0));
      histoV2BkgPtInt->SetBinContent(1, v2BkgFunction[pt]->Eval(mean[pt]));
      histoV2BkgPtInt->SetBinError(1, PzBkgError);
      histoMeanPtInt->SetBinContent(1, mean[pt]);
      histoMeanPtInt->SetBinError(1, errmean[pt]);
      histoMean1PtInt->SetBinContent(1, functions1[pt]->GetParameter(1));
      histoMean1PtInt->SetBinError(1, total[pt]->GetParError(1));
      histoMean2PtInt->SetBinContent(1, functions2[pt]->GetParameter(1));
      histoMean2PtInt->SetBinError(1, total[pt]->GetParError(4));
      histoSigmaPtInt->SetBinContent(1, sigma[pt]);
      histoSigmaPtInt->SetBinError(1, errsigma[pt]);
      histoSigma1PtInt->SetBinContent(1, functions1[pt]->GetParameter(2));
      histoSigma1PtInt->SetBinError(1, total[pt]->GetParError(2));
      histoSigma2PtInt->SetBinContent(1, functions2[pt]->GetParameter(2));
      histoSigma2PtInt->SetBinError(1, total[pt]->GetParError(5));
      histoSigmaPtIntWeighted->SetBinContent(1, sigmaw[pt]);
      histoSigmaPtIntWeighted->SetBinError(1, errsigmaw[pt]);
      histoPurityPtInt->SetBinContent(1, SSB[pt]);
      histoPurityPtInt->SetBinError(1, errSSB[pt]);
      histoYieldPtInt->SetBinContent(1, Yield[pt] / NEvents / histoYieldPtInt->GetBinWidth(1));
      histoYieldPtInt->SetBinError(1, ErrYield[pt] / NEvents / histoYieldPtInt->GetBinWidth(1));
      histoYieldFractionPtInt->SetBinContent(1, YieldFromFit[pt] / FitIntegral[pt]);
      histoYieldFractionPtInt->SetBinError(1, 0);
      histoSignificancePtInt->SetBinContent(1, Yield[pt] / ErrYield[pt]);
      histoSignificancePtInt->SetBinError(1, 0);
      histoTotPtInt->SetBinContent(1, entries_range[pt] / NEvents / histoYieldPtInt->GetBinWidth(1));
      histoTotPtInt->SetBinError(1, sqrt(entries_range[pt]) / NEvents / histoYieldPtInt->GetBinWidth(1));
      histoBPtInt->SetBinContent(1, b[pt] / NEvents / histoYieldPtInt->GetBinWidth(1));
      histoBPtInt->SetBinError(1, errb[pt] / NEvents / histoYieldPtInt->GetBinWidth(1));
    }

    hV2MassIntegrated[pt] = (TH1F *)hmassVsV2C[pt]->ProjectionX(Form("V2CvsMass_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt), hmassVsV2C[pt]->GetYaxis()->FindBin(LowLimit[pt] + 0.00001), hmassVsV2C[pt]->GetYaxis()->FindBin(UpLimit[pt] - 0.00001));
    if (isV2)
      StyleHisto(hV2MassIntegrated[pt], 0, 1.2 * hV2MassIntegrated[pt]->GetBinContent(hV2MassIntegrated[pt]->GetMaximumBin()), 1, 20, "v_{2}", "Counts", SPt[pt] + " GeV/#it{c}", 1, -1, 1, 1.4, 1.6, 0.7);
    // else StyleHisto(hV2MassIntegrated[pt], 0, 1.2 * hV2MassIntegrated[pt]->GetBinContent(hV2MassIntegrated[pt]->GetMaximumBin()), 1, 20, "v_{2}", "Counts", SPt[pt] + " GeV/#it{c}", 1, -100, 100, 1.4, 1.6, 0.7);
    //  hV2MassIntegrated[pt]->Rebin(2);

    hCos2ThetaMassIntegrated[pt] = (TH1F *)hmassVsCos2Theta[pt]->ProjectionX(Form("Cos2ThetavsMass_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt), hmassVsCos2Theta[pt]->GetYaxis()->FindBin(LowLimit[pt]), hmassVsCos2Theta[pt]->GetYaxis()->FindBin(UpLimit[pt]));
    StyleHisto(hCos2ThetaMassIntegrated[pt], 0, 1.2 * hCos2ThetaMassIntegrated[pt]->GetBinContent(hCos2ThetaMassIntegrated[pt]->GetMaximumBin()), 1, 20, "cos(2#theta)", "Counts", SPt[pt] + " GeV/#it{c}", 1, -1, 1, 1.4, 1.6, 0.7);
    // hCos2ThetaMassIntegrated[pt]->Rebin(2);

    // if (SSB[pt] > LimitForV2woFit || PtBins[pt] > 2.)
    if (SSB[pt] > LimitForV2woFit || (PtBins[pt] > 2. && CentFT0CMin >= 40 && pt != numPtBinsVar))
      isV2FromFit[pt] = 0;
    else
      isV2FromFit[pt] = 1;
    if (isTightMassCut)
      isV2FromFit[pt] = 0;
    if (isV2FromFit[pt])
    {
      hV2[pt]->Draw("e");
      v2FitFunction[pt]->Draw("same");
    }
    if (!isV2FromFit[pt])
      hV2MassIntegrated[pt]->Draw("");

    TLegend *legendV2 = new TLegend(0.2, 0.8, 0.4, 0.9);
    legendV2->SetTextSize(0.055);
    legendV2->AddEntry("", Form("Mean = %.3f +- %.3f", hV2MassIntegrated[pt]->GetMean(), hV2MassIntegrated[pt]->GetMeanError()), "");
    legendV2->Draw("same");

    // canvas for cos2theta
    if (pt < 4)
      canvasCos2Theta[0]->cd(pt + 4 + 1);
    else if (pt < 8)
      canvasCos2Theta[1]->cd(pt + 4 + 1 - 4);
    else if (pt < 12)
      canvasCos2Theta[2]->cd(pt + 4 + 1 - 8);
    else if (pt < 16)
      canvasCos2Theta[3]->cd(pt + 4 + 1 - 12);
    Cos2ThetaFitFunction[pt] = new TF1(Form("cosfunction%i", pt), v2fitarray[pt], liminf[ChosenPart], limsup[ChosenPart], 3);
    Cos2ThetaFitFunction[pt]->SetLineColor(kRed + 1);

    cout << "\n\n Cos2Theta fit function " << endl;
    hCos2Theta[pt]->Fit(Cos2ThetaFitFunction[pt], "R0");
    cout << "Cos2 signal: " << Cos2ThetaFitFunction[pt]->GetParameter(0) << " +- " << Cos2ThetaFitFunction[pt]->GetParError(0) << endl;
    cout << "Cos2 bkg at mass peak: " << Cos2ThetaFitFunction[pt]->GetParameter(1) + Cos2ThetaFitFunction[pt]->GetParameter(2) * mean[pt] << endl;

    Cos2ThetaBkgFunction[pt] = new TF1(Form("cosbkgfunction%i", pt), v2bkgfit, liminf[ChosenPart], limsup[ChosenPart], 2);
    Cos2ThetaBkgFunction[pt]->FixParameter(0, Cos2ThetaFitFunction[pt]->GetParameter(1));
    Cos2ThetaBkgFunction[pt]->FixParameter(1, Cos2ThetaFitFunction[pt]->GetParameter(2));
    Cos2ThetaBkgFunction[pt]->SetLineColor(kBlack);
    Cos2ThetaBkgFunction[pt]->SetLineStyle(8);

    if (isV2FromFit[pt])
    {
      hCos2Theta[pt]->Draw("e");
      Cos2ThetaFitFunction[pt]->Draw("same");
    }
    // if (!isV2FromFit[pt]) hCos2ThetaMassIntegrated[pt]->Draw("");
    hCos2Theta[pt]->Draw("e");
    Cos2ThetaFitFunction[pt]->Draw("same");
    Float_t BinMax = hInvMass[pt]->GetMaximumBin();
    // Float_t BinMax = hCos2Theta[pt]->FindBin(mean[pt]);

    if (pt < numPtBinsVar)
    {
      histoCos2Theta->SetBinContent(pt + 1, Cos2ThetaFitFunction[pt]->GetParameter(0));
      histoCos2Theta->SetBinError(pt + 1, Cos2ThetaFitFunction[pt]->GetParError(0));
      // histoCos2ThetaPeakPos->SetBinContent(pt + 1, hCos2Theta[pt]->GetBinContent(BinMax));
      // histoCos2ThetaPeakPos->SetBinError(pt + 1, hCos2Theta[pt]->GetBinError(BinMax));
      histoCos2ThetaPeakPos->SetBinContent(pt + 1, hCos2Theta[pt]->GetBinContent(hCos2Theta[pt]->FindBin(mean[pt])));
      histoCos2ThetaPeakPos->SetBinError(pt + 1, hCos2Theta[pt]->GetBinError(hCos2Theta[pt]->FindBin(mean[pt])));
      // histoCos2ThetaPeakPos->SetBinContent(pt + 1, hCos2Theta[pt]->GetBinContent(hCos2Theta[pt]->FindBin(1.115)));
      // histoCos2ThetaPeakPos->SetBinError(pt + 1, hCos2Theta[pt]->GetBinError(hCos2Theta[pt]->FindBin(1.115)));
    }
    else
    {
      histoCos2ThetaPtInt->SetBinContent(1, Cos2ThetaFitFunction[pt]->GetParameter(0));
      histoCos2ThetaPtInt->SetBinError(1, Cos2ThetaFitFunction[pt]->GetParError(0));
      // histoCos2ThetaPtIntPeakPos->SetBinContent(1, hCos2Theta[pt]->GetBinContent(BinMax));
      // histoCos2ThetaPtIntPeakPos->SetBinError(1, hCos2Theta[pt]->GetBinError(BinMax));
      histoCos2ThetaPtIntPeakPos->SetBinContent(1, hCos2Theta[pt]->GetBinContent(hCos2Theta[pt]->FindBin(mean[pt])));
      histoCos2ThetaPtIntPeakPos->SetBinError(1, hCos2Theta[pt]->GetBinError(hCos2Theta[pt]->FindBin(mean[pt])));
      // histoCos2ThetaPtIntPeakPos->SetBinContent(1, hCos2Theta[pt]->GetBinContent(hCos2Theta[pt]->FindBin(1.115)));
      // histoCos2ThetaPtIntPeakPos->SetBinError(1, hCos2Theta[pt]->GetBinError(hCos2Theta[pt]->FindBin(1.115)));
    }
    TLegend *legendCos2Theta = new TLegend(0.2, 0.8, 0.4, 0.9);
    legendCos2Theta->SetTextSize(0.055);
    legendCos2Theta->AddEntry("", Form("Mean = %.3f +- %.3f", hCos2ThetaMassIntegrated[pt]->GetMean(), hCos2ThetaMassIntegrated[pt]->GetMeanError()), "");
    legendCos2Theta->Draw("same");

    if (isV2)
      cout << "\nv2 (no fit): " << hV2MassIntegrated[pt]->GetMean() << " +- " << hV2MassIntegrated[pt]->GetMeanError() << endl;
    else
      cout << "\nPzs2 (no fit): " << hV2MassIntegrated[pt]->GetMean() << " +- " << hV2MassIntegrated[pt]->GetMeanError() << endl;
    hV2MassIntegrated[pt]->ResetStats();
    hCos2ThetaMassIntegrated[pt]->ResetStats();

    if (pt < numPtBinsVar)
    {
      histoV2NoFit->SetBinContent(pt + 1, hV2MassIntegrated[pt]->GetMean());
      histoV2NoFit->SetBinError(pt + 1, hV2MassIntegrated[pt]->GetMeanError());
      histoCos2ThetaNoFit->SetBinContent(pt + 1, hCos2ThetaMassIntegrated[pt]->GetMean());
      histoCos2ThetaNoFit->SetBinError(pt + 1, hCos2ThetaMassIntegrated[pt]->GetMeanError());
      for (Int_t eta = 0; eta < numBinsEta; eta++)
      {
        Int_t bin2D = histoCos2ThetaNoFit2D->GetBin(pt + 1, eta + 1);
        histoCos2ThetaNoFit2D->SetBinContent(bin2D, hCos2ThetaMassIntegrated[pt]->GetMean());
        histoCos2ThetaNoFit2D->SetBinError(bin2D, hCos2ThetaMassIntegrated[pt]->GetMeanError());
      }
    }
    else
    {
      histoV2PtIntNoFit->SetBinContent(1, hV2MassIntegrated[pt]->GetMean());
      histoV2PtIntNoFit->SetBinError(1, hV2MassIntegrated[pt]->GetMeanError());
      histoCos2ThetaPtIntNoFit->SetBinContent(1, hCos2ThetaMassIntegrated[pt]->GetMean());
      histoCos2ThetaPtIntNoFit->SetBinError(1, hCos2ThetaMassIntegrated[pt]->GetMeanError());
    }
    if (v2type == 1)
    {
      fitV2SP[pt] = new TF1(Form("fitV2SP%i", pt), "gaus", -5, 5);
      hV2MassIntegrated[pt]->Fit(fitV2SP[pt], "R+");
      // histoV2NoFit->SetBinContent(pt + 1, fitV2SP[pt]->GetParameter(1));
      // histoV2NoFit->SetBinError(pt + 1, fitV2SP[pt]->GetParError(1));
    }
    // cout << histoV2NoFit->GetBinCenter(pt + 1) << " bin c: " << histoV2NoFit->GetBinContent(pt + 1) << endl;
    if (pt < numPtBinsVar)
    {
      if (!isV2FromFit[pt])
      {
        histoV2Mixed->SetBinContent(pt + 1, hV2MassIntegrated[pt]->GetMean());
        histoV2Mixed->SetBinError(pt + 1, hV2MassIntegrated[pt]->GetMeanError());
        histoCos2ThetaMixed->SetBinContent(pt + 1, hCos2ThetaMassIntegrated[pt]->GetMean());
        histoCos2ThetaMixed->SetBinError(pt + 1, hCos2ThetaMassIntegrated[pt]->GetMeanError());
      }
      else
      {
        histoV2Mixed->SetBinContent(pt + 1, v2FitFunction[pt]->GetParameter(0));
        histoV2Mixed->SetBinError(pt + 1, v2FitFunction[pt]->GetParError(0));
        histoCos2ThetaMixed->SetBinContent(pt + 1, Cos2ThetaFitFunction[pt]->GetParameter(0));
        histoCos2ThetaMixed->SetBinError(pt + 1, Cos2ThetaFitFunction[pt]->GetParError(0));
      }
    }
    else
    {
      if (!isV2FromFit[pt])
      {
        histoV2PtIntMixed->SetBinContent(1, hV2MassIntegrated[pt]->GetMean());
        histoV2PtIntMixed->SetBinError(1, hV2MassIntegrated[pt]->GetMeanError());
        histoCos2ThetaPtIntMixed->SetBinContent(1, hCos2ThetaMassIntegrated[pt]->GetMean());
        histoCos2ThetaPtIntMixed->SetBinError(1, hCos2ThetaMassIntegrated[pt]->GetMeanError());
      }
      else
      {
        histoV2PtIntMixed->SetBinContent(1, v2FitFunction[pt]->GetParameter(0));
        histoV2PtIntMixed->SetBinError(1, v2FitFunction[pt]->GetParError(0));
        histoCos2ThetaPtIntMixed->SetBinContent(1, Cos2ThetaFitFunction[pt]->GetParameter(0));
        histoCos2ThetaPtIntMixed->SetBinError(1, Cos2ThetaFitFunction[pt]->GetParError(0));
      }
    }
  }

  TCanvas *canvasMass = new TCanvas("canvasMass", "canvasMass", 800, 1800);
  canvasMass->Divide(2, 3);
  StyleCanvas(canvasMass, 0.15, 0.05, 0.05, 0.15);

  Int_t index = 0;
  for (Int_t pt = 0; pt < numPtBinsVar + 1; pt++)
  {
    if (!isPtAnalysis)
    {
      if (pt == numPtBinsVar)
        continue; // skip the integrated
    }
    if (ParticleType == 0 && pt == 0)
      continue;

    if (pt == 0)
      index = 1;
    else if (pt == 1)
      index = 2;
    else if (pt == 2)
      index = 3;
    else if (pt == 3)
      index = 4;
    else if (pt == 4)
      index = 5;
    else if (pt == 7) // 7
      index = 6;
    else
      continue;

    canvasMass->cd(index);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.18);
    lineP3Sigma[pt] = new TLine(UpLimit[pt], 0, UpLimit[pt], 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
    lineM3Sigma[pt] = new TLine(LowLimit[pt], 0, LowLimit[pt], 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
    hInvMassDraw[pt] = (TH1F *)hInvMass[pt]->Clone(Form("hInvMassDraw%i", pt));
    hInvMassDraw[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[ChosenPart], histoMassRangeUp[ChosenPart]);
    hInvMassDraw[pt]->GetYaxis()->SetRangeUser(1, 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
    hInvMassDraw[pt]->Draw("");
    functions1[pt]->Draw("same");
    functions2[pt]->Draw("same");
    totalSignal[pt]->Draw("same");
    if (BkgType == 0)
    {
      bkg1[pt]->SetLineColor(1);
      bkg1[pt]->SetLineStyle(2);
      bkg1[pt]->Draw("same");
    }
    else if (BkgType == 1)
    {
      bkg2[pt]->SetLineColor(1);
      bkg2[pt]->SetLineStyle(2);
      bkg2[pt]->Draw("same");
    }
    else if (BkgType == 2)
      bkg3[pt]->Draw("same");
    else
      bkg4[pt]->Draw("same");

    lineP3Sigma[pt]->Draw("same");
    lineM3Sigma[pt]->Draw("same");
  }
  TCanvas *canvasSummary = new TCanvas("canvasSummary", "canvasSummary", 1900, 1200);
  canvasSummary->Divide(5, 3);

  TString titleX = titlePt;
  if (!isPtAnalysis)
    titleX = "2(#varphi-#Psi_{EP})";
  canvasSummary->cd(1);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoMean, gaussDisplayRangeLow[ChosenPart], gaussDisplayRangeUp[ChosenPart], 1, 1, titleX, "#mu (GeV/c^{2})", "histoMean", 0, 0, 0, 1.4, 1.4, 1.2);
  histoMean->Draw("");
  histoMeanPtInt->Draw("same");
  canvasSummary->cd(2);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoSigma, 0, 0.010, 1, 1, titleX, "#sigma (GeV/c^{2})", "histoSigma", 0, 0, 0, 1.4, 1.4, 1.2);
  histoSigma->Draw("");
  histoSigmaPtInt->Draw("same");
  canvasSummary->cd(3);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoSigmaWeighted, 0, 0.010, 1, 1, titleX, "#sigma_{w} (GeV/c^{2})", "histoSigmaWeighted", 0, 0, 0, 1.4, 1.4, 1.2);
  histoSigmaWeighted->Draw("");
  histoSigmaPtIntWeighted->Draw("same");
  canvasSummary->cd(4);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoPurity, 0, 1, 1, 1, titleX, "S / (S+B)", "histoPurity", 0, 0, 0, 1.4, 1.4, 1.2);
  histoPurity->Draw("");
  histoPurityPtInt->Draw("same");
  canvasSummary->cd(5);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoYield, 0, 1.2 * histoYield->GetBinContent(histoYield->GetMaximumBin()), 1, 1, titleX, titleYield, "histoYield", 0, 0, 0, 1.4, 1.4, 1.2);
  StyleHisto(histoYieldNN, 0, 1.2 * histoYieldNN->GetBinContent(histoYieldNN->GetMaximumBin()), 1, 1, titleX, titleYieldNN, "histoYieldNotNormByEvts", 0, 0, 0, 1.4, 1.4, 1.2);
  histoYield->Draw("");
  histoYieldPtInt->Draw("same");
  canvasSummary->cd(6);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoTot, 0, 1.2 * histoTot->GetBinContent(histoTot->GetMaximumBin()), 1, 1, titleX, "S+B", "histoTot", 0, 0, 0, 1.4, 1.4, 1.2);
  StyleHisto(histoRelErrYield, 0, 1.2 * histoRelErrYield->GetBinContent(histoRelErrYield->GetMaximumBin()), 1, 1, titleX, "#sigma_{Y}/Y", "histoRelErrorYield", 0, 0, 0, 1.4, 1.4, 1.2);
  // histoTot->Draw("same");
  histoRelErrYield->Draw("same");
  canvasSummary->cd(7);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoB, 0, 1.2 * histoB->GetBinContent(histoB->GetMaximumBin()), 1, 1, titleX, "S+B", "histoB", 0, 0, 0, 1.4, 1.4, 1.2);
  // histoB->Draw("same");
  if (isV2)
  {
    histoV2NoFit->SetTitle("v2 without fit");
    histoV2NoFitErr->SetTitle("Error of v2 without fit");
    histoV2MixedErr->SetTitle("Error of v2 mixed");
    histoV2Err->SetTitle("Error of v2");
  }
  else
  {
    histoV2NoFit->SetTitle("Pz,s2 without fit");
    histoV2NoFitErr->SetTitle("Error of Pz,s2 without fit");
    histoV2MixedErr->SetTitle("Error of Pz,s2 mixed");
    histoV2Err->SetTitle("Error of Pz,s2");
    if (!isPtAnalysis)
    {
      histoV2NoFit->SetTitle("Pz without fit");
      histoV2NoFitErr->SetTitle("Error of Pz without fit");
      histoV2MixedErr->SetTitle("Error of Pz mixed");
      histoV2Err->SetTitle("Error of Pz");
    }
  }
  // subtract baseline
  Float_t baseline = 0;
  Float_t baselineNoFit = 0;
  Float_t baselineMixed = 0;
  TH1F *histobaseline = (TH1F *)histoV2->Clone("histobaseline");
  TH1F *histobaselineNoFit = (TH1F *)histoV2NoFit->Clone("histobaselineNoFit");
  TH1F *histobaselineMixed = (TH1F *)histoV2Mixed->Clone("histobaselineMixed");
  histobaseline->Reset();
  histobaselineNoFit->Reset();
  histobaselineMixed->Reset();
  if (!isPtAnalysis)
  {
    for (Int_t pt = 1; pt <= histoV2->GetNbinsX(); pt++)
    {
      baseline += histoV2->GetBinContent(pt);
      baselineNoFit += histoV2NoFit->GetBinContent(pt);
      baselineMixed += histoV2Mixed->GetBinContent(pt);
    }
    baseline /= histoV2->GetNbinsX();
    baselineNoFit /= histoV2->GetNbinsX();
    baselineMixed /= histoV2->GetNbinsX();
    for (Int_t pt = 1; pt <= histoV2->GetNbinsX(); pt++)
    {
      histobaseline->SetBinContent(pt, baseline);
      histobaselineNoFit->SetBinContent(pt, baselineNoFit);
      histobaselineMixed->SetBinContent(pt, baselineMixed);
    }
    histoV2->Add(histobaseline, -1);
    histoV2NoFit->Add(histobaselineNoFit, -1);
    histoV2Mixed->Add(histobaselineMixed, -1);
  }

  // acceptance correction for polarization
  TFile *fileAcceptance = new TFile(SAcceptanceFile, "READ");
  if (!fileAcceptance)
  {
    cout << "Acceptance file not found" << endl;
    return;
  }
  if (isAcceptanceFromExternalFile)
  {
    TList *dir = (TList *)fileAcceptance->Get("ccdb_object");
    if (!dir)
    {
      cout << "directory not found" << endl;
      return;
    }
    histoCos2Theta = (TH1F *)dir->FindObject("histoCos2ThetaNoFit");
    if (!histoCos2Theta)
    {
      cout << "histogram not found" << endl;
      return;
    }
    histoCos2ThetaPtInt = (TH1F *)dir->FindObject("histoCos2ThetaPtIntNoFit");
    histoCos2ThetaNoFit = (TH1F *)dir->FindObject("histoCos2ThetaNoFit");
    histoCos2ThetaPtIntNoFit = (TH1F *)dir->FindObject("histoCos2ThetaPtIntNoFit");
    histoCos2ThetaMixed = (TH1F *)dir->FindObject("histoCos2ThetaNoFit");
    histoCos2ThetaPtIntMixed = (TH1F *)dir->FindObject("histoCos2ThetaPtIntNoFit");
  }
  if (!isV2 && isApplyAcceptanceCorrection)
  {
    histoV2->Divide(histoCos2Theta);
    histoV2NoFit->Divide(histoCos2ThetaNoFit);
    histoV2Mixed->Divide(histoCos2ThetaMixed);
    histoV2PtInt->Divide(histoCos2ThetaPtInt);
    histoV2PtIntNoFit->Divide(histoCos2ThetaPtIntNoFit);
    histoV2PtIntMixed->Divide(histoCos2ThetaPtIntMixed);
  }

  if (!isV2)
  {
    if (isPolFromLambda)
    {
      // histoV2->Scale(1. / AlphaLambda[ChosenPart] / CXiToLambda);
      // histoV2NoFit->Scale(1. / AlphaLambda[ChosenPart] / CXiToLambda);
      // histoV2Mixed->Scale(1. / AlphaLambda[ChosenPart] / CXiToLambda);
      // histoV2PtInt->Scale(1. / AlphaLambda[ChosenPart] / CXiToLambda);
      // histoV2PtIntNoFit->Scale(1. / AlphaLambda[ChosenPart] / CXiToLambda);
      // histoV2PtIntMixed->Scale(1. / AlphaLambda[ChosenPart] / CXiToLambda);
      // NOTE: Alpha already integrated in ProcessTree.C
      histoV2->Scale(1. / CXiToLambda);
      histoV2NoFit->Scale(1. / CXiToLambda);
      histoV2Mixed->Scale(1. / CXiToLambda);
      histoV2PtInt->Scale(1. / CXiToLambda);
      histoV2PtIntNoFit->Scale(1. / CXiToLambda);
      histoV2PtIntMixed->Scale(1. / CXiToLambda);
    }
    else
    {
      // NOTE: Alpha already integrated in ProcessTree.C
      // histoV2->Scale(1. / AlphaH[ChosenPart]);
      // histoV2NoFit->Scale(1. / AlphaH[ChosenPart]);
      // histoV2Mixed->Scale(1. / AlphaH[ChosenPart]);
      // histoV2PtInt->Scale(1. / AlphaH[ChosenPart]);
      // histoV2PtIntNoFit->Scale(1. / AlphaH[ChosenPart]);
      // histoV2PtIntMixed->Scale(1. / AlphaH[ChosenPart]);
    }
  }

  // scaling by resolution
  if (!ExtrisApplyResoOnTheFly && isV2 == 0) // reso applied on the fly only for polarization
  {
    histoV2->Scale(1. / ftcReso[mul]);
    histoV2NoFit->Scale(1. / ftcReso[mul]);
    histoV2Mixed->Scale(1. / ftcReso[mul]);
    histoV2PtInt->Scale(1. / ftcReso[mul]);
    histoV2PtIntNoFit->Scale(1. / ftcReso[mul]);
    histoV2PtIntMixed->Scale(1. / ftcReso[mul]);
  }
  // Histograms with errors
  for (Int_t pt = 0; pt < histoV2->GetNbinsX(); pt++)
  {
    if (isV2)
      histoV2Err->SetBinContent(pt + 1, histoV2->GetBinError(pt + 1) / histoV2->GetBinContent(pt + 1));
    else
      histoV2Err->SetBinContent(pt + 1, histoV2->GetBinError(pt + 1));
    histoV2Err->SetBinError(pt + 1, 0);
    if (isV2)
      histoV2NoFitErr->SetBinContent(pt + 1, histoV2NoFit->GetBinError(pt + 1) / histoV2NoFit->GetBinContent(pt + 1));
    else
      histoV2NoFitErr->SetBinContent(pt + 1, histoV2NoFit->GetBinError(pt + 1));
    histoV2NoFitErr->SetBinError(pt + 1, 0);
    if (isV2)
      histoV2MixedErr->SetBinContent(pt + 1, histoV2Mixed->GetBinError(pt + 1) / histoV2Mixed->GetBinContent(pt + 1));
    else
      histoV2MixedErr->SetBinContent(pt + 1, histoV2Mixed->GetBinError(pt + 1));
    histoV2MixedErr->SetBinError(pt + 1, 0);
  }
  if (isV2)
    histoV2PtIntErr->SetBinContent(1, histoV2PtInt->GetBinError(1) / histoV2PtInt->GetBinContent(1));
  else
    histoV2PtIntErr->SetBinContent(1, histoV2PtInt->GetBinError(1));
  histoV2PtIntErr->SetBinError(1, 0);
  if (isV2)
    histoV2PtIntNoFitErr->SetBinContent(1, histoV2PtIntNoFit->GetBinError(1) / histoV2PtIntNoFit->GetBinContent(1));
  else
    histoV2PtIntNoFitErr->SetBinContent(1, histoV2PtIntNoFit->GetBinError(1));
  histoV2PtIntNoFitErr->SetBinError(1, 0);
  if (isV2)
    histoV2PtIntMixedErr->SetBinContent(1, histoV2PtIntMixed->GetBinError(1) / histoV2PtIntMixed->GetBinContent(1));
  else
    histoV2PtIntMixedErr->SetBinContent(1, histoV2PtIntMixed->GetBinError(1));
  histoV2PtIntMixedErr->SetBinError(1, 0);
  if (isV2)
  {
    histoV2NoFitErr->SetTitle("Rel. error of v2 without fit");
    histoV2MixedErr->SetTitle("Rel. error of v2 mixed");
    histoV2Err->SetTitle("Rel error of v2");
  }
  else
  {
    histoV2NoFitErr->SetTitle("Error of Pz,s2 without fit");
    histoV2MixedErr->SetTitle("Error of Pz,s2 mixed");
    histoV2Err->SetTitle("Error of Pz,s2");
    if (!isPtAnalysis)
    {
      histoV2NoFitErr->SetTitle("Error of Pz without fit");
      histoV2MixedErr->SetTitle("Error of Pz mixed");
      histoV2Err->SetTitle("Error of Pz");
    }
  }
  histoV2PtIntErr->SetLineColor(kMagenta);
  histoV2PtIntNoFitErr->SetLineColor(kMagenta);
  histoV2PtIntMixedErr->SetLineColor(kMagenta);

  // v2 corrected by weights a posteriori
  histoV2MixedCorr = (TH1F *)histoV2Mixed->Clone("histoV2MixedCorr");
  TFile *fileV2Correction = new TFile("../V2Corr.root", "READ");
  TH1F *histoV2Corr = (TH1F *)fileV2Correction->Get(Form("v2CorrCent%i", mul));
  // this histogram is already corrected by resolution
  if (!histoV2Corr && (mul != numCentMax))
  {
    cout << "Error: histoV2Corr not found" << endl;
    // return;
  }
  if (ExtrisApplyEffWeights && (mul != numCentMax))
    histoV2MixedCorr->Add(histoV2Corr, 1);
  for (Int_t pt = 0; pt < numPtBinsVar; pt++)
  {
    histoV2MixedCorr->SetBinError(pt + 1, histoV2Mixed->GetBinError(pt + 1));
  }

  histoV2NoFit->Draw();
  histoV2PtIntNoFit->Draw("same");

  canvasSummary->cd(8);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  if (isV2)
    histoV2->SetTitle("v2 from fit");
  else
  {
    histoV2->SetTitle("Pz,s2 from fit");
    if (!isPtAnalysis)
      histoV2->SetTitle("Pz from fit");
  }
  histoV2->Draw();
  histoV2PtInt->Draw("same");

  cout << "************ Rel error of histoV2 ********" << endl;
  for (Int_t b = 1; b <= histoV2->GetNbinsX(); b++)
  {
    cout << "Centre of the bin: " << histoV2->GetBinCenter(b) << " rel. error: " << histoV2->GetBinError(b) / histoV2->GetBinContent(b) << endl;
  }

  canvasSummary->cd(9);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoCos2ThetaNoFit, 0, 0.4, 1, 1, titleX, TitleCos2Theta, "histoCos2ThetaNoFit", 0, 0, 0, 1.4, 1.4, 1.2);
  if (!isV2)
  {
    histoCos2ThetaNoFit->SetTitle("cos^{2}(#theta*) without fit");
    histoCos2ThetaNoFit->Draw();
    histoCos2ThetaPtIntNoFit->Draw("same");
  }

  canvasSummary->cd(10);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoCos2Theta, 0, 0.4, 1, 1, titleX, TitleCos2Theta, "histoCos2Theta", 0, 0, 0, 1.4, 1.4, 1.2);
  if (!isV2)
  {
    histoCos2Theta->SetTitle("cos^{2}(#theta*) from fit");
    histoCos2Theta->Draw();
    histoCos2ThetaPtInt->Draw("same");
  }

  canvasSummary->cd(11);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  histoV2PtIntNoFitErr->GetYaxis()->SetRangeUser(0, histoV2Err->GetBinContent(histoV2Err->GetMaximumBin()) * 1.2);
  histoV2PtIntNoFitErr->Draw("");
  histoV2NoFitErr->Draw("same");

  canvasSummary->cd(12);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  histoV2PtIntErr->GetYaxis()->SetRangeUser(0, histoV2Err->GetBinContent(histoV2Err->GetMaximumBin()) * 1.2);
  histoV2PtIntErr->Draw("");
  histoV2Err->Draw("same");

  canvasSummary->cd(13);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  histoYieldFraction->Draw("");
  histoYieldFractionPtInt->Draw("same");

  TString Soutputfile;
  TString SoutputfileAcceptance;
  Soutputfile = "../OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_" + inputFileName + "_" + ParticleName[ChosenPart];
  SoutputfileAcceptance = "../AcceptancePlots/Acceptance_" + inputFileName + "_" + ParticleName[ChosenPart];
  Soutputfile += IsOneOrTwoGauss[UseTwoGauss];
  Soutputfile += SIsBkgParab[BkgType];
  Soutputfile += Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);
  SoutputfileAcceptance += Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);
  if (isApplyWeights)
    Soutputfile += "_Weighted";
  if (isApplyCentWeight)
    Soutputfile += "_CentWeighted";
  if (v2type == 1)
    Soutputfile += "_SP";
  if (!useCommonBDTValue)
    Soutputfile += "_BDTCentDep";
  if (isRun2Binning)
    Soutputfile += "_Run2Binning";
  if (!isPtAnalysis)
  {
    Soutputfile += "_vsPsi";
    SoutputfileAcceptance += "_vsPsi";
  }
  if (!isV2 && isPolFromLambda)
  {
    Soutputfile += "_PolFromLambda";
    SoutputfileAcceptance += "_PolFromLambda";
  }
  if (ExtrisApplyEffWeights)
  {
    Soutputfile += "_EffW";
  }
  if (isSysMultTrial && ChosenPart < 6)
    Soutputfile += SBDT;

  Soutputfile += SEtaSysChoice[EtaSysChoice];
  SoutputfileAcceptance += SEtaSysChoice[EtaSysChoice];
  if (!isRapiditySel)
    Soutputfile += "_Eta08";
  Soutputfile += STHN[ExtrisFromTHN];
  if (useMixedBDTValueInFitMacro)
    Soutputfile += "_MixedBDT";
  if (isProducedAcceptancePlots)
    Soutputfile += "_AcceptancePlots";
  if (isTightMassCut)
  {
    if (ExtrisSysMassCut == 0)
      Soutputfile += Form("_TightMassCut%.1f", sigmacentral);
    else
      Soutputfile += Form("_TightMassCutSyst%i", indexMassCut);
  }
  if (isReducedPtBins)
    Soutputfile += "_ReducedPtBins";
  if (ChosenPart >= 6 && ExtrisSysLambdaMultTrial)
  {
    if (isLoosest)
      Soutputfile += "_isLoosest";
    else if (isTightest)
      Soutputfile += "_isTightest";
    else
      Soutputfile += Form("_SysMultTrial_%i", indexMultTrial);
    Soutputfile += "_isSysLambdaMultTrial";
  }
  if (ExtrisApplyResoOnTheFly)
    Soutputfile += "_ResoOnTheFly";
  // if (ChosenPart >= 6)
  // Soutputfile += "_CorrectReso_TestLeassPtBins";
  if (ChosenPart == 0)
    Soutputfile += "_EPReso";
  if (isBkgPol == 0)
    Soutputfile += "_isBkgPol0";
  if (isTightMassForAcceptancePurity && ChosenPart >= 6)
    Soutputfile += "_isTightMassForAcceptancePurity";
  if (isTighterPzFitRange)
    Soutputfile += "_TighterPzFitRange";
  // Soutputfile += "_SystReso";
  if (ChosenPart >= 6 && !isMassCutForAcceptance && isProducedAcceptancePlots)
  {
    Soutputfile += "_NoMassCutForAcceptance";
    SoutputfileAcceptance += "_NoMassCutForAcceptance";
  }
  // Soutputfile += "_TestMoreBins";

  // save canvases
  canvas[0]->SaveAs(Soutputfile + ".pdf(");
  canvas[1]->SaveAs(Soutputfile + ".pdf");
  canvas[2]->SaveAs(Soutputfile + ".pdf");
  canvas[3]->SaveAs(Soutputfile + ".pdf");
  canvasSummary->SaveAs(Soutputfile + ".pdf)");
  canvasMass->SaveAs(Soutputfile + "_MassPlot.pdf");
  canvasMass->SaveAs(Soutputfile + "_MassPlot.png");

  TFile *outputfile = new TFile(Soutputfile + ".root", "RECREATE");
  for (Int_t i = 0; i < numCanvas; i++)
  {
    outputfile->WriteTObject(canvas[i]);
  }
  for (Int_t pt = 0; pt <= numPtBinsVar; pt++)
  {
    // outputfile->WriteTObject(hInvMass[pt]);
    // outputfile->WriteTObject(hV2[pt]);
    // outputfile->WriteTObject(hCos2Theta[pt]);
  }
  outputfile->WriteTObject(histoAppliedBDT);
  outputfile->WriteTObject(histoAppliedBDTPtInt);
  outputfile->WriteTObject(histoYield);
  outputfile->WriteTObject(histoYieldPtInt);
  outputfile->WriteTObject(histoYieldNN);
  outputfile->WriteTObject(histoTot);
  outputfile->WriteTObject(histoTotPtInt);
  outputfile->WriteTObject(histoB);
  outputfile->WriteTObject(histoBPtInt);
  outputfile->WriteTObject(histoMean);
  outputfile->WriteTObject(histoMeanPtInt);
  outputfile->WriteTObject(histoMean1);
  outputfile->WriteTObject(histoMean1PtInt);
  outputfile->WriteTObject(histoMean2);
  outputfile->WriteTObject(histoMean2PtInt);
  outputfile->WriteTObject(histoSigma);
  outputfile->WriteTObject(histoSigmaPtInt);
  outputfile->WriteTObject(histoSigma1);
  outputfile->WriteTObject(histoSigma1PtInt);
  outputfile->WriteTObject(histoSigma2);
  outputfile->WriteTObject(histoSigma2PtInt);
  outputfile->WriteTObject(histoSigmaWeighted);
  outputfile->WriteTObject(histoSigmaPtIntWeighted);
  outputfile->WriteTObject(histoPurity);
  outputfile->WriteTObject(histoPurityPtInt);
  outputfile->WriteTObject(histoSignificance);
  outputfile->WriteTObject(histoSignificancePtInt);

  outputfile->WriteTObject(histoV2);
  outputfile->WriteTObject(histoV2NoFit);
  outputfile->WriteTObject(histoV2Mixed);
  outputfile->WriteTObject(histoV2MixedCorr);
  outputfile->WriteTObject(histoV2Err);
  outputfile->WriteTObject(histoV2NoFitErr);
  outputfile->WriteTObject(histoV2MixedErr);
  outputfile->WriteTObject(histoV2Bkg);

  outputfile->WriteTObject(histoV2PtInt);
  outputfile->WriteTObject(histoV2PtIntNoFit);
  outputfile->WriteTObject(histoV2PtIntMixed);
  outputfile->WriteTObject(histoV2PtIntErr);
  outputfile->WriteTObject(histoV2PtIntNoFitErr);
  outputfile->WriteTObject(histoV2PtIntMixedErr);
  outputfile->WriteTObject(histoV2BkgPtInt);

  outputfile->WriteTObject(histoCos2Theta);
  outputfile->WriteTObject(histoCos2ThetaNoFit);
  outputfile->WriteTObject(histoCos2ThetaMixed);
  outputfile->WriteTObject(histoCos2ThetaPeakPos);
  outputfile->WriteTObject(histoCos2ThetaPtInt);
  outputfile->WriteTObject(histoCos2ThetaPtIntNoFit);
  outputfile->WriteTObject(histoCos2ThetaPtIntMixed);
  outputfile->WriteTObject(histoCos2ThetaPtIntPeakPos);
  outputfile->Close();

  TFile *outputfile2;

  if (isProducedAcceptancePlots)
  {
    outputfile2 = new TFile(SoutputfileAcceptance + ".root", "RECREATE");
    TList *listAcceptance = new TList();
    listAcceptance->Add(histoCos2ThetaNoFit);
    listAcceptance->Add(histoCos2ThetaNoFit2D);
    listAcceptance->Add(histoCos2ThetaPtIntNoFit);
    listAcceptance->Write("ccdb_object", TObject::kSingleKey);
    outputfile2->Close();
    cout << "I stored the acceptance plots in the file: " << SoutputfileAcceptance << ".root" << endl;
  }
  // Acceptance plots - comparison
  TCanvas *canvasAcc = new TCanvas("canvasAcc", "canvasAcc", 800, 1200);
  Float_t LLUpperPad = 0.44;
  Float_t ULLowerPad = 0.44;
  TPad *padA = new TPad("padA", "padA", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padAL = new TPad("padAL", "padAL", 0, 0, 1, ULLowerPad);

  StylePad(padA, 0.18, 0.01, 0.03, 0.);           // L, R, T, B
  StylePad(padAL, 0.18, 0.01, 0.02, 0.3);         // L, R, T, B
  StyleCanvas(canvasAcc, 0.15, 0.03, 0.02, 0.14); // L, R, T, B
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  histoCos2ThetaNoFit->SetLineColor(kRed);
  histoCos2ThetaNoFit->SetMarkerColor(kRed);
  histoCos2ThetaNoFit->SetMarkerStyle(20);
  histoCos2ThetaPtIntNoFit->SetLineColor(kRed);
  histoCos2ThetaPtIntNoFit->SetMarkerColor(kRed);
  histoCos2ThetaPtIntNoFit->SetMarkerStyle(20);
  histoCos2Theta->SetLineColor(kBlue);
  histoCos2Theta->SetMarkerColor(kBlue);
  histoCos2Theta->SetMarkerStyle(33);
  histoCos2ThetaPtInt->SetLineColor(kBlue);
  histoCos2ThetaPtInt->SetMarkerColor(kBlue);
  histoCos2ThetaPtInt->SetMarkerStyle(33);
  histoCos2ThetaPeakPos->SetLineColor(kGreen + 2);
  histoCos2ThetaPeakPos->SetMarkerColor(kGreen + 2);
  histoCos2ThetaPeakPos->SetMarkerStyle(21);
  histoCos2ThetaPtIntPeakPos->SetLineColor(kGreen + 2);
  histoCos2ThetaPtIntPeakPos->SetMarkerColor(kGreen + 2);
  histoCos2ThetaPtIntPeakPos->SetMarkerStyle(21);
  TLegend *legendAcc = new TLegend(0.4, 0.2, 0.7, 0.5);
  legendAcc->SetFillStyle(0);
  legendAcc->SetTextSize(0.037);
  legendAcc->SetTextAlign(12);
  legendAcc->AddEntry(histoCos2Theta, "cos^{2}_{sig, fit} from fit", "pl");
  legendAcc->AddEntry(histoCos2ThetaNoFit, "cos^{2} assuming P = 1", "pl");
  legendAcc->AddEntry(histoCos2ThetaPeakPos, "cos^{2} at peak position", "pl");
  histoCos2ThetaNoFit->SetTitle("");
  padA->Draw();
  padA->cd();

  if (!isV2)
  {
    histoCos2ThetaNoFit->Draw();
    histoCos2ThetaPtIntNoFit->Draw("same");
    histoCos2Theta->Draw("same");
    histoCos2ThetaPtInt->Draw("same");
    histoCos2ThetaPeakPos->Draw("same");
    histoCos2ThetaPtIntPeakPos->Draw("same");
  }
  legendAcc->Draw();
  canvasAcc->cd();
  padAL->Draw();
  padAL->cd();
  TH1F *hRatio1 = (TH1F *)histoCos2ThetaNoFit->Clone("hRatio1");
  TH1F *hRatio2 = (TH1F *)histoCos2ThetaPeakPos->Clone("hRatio2");
  hRatio1->Divide(histoCos2Theta);
  hRatio2->Divide(histoCos2Theta);
  hRatio1->GetYaxis()->SetRangeUser(0.7, 1.05);
  hRatio1->GetYaxis()->SetTitle("Ratio to cos^{2}_{sig, fit}");
  TF1 *lineat1 = new TF1("lineat1", "1", 0, 8);
  lineat1->SetLineColor(kBlack);
  lineat1->SetLineStyle(2);
  if (!isV2)
  {
    hRatio1->Draw();
    hRatio2->Draw("same");
    lineat1->Draw("same");
  }
  canvasAcc->SaveAs(Soutputfile + "_AccPlot.pdf");
  canvasAcc->SaveAs(Soutputfile + "_AccPlot.png");
  canvasAcc->SaveAs(Form("../AcceptanceComparison_%i-%i.pdf", CentFT0C[mul], CentFT0C[mul + 1]));
  canvasAcc->SaveAs(Form("../AcceptanceComparison_%i-%i.png", CentFT0C[mul], CentFT0C[mul + 1]));

  // Performance plot
  Int_t ChosenPt = 8; // 8
  if (ParticleType == 1 || ParticleType == 2)
    ChosenPt = numPtBinsVar;
  Float_t LowLimitMass[numPart] = {1.29, 1.65, 1.29, 1.29, 1.65, 1.65, 1.1, 1.1, 1.1};
  Float_t UpLimitMass[numPart] = {1.35, 1.7, 1.35, 1.35, 1.7, 1.7, 1.13, 1.13, 1.13};
  Float_t UpperCutHisto = 1.7;
  if (ParticleType == 1)
    UpperCutHisto = 1.8;
  else if (ParticleType == 2)
    UpperCutHisto = 2.1;

  TString TitleXMass = "#it{m}_{#Lambda#pi} (GeV/#it{c}^{2})";
  if (ParticleType == 0)
    TitleXMass = "#it{m}_{#LambdaK} (GeV/#it{c}^{2})";
  else if (ParticleType == 2)
    TitleXMass = "#it{m}_{p#pi} (GeV/#it{c}^{2})";

  TCanvas *canvasMassP = new TCanvas("canvasMassP", "canvasMassP", 800, 800);
  StyleCanvas(canvasMassP, 0.15, 0.03, 0.02, 0.14); // L, R, T, B

  TLegend *legend = new TLegend(0.2, 0.68, 0.71, 0.95);
  legend->SetFillStyle(0);
  legend->SetMargin(0);
  legend->SetTextSize(0.037);
  legend->SetTextAlign(12);
  legend->AddEntry("", "#bf{ALICE Performance}", "");
  legend->AddEntry("", Form("Run 3 Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, %i-%i%s", CentFT0CMin, CentFT0CMax, "%"), "");
  if (ParticleType == 1)
    legend->AddEntry("", "#Xi^{#minus} #rightarrow #Lambda #pi^{#minus} #rightarrow p #pi^{#minus} #pi^{#minus} + c.c.", "");
  else if (ParticleType == 0)
    legend->AddEntry("", "#Omega^{#minus} #rightarrow #Lambda K^{#minus} #rightarrow p #pi^{#minus} K^{#minus} + c.c.", "");
  else if (ParticleType == 2)
    legend->AddEntry("", "#Lambda #rightarrow p #pi^{#minus} + c.c.", "");
  if (ChosenPt == numPtBinsVar)
    legend->AddEntry("", Form("|#it{#eta}| < 0.8, %.1f < #it{p}_{T} < %.1f GeV/#it{c}", PtBins[0], PtBins[numPtBinsVar]), "");
  else
    legend->AddEntry("", Form("|#it{#eta}| < 0.8, %.1f < #it{p}_{T} < %.1f GeV/#it{c}", PtBins[ChosenPt], PtBins[ChosenPt + 1]), "");
  // legend->AddEntry("", Form("BDT, Signif.(4#sigma) = %.0f #pm %.0f", Signif[ChosenPt], errSignif[ChosenPt]), "");
  if (ChosenPart >= 6)
    legend->AddEntry("", Form("S/(S+B)(2#sigma) = %.4f #pm %.4f", SSB[ChosenPt], errSSB[ChosenPt]), "");
  else
    legend->AddEntry("", Form("BDT selected, S/(S+B)(2#sigma) = %.3f #pm %.3f", SSB[ChosenPt], errSSB[ChosenPt]), "");
  // legend->AddEntry("", Form("BDT selected, Purity(2#sigma) = %.3f #pm %.3f", SSB[ChosenPt], errSSB[ChosenPt]), "");

  TLegend *legendfit = new TLegend(0.2, 0.57, 0.71, 0.65);
  legendfit->SetFillStyle(0);
  legendfit->SetMargin(0.1);
  legendfit->SetTextSize(0.031);
  legendfit->SetTextAlign(12);

  TLegend *legendfit2 = new TLegend(0.24, 0.54, 0.71, 0.65);
  legendfit2->SetFillStyle(0);
  legendfit2->SetMargin(0.1);
  legendfit2->SetTextSize(0.033);
  legendfit2->SetTextAlign(12);

  TH1F *histo = hInvMass[ChosenPt];
  Float_t histoIntegral = histo->Integral("width");
  histo->Scale(1. / histoIntegral);
  TString titleyNorm = Form("Normalized counts/(%.1f MeV/#it{c}^{2})", histo->GetBinWidth(1) * 1000);
  if (ParticleType == 1)
    titleyNorm = Form("Normalized counts/(%.1f MeV/#it{c}^{2})", (float)(histo->GetBinWidth(1) * 1000));
  StyleHisto(histo, 0.0001, UpperCutHisto * histo->GetBinContent(histo->GetMaximumBin()), 1, 20,
             TitleXMass, titleyNorm, "", 1, LowLimitMass[ChosenPart] + 0.001, UpLimitMass[ChosenPart] - 0.001, 1.2, 1.8, 1.2);
  histo->GetXaxis()->SetRangeUser(XRangeMin[ChosenPart], XRangeMax[ChosenPart]);
  histo->GetYaxis()->SetRangeUser(0.0001, 299);
  if (ChosenPart >= 6)
    histo->GetYaxis()->SetRangeUser(0.0001, 370);
  histo->GetXaxis()->SetLabelSize(0.043);
  histo->GetXaxis()->SetTitleSize(0.045);
  histo->GetYaxis()->SetLabelSize(0.043);
  histo->GetYaxis()->SetTitleSize(0.045);
  histo->GetYaxis()->SetTitleOffset(1.6);
  histo->DrawClone("pe");
  lineP3SigmaNorm[ChosenPt] = new TLine(UpLimit[ChosenPt], 0, UpLimit[ChosenPt], histo->GetMaximum());
  lineM3SigmaNorm[ChosenPt] = new TLine(LowLimit[ChosenPt], 0, LowLimit[ChosenPt], histo->GetMaximum());
  lineP3SigmaNorm[ChosenPt]->SetLineStyle(2);
  lineM3SigmaNorm[ChosenPt]->SetLineStyle(2);
  // lineP3SigmaNorm[ChosenPt]->Draw("same");
  // lineM3SigmaNorm[ChosenPt]->Draw("same");

  TF1 *totalPNorm = new TF1("totalP", "gaus(0)+gaus(3)+pol2(6)", XRangeMin[ChosenPart], XRangeMax[ChosenPart]);
  totalPNorm->SetParameter(0, total[ChosenPt]->GetParameter(0) / histoIntegral);
  totalPNorm->SetParameter(1, total[ChosenPt]->GetParameter(1));
  totalPNorm->SetParameter(2, total[ChosenPt]->GetParameter(2));
  totalPNorm->SetParameter(3, total[ChosenPt]->GetParameter(3) / histoIntegral);
  totalPNorm->SetParameter(4, total[ChosenPt]->GetParameter(4));
  totalPNorm->SetParameter(5, total[ChosenPt]->GetParameter(5));
  totalPNorm->SetParameter(6, total[ChosenPt]->GetParameter(6) / histoIntegral);
  totalPNorm->SetParameter(7, total[ChosenPt]->GetParameter(7) / histoIntegral);
  totalPNorm->SetParameter(8, total[ChosenPt]->GetParameter(8) / histoIntegral);

  legendfit->AddEntry(totalPNorm, "Gaussian fits + bkg.", "l");
  legendfit2->AddEntry(totalPNorm, "Gaussian fits + bkg.", "l");
  TF1 *bkg;
  if (BkgType == 0)
    bkg = bkg1[ChosenPt];
  else if (BkgType == 1)
    bkg = bkg2[ChosenPt];
  else if (BkgType == 2)
    bkg = bkg3[ChosenPt];
  else
    bkg = bkg4[ChosenPt];
  bkg->SetParameter(0, bkg->GetParameter(0) / histoIntegral);
  bkg->SetParameter(1, bkg->GetParameter(1) / histoIntegral);
  if (BkgType == 1)
    bkg->SetParameter(2, bkg->GetParameter(2) / histoIntegral);
  else if (BkgType == 2)
  {
    bkg->SetParameter(2, bkg->GetParameter(2) / histoIntegral);
    bkg->SetParameter(3, bkg->GetParameter(3) / histoIntegral);
  }
  bkg->SetLineColor(kBlack);
  bkg->SetLineStyle(8);
  bkg->Draw("same");
  legendfit->AddEntry(bkg, "bkg.", "l");
  legendfit2->AddEntry(bkg, "bkg.", "l");
  if (ExtrisApplyResoOnTheFly)
    legendfit2->AddEntry("", Form("Pz = %.5f + %.5f", v2FitFunction[ChosenPt]->GetParameter(0), v2FitFunction[ChosenPt]->GetParError(0)));
  else
    legendfit2->AddEntry("", Form("Pz (no reso) = %.5f + %.5f", v2FitFunction[ChosenPt]->GetParameter(0), v2FitFunction[ChosenPt]->GetParError(0)));
  totalPNorm->SetRange(LowLimitMass[ChosenPart], UpLimitMass[ChosenPart]);
  totalPNorm->SetLineColor(kRed + 1);
  totalPNorm->Draw("same");
  legend->Draw("");
  legendfit->Draw("");
  canvasMassP->SaveAs("../PerformancePlots/MassFit" + ParticleName[ChosenPart] + Form("_Cent%i-%i_Pt%i.pdf", CentFT0CMin, CentFT0CMax, ChosenPt));
  canvasMassP->SaveAs("../PerformancePlots/MassFit" + ParticleName[ChosenPart] + Form("_Cent%i-%i_Pt%i.png", CentFT0CMin, CentFT0CMax, ChosenPt));
  canvasMassP->SaveAs("../PerformancePlots/MassFit" + ParticleName[ChosenPart] + Form("_Cent%i-%i_Pt%i.eps", CentFT0CMin, CentFT0CMax, ChosenPt));

  TCanvas *canvasCos2P = new TCanvas("canvasCos2P", "canvasCos2P", 800, 800);
  StyleCanvas(canvasCos2P, 0.15, 0.03, 0.02, 0.14); // L, R, T, B
  TH1F *histoCos2 = hCos2ThetaMassIntegrated[ChosenPt];
  histoCos2->Scale(1. / histoCos2->Integral(""));
  TString titleyNormCos2 = "Normalized counts";
  TString TitleCos2 = "cos^{2}(#theta_{#Lambda}*)";
  if (isPolFromLambda)
    TitleCos2 = "cos^{2}(#theta_{p}*)";

  StyleHisto(histoCos2, 0.0001, 1.2 * histoCos2->GetBinContent(histoCos2->GetMaximumBin()), 1, 20,
             TitleCos2, titleyNormCos2, "", 1, 0, 1, 1.2, 1.4, 1.2);
  histoCos2->Draw("pe");
  TLegend *legendCos2P = new TLegend(0.3, 0.85, 0.5, 0.95);
  legendCos2P->SetTextSize(0.035);
  legendCos2P->AddEntry("", Form("Mean = %.3f", histoCos2->GetMean()), "");
  legendCos2P->Draw("same");
  canvasCos2P->SaveAs("../PerformancePlots/Cos2Theta" + ParticleName[ChosenPart] + Form("_Cent%i-%i_Pt%i.pdf", CentFT0CMin, CentFT0CMax, ChosenPt));
  canvasCos2P->SaveAs("../PerformancePlots/Cos2Theta" + ParticleName[ChosenPart] + Form("_Cent%i-%i_Pt%i.png", CentFT0CMin, CentFT0CMax, ChosenPt));

  TCanvas *canvasCosSinP = new TCanvas("canvasCosSinP", "canvasCosSinP", 800, 800);
  StyleCanvas(canvasCosSinP, 0.2, 0.03, 0.02, 0.14); // L, R, T, B
  TH1F *histoCosSin = (TH1F *)hV2MassIntegrated[ChosenPt]->Clone("histoCosSin");
  histoCosSin->Scale(1. / histoCosSin->Integral(""));
  TString titleyNormCosSin = "Normalized counts";
  TString TitleCosSin = "1/#alpha_{#Xi} cos(#theta_{#Lambda}*) sin(2(#varphi_{#Xi}-#Psi_{2}))";
  if (!isApplyAcceptanceCorrection)
  {
    TitleCosSin = "1/#LTcos^{2}(#theta_{#Lambda}*)#GT 1/#alpha_{#Xi} cos(#theta_{#Lambda}*) sin(2(#varphi_{#Xi}-#Psi_{2}))";
  }
  if (isPolFromLambda)
  {
    TitleCosSin = "1/#alpha_{#Lambda} cos(#theta_{p}*) sin(2(#varphi_{#Xi}-#Psi_{2}))";
    if (!isApplyAcceptanceCorrection)
    {
      TitleCosSin = "1/#LTcos^{2}(#theta_{p}*)#GT 1/#alpha_{#Lambda} cos(#theta_{p}*) sin(2(#varphi_{#Xi}-#Psi_{2}))";
    }
  }
  StyleHisto(histoCosSin, 0.0001, 1.2 * histoCosSin->GetBinContent(histoCosSin->GetMaximumBin()), 1, 20,
             TitleCosSin, titleyNormCosSin, "", 1, -15, 15, 1.2, 1.8, 1.2);
  histoCosSin->Draw("pe");
  TLegend *legendCosSin = new TLegend(0.3, 0.85, 0.5, 0.95);
  legendCosSin->SetTextSize(0.035);
  legendCosSin->AddEntry("", Form("Mean = %.5f +- %.5f", histoCosSin->GetMean(), histoCosSin->GetMeanError()), "");
  legendCosSin->Draw("same");
  canvasCosSinP->SaveAs("../PerformancePlots/CosSinTheta" + ParticleName[ChosenPart] + Form("_Cent%i-%i_Pt%i.pdf", CentFT0CMin, CentFT0CMax, ChosenPt));
  canvasCosSinP->SaveAs("../PerformancePlots/CosSinTheta" + ParticleName[ChosenPart] + Form("_Cent%i-%i_Pt%i.png", CentFT0CMin, CentFT0CMax, ChosenPt));

  TCanvas *canvasP = new TCanvas("canvasP", "canvasP", 800, 1100);
  TCanvas *canvasCos2 = new TCanvas("canvasCos2", "canvasCos2", 800, 1100);
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);
  TPad *pad2 = new TPad("pad2", "pad2", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL2 = new TPad("padL2", "padL2", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B
  StylePad(pad2, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL2, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

  TLegend *LegendTitle = new TLegend(0.24, 0.65, 0.75, 0.95);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetMargin(0);
  LegendTitle->SetTextSize(0.038);
  LegendTitle->SetTextAlign(12);
  LegendTitle->AddEntry("", "#bf{ALICE Performance}", "");
  if (isOOCentrality)
    LegendTitle->AddEntry("", Form("Run 3 OO #sqrt{#it{s}_{NN}} = 5.36 TeV, %i-%i%s", CentFT0CMin, CentFT0CMax, "%"), "");
  else
    LegendTitle->AddEntry("", Form("Run 3 Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, %i-%i%s", CentFT0CMin, CentFT0CMax, "%"), "");
  if (ParticleType == 1)
    LegendTitle->AddEntry("", "#Xi^{#minus} #rightarrow #Lambda #pi^{#minus} #rightarrow p #pi^{#minus} #pi^{#minus} + c.c.", "");
  else if (ParticleType == 0)
    LegendTitle->AddEntry("", "#Omega^{#minus} #rightarrow #Lambda K^{#minus} #rightarrow p #pi^{#minus} K^{#minus} + c.c.", "");
  else if (ParticleType == 2)
    LegendTitle->AddEntry("", "#Lambda #rightarrow p #pi^{#minus} + c.c.", "");

  if (ChosenPt == numPtBinsVar)
    LegendTitle->AddEntry("", Form("|#it{#eta}| < 0.8, %.1f < #it{p}_{T} < %.1f GeV/#it{c}", PtBins[0], PtBins[numPtBinsVar]), "");
  else
    LegendTitle->AddEntry("", Form("|#it{#eta}| < 0.8, %.1f < #it{p}_{T} < %.1f GeV/#it{c}", PtBins[ChosenPt], PtBins[ChosenPt + 1]), "");
  if (ChosenPart >= 6)
    LegendTitle->AddEntry("", Form("Signif.(2#sigma) = %.0f #pm %.0f", Signif[ChosenPt], errSignif[ChosenPt]), "");
  else
    LegendTitle->AddEntry("", Form("BDT, Signif.(2#sigma) = %.0f #pm %.0f", Signif[ChosenPt], errSignif[ChosenPt]), "");

  TLegend *legendCos2 = new TLegend(0.6, 0.37, 0.85, 0.52);
  legendCos2->SetFillStyle(0);
  legendCos2->SetMargin(0);
  legendCos2->SetTextSize(0.05);
  legendCos2->SetTextAlign(12);
  Float_t Cos2Signal = 0;
  Float_t ErrCos2Signal = 0;
  if (ChosenPt == numPtBinsVar)
  {
    Cos2Signal = histoCos2ThetaPtInt->GetBinContent(1);
    ErrCos2Signal = histoCos2ThetaPtInt->GetBinError(1);
  }
  else
  {
    Cos2Signal = histoCos2Theta->GetBinContent(ChosenPt + 1);
    ErrCos2Signal = histoCos2Theta->GetBinError(ChosenPt + 1);
  }
  legendCos2->AddEntry("", TitleCos2Theta_sig + Form(" = %.3f", Cos2Signal), "");

  Float_t LimSupSpectra = 9.99;
  Float_t LimInfSpectra = 0.2 * 1e-5;
  Float_t xTitle = 15;
  Float_t xOffset = 4;
  Float_t yTitle = 30;
  Float_t yOffset = 2.4;

  Float_t xLabel = 30;
  Float_t yLabel = 30;
  Float_t xLabelOffset = 0.05;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.03;
  Float_t tickY = 0.042;

  TLegend *legendChi2 = new TLegend(0.4, 0.82, 0.65, 0.92);
  legendChi2->SetFillStyle(0);
  legendChi2->SetMargin(0);
  legendChi2->SetTextSize(0.05);
  legendChi2->SetTextAlign(12);
  legendChi2->AddEntry("", Form("#chi^{2}/NDF = %.2f/%i", v2FitFunction[ChosenPt]->GetChisquare(), v2FitFunction[ChosenPt]->GetNDF()), "");

  TLegend *legendMassChi2 = new TLegend(0.4, 0.82, 0.65, 0.92);
  legendMassChi2->SetFillStyle(0);
  legendMassChi2->SetMargin(0);
  legendMassChi2->SetTextSize(0.05);
  legendMassChi2->SetTextAlign(12);
  legendMassChi2->AddEntry("", Form("#chi^{2}/NDF = %.2f/%i", total[ChosenPt]->GetChisquare(), total[ChosenPt]->GetNDF()), "");

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, hInvMass[ChosenPt]->GetXaxis()->GetXmin(), hInvMass[ChosenPt]->GetXaxis()->GetXmax());
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasP->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, 1e-3, 1.2 * hInvMass[ChosenPt]->GetMaximum(), 1, 1, TitleXMass, titleyNorm, "", 1, 1.15, 1.6);
  if (ChosenPart >= 6)
    StyleHistoYield(hDummy, 1e-3, hInvMass[ChosenPt]->GetMaximum(), 1, 1, TitleXMass, titleyNorm, "", 1, 1.15, 1.8);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(XRangeMin[ChosenPart], XRangeMax[ChosenPart]);
  pad1->Draw();
  pad1->cd();
  // gPad->SetLogy();
  hDummy->Draw("same");
  hInvMass[ChosenPt]->Draw("hist same pe");
  totalPNorm->Draw("same");
  LegendTitle->Draw("");
  legendfit2->Draw("");
  bkg->Draw("same");
  // legendMassChi2->Draw("");

  canvasCos2->cd();
  pad2->Draw();
  pad2->cd();
  hDummy->Draw("same");
  hInvMass[ChosenPt]->Draw("hist same pe");
  totalPNorm->Draw("same");
  LegendTitle->Draw("");
  legendfit2->Draw("");
  bkg->Draw("same");

  Float_t LimSupMultRatio = 5.1;
  Float_t LimInfMultRatio = 1e-2;
  Float_t YoffsetSpectraRatio = 1.1;
  Float_t xTitleR = 30;
  Float_t xOffsetR = 1.5;
  Float_t yTitleR = 30;
  Float_t yOffsetR = 2.4;

  Float_t xLabelR = 30;
  Float_t yLabelR = 30;
  Float_t xLabelOffsetR = 0.02;
  Float_t yLabelOffsetR = 0.014;

  Float_t tickXR = 0.035;
  Float_t tickYR = 0.042;

  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, XRangeMin[ChosenPart], XRangeMax[ChosenPart]);
  for (Int_t i = 1; i <= hDummyRatio->GetNbinsX(); i++)
    hDummyRatio->SetBinContent(i, 1e-12);

  TString TitleDummyRatio = "v_{2}";
  if (!isV2)
  {
    // TitleDummyRatio = "P_{z,s2}";
    TitleDummyRatio = "#LT 1/#alpha_{#Xi} cos(#theta_{#Lambda}*) sin(2(#varphi_{#Xi}-#Psi_{2})) #GT";
    if (ChosenPart >= 6)
      TitleDummyRatio = "#LT 1/#alpha_{#Lambda} cos(#theta_{p}*) sin(2(#varphi_{#Lambda}-#Psi_{2})) #GT";
    if (isPolFromLambda)
      TitleDummyRatio = "#LT 1/#alpha_{#Lambda} cos(#theta_{p}*) sin(2(#varphi_{#Xi}-#Psi_{2})) #GT";
  }
  StyleHistoYield(hDummyRatio, 1e-5, 0.15 - 1e-5, 1, 1, TitleXMass, TitleDummyRatio, "", 1, 1.15, YoffsetSpectraRatio);
  hDummyRatio->GetXaxis()->SetRangeUser(XRangeMin[ChosenPart], XRangeMax[ChosenPart]);
  SetFont(hDummyRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetHistoTextSize(hV2[ChosenPt], xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  hDummyRatio->GetYaxis()->SetNdivisions(704);
  SetTickLength(hDummyRatio, tickXR, tickYR);

  canvasP->cd();
  padL1->Draw();
  padL1->cd();
  if (!isV2)
  {
    hDummyRatio->GetYaxis()->SetRangeUser(-0.015, 0.015);
    if (mul > 3)
      hDummyRatio->GetYaxis()->SetRangeUser(-0.02, 0.02);
    if (mul > 6)
      hDummyRatio->GetYaxis()->SetRangeUser(-0.04, 0.04);
    if (mul > 7)
      hDummyRatio->GetYaxis()->SetRangeUser(-0.005, 0.005);
    if (ChosenPart >= 6)
    {
      hDummyRatio->GetYaxis()->SetRangeUser(-0.02, 0.02);
      if (mul > 0)
        hDummyRatio->GetYaxis()->SetRangeUser(-0.03, 0.03);
      if (mul > 2)
        hDummyRatio->GetYaxis()->SetRangeUser(-0.05, 0.05);
      if (mul > 4)
        hDummyRatio->GetYaxis()->SetRangeUser(-0.08, 0.08);
      if (mul == 6)
        hDummyRatio->GetYaxis()->SetRangeUser(-0.1, 0.1);
      if (mul == 7)
        hDummyRatio->GetYaxis()->SetRangeUser(-0.2, 0.2);
      if (mul == 8)
        hDummyRatio->GetYaxis()->SetRangeUser(-0.5, 0.5);
    }
  }
  hDummyRatio->Draw("same");
  hV2[ChosenPt]->GetXaxis()->SetRangeUser(XRangeMin[ChosenPart], XRangeMax[ChosenPart]);
  hV2[ChosenPt]->GetYaxis()->SetRangeUser(0, 0.15);
  hV2[ChosenPt]->SetTitle("");
  hV2[ChosenPt]->Draw("same e");
  v2FitFunction[ChosenPt]->Draw("same");
  v2BkgFunction[ChosenPt]->Draw("same");
  //  hV2MassIntegrated[ChosenPt]->Draw("");
  legendChi2->Draw("");
  TString SIsPolFromLambda[2] = {"", "_isPolFromLambda"};
  canvasP->SaveAs("../PerformancePlots/MassAnd" + NameAnalysis[!isV2] + ParticleName[ChosenPart] + SIsPolFromLambda[isPolFromLambda] + Form("_Cent%i-%i_Pt%i.pdf", CentFT0CMin, CentFT0CMax, ChosenPt));
  canvasP->SaveAs("../PerformancePlots/MassAnd" + NameAnalysis[!isV2] + ParticleName[ChosenPart] + SIsPolFromLambda[isPolFromLambda] + Form("_Cent%i-%i_Pt%i.png", CentFT0CMin, CentFT0CMax, ChosenPt));
  canvasP->SaveAs("../PerformancePlots/MassAnd" + NameAnalysis[!isV2] + ParticleName[ChosenPart] + SIsPolFromLambda[isPolFromLambda] + Form("_Cent%i-%i_Pt%i.eps", CentFT0CMin, CentFT0CMax, ChosenPt));

  canvasCos2->cd();
  padL2->Draw();
  padL2->cd();
  TH1F *hDummyRatioClone = (TH1F *)hDummyRatio->Clone("hDummyRatioClone");
  hDummyRatioClone->GetYaxis()->SetRangeUser(0, 0.5);
  hDummyRatioClone->GetYaxis()->SetTitle(TitleCos2Theta);
  hDummyRatioClone->Draw("same");
  hCos2Theta[ChosenPt]->GetXaxis()->SetRangeUser(XRangeMin[ChosenPart], XRangeMax[ChosenPart]);
  hCos2Theta[ChosenPt]->SetTitle("");
  hCos2Theta[ChosenPt]->Draw("same e");
  Cos2ThetaFitFunction[ChosenPt]->Draw("same");
  Cos2ThetaBkgFunction[ChosenPt]->Draw("same");
  legendCos2->Draw("");
  //  hV2MassIntegrated[ChosenPt]->Draw("");
  canvasCos2->SaveAs("../PerformancePlots/MassAndCosTheta" + ParticleName[ChosenPart] + SIsPolFromLambda[isPolFromLambda] + Form("_Cent%i-%i_Pt%i.pdf", CentFT0CMin, CentFT0CMax, ChosenPt));
  canvasCos2->SaveAs("../PerformancePlots/MassAndCosTheta" + ParticleName[ChosenPart] + SIsPolFromLambda[isPolFromLambda] + Form("_Cent%i-%i_Pt%i.png", CentFT0CMin, CentFT0CMax, ChosenPt));
  canvasCos2->SaveAs("../PerformancePlots/MassAndCosTheta" + ParticleName[ChosenPart] + SIsPolFromLambda[isPolFromLambda] + Form("_Cent%i-%i_Pt%i.eps", CentFT0CMin, CentFT0CMax, ChosenPt));

  // cout << "\nA partire dal file:\n"
  //      << SPathIn << endl;
  cout << "\nHo creato il file: " << Soutputfile << ".root" << endl;

  cout << "\n\nFor the result intergated in pt these are important info:" << endl;
  cout << "BDT score > " << BDTscoreCutPtInt_checkValue << endl;
  cout << "Input file: " << SPathInPtInt << endl;
  cout << "Is the final V2 from fit? " << isV2FromFit[numPtBinsVar] << endl;
  if (ExtrisApplyResoOnTheFly)
    cout << "The resolution was applied on the fly" << endl;
  else
    cout << "The resolution is: " << ftcReso[mul] << endl;
  cout << "The acceptance correction was applied? " << isApplyAcceptanceCorrection << endl;
  cout << "The purity of the pt integrated sample is: " << histoPurityPtInt->GetBinContent(1) << endl;
  cout << "\nResult (pt integrated measurement, no fit): " << histoV2PtIntNoFit->GetBinContent(1) << endl;
  cout << "Result (pt integrated measurement, fit): " << histoV2PtInt->GetBinContent(1) << endl;
  // cout << "Result before application of the resolution: " << hV2MassIntegrated[numPtBins]->GetMean() << " +- " << hV2MassIntegrated[numPtBins]->GetMeanError() << endl;
  cout << "Result before application of the resolution: " << hV2MassIntegrated[ChosenPt]->GetMean() << " +- " << hV2MassIntegrated[ChosenPt]->GetMeanError() << endl;
  cout << "\nError of pt integrated measurement (no fit): " << histoV2PtIntNoFitErr->GetBinContent(1) << endl;
  cout << "Error of pt integrated measurement (fit): " << histoV2PtIntErr->GetBinContent(1) << endl;
  cout << "\nPz, bkg in correspondence of mass peak " << histoV2BkgPtInt->GetBinContent(1) << " +- " << histoV2BkgPtInt->GetBinError(1) << " nsigma from zero = " << abs(histoV2BkgPtInt->GetBinContent(1)) / histoV2BkgPtInt->GetBinError(1) << endl;

  if (!ExtrisSysMassCut)
    cout << "Purity, significance and yields computed in mass interval of: " << sigmacentral << " sigmas " << endl;
  cout << "This interval is: " << LowLimit[ChosenPt] << " - " << UpLimit[ChosenPt] << " GeV/c^2" << endl;
  cout << "In this interval, the integral of the signal function is:" << endl;
  cout << histoYieldFractionPtInt->GetBinContent(1) << " of the total integral" << endl;

  if (isProducedAcceptancePlots)
  {
    cout << "I stored the acceptance plots in the file: " << SoutputfileAcceptance << ".root" << endl;
  }
  if (isAcceptanceFromExternalFile && isApplyAcceptanceCorrection)
  {
    cout << "I took the acceptance from an external file: " << SAcceptanceFile << endl;
    cout << "The acceptance is computed without invariant mass fit " << endl;
    cout << "The acceptance value is : " << histoCos2ThetaPtIntNoFit->GetBinContent(1) << endl;
  }
  cout << "\nSignificance of the Pz,s2 measurement: " << histoV2PtInt->GetBinContent(1) / histoV2PtIntErr->GetBinContent(1) << endl;
}
