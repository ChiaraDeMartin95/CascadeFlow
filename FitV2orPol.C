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
#include "CommonVar.h"

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

const Float_t UpperLimitLSBOmega = 1.655; // upper limit of fit of left sidebands for omega
const Float_t LowerLimitRSBOmega = 1.689; // lower limit of fit of right sidebands for omega
const Float_t UpperLimitLSBXi = 1.302;    // upper limit of fit of left sidebands for Xi
const Float_t LowerLimitRSBXi = 1.34;     // lower limit of fit of right sidebands for Xi

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
    bkgfraction = TH1D(Form("fraction%s_%s", bkg->GetName(), total->GetName()), "", 1600, min, max);
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
  if (reject && x[0] > LimInf && x[0] < LimSup)
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0];
}

TString titleYield = "dN/dp_{T}";
TString titlePt = "p_{T} (GeV/c)";
TString TitleInvMass[numPart] = {"(#Lambda, #pi)", "(#Lambda, K)"};
TString SInvMass = "invariant mass (GeV/c^{2})";

// fit ranges
Float_t min_range_signal[numPart] = {1.3, 1.65}; // gauss fit range
Float_t max_range_signal[numPart] = {1.335, 1.69};
Float_t liminf[numPart] = {1.29, 1.63}; // bkg and total fit range
Float_t limsup[numPart] = {1.352, 1.71};

// visualisation ranges
Float_t LowMassRange[numPart] = {1.31, 1.655}; // range to compute approximate yield (signal + bkg)
Float_t UpMassRange[numPart] = {1.33, 1.685};
Float_t gaussDisplayRangeLow[numPart] = {1.29, 1.63}; // display range of gauss functions (from total fit)
Float_t gaussDisplayRangeUp[numPart] = {1.35, 1.71};
Float_t bkgDisplayRangeLow[numPart] = {1.29, 1.626}; // display range of bkg function (from total fit)
Float_t bkgDisplayRangeUp[numPart] = {1.35, 1.72};
Float_t histoMassRangeLow[numPart] = {1.29, 1.626}; // display range of mass histograms
Float_t histoMassRangeUp[numPart] = {1.35, 1.72};

// Event plane resolution
Float_t ftcReso[numCent + 1] = {0};

void FitV2orPol(
    Int_t indexMultTrial = 0,
    Int_t mul = 0,
    Bool_t isXi = ChosenParticleXi,
    TString inputFileName = SinputFileName,
    Int_t EtaSysChoice = ExtrEtaSysChoice,
    Int_t BkgType = ExtrBkgType,
    Bool_t isLogy = 1,
    Int_t part = ExtrParticle,
    Bool_t isYAxisMassZoomed = 0,
    Bool_t UseTwoGauss = ExtrUseTwoGauss,
    Bool_t isMeanFixedPDG = 0,
    Float_t sigmacentral = 4.2,
    Bool_t isSysMultTrial = ExtrisSysMultTrial)
{
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

  cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << endl;
  TString histoNameMassvsV2 = "";
  TString ProfileV2 = "";

  if (isSysMultTrial)
    inputFileName = SinputFileNameSyst;
  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trialsBDT)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;
  TString SBDT = "";
  if (BDTscoreCut != DefaultBDTscoreCut)
    SBDT = Form("_BDT%.3f", BDTscoreCut);

  Int_t NEvents = 0;
  TString PathInEvents = "TreeForAnalysis/AnalysisResults_" + inputFileName + ".root";
  TFile *fileEvt = new TFile(PathInEvents, "");
  if (!fileEvt)
  {
    cout << "File Evt does not exist" << endl;
    return;
  }
  TDirectoryFile *dirEvt = (TDirectoryFile *)fileEvt->Get("lf-cascade-flow");
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
      fileResoName = ResoFileName_SPCFW;
    else
      fileResoName = ResoFileName_SPLF;
  }
  else
  {
    if (SinputFileName.Index("CFW") != -1)
      fileResoName = ResoFileName_EPCFW;
    else
      fileResoName = ResoFileName_EPLF;
  }
  fileResoName += ".root";
  TFile *fileResoEP = new TFile(fileResoName, "");
  TH1F *hReso = (TH1F *)fileResoEP->Get("hReso");
  TH1F *hReso080 = (TH1F *)fileResoEP->Get("hReso080");
  cout << "Reso name: " << fileResoName << endl;
  if (mul == numCent)
    ftcReso[mul] = hReso080->GetBinContent(1);
  else
    ftcReso[mul] = hReso->GetBinContent(hReso->FindBin(CentFT0C[mul] + 0.001));
  cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << endl;
  cout << "Resolution: " << ftcReso[mul] << endl;

  Float_t UpperLimitLSB = 0;
  Float_t LowerLimitRSB = 0;
  if (part == 0)
  {
    UpperLimitLSB = UpperLimitLSBXi;
    LowerLimitRSB = LowerLimitRSBXi;
  }
  else if (part == 1)
  {
    UpperLimitLSB = UpperLimitLSBOmega;
    LowerLimitRSB = LowerLimitRSBOmega;
  }

  if (mul > numCent + 1)
  {
    cout << "Multiplicity out of range" << endl;
    return;
  }

  TString SPathIn = "OutputAnalysis/V2_" + inputFileName + "_" + ParticleName[!isXi] + ChargeName[ExtrCharge + 1] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (isApplyWeights)
    SPathIn += "_Weighted";
  if (v2type == 1)
    SPathIn += "_SP";
  if (!useCommonBDTValue)
    SPathIn += "_BDTCentDep";
  if (isRun2Binning)
    SPathIn += "_Run2Binning";
  SPathIn += ".root";

  TFile *filein = new TFile(SPathIn, "");
  if (!filein)
  {
    cout << "FileIn not available" << endl;
    return;
  }

  TString SPt[numPtBins] = {""};
  TH1F *hInvMass[numPtBins];
  TH1F *hInvMassDraw[numPtBins];
  TH1F *hV2[numPtBins];
  TH2F *hmassVsV2C[numPtBins];
  TH1F *hV2MassIntegrated[numPtBins];
  TF1 *fitV2SP[numPtBins];

  Int_t numCanvas = 4;
  TCanvas *canvas[numCanvas];
  for (Int_t c = 0; c < numCanvas; c++)
  {
    canvas[c] = new TCanvas(Form("canvas_%i", c), Form("canvas%i", c), 1800, 1400);
    canvas[c]->Divide(4, 2);
    StyleCanvas(canvas[c], 0.15, 0.05, 0.05, 0.15);
  }

  TH1F *histoMean = new TH1F("histoMean", "histoMean", numPtBins, PtBins);
  TH1F *histoSigma = new TH1F("histoSigma", "histoSigma", numPtBins, PtBins);
  TH1F *histoSigmaWeighted = new TH1F("histoSigmaWeighted", "histoSigmaWeighted", numPtBins, PtBins);
  TH1F *histoPurity = new TH1F("histoPurity", "histoPurity", numPtBins, PtBins);
  TH1F *histoSignificance = new TH1F("histoSignificance", "histoSignificance", numPtBins, PtBins);
  TH1F *histoYield = new TH1F("histoYield", "histoYield", numPtBins, PtBins);
  TH1F *histoRelErrYield = new TH1F("histoRelErrYield", "histoRelErrYield", numPtBins, PtBins);
  TH1F *histoTot = new TH1F("histoTot", "histoTot", numPtBins, PtBins);
  TH1F *histoB = new TH1F("histoB", "histoB", numPtBins, PtBins);
  TString ShistoV2 = "histoV2";
  if (!isV2) // polarization
    ShistoV2 = "histoPzs2";
  TH1F *histoV2 = new TH1F(ShistoV2, ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", numPtBins, PtBins);
  TH1F *histoV2NoFit = new TH1F(ShistoV2 + "NoFit", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", numPtBins, PtBins);
  TH1F *histoV2Mixed = new TH1F(ShistoV2 + "Mixed", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", numPtBins, PtBins);

  for (Int_t pt = 0; pt < numPtBins; pt++)
  {
    if (!isXi && pt == 0)
      continue;
    SPt[pt] = Form("%.2f < p_{T} < %.2f", PtBins[pt], PtBins[pt + 1]);
    cout << "Analysed pt interval: " << PtBins[pt] << "-" << PtBins[pt + 1] << endl;
    cout << PtBins[pt] << endl;

    hInvMass[pt] = (TH1F *)filein->Get(Form("mass_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt));
    if (!hInvMass[pt])
    {
      cout << "Histogram inv. mass not available" << endl;
      return;
    }

    StyleHisto(hInvMass[pt], 0, 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), 1, 20, TitleInvMass[part] + " " + SInvMass, "Counts", SPt[pt] + " GeV/#it{c}", 1, histoMassRangeLow[part], histoMassRangeUp[part], 1.4, 1.6, 0.7);

    if (isV2)
    {
      histoNameMassvsV2 = Form("MassvsV2C_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      ProfileV2 = Form("V2C_cent%i-%i_pt%i_Profile", CentFT0CMin, CentFT0CMax, pt);
    }
    else
    { // polarization
      histoNameMassvsV2 = Form("MassvsPzs2_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      ProfileV2 = Form("Pzs2_cent%i-%i_pt%i_Profile", CentFT0CMin, CentFT0CMax, pt);
    }

    hmassVsV2C[pt] = (TH2F *)filein->Get(histoNameMassvsV2);
    if (!hmassVsV2C[pt])
    {
      cout << "Histogram hmassVsV2C not available" << endl;
      return;
    }

    // hV2[pt] = (TH1F *)filein->Get(Form("V2C_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt));
    hV2[pt] = (TH1F *)filein->Get(ProfileV2);
    if (!hV2[pt])
    {
      cout << "Histogram hV2 not available" << endl;
      return;
    }

    StyleHisto(hV2[pt], -0.2, 0.2, 1, 20, titlePt, "v_{2}", TitleInvMass[part] + " " + SInvMass, 1, 0, 100, 1.4, 1.6, 0.7);
    hV2[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[part], histoMassRangeUp[part]);

    if (isYAxisMassZoomed)
    {
      if (part == 0)
        hInvMass[pt]->GetYaxis()->SetRangeUser(0, 2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(1.65)));
      else
        hInvMass[pt]->GetYaxis()->SetRangeUser(0, 2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(1.29)));
    }
    if (isLogy)
    {
      if (hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(histoMassRangeUp[part])) > hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(histoMassRangeLow[part])))
      {
        hInvMass[pt]->SetMinimum(0.8 * hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(histoMassRangeLow[part])));
      }
      else
      {
        hInvMass[pt]->SetMinimum(0.8 * hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(histoMassRangeUp[part])));
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
    // hInvMass[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[part], histoMassRangeUp[part]);
    hInvMass[pt]->Draw("e same");
  }

  // fits

  TF1 **functionsFirst = new TF1 *[numPtBins];
  TF1 **functionsSecond = new TF1 *[numPtBins];
  TF1 **functions1 = new TF1 *[numPtBins];
  TF1 **functions2 = new TF1 *[numPtBins];
  TF1 **bkg1 = new TF1 *[numPtBins];
  TF1 **bkg2 = new TF1 *[numPtBins];
  TF1 **bkg3 = new TF1 *[numPtBins];
  TF1 **bkg4 = new TF1 *[numPtBins];
  TF1 **bkgretta = new TF1 *[numPtBins]; // initial bkg fit
  TF1 **bkgparab = new TF1 *[numPtBins]; // initial bkg fit
  TF1 **bkgpol3 = new TF1 *[numPtBins];  // initial bkg fit
  TF1 **bkgexpo = new TF1 *[numPtBins];  // initial bkg fit
  TF1 **total = new TF1 *[numPtBins];
  TF1 **totalbis = new TF1 *[numPtBins];

  v2fit v2fitarray[numPtBins];
  TF1 **v2FitFunction = new TF1 *[numPtBins];
  TF1 **v2BkgFunction = new TF1 *[numPtBins];

  Double_t parTwoGaussRetta[numPtBins + 1][8];
  Double_t parTwoGaussParab[numPtBins + 1][9];
  Double_t parTwoGaussPol3[numPtBins + 1][10];
  Double_t parTwoGaussExpo[numPtBins + 1][8];
  Double_t parOneGaussRetta[numPtBins + 1][5];
  Double_t parOneGaussParab[numPtBins + 1][6];
  Double_t parOneGaussPol3[numPtBins + 1][7];
  Double_t parOneGaussExpo[numPtBins + 1][5];

  TFitResultPtr fFitResultPtr0[numPtBins];
  TFitResultPtr fFitResultPtr1[numPtBins];

  Float_t mean[numPtBins] = {0};
  Float_t errmean[numPtBins] = {0};
  Float_t sigma[numPtBins] = {0};
  Float_t errsigma[numPtBins] = {0};
  Float_t sigmaw[numPtBins] = {0};
  Float_t errsigmaw[numPtBins] = {0};
  Float_t w1[numPtBins] = {0};
  Float_t w2[numPtBins] = {0};
  Float_t I12[numPtBins] = {0};
  Float_t Yield[numPtBins] = {0};
  Float_t ErrYield[numPtBins] = {0};
  Float_t LowLimit[numPtBins] = {0};
  Float_t UpLimit[numPtBins] = {0};
  Float_t b[numPtBins] = {0};
  Float_t errb[numPtBins] = {0};
  Float_t SSB[numPtBins] = {0};
  Float_t errSSB[numPtBins] = {0};
  Float_t Signif[numPtBins] = {0};
  Float_t errSignif[numPtBins] = {0};
  Float_t entries_range[numPtBins] = {0};
  Float_t TotYield = 0;
  Float_t TotSigBkg = 0;

  TLine *lineP3Sigma[numPtBins];
  TLine *lineM3Sigma[numPtBins];

  Float_t bTest[numPtBins] = {0};
  Float_t errbTest[numPtBins] = {0};
  Float_t SignalTest[numPtBins] = {0};
  Float_t errSignalTest[numPtBins] = {0};
  Float_t YieldTest[numPtBins] = {0};
  Float_t ErrYieldTest[numPtBins] = {0};
  TH1F *hYieldTest[numPtBins];
  TH1F *hYieldRelErrorTest[numPtBins];
  TH1F *hYieldRelErrorTestRelative[numPtBins];
  Float_t LowLimitTest[numPtBins] = {0};
  Float_t UpLimitTest[numPtBins] = {0};
  Float_t LowBin0[numPtBins] = {0};
  Float_t UpBin0[numPtBins] = {0};
  Float_t LowBin[numPtBins] = {0};
  Float_t UpBin[numPtBins] = {0};
  Int_t numMassInt = 10;
  TF1 *LineAt1 = new TF1("LineAt1", "[0]+[1]*x", 0, numMassInt);
  LineAt1->SetParameter(0, 1);
  TF1 *LineAt995 = new TF1("LineAt995", "[0]+[1]*x", 0, numMassInt);
  LineAt995->SetParameter(0, 0.995);

  TLine *lineBkgLimitA[numPtBins];
  TLine *lineBkgLimitB[numPtBins];
  TLine *lineBkgLimitC[numPtBins];
  TLine *lineBkgLimitD[numPtBins];

  Bool_t isV2FromFit[numPtBins] = {0};

  for (Int_t pt = 0; pt < numPtBins; pt++)
  {

    if (!isXi && pt == 0)
      continue;
    if (pt < 4)
      canvas[0]->cd(pt + 1);
    else if (pt < 8)
      canvas[1]->cd(pt + 1 - 4);
    else if (pt < 12)
      canvas[2]->cd(pt + 1 - 8);
    else if (pt < 16)
      canvas[3]->cd(pt + 1 - 12);

    functionsFirst[pt] = new TF1(Form("1f_%i", pt), "gaus", min_range_signal[part], max_range_signal[part]);
    functionsFirst[pt]->SetLineColor(881);
    functionsFirst[pt]->SetParameter(1, ParticleMassPDG[part]);
    functionsFirst[pt]->SetParName(0, "norm");
    functionsFirst[pt]->SetParName(1, "mean");
    functionsFirst[pt]->SetParName(2, "sigma");
    functionsFirst[pt]->SetParLimits(1, min_range_signal[part], max_range_signal[part]);
    functionsFirst[pt]->SetParLimits(2, 0.001, 0.1);
    functionsFirst[pt]->SetParLimits(0, 0, 1.1 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));

    functionsSecond[pt] = new TF1(Form("2f_%i", pt), "gaus", min_range_signal[part], max_range_signal[part]);
    functionsSecond[pt]->SetLineColor(867);
    functionsSecond[pt]->SetParameter(1, ParticleMassPDG[part]);
    functionsSecond[pt]->SetParName(0, "norm");
    functionsSecond[pt]->SetParName(1, "mean");
    functionsSecond[pt]->SetParName(2, "sigma");
    functionsSecond[pt]->SetParLimits(1, min_range_signal[part], max_range_signal[part]);
    functionsSecond[pt]->SetParLimits(2, 0.001, 0.15);
    functionsSecond[pt]->SetParLimits(0, 0, 1.1 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));

    functions1[pt] = new TF1(Form("1f_%i_final", pt), "gaus", gaussDisplayRangeLow[part], gaussDisplayRangeUp[part]);
    functions1[pt]->SetLineColor(kRed); // 867
    functions1[pt]->SetParName(0, "norm");
    functions1[pt]->SetParName(1, "mean");
    functions1[pt]->SetParName(2, "sigma");

    functions2[pt] = new TF1(Form("2f_%i_final", pt), "gaus", gaussDisplayRangeLow[part], gaussDisplayRangeUp[part]);
    functions2[pt]->SetLineColor(kMagenta); // 891
    functions2[pt]->SetParName(0, "norm");
    functions2[pt]->SetParName(1, "mean");
    functions2[pt]->SetParName(2, "sigma");

    bkg1[pt] = new TF1(Form("bkg1%i", pt), "pol1", bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
    bkg1[pt]->SetLineColor(418);
    bkg1[pt]->SetLineStyle(2);

    bkg2[pt] = new TF1(Form("bkg2%i", pt), "pol2", bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
    bkg2[pt]->SetLineColor(1);
    bkg2[pt]->SetLineStyle(2);

    bkg3[pt] = new TF1(Form("bkg3%i", pt), "pol3", bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
    bkg3[pt]->SetLineColor(kOrange + 7);
    bkg3[pt]->SetLineStyle(2);

    bkg4[pt] = new TF1(Form("bkg4%i", pt), "expo", bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
    bkg4[pt]->SetLineColor(kOrange + 7);
    bkg4[pt]->SetLineStyle(2);

    bkgretta[pt] = new TF1(Form("retta%i", pt), fretta, liminf[part], limsup[part], 3);
    bkgretta[pt]->SetLineColor(kGreen + 3);
    bkgretta[pt]->FixParameter(2, part);

    bkgparab[pt] = new TF1(Form("parab%i", pt), fparab, liminf[part], limsup[part], 4);
    bkgparab[pt]->SetLineColor(kAzure + 7);
    bkgparab[pt]->FixParameter(3, part);

    bkgpol3[pt] = new TF1(Form("pol3%i", pt), fpol3, liminf[part], limsup[part], 5);
    bkgpol3[pt]->SetLineColor(kRed + 7);
    bkgpol3[pt]->FixParameter(4, part);

    bkgexpo[pt] = new TF1(Form("expo%i", pt), fexpo, liminf[part], limsup[part], 3);
    bkgexpo[pt]->SetLineColor(kGreen + 2);
    bkgexpo[pt]->FixParameter(2, part);

    Bool_t UseTwoGaussUpdated = 1;
    TF1 *totalFunction = nullptr, *bkgFunction = nullptr;
    if (UseTwoGauss)
    {
      cout << "\n\e[35mFit with two gauss \e[39m"
           << " Pt: " << PtBins[pt] << "-" << PtBins[pt + 1] << endl;

      if (BkgType == 0)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+gaus(3)+pol1(6)", liminf[part], limsup[part]);
      else if (BkgType == 1)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+gaus(3)+pol2(6)", liminf[part], limsup[part]);
      else if (BkgType == 2)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+gaus(3)+pol3(6)", liminf[part], limsup[part]);
      else
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+gaus(3)+expo(6)", liminf[part], limsup[part]);
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

      bkg1[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkg2[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkg3[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkg4[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkgparab[pt]->SetRange(liminf[part], limsup[part]);
      bkgretta[pt]->SetRange(liminf[part], limsup[part]);
      bkgexpo[pt]->SetRange(liminf[part], limsup[part]);
      bkgpol3[pt]->SetRange(liminf[part], limsup[part]);
      total[pt]->SetRange(liminf[part], limsup[part]);

      cout << "\n\n fit bkg " << endl;
      if (BkgType == 0)
        hInvMass[pt]->Fit(bkgretta[pt], "RB0");
      if (BkgType == 1)
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
      if (BkgType == 1)
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
      if (isXi)
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.31, 1.335);
        total[pt]->SetParLimits(2, 0.0012, 0.010);
        total[pt]->SetParLimits(3, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin())); // maximum was wothout 0.3
        total[pt]->SetParLimits(4, 1.31, 1.335);
        total[pt]->SetParLimits(5, 0.001, 0.01);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, ParticleMassPDG[part]);
          total[pt]->FixParameter(4, ParticleMassPDG[part]);
        }
      }
      else
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.66, 1.68);
        total[pt]->SetParLimits(2, 0.002, 0.01);
        total[pt]->SetParLimits(3, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin())); // maximum was wothout 0.3
        total[pt]->SetParLimits(4, 1.66, 1.68);
        total[pt]->SetParLimits(5, 0.001, 0.01);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, ParticleMassPDG[part]);
          total[pt]->FixParameter(4, ParticleMassPDG[part]);
        }
      }

      fFitResultPtr0[pt] = hInvMass[pt]->Fit(total[pt], "SRB0"); // per errore gaussiana, S indica che il risultato del fit e' accessibile da fFitResultPtr0
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

      totalbis[pt] = (TF1 *)total[pt]->Clone();
      fFitResultPtr1[pt] = fFitResultPtr0[pt];

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

      totalFunction = total[pt];
      if (BkgType == 0)
      {
        bkg1[pt]->FixParameter(0, total[pt]->GetParameter(6));
        bkg1[pt]->FixParameter(1, total[pt]->GetParameter(7));
        bkgFunction = bkg1[pt];
      }
      else if (BkgType == 1)
      {
        bkg2[pt]->FixParameter(0, total[pt]->GetParameter(6));
        bkg2[pt]->FixParameter(1, total[pt]->GetParameter(7));
        bkg2[pt]->FixParameter(2, total[pt]->GetParameter(8));
        bkgFunction = bkg2[pt];
      }
      else if (BkgType == 2)
      {
        bkg3[pt]->FixParameter(0, total[pt]->GetParameter(6));
        bkg3[pt]->FixParameter(1, total[pt]->GetParameter(7));
        bkg3[pt]->FixParameter(2, total[pt]->GetParameter(8));
        bkg3[pt]->FixParameter(3, total[pt]->GetParameter(9));
        bkgFunction = bkg3[pt];
      }
      else
      {
        bkg4[pt]->FixParameter(0, total[pt]->GetParameter(6));
        bkg4[pt]->FixParameter(1, total[pt]->GetParameter(7));
        // bkg4[pt]->FixParameter(2, total[pt]->GetParameter(8));
        bkgFunction = bkg4[pt];
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

        hInvMass[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[part], histoMassRangeUp[part]);
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
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+pol1(3)", liminf[part], limsup[part]);
      if (BkgType == 1)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+pol2(3)", liminf[part], limsup[part]);
      else if (BkgType == 2)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+pol3(3)", liminf[part], limsup[part]);
      else
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+expo(3)", liminf[part], limsup[part]);

      total[pt]->SetLineColor(7);
      total[pt]->SetParName(0, "norm");
      total[pt]->SetParName(1, "mean");
      total[pt]->SetParName(2, "sigma");

      cout << "\n\n fit gauss " << endl;
      hInvMass[pt]->Fit(functionsFirst[pt], "RB");

      bkg1[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkg2[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkg3[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkg4[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkgparab[pt]->SetRange(liminf[part], limsup[part]);
      bkgretta[pt]->SetRange(liminf[part], limsup[part]);
      bkgpol3[pt]->SetRange(liminf[part], limsup[part]);
      bkgexpo[pt]->SetRange(liminf[part], limsup[part]);
      total[pt]->SetRange(liminf[part], limsup[part]);

      cout << "\n\n fit bkg " << endl;
      if (BkgType == 0)
        hInvMass[pt]->Fit(bkgretta[pt], "RB0");
      if (BkgType == 1)
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
      if (BkgType == 1)
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
      if (isXi)
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.31, 1.335);
        total[pt]->SetParLimits(2, 0.0012, 0.010);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, ParticleMassPDG[part]);
        }
      }
      else
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.66, 1.68);
        total[pt]->SetParLimits(2, 0.001, 0.02);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, ParticleMassPDG[part]);
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

      totalFunction = total[pt];
      if (BkgType == 0)
      {
        bkg1[pt]->FixParameter(0, total[pt]->GetParameter(3));
        bkg1[pt]->FixParameter(1, total[pt]->GetParameter(4));
        bkgFunction = bkg1[pt];
      }
      else if (BkgType == 1)
      {
        bkg2[pt]->FixParameter(0, total[pt]->GetParameter(3));
        bkg2[pt]->FixParameter(1, total[pt]->GetParameter(4));
        bkg2[pt]->FixParameter(2, total[pt]->GetParameter(5));
        bkgFunction = bkg2[pt];
      }
      else if (BkgType == 2)
      {
        bkg3[pt]->FixParameter(0, total[pt]->GetParameter(3));
        bkg3[pt]->FixParameter(1, total[pt]->GetParameter(4));
        bkg3[pt]->FixParameter(2, total[pt]->GetParameter(5));
        bkg3[pt]->FixParameter(3, total[pt]->GetParameter(6));
        bkgFunction = bkg3[pt];
      }
      else
      {
        bkg4[pt]->FixParameter(0, total[pt]->GetParameter(3));
        bkg4[pt]->FixParameter(1, total[pt]->GetParameter(4));
        // bkg4[pt]->FixParameter(2, total[pt]->GetParameter(5));
        bkgFunction = bkg4[pt];
      }
      if (pt < 4)
        canvas[0]->cd(pt + 1);
      else if (pt < 8)
        canvas[1]->cd(pt + 1 - 4);
      else if (pt < 12)
        canvas[2]->cd(pt + 1 - 8);
      else if (pt < 16)
        canvas[3]->cd(pt + 1 - 12);

      hInvMass[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[part], histoMassRangeUp[part]);
      hInvMass[pt]->Draw("same e");
      mean[pt] = total[pt]->GetParameter(1);
      errmean[pt] = total[pt]->GetParError(1);
      sigma[pt] = total[pt]->GetParameter(2);
      errsigma[pt] = total[pt]->GetParError(2);
    }

    cout << "\nMean: " << mean[pt] << " +/- " << errmean[pt] << endl;
    cout << "Sigma: " << sigma[pt] << " +/- " << errsigma[pt] << endl;

    TLine *linebkgFitLL = new TLine(liminf[part], 0, liminf[part], hInvMass[pt]->GetMaximum()); // low limit of left SB
    TLine *linebkgFitRR = new TLine(limsup[part], 0, limsup[part], hInvMass[pt]->GetMaximum()); // upper limit of right SB
    linebkgFitLL->SetLineColor(kBlue);
    linebkgFitRR->SetLineColor(kBlue);
    TLine *linebkgFitLR = new TLine(UpperLimitLSB, 0, UpperLimitLSB, hInvMass[pt]->GetMaximum()); // upper limit of left SB
    TLine *linebkgFitRL = new TLine(LowerLimitRSB, 0, LowerLimitRSB, hInvMass[pt]->GetMaximum()); // lower limit of right SB
    linebkgFitLR->SetLineColor(kBlue);
    linebkgFitRL->SetLineColor(kBlue);
    lineBkgLimitA[pt] = new TLine(liminf[part], 0, liminf[part], hInvMass[pt]->GetMaximum());
    lineBkgLimitB[pt] = new TLine(UpperLimitLSB, 0, UpperLimitLSB, hInvMass[pt]->GetMaximum());
    lineBkgLimitC[pt] = new TLine(limsup[part], 0, limsup[part], hInvMass[pt]->GetMaximum());
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

    LowLimit[pt] = hInvMass[pt]->GetXaxis()->GetBinLowEdge(hInvMass[pt]->GetXaxis()->FindBin(mean[pt] - sigmacentral * sigma[pt]));
    UpLimit[pt] = hInvMass[pt]->GetXaxis()->GetBinUpEdge(hInvMass[pt]->GetXaxis()->FindBin(mean[pt] + sigmacentral * sigma[pt]));

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
    if (BkgType == 1)
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
    for (Int_t l = hInvMass[pt]->GetXaxis()->FindBin(mean[pt] - sigmacentral * sigma[pt]); l <= hInvMass[pt]->GetXaxis()->FindBin(mean[pt] + sigmacentral * sigma[pt]); l++)
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

    //*********************************************

    histoYield->SetBinContent(pt + 1, Yield[pt] / NEvents / histoYield->GetBinWidth(pt + 1));
    histoYield->SetBinError(pt + 1, ErrYield[pt] / NEvents / histoYield->GetBinWidth(pt + 1));

    histoRelErrYield->SetBinContent(pt + 1, ErrYield[pt] / Yield[pt]);
    histoRelErrYield->SetBinError(pt + 1, 0);

    histoTot->SetBinContent(pt + 1, entries_range[pt] / NEvents / histoYield->GetBinWidth(pt + 1));
    histoTot->SetBinError(pt + 1, sqrt(entries_range[pt]) / NEvents / histoYield->GetBinWidth(pt + 1));

    histoB->SetBinContent(pt + 1, b[pt] / NEvents / histoYield->GetBinWidth(pt + 1));
    histoB->SetBinError(pt + 1, errb[pt] / NEvents / histoYield->GetBinWidth(pt + 1));

    histoMean->SetBinContent(pt + 1, mean[pt]);
    histoMean->SetBinError(pt + 1, errmean[pt]);

    histoSigma->SetBinContent(pt + 1, sigma[pt]);
    histoSigma->SetBinError(pt + 1, errsigma[pt]);

    histoSigmaWeighted->SetBinContent(pt + 1, sigmaw[pt]);
    histoSigmaWeighted->SetBinError(pt + 1, errsigmaw[pt]);

    histoPurity->SetBinContent(pt + 1, SSB[pt]);
    histoPurity->SetBinError(pt + 1, errSSB[pt]);

    histoSignificance->SetBinContent(pt + 1, Yield[pt] / ErrYield[pt]);
    histoSignificance->SetBinError(pt + 1, 0);

    v2fitarray[pt].setBkgFraction(bkgFunction, totalFunction, liminf[part], limsup[part]);
    v2FitFunction[pt] = new TF1(Form("v2function%i", pt), v2fitarray[pt], liminf[part], limsup[part], 3);
    v2FitFunction[pt]->SetLineColor(kRed + 1);
    if (pt < 4)
      canvas[0]->cd(pt + 4 + 1);
    else if (pt < 8)
      canvas[1]->cd(pt + 4 + 1 - 4);
    else if (pt < 12)
      canvas[2]->cd(pt + 4 + 1 - 8);
    else if (pt < 16)
      canvas[3]->cd(pt + 4 + 1 - 12);
    hV2[pt]->Fit(v2FitFunction[pt], "R0");

    v2BkgFunction[pt] = new TF1(Form("v2bkgfunction%i", pt), v2bkgfit, liminf[part], limsup[part], 2);
    v2BkgFunction[pt]->FixParameter(0, v2FitFunction[pt]->GetParameter(1));
    v2BkgFunction[pt]->FixParameter(1, v2FitFunction[pt]->GetParameter(2));
    v2BkgFunction[pt]->SetLineColor(kBlack);
    v2BkgFunction[pt]->SetLineStyle(8);

    hV2[pt]->GetXaxis()->SetTitle(TitleInvMass[part] + " " + SInvMass);
    hV2[pt]->SetTitle(SPt[pt] + " GeV/#it{c}");
    histoV2->SetBinContent(pt + 1, v2FitFunction[pt]->GetParameter(0));
    histoV2->SetBinError(pt + 1, v2FitFunction[pt]->GetParError(0));

    hV2MassIntegrated[pt] = (TH1F *)hmassVsV2C[pt]->ProjectionX(Form("V2CvsMass_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt), hmassVsV2C[pt]->GetYaxis()->FindBin(LowLimit[pt]), hmassVsV2C[pt]->GetYaxis()->FindBin(UpLimit[pt]));
    StyleHisto(hV2MassIntegrated[pt], 0, 1.2 * hV2MassIntegrated[pt]->GetBinContent(hV2MassIntegrated[pt]->GetMaximumBin()), 1, 20, "v_{2}", "Counts", SPt[pt] + " GeV/#it{c}", 1, -1, 1, 1.4, 1.6, 0.7);
    hV2MassIntegrated[pt]->Rebin(2);

    if (SSB[pt] > LimitForV2woFit || PtBins[pt] > 2.)
      isV2FromFit[pt] = 0;
    else
      isV2FromFit[pt] = 1;
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
    if (isV2)
      cout << "\nv2 (no fit): " << hV2MassIntegrated[pt]->GetMean() << " +- " << hV2MassIntegrated[pt]->GetMeanError() << endl;
    else
      cout << "\nPzs2 (no fit): " << hV2MassIntegrated[pt]->GetMean() << " +- " << hV2MassIntegrated[pt]->GetMeanError() << endl;
    hV2MassIntegrated[pt]->ResetStats();
    histoV2NoFit->SetBinContent(pt + 1, hV2MassIntegrated[pt]->GetMean());
    histoV2NoFit->SetBinError(pt + 1, hV2MassIntegrated[pt]->GetMeanError());
    if (v2type == 1)
    {
      fitV2SP[pt] = new TF1(Form("fitV2SP%i", pt), "gaus", -5, 5);
      hV2MassIntegrated[pt]->Fit(fitV2SP[pt], "R+");
      // histoV2NoFit->SetBinContent(pt + 1, fitV2SP[pt]->GetParameter(1));
      // histoV2NoFit->SetBinError(pt + 1, fitV2SP[pt]->GetParError(1));
    }
    cout << histoV2NoFit->GetBinCenter(pt + 1) << " bin c: " << histoV2NoFit->GetBinContent(pt + 1) << endl;

    if (!isV2FromFit[pt])
    {
      histoV2Mixed->SetBinContent(pt + 1, hV2MassIntegrated[pt]->GetMean());
      histoV2Mixed->SetBinError(pt + 1, hV2MassIntegrated[pt]->GetMeanError());
    }
    else
    {
      histoV2Mixed->SetBinContent(pt + 1, v2FitFunction[pt]->GetParameter(0));
      histoV2Mixed->SetBinError(pt + 1, v2FitFunction[pt]->GetParError(0));
    }
  }

  TCanvas *canvasMass = new TCanvas("canvasMass", "canvasMass", 800, 1800);
  canvasMass->Divide(2, 3);
  StyleCanvas(canvasMass, 0.15, 0.05, 0.05, 0.15);

  Int_t index = 0;
  for (Int_t pt = 0; pt < numPtBins; pt++)
  {
    if (!isXi && pt == 0)
      continue;
    if (pt == 2)
      index = 1;
    else if (pt == 3)
      index = 2;
    else if (pt == 5)
      index = 3;
    else if (pt == 7)
      index = 4;
    else if (pt == 9)
      index = 5;
    else if (pt == 11)
      index = 6;
    else
      continue;

    canvasMass->cd(index);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.18);

    hInvMassDraw[pt] = (TH1F *)hInvMass[pt]->Clone(Form("hInvMassDraw%i", pt));
    hInvMassDraw[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[part], histoMassRangeUp[part]);
    hInvMassDraw[pt]->GetYaxis()->SetRangeUser(1, 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
    hInvMassDraw[pt]->Draw("");
    functions1[pt]->Draw("same");
    functions2[pt]->Draw("same");
    bkg2[pt]->SetLineColor(1);
    bkg2[pt]->SetLineStyle(2);
    bkg1[pt]->SetLineColor(1);
    bkg1[pt]->SetLineStyle(2);
    if (BkgType == 0)
      bkg1[pt]->Draw("same");
    else if (BkgType == 1)
      bkg2[pt]->Draw("same");
    else
      bkg3[pt]->Draw("same");
    lineP3Sigma[pt] = new TLine(UpLimit[pt], 0, UpLimit[pt], 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
    lineM3Sigma[pt] = new TLine(LowLimit[pt], 0, LowLimit[pt], 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));

    lineP3Sigma[pt]->Draw("same");
    lineM3Sigma[pt]->Draw("same");
  }

  TCanvas *canvasSummary = new TCanvas("canvasSummary", "canvasSummary", 1700, 1000);
  canvasSummary->Divide(4, 2);

  canvasSummary->cd(1);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoMean, gaussDisplayRangeLow[part], gaussDisplayRangeUp[part], 1, 1, titlePt, "#mu (GeV/c^{2})", "histoMean", 0, 0, 0, 1.4, 1.4, 1.2);
  histoMean->Draw("");
  canvasSummary->cd(2);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoSigma, 0, 0.010, 1, 1, titlePt, "#sigma (GeV/c^{2})", "histoSigma", 0, 0, 0, 1.4, 1.4, 1.2);
  histoSigma->Draw("");
  canvasSummary->cd(3);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoSigmaWeighted, 0, 0.010, 1, 1, titlePt, "#sigma_{w} (GeV/c^{2})", "histoSigmaWeighted", 0, 0, 0, 1.4, 1.4, 1.2);
  histoSigmaWeighted->Draw("");
  canvasSummary->cd(4);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoPurity, 0, 1, 1, 1, titlePt, "S / (S+B)", "histoPurity", 0, 0, 0, 1.4, 1.4, 1.2);
  histoPurity->Draw("");
  canvasSummary->cd(5);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoYield, 0, 1.2 * histoYield->GetBinContent(histoYield->GetMaximumBin()), 1, 1, titlePt, titleYield, "histoYield", 0, 0, 0, 1.4, 1.4, 1.2);
  histoYield->Draw("same");
  canvasSummary->cd(6);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoTot, 0, 1.2 * histoTot->GetBinContent(histoTot->GetMaximumBin()), 1, 1, titlePt, "S+B", "histoTot", 0, 0, 0, 1.4, 1.4, 1.2);
  StyleHisto(histoRelErrYield, 0, 1.2 * histoRelErrYield->GetBinContent(histoRelErrYield->GetMaximumBin()), 1, 1, titlePt, "#sigma_{Y}/Y", "histoRelErrorYield", 0, 0, 0, 1.4, 1.4, 1.2);
  // histoTot->Draw("same");
  histoRelErrYield->Draw("same");
  canvasSummary->cd(7);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoB, 0, 1.2 * histoB->GetBinContent(histoB->GetMaximumBin()), 1, 1, titlePt, "S+B", "histoB", 0, 0, 0, 1.4, 1.4, 1.2);
  // histoB->Draw("same");
  if (isV2)
    histoV2NoFit->SetTitle("v2 without fit");
  else
    histoV2NoFit->SetTitle("Pz,s2 without fit");

  histoV2NoFit->Scale(1. / ftcReso[mul]);
  histoV2Mixed->Scale(1. / ftcReso[mul]);
  if (!isV2)
  {
    histoV2NoFit->Scale(1. / AlphaH[!isXi]);
    histoV2Mixed->Scale(1. / AlphaH[!isXi]);
  }
  histoV2NoFit->Draw();

  canvasSummary->cd(8);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  if (isV2)
    histoV2->SetTitle("v2 from fit");
  else
    histoV2->SetTitle("Pz,s2 from fit");
  histoV2->Scale(1. / ftcReso[mul]);
  if (!isV2)
    histoV2->Scale(1. / AlphaH[!isXi]);
  histoV2->Draw();

  TString Soutputfile;
  Soutputfile = "OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_" + inputFileName + "_" + ParticleName[!isXi] + ChargeName[ExtrCharge + 1];
  Soutputfile += IsOneOrTwoGauss[UseTwoGauss];
  Soutputfile += SIsBkgParab[BkgType];
  Soutputfile += Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);
  Soutputfile += SEtaSysChoice[EtaSysChoice];
  Soutputfile += SBDT;
  if (isApplyWeights)
    Soutputfile += "_Weighted";
  if (v2type == 1)
    Soutputfile += "_SP";
  if (!useCommonBDTValue)
    Soutputfile += "_BDTCentDep";
  if (isRun2Binning)
    Soutputfile += "_Run2Binning";
  // Soutputfile += "_Test";

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
  outputfile->WriteTObject(histoYield);
  outputfile->WriteTObject(histoTot);
  outputfile->WriteTObject(histoB);
  outputfile->WriteTObject(histoMean);
  outputfile->WriteTObject(histoSigma);
  outputfile->WriteTObject(histoSigmaWeighted);
  outputfile->WriteTObject(histoPurity);
  outputfile->WriteTObject(histoSignificance);
  outputfile->WriteTObject(histoV2);
  outputfile->WriteTObject(histoV2NoFit);
  outputfile->WriteTObject(histoV2Mixed);
  outputfile->Close();

  // Performance plot
  Int_t ChosenPt = 8;
  if (isXi)
    ChosenPt = 2;
  Float_t LowLimitMass[numPart] = {1.29, 1.65};
  Float_t UpLimitMass[numPart] = {1.35, 1.7};
  Float_t UpperCutHisto = 1.7;
  if (isXi)
    UpperCutHisto = 1.8;
  Float_t XRangeMin[numPart] = {1.3, 1.656};
  Float_t XRangeMax[numPart] = {1.343, 1.688};

  TString TitleXMass = "#it{m}_{#Lambda#pi} (GeV/#it{c}^{2})";
  if (!isXi)
    TitleXMass = "#it{m}_{#LambdaK} (GeV/#it{c}^{2})";

  TCanvas *canvasMassP = new TCanvas("canvasMassP", "canvasMassP", 800, 800);
  StyleCanvas(canvasMassP, 0.15, 0.03, 0.02, 0.14); // L, R, T, B

  TLegend *legend = new TLegend(0.2, 0.68, 0.71, 0.95);
  legend->SetFillStyle(0);
  legend->SetMargin(0);
  legend->SetTextSize(0.037);
  legend->SetTextAlign(12);
  legend->AddEntry("", "#bf{ALICE Performance}", "");
  legend->AddEntry("", Form("Run 3 Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, %i-%i%s", CentFT0CMin, CentFT0CMax, "%"), "");
  if (isXi)
    legend->AddEntry("", "#Xi^{#minus} #rightarrow #Lambda #pi^{#minus} #rightarrow p #pi^{#minus} #pi^{#minus} + c.c.", "");
  else
    legend->AddEntry("", "#Omega^{#minus} #rightarrow #Lambda K^{#minus} #rightarrow p #pi^{#minus} K^{#minus} + c.c.", "");
  legend->AddEntry("", Form("|#it{#eta}| < 0.8, %.1f < #it{p}_{T} < %.1f GeV/#it{c}", PtBins[ChosenPt], PtBins[ChosenPt + 1]), "");
  legend->AddEntry("", Form("BDT, Signif.(4#sigma) = %.0f #pm %.0f", Signif[ChosenPt], errSignif[ChosenPt]), "");

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
  if (isXi)
    titleyNorm = Form("Normalized counts/(%i MeV/#it{c}^{2})", (int)(histo->GetBinWidth(1) * 1000));
  StyleHisto(histo, 0.0001, UpperCutHisto * histo->GetBinContent(histo->GetMaximumBin()), 1, 20,
             TitleXMass, titleyNorm, "", 1, LowLimitMass[part] + 0.001, UpLimitMass[part] - 0.001, 1.2, 1.8, 1.2);
  histo->GetXaxis()->SetRangeUser(XRangeMin[!isXi], XRangeMax[!isXi]);
  histo->GetXaxis()->SetLabelSize(0.043);
  histo->GetXaxis()->SetTitleSize(0.045);
  histo->GetYaxis()->SetLabelSize(0.043);
  histo->GetYaxis()->SetTitleSize(0.045);
  histo->GetYaxis()->SetTitleOffset(1.6);

  histo->DrawClone("pe");

  TF1 *totalPNorm = new TF1("totalP", "gaus(0)+gaus(3)+pol2(6)", XRangeMin[!isXi], XRangeMax[!isXi]);
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
  bkg->SetParameter(2, bkg->GetParameter(2) / histoIntegral);
  bkg->SetLineColor(kBlack);
  bkg->SetLineStyle(8);
  bkg->Draw("same");
  legendfit->AddEntry(bkg, "bkg.", "l");
  legendfit2->AddEntry(bkg, "bkg.", "l");
  totalPNorm->SetRange(LowLimitMass[part], UpLimitMass[part]);
  totalPNorm->SetLineColor(kRed + 1);
  totalPNorm->Draw("same");
  legend->Draw("");
  legendfit->Draw("");
  canvasMassP->SaveAs("PerformancePlots/MassFit" + ParticleName[!isXi] + Form("_Pt%i.pdf", ChosenPt));
  canvasMassP->SaveAs("PerformancePlots/MassFit" + ParticleName[!isXi] + Form("_Pt%i.png", ChosenPt));
  canvasMassP->SaveAs("PerformancePlots/MassFit" + ParticleName[!isXi] + Form("_Pt%i.eps", ChosenPt));

  TCanvas *canvasP = new TCanvas("canvasP", "canvasP", 800, 1100);
  Float_t LLUpperPad = 0.44;
  Float_t ULLowerPad = 0.44;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

  TLegend *LegendTitle = new TLegend(0.24, 0.65, 0.75, 0.95);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetMargin(0);
  LegendTitle->SetTextSize(0.038);
  LegendTitle->SetTextAlign(12);
  LegendTitle->AddEntry("", "#bf{ALICE Performance}", "");
  LegendTitle->AddEntry("", Form("Run 3 Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, %i-%i%s", CentFT0CMin, CentFT0CMax, "%"), "");
  if (isXi)
    LegendTitle->AddEntry("", "#Xi^{#minus} #rightarrow #Lambda #pi^{#minus} #rightarrow p #pi^{#minus} #pi^{#minus} + c.c.", "");
  else
    LegendTitle->AddEntry("", "#Omega^{#minus} #rightarrow #Lambda K^{#minus} #rightarrow p #pi^{#minus} K^{#minus} + c.c.", "");
  LegendTitle->AddEntry("", Form("|#it{#eta}| < 0.8, %.1f < #it{p}_{T} < %.1f GeV/#it{c}", PtBins[ChosenPt], PtBins[ChosenPt + 1]), "");
  LegendTitle->AddEntry("", Form("BDT, Signif.(4#sigma) = %.0f #pm %.0f", Signif[ChosenPt], errSignif[ChosenPt]), "");

  Float_t LimSupSpectra = 9.99;
  Float_t LimInfSpectra = 0.2 * 1e-5;
  Float_t xTitle = 15;
  Float_t xOffset = 4;
  Float_t yTitle = 30;
  Float_t yOffset = 2;

  Float_t xLabel = 30;
  Float_t yLabel = 30;
  Float_t xLabelOffset = 0.05;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.03;
  Float_t tickY = 0.042;

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, hInvMass[ChosenPt]->GetXaxis()->GetXmin(), hInvMass[ChosenPt]->GetXaxis()->GetXmax());
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasP->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, 1e-3, hInvMass[ChosenPt]->GetMaximum(), 1, 1, TitleXMass, titleyNorm, "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(XRangeMin[!isXi], XRangeMax[!isXi]);
  pad1->Draw();
  pad1->cd();
  // gPad->SetLogy();
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
  Float_t yOffsetR = 2;

  Float_t xLabelR = 30;
  Float_t yLabelR = 30;
  Float_t xLabelOffsetR = 0.02;
  Float_t yLabelOffsetR = 0.02;

  Float_t tickXR = 0.035;
  Float_t tickYR = 0.042;

  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, XRangeMin[!isXi], XRangeMax[!isXi]);
  for (Int_t i = 1; i <= hDummyRatio->GetNbinsX(); i++)
    hDummyRatio->SetBinContent(i, 1e-12);

  StyleHistoYield(hDummyRatio, 1e-5, 0.15 - 1e-5, 1, 1, TitleXMass, "v_{2}", "", 1, 1.15, YoffsetSpectraRatio);
  hDummyRatio->GetXaxis()->SetRangeUser(XRangeMin[!isXi], XRangeMax[!isXi]);
  SetFont(hDummyRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetHistoTextSize(hV2[ChosenPt], xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  hDummyRatio->GetYaxis()->SetNdivisions(704);
  SetTickLength(hDummyRatio, tickXR, tickYR);

  canvasP->cd();
  padL1->Draw();
  padL1->cd();
  hDummyRatio->Draw("same");
  hV2[ChosenPt]->GetXaxis()->SetRangeUser(XRangeMin[!isXi], XRangeMax[!isXi]);
  hV2[ChosenPt]->GetYaxis()->SetRangeUser(0, 0.15);
  hV2[ChosenPt]->SetTitle("");
  hV2[ChosenPt]->Draw("same e");
  v2FitFunction[ChosenPt]->Draw("same");
  v2BkgFunction[ChosenPt]->Draw("same");
  //  hV2MassIntegrated[ChosenPt]->Draw("");
  canvasP->SaveAs("PerformancePlots/MassAnd" + NameAnalysis[!isV2] + ParticleName[!isXi] + Form("_Pt%i.pdf", ChosenPt));
  canvasP->SaveAs("PerformancePlots/MassAnd" + NameAnalysis[!isV2] + ParticleName[!isXi] + Form("_Pt%i.png", ChosenPt));
  canvasP->SaveAs("PerformancePlots/MassAnd" + NameAnalysis[!isV2] + ParticleName[!isXi] + Form("_Pt%i.eps", ChosenPt));

  cout << "\nA partire dal file:\n"
       << SPathIn << endl;
  cout << "\nHo creato il file: " << Soutputfile << ".root" << endl;
}
