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

Float_t DefineMixedBDTValue(Int_t mul = 0, Int_t pt = 0)
{
  if (CentFT0C[mul] >= 40)
  {
    return 0.4;
  }
  else
  {
    if (PtBins[pt] < 1)
    {
      return 0.96;
    }
    else if (PtBins[pt] < 1.2)
    {
      return 0.8;
    }
    else if (PtBins[pt] < 2)
    {
      return 0.64;
    }
    else
    {
      return 0.4;
    }
  }
}

TString titleYield = "1/N_{evt} dN/dp_{T}";
TString titleYieldNN = "dN/dp_{T}";
TString titlePt = "p_{T} (GeV/c)";
TString TitleInvMass[numPart] = {"(#Lambda, #pi)", "(#Lambda, K)", "(#Lambda, #pi^{-})", "(#overline{#Lambda}, #pi^{+})", "(#Lambda, K^{-})", "(#overline{#Lambda}, K^{+})"};
TString SInvMass = "invariant mass (GeV/c^{2})";

// fit ranges
Float_t min_range_signal[numPart] = {1.3, 1.65, 1.3, 1.3, 1.65, 1.65}; // gauss fit range
Float_t max_range_signal[numPart] = {1.335, 1.69, 1.335, 1.335, 1.69, 1.69};
Float_t liminf[numPart] = {1.29, 1.63, 1.29, 1.29, 1.63, 1.63}; // bkg and total fit range
Float_t limsup[numPart] = {1.352, 1.71, 1.352, 1.352, 1.71, 1.71};
Float_t XRangeMin[numPart] = {1.301, 1.656, 1.3, 1.3, 1.656, 1.656};
Float_t XRangeMax[numPart] = {1.344, 1.688, 1.343, 1.343, 1.688, 1.688};

// visualisation ranges
Float_t LowMassRange[numPart] = {1.31, 1.655, 1.31, 1.31, 1.655, 1.655}; // range to compute approximate yield (signal + bkg)
Float_t UpMassRange[numPart] = {1.33, 1.685, 1.33, 1.33, 1.685, 1.685};
Float_t gaussDisplayRangeLow[numPart] = {1.29, 1.63, 1.29, 1.29, 1.63, 1.63}; // display range of gauss functions (from total fit)
Float_t gaussDisplayRangeUp[numPart] = {1.35, 1.71, 1.35, 1.35, 1.71, 1.71};
Float_t bkgDisplayRangeLow[numPart] = {1.29, 1.626, 1.29, 1.29, 1.626, 1.626}; // display range of bkg function (from total fit)
Float_t bkgDisplayRangeUp[numPart] = {1.35, 1.72, 1.35, 1.35, 1.72, 1.72};
Float_t histoMassRangeLow[numPart] = {1.301, 1.626, 1.301, 1.301, 1.626, 1.626}; // display range of mass histograms
Float_t histoMassRangeUp[numPart] = {1.344, 1.72, 1.344, 1.344, 1.72, 1.72};

// Event plane resolution
Float_t ftcReso[numCent + 1] = {0};

void Acceptance(Int_t indexMultTrial = 0,
                Int_t ChosenPart = ChosenParticle,
                Bool_t isRapiditySel = ExtrisRapiditySel,
                TString inputFileName = SinputFileName,
                Int_t EtaSysChoice = ExtrEtaSysChoice,
                Bool_t isSysMultTrial = ExtrisSysMultTrial)
{

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

  TString SinputFile = "OutputAnalysis/Output" + STHN[ExtrisFromTHN] + "_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice]; // + SBDT;
  if (isApplyWeights)
    SinputFile += "_Weighted";
  if (v2type == 1)
    SinputFile += "_SP";
  if (!useCommonBDTValue)
    SinputFile += "_BDTCentDep";
  if (isRun2Binning)
    SinputFile += "_Run2Binning";
  if (!isRapiditySel)
    SinputFile += "_Eta08";
  SinputFile += SBDT;
  SinputFile += ".root";
  cout << "Input file: " << SinputFile << endl;
  TFile *inputFile = new TFile(SinputFile);
  if (!inputFile)
  {
    cout << "File not found" << endl;
    return;
  }

  TH3D *hEtaVsPtVsCos2ThetaLambdaFromC[numCent + 1];
  TString hNameCos2ThetaLambdaFromC_Eta3D[numCent + 1] = {""};
  
  TH2F *hEtaVsCos2ThetaLambdaFromC[numCent + 1][numPtBins + 1];
  TString hNameEtaCos2ThetaLambdaFromC[numCent + 1][numPtBins + 1] = {""};
  TH1F *hCos2ThetaLambdaFromCVsEta[numCent + 1][numPtBins + 1];
  TString hNameCos2ThetaLambdaFromCVsEta[numCent + 1][numPtBins + 1] = {""};
  TH1F *hEta[numCent + 1][numPtBins + 1];
  TString hNameEta[numCent + 1][numPtBins + 1] = {""};
  
  TH2F *hPtVsCos2ThetaLambdaFromC[numCent + 1][numEtaBins + 1];
  TString hNamePtVsCos2ThetaLambdaFromC[numCent + 1][numEtaBins + 1] = {""};
  TH1F *hCos2ThetaLambdaFromCVsPt[numCent + 1][numEtaBins + 1];
  TString hNameCos2ThetaLambdaFromCVsPt[numCent + 1][numEtaBins + 1] = {""};
  TH1F *hPtLambda[numCent + 1][numEtaBins + 1];
  TString hNamePtLambda[numCent + 1][numEtaBins + 1] = {""};

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;

  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    if (cent == numCent)
    { // 0-80%
      CentFT0CMin = 0;
      CentFT0CMax = 80;
    }
    else
    {
      CentFT0CMin = CentFT0C[cent];
      CentFT0CMax = CentFT0C[cent + 1];
    }
    hNameCos2ThetaLambdaFromC_Eta3D[cent] = Form("etaVsPtVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hEtaVsPtVsCos2ThetaLambdaFromC[cent] = (TH3D *)inputFile->Get(hNameCos2ThetaLambdaFromC_Eta3D[cent]);
    if (!hEtaVsPtVsCos2ThetaLambdaFromC[cent])
    {
      cout << "Histogram hEtaVsPtVsCos2ThetaLambdaFromC not available" << endl;
    }

    for (Int_t pt = 0; pt < numPtBins + 1; pt++)
    {
      hNameEtaCos2ThetaLambdaFromC[cent][pt] = Form("EtaVsCos2ThetaLambdaFromC_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameCos2ThetaLambdaFromCVsEta[cent][pt] = Form("Cos2ThetaLambdaFromCVsEta_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameEta[cent][pt] = Form("Eta_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      if (pt == numPtBins)
      {
        hEtaVsPtVsCos2ThetaLambdaFromC[cent]->GetYaxis()->SetRangeUser(PtBins[0] + 0.0001, PtBins[numPtBins] - 0.0001);
      }
      else
      {
        hEtaVsPtVsCos2ThetaLambdaFromC[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
      }
      hEtaVsCos2ThetaLambdaFromC[cent][pt] = (TH2F *)hEtaVsPtVsCos2ThetaLambdaFromC[cent]->Project3D("xze"); // eta vs cos2thetaLambdaFromC
      hEtaVsCos2ThetaLambdaFromC[cent][pt]->SetName(hNameEtaCos2ThetaLambdaFromC[cent][pt]);

      hEta[cent][pt] = (TH1F *)hEtaVsPtVsCos2ThetaLambdaFromC[cent]->Project3D("xe");
      hEta[cent][pt]->SetName(hNameEta[cent][pt]);

      hCos2ThetaLambdaFromCVsEta[cent][pt] = (TH1F *)hEta[cent][pt]->Clone(hNameCos2ThetaLambdaFromCVsEta[cent][pt]);
      hCos2ThetaLambdaFromCVsEta[cent][pt]->Reset();
      for (Int_t bin = 0; bin < hEta[cent][pt]->GetNbinsX(); bin++)
      {
        TH1D *htemp = (TH1D *)hEtaVsCos2ThetaLambdaFromC[cent][pt]->ProjectionX(Form("_htemp_%i", bin), bin + 1, bin + 1);
        hCos2ThetaLambdaFromCVsEta[cent][pt]->SetBinContent(bin + 1, htemp->GetMean());
        hCos2ThetaLambdaFromCVsEta[cent][pt]->SetBinError(bin + 1, htemp->GetMeanError());
      }
    }
    for (Int_t eta = 0; eta < numEtaBins + 1; eta++)
    {
      
      hNameCos2ThetaLambdaFromCVsPt[cent][eta] = Form("Cos2ThetaLambdaFromCVsPt_cent%i-%i_eta%i", CentFT0CMin, CentFT0CMax, eta);
      hNamePtLambda[cent][eta] = Form("PtLambda_cent%i-%i_eta%i", CentFT0CMin, CentFT0CMax, eta);

      if (eta == numEtaBins)
      {
        hEtaVsPtVsCos2ThetaLambdaFromC[cent]->GetXaxis()->SetRangeUser(EtaBins[0] + 0.0001, EtaBins[numEtaBins] - 0.0001);
      }
      else
      {
        hEtaVsPtVsCos2ThetaLambdaFromC[cent]->GetYaxis()->SetRangeUser(EtaBins[eta] + 0.0001, EtaBins[eta + 1] - 0.0001);
      }

      hPtVsCos2ThetaLambdaFromC[cent][eta] = (TH2F *)hEtaVsPtVsCos2ThetaLambdaFromC[cent]->Project3D("yze");
      hPtVsCos2ThetaLambdaFromC[cent][eta]->SetName(hNamePtVsCos2ThetaLambdaFromC[cent][eta]);
      
      hPtLambda[cent][eta] = (TH1F *)hEtaVsPtVsCos2ThetaLambdaFromC[cent]->Project3D("ye");
      hPtLambda[cent][eta]->SetName(hNamePtLambda[cent][eta]);

      hCos2ThetaLambdaFromCVsPt[cent][eta] = (TH1F *)hPtLambda[cent][eta]->Clone(hNameCos2ThetaLambdaFromCVsPt[cent][eta]);
      hCos2ThetaLambdaFromCVsPt[cent][eta]->Reset();

      for (Int_t bin = 0; bin < hPtLambda[cent][eta]->GetNbinsX(); bin++)
      {
        TH1D *htemp = (TH1D *)hPtVsCos2ThetaLambdaFromC[cent][eta]->ProjectionX(Form("_htemp_%i", bin), bin + 1, bin + 1);
        hCos2ThetaLambdaFromCVsPt[cent][eta]->SetBinContent(bin + 1, htemp->GetMean());
        hCos2ThetaLambdaFromCVsPt[cent][eta]->SetBinError(bin + 1, htemp->GetMeanError());
      }
    }
  }

  TCanvas *canvasAcceptance = new TCanvas("canvasAcceptance", "canvasAcceptance", 1200, 800);
  StyleCanvas(canvasAcceptance, 0.15, 0.1, 0.1, 0.1);
  canvasAcceptance->Divide(2, 2);
  // cos2 vs eta for different pt bins
  canvasAcceptance->cd(1);
  for (Int_t pt = 0; pt < numPtBins + 1; pt++)
  {
    Int_t cent = 7;
    hCos2ThetaLambdaFromCVsEta[cent][pt]->GetYaxis()->SetRangeUser(0, 0.5);
    hCos2ThetaLambdaFromCVsEta[cent][pt]->SetMarkerStyle(20);
    hCos2ThetaLambdaFromCVsEta[cent][pt]->SetMarkerSize(0.5);
    hCos2ThetaLambdaFromCVsEta[cent][pt]->SetMarkerColor(ColorMult[pt]);
    hCos2ThetaLambdaFromCVsEta[cent][pt]->SetLineColor(ColorMult[pt]);
    hCos2ThetaLambdaFromCVsEta[cent][pt]->Draw("same");
  }

  TString SOutputFile = "AcceptancePlots/Acceptance_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (isApplyWeights)
    SOutputFile += "_Weighted";
  if (v2type == 1)
    SOutputFile += "_SP";
  if (!useCommonBDTValue)
    SOutputFile += "_BDTCentDep";
  if (isRun2Binning)
    SOutputFile += "_Run2Binning";
  if (ExtrisApplyEffWeights)
  {
    SOutputFile += "_EffW";
  }
  SOutputFile += "_WithAlpha";
  if (!isRapiditySel || ExtrisFromTHN)
    SOutputFile += "_Eta08";
  SOutputFile += STHN[ExtrisFromTHN] + ".root";
  TFile *file = new TFile(SOutputFile, "RECREATE");
  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    hEtaVsPtVsCos2ThetaLambdaFromC[cent]->Write();
    for (Int_t pt = 0; pt < numPtBins + 1; pt++)
    {
      hEtaVsCos2ThetaLambdaFromC[cent][pt]->Write();
      hEta[cent][pt]->Write();
      hCos2ThetaLambdaFromCVsEta[cent][pt]->Write();
    }
    for (Int_t eta = 0; eta < numEtaBins + 1; eta++)
    {
      hPtVsCos2ThetaLambdaFromC[cent][eta]->Write();
      hPtLambda[cent][eta]->Write();
      hCos2ThetaLambdaFromCVsPt[cent][eta]->Write();
    }
  }
  file->Close();
  cout << "I produced the file: " << SOutputFile << endl;
}
