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
#include "CommonVarPub.h"
// #include "CommonVar_v2.h"
// #include "CommonVarOmega.h"
#include "CommonVarXi.h"
// #include "CommonVarLambda.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TF1Convolution.h"

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

const Float_t UpperLimitLSBOmega = 1.663;  // upper limit of fit of left sidebands for omega
const Float_t LowerLimitRSBOmega = 1.678;  // lower limit of fit of right sidebands for omega
const Float_t UpperLimitLSBXi = 1.308;     // upper limit of fit of left sidebands for Xi
const Float_t LowerLimitRSBXi = 1.335;     // lower limit of fit of right sidebands for Xi
const Float_t UpperLimitLSBLambda = 1.107; // upper limit of fit of left sidebands for Lambda
const Float_t LowerLimitRSBLambda = 1.124; // lower limit of fit of right sidebands for Lambda
const Int_t numBinsEta = 8;
const Bool_t reject = 1;

// fit ranges
Float_t min_range_signal[numPart] = {1.31, 1.65, 1.31, 1.31, 1.65, 1.65, 1.11, 1.11, 1.11}; // gauss fit range
Float_t max_range_signal[numPart] = {1.335, 1.69, 1.335, 1.335, 1.69, 1.69, 1.12, 1.12, 1.12};
Float_t liminf[numPart] = {1.3, 1.656, 1.3, 1.3, 1.656, 1.656, 1.1, 1.1, 1.1}; // bkg and total fit range
Float_t limsup[numPart] = {1.345, 1.686, 1.345, 1.345, 1.686, 1.686, 1.13, 1.13, 1.13};
Float_t liminfBkg[numPart] = {1.3, 1.656, 1.3, 1.3, 1.656, 1.656, 1.1, 1.1, 1.1}; // bkg and total fit range
Float_t limsupBkg[numPart] = {1.345, 1.686, 1.345, 1.345, 1.686, 1.686, 1.13, 1.13, 1.13};
Float_t liminfV2[numPart] = {1.304, 1.656, 1.29, 1.29, 1.656, 1.656, 1.1, 1.1, 1.1}; // v2 fit range
Float_t limsupV2[numPart] = {1.339, 1.686, 1.352, 1.352, 1.686, 1.686, 1.3, 1.3, 1.3};
Float_t XRangeMin[numPart] = {1.303, 1.656, 1.3, 1.3, 1.656, 1.656, 1.1, 1.1, 1.1};
Float_t XRangeMax[numPart] = {1.340, 1.688, 1.343, 1.343, 1.688, 1.688, 1.13, 1.13, 1.13};

// visualisation ranges
Float_t LowMassRange[numPart] = {1.31, 1.656, 1.31, 1.31, 1.656, 1.656, 1.1, 1.1, 1.1}; // range to compute approximate yield (signal + bkg)
Float_t UpMassRange[numPart] = {1.33, 1.686, 1.33, 1.33, 1.686, 1.686, 1.13, 1.13, 1.13};
Float_t gaussDisplayRangeLow[numPart] = {1.29, 1.656, 1.29, 1.29, 1.656, 1.656, 1.1, 1.1, 1.1}; // display range of gauss functions (from total fit)
Float_t gaussDisplayRangeUp[numPart] = {1.35, 1.69, 1.35, 1.35, 1.69, 1.69, 1.13, 1.13, 1.13};
Float_t bkgDisplayRangeLow[numPart] = {1.29, 1.656, 1.29, 1.29, 1.656, 1.656, 1.1, 1.1, 1.1}; // display range of bkg function (from total fit)
Float_t bkgDisplayRangeUp[numPart] = {1.35, 1.686, 1.35, 1.35, 1.686, 1.686, 1.13, 1.13, 1.13};
Float_t histoMassRangeLow[numPart] = {1.301, 1.656, 1.301, 1.301, 1.656, 1.656, 1.1, 1.1, 1.1}; // display range of mass histograms
Float_t histoMassRangeUp[numPart] = {1.344, 1.686, 1.344, 1.344, 1.686, 1.686, 1.13, 1.13, 1.13};

TString titleYield = "1/N_{evt} dN/dp_{T}";
TString titleYieldNN = "dN/dp_{T}";
TString titlePt = "p_{T} (GeV/c)";
TString TitleInvMass[numPart] = {"(#Lambda, #pi)", "(#Lambda, K)", "(#Lambda, #pi^{-})", "(#overline{#Lambda}, #pi^{+})", "(#Lambda, K^{-})", "(#overline{#Lambda}, K^{+})", "(p, #pi)", "(p, #pi^{-})", "(#overline{p}, #pi^{+})"};
TString SInvMass = "invariant mass (GeV/c^{2})";
TString TitleXMass = "#it{M}_{#Lambda#pi} (GeV/#it{c}^{2})";

void PlotXiPzs2VsMass()
{

  Int_t ChosenPart = 0;
  Int_t CentFT0CMin1 = 0;
  Int_t CentFT0CMax1 = 10;
  Int_t CentFT0CMin2 = 60;
  Int_t CentFT0CMax2 = 70;
  TString inputFile1 = Form("../OutputAnalysis/FitPzs2_LHC23_PbPb_pass5_Train540301_Xi_DSCB_BkgCheb_GaussConv_Cent%i-%i_PolFromLambda_Eta08_FromTHN_MixedBDT_TightMassCut2.1_EPReso_isBkgPol0.root", CentFT0CMin1, CentFT0CMax1);
  TString inputFile2 = Form("../OutputAnalysis/FitPzs2_LHC23_PbPb_pass5_Train540301_Xi_DSCB_BkgCheb_GaussConv_Cent%i-%i_PolFromLambda_Eta08_FromTHN_MixedBDT_TightMassCut2.1_EPReso_isBkgPol0.root", CentFT0CMin2, CentFT0CMax2);

  TFile *file0 = TFile::Open(inputFile1);
  if (!file0)
  {
    std::cerr << "Error: Could not open file " << inputFile1 << std::endl;
    return;
  }
  TFile *file1 = TFile::Open(inputFile2);
  if (!file1)
  {
    std::cerr << "Error: Could not open file " << inputFile2 << std::endl;
    return;
  }

  TH1F *hInvMass1 = (TH1F *)file0->Get(Form("mass_cent%i-%i_pt7", CentFT0CMin1, CentFT0CMax1));
  TH1F *hInvMass2 = (TH1F *)file1->Get(Form("mass_cent%i-%i_pt7", CentFT0CMin2, CentFT0CMax2));
  TH1F *hV2P1 = (TH1F *)file0->Get(Form("Pzs2LambdaFromC_cent%i-%i_pt7_Profile", CentFT0CMin1, CentFT0CMax1));
  TH1F *hV2P2 = (TH1F *)file1->Get(Form("Pzs2LambdaFromC_cent%i-%i_pt7_Profile", CentFT0CMin2, CentFT0CMax2));
  TF1 *ffitPzs1 = (TF1 *)file0->Get("v2functionExt7");
  ffitPzs1->SetName(Form("v2functionExt7_cent%i-%i", CentFT0CMin1, CentFT0CMax1));
  TF1 *ffitBkg1 = (TF1 *)file0->Get("v2bkgfunction7");
  ffitBkg1->SetName(Form("v2bkgfunction7_cent%i-%i", CentFT0CMin1, CentFT0CMax1));
  TF1 *ffitPzs2 = (TF1 *)file1->Get("v2functionExt7");
  ffitPzs2->SetName(Form("v2functionExt7_cent%i-%i", CentFT0CMin2, CentFT0CMax2));
  TF1 *ffitBkg2 = (TF1 *)file1->Get("v2bkgfunction7");
  ffitBkg2->SetName(Form("v2bkgfunction7_cent%i-%i", CentFT0CMin2, CentFT0CMax2));
  // invariant-mass fit functions (normalised to the integral of the mass histogram)
  TF1 *ftotal1 = (TF1 *)file0->Get("totalP");
  TF1 *ftotal2 = (TF1 *)file1->Get("totalP");
  TF1 *fbkg1 = (TF1 *)file0->Get("bkg57");
  TF1 *fbkg2 = (TF1 *)file1->Get("bkg57");
  if (!ftotal1 || !fbkg1)
    std::cerr << "Warning: totalP or bkg57 not found in " << inputFile1 << std::endl;
  if (!ftotal2 || !fbkg2)
    std::cerr << "Warning: totalP or bkg57 not found in " << inputFile2 << std::endl;
  if (ftotal1)
    ftotal1->SetName(Form("totalP_cent%i-%i", CentFT0CMin1, CentFT0CMax1));
  if (ftotal2)
    ftotal2->SetName(Form("totalP_cent%i-%i", CentFT0CMin2, CentFT0CMax2));
  if (fbkg1)
    fbkg1->SetName(Form("bkg57_cent%i-%i", CentFT0CMin1, CentFT0CMax1));
  if (fbkg2)
    fbkg2->SetName(Form("bkg57_cent%i-%i", CentFT0CMin2, CentFT0CMax2));
  // the mass histograms are stored before normalisation: normalise them as done for the fit functions
  hInvMass1->Scale(1. / hInvMass1->Integral("width"));
  hInvMass2->Scale(1. / hInvMass2->Integral("width"));

  const Int_t numCent = 2;
  Int_t CentMin[numCent] = {CentFT0CMin1, CentFT0CMin2};
  Int_t CentMax[numCent] = {CentFT0CMax1, CentFT0CMax2};
  TH1F *hInvMass[numCent] = {hInvMass1, hInvMass2};
  TH1F *hV2P[numCent] = {hV2P1, hV2P2};
  TF1 *ffitPzs[numCent] = {ffitPzs1, ffitPzs2};
  TF1 *ffitBkg[numCent] = {ffitBkg1, ffitBkg2};
  TF1 *ftotal[numCent] = {ftotal1, ftotal2};
  TF1 *fbkg[numCent] = {fbkg1, fbkg2};

  Float_t LLUpperPad = 0.44;
  Float_t ULLowerPad = 0.44;
  TCanvas *canvasP = new TCanvas("canvasP", "canvasP", 1600, 1100);
  StyleCanvas(canvasP, 0., 0., 0., 0.); // removes stat boxes and legend borders
  TPad *pad1[numCent];
  TPad *padL1[numCent];

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

  Float_t tickX = 0.035; // 0.03
  Float_t tickY = 0.025; // 0.042

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
  Float_t yLabelOffsetR = 0.014;

  Float_t tickXR = 0.035; // 0.035
  Float_t tickYR = 0.025; // 0.042

  Float_t histoIntegral1 = hInvMass1->Integral("width");
  Float_t histoIntegral2 = hInvMass2->Integral("width");
  hInvMass1->Scale(1. / histoIntegral1);
  hInvMass2->Scale(1. / histoIntegral2);

  Float_t MaxInvMass = TMath::Max(hInvMass1->GetMaximum(), hInvMass2->GetMaximum()); // common y range for the mass pads
  Float_t MinInvMass = MaxInvMass;
  for (Int_t ic = 0; ic < 2; ic++)
  {
    TH1F *h = (ic == 0) ? hInvMass1 : hInvMass2;
    for (Int_t i = h->FindBin(XRangeMin[ChosenPart]); i <= h->FindBin(XRangeMax[ChosenPart]); i++)
      if (h->GetBinContent(i) > 0 && h->GetBinContent(i) < MinInvMass)
        MinInvMass = h->GetBinContent(i);
  }
  TString titleyNorm = Form("Normalized counts/(%.1f MeV/#it{c}^{2})", hInvMass1->GetBinWidth(1) * 1000);
  TString TitleDummyRatio = "#it{P}_{z, s2}";

  TH1F *hDummy[numCent];
  TH1F *hDummyRatio[numCent];
  TLegend *legendfit2[numCent];
  TLegend *legendChi2[numCent];
  TLegend *LegendTitle[numCent];
  TLegend *LegendCent[numCent];

  for (Int_t ic = 0; ic < numCent; ic++)
  {
    // pad widths chosen so that the two frames have the same width: W0 * (1 - 0.18) = W1 * (1 - 0.04)
    Float_t xSplit = (1 - 0.04) / ((1 - 0.18) + (1 - 0.04));
    Float_t xLow = (ic == 0) ? 0. : xSplit;
    Float_t xUp = (ic == 0) ? xSplit : 1.;
    pad1[ic] = new TPad(Form("pad1_%i", ic), Form("pad1_%i", ic), xLow, LLUpperPad, xUp, 1); // xlow, ylow, xup, yup
    padL1[ic] = new TPad(Form("padL1_%i", ic), Form("padL1_%i", ic), xLow, 0, xUp, ULLowerPad);
    // no space between left and right pads: left pad has only the left margin, right pad only the right margin
    Float_t LMargin = (ic == 0) ? 0.18 : 0.;
    Float_t RMargin = (ic == 0) ? 0. : 0.04;
    Float_t xShift = (ic == 0) ? 0. : -0.18;          // shift of legends to keep the same position w.r.t. the frame
    StylePad(pad1[ic], LMargin, RMargin, 0.07, 0.);   // L, R, T, B
    StylePad(padL1[ic], LMargin, RMargin, 0.02, 0.3); // L, R, T, B

    legendfit2[ic] = new TLegend(0.255 + xShift, 0.69, 0.72 + xShift, 0.90);
    legendfit2[ic]->SetFillStyle(0);
    legendfit2[ic]->SetMargin(0.1);
    legendfit2[ic]->SetTextSize(0.042);
    legendfit2[ic]->SetTextAlign(12);

    legendChi2[ic] = new TLegend(0.4 + xShift, 0.82, 0.65 + xShift, 0.92);
    legendChi2[ic]->SetFillStyle(0);
    legendChi2[ic]->SetMargin(0);
    legendChi2[ic]->SetTextSize(0.05);
    legendChi2[ic]->SetTextAlign(12);

    // general info (left pad only)
    LegendTitle[ic] = new TLegend(0.24 + xShift, 0.69, 0.75 + xShift, 0.90);
    LegendTitle[ic]->SetFillStyle(0);
    LegendTitle[ic]->SetMargin(0);
    LegendTitle[ic]->SetTextSize(0.05);
    LegendTitle[ic]->SetTextAlign(12);
    LegendTitle[ic]->AddEntry("", "ALICE, Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
    LegendTitle[ic]->AddEntry("", "#Xi^{#minus} #rightarrow #Lambda #pi^{#minus} #rightarrow p #pi^{#minus} #pi^{#minus} + c.c.", "");
    LegendTitle[ic]->AddEntry("", Form("|#it{#eta}| < 0.8, #it{p}_{T} > %.1f GeV/#it{c}", PtBins[0]), "");

    // centrality (both pads, same height)
    LegendCent[ic] = new TLegend(0.24 + xShift, 0.58, 0.75 + xShift, 0.65);
    LegendCent[ic]->SetFillStyle(0);
    LegendCent[ic]->SetMargin(0);
    LegendCent[ic]->SetTextSize(0.05);
    LegendCent[ic]->SetTextAlign(12);
    LegendCent[ic]->AddEntry("", Form("FT0C centrality: %i#minus%i%s", CentMin[ic], CentMax[ic], "%"), "");

    // upper pad: invariant mass
    hDummy[ic] = new TH1F(Form("hDummy_%i", ic), Form("hDummy_%i", ic), 10000, hInvMass[ic]->GetXaxis()->GetXmin(), hInvMass[ic]->GetXaxis()->GetXmax());
    for (Int_t i = 1; i <= hDummy[ic]->GetNbinsX(); i++)
      hDummy[ic]->SetBinContent(i, 1e-12);
    SetFont(hDummy[ic]);
    StyleHistoYield(hDummy[ic], 0.5 * MinInvMass, 50 * MaxInvMass, 1, 1, TitleXMass, titleyNorm, "", 1, 1.15, 1.8);
    SetHistoTextSize(hDummy[ic], xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
    SetTickLength(hDummy[ic], tickX, tickY);
    hDummy[ic]->GetXaxis()->SetRangeUser(XRangeMin[ChosenPart], XRangeMax[ChosenPart]);
    canvasP->cd();
    pad1[ic]->Draw();
    pad1[ic]->cd();
    gPad->SetLogy();
    if (ic > 0)
    { // same y range as left pad: no title and no labels
      hDummy[ic]->GetYaxis()->SetTitle("");
      hDummy[ic]->GetYaxis()->SetLabelSize(0);
    }
    hDummy[ic]->Draw("same");
    hInvMass[ic]->SetMarkerSize(1.0);
    hInvMass[ic]->Draw("hist same pe");
    legendfit2[ic]->AddEntry(hInvMass[ic], "Data", "pe");
    if (ftotal[ic])
    {
      ftotal[ic]->SetLineColor(kMagenta + 1);
      ftotal[ic]->Draw("same");
      legendfit2[ic]->AddEntry(ftotal[ic], "Total fit", "l");
    }
    if (fbkg[ic])
    {
      fbkg[ic]->SetLineColor(kBlue + 2);
      fbkg[ic]->SetLineStyle(2);
      fbkg[ic]->Draw("same");
      legendfit2[ic]->AddEntry(fbkg[ic], "Background", "l");
    }
    // totalPNorm->Draw("same");
    if (ic == 0)
      LegendTitle[ic]->Draw("");
    LegendCent[ic]->Draw("");
    if (ic == numCent - 1)
      legendfit2[ic]->Draw(""); // fit legend only in the right pad
    // bkg->Draw("same");

    // lower pad: Pzs2 vs mass
    hDummyRatio[ic] = new TH1F(Form("hDummyRatio_%i", ic), Form("hDummyRatio_%i", ic), 10000, XRangeMin[ChosenPart], XRangeMax[ChosenPart]);
    for (Int_t i = 1; i <= hDummyRatio[ic]->GetNbinsX(); i++)
      hDummyRatio[ic]->SetBinContent(i, -9999);
    StyleHistoYield(hDummyRatio[ic], 1e-5, 0.15 - 1e-5, 1, 1, TitleXMass, TitleDummyRatio, "", 1, 1.15, YoffsetSpectraRatio);
    hDummyRatio[ic]->GetXaxis()->SetRangeUser(XRangeMin[ChosenPart], XRangeMax[ChosenPart]);
    SetFont(hDummyRatio[ic]);
    SetHistoTextSize(hDummyRatio[ic], xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
    hDummyRatio[ic]->GetYaxis()->SetNdivisions(704);
    SetTickLength(hDummyRatio[ic], tickXR, tickYR);
    canvasP->cd();
    padL1[ic]->Draw();
    padL1[ic]->cd();
    hDummyRatio[ic]->GetYaxis()->SetRangeUser(-0.0036, 0.0036);
    hDummyRatio[ic]->GetYaxis()->SetRangeUser(-0.01, 0.01);
    if (ic > 0)
    { // same y range as left pad: no title and no labels
      hDummyRatio[ic]->GetYaxis()->SetTitle("");
      hDummyRatio[ic]->GetYaxis()->SetLabelSize(0);
    }
    else
      hDummyRatio[ic]->GetXaxis()->ChangeLabel(-1, -1, 0); // hide last x label, overlapping with the first one of the right pad
    hDummyRatio[ic]->Draw("same");
    hV2P[ic]->SetMarkerSize(0.9);
    hV2P[ic]->SetTitle("");
    hV2P[ic]->Draw("same e0"); // e0: error bars drawn also for points outside the y range
    ffitPzs[ic]->SetLineColor(kMagenta + 1);
    ffitPzs[ic]->Draw("same");
    ffitBkg[ic]->SetLineColor(kBlue + 2);
    ffitBkg[ic]->SetLineStyle(2);
    ffitBkg[ic]->Draw("same");
    legendChi2[ic]->Draw("");
  }
  canvasP->Modified();
  canvasP->Update();

  canvasP->SaveAs("../XiPzs2VsMass.pdf");
  canvasP->SaveAs("../XiPzs2VsMass.png");
}
