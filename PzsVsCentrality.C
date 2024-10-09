#include "Riostream.h"
#include "TTimer.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TH3F.h>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TProfile.h"
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TFile.h>
#include <TLine.h>
#include <TSpline.h>
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
#include "CommonVar.h"
#include "ErrRatioCorr.C"

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(1.5);
  histo->GetXaxis()->SetTitle(TitleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(TitleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.3);
  histo->SetTitle(title);
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

void StyleCanvas(TCanvas *canvas, Float_t TopMargin, Float_t BottomMargin, Float_t LeftMargin, Float_t RightMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  gPad->SetTopMargin(TopMargin);
  gPad->SetLeftMargin(LeftMargin);
  gPad->SetBottomMargin(BottomMargin);
  gPad->SetRightMargin(RightMargin);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
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

Float_t YLow[numPart] = {-0.002};
Float_t YUp[numPart] = {0.002};

void PzsVsCentrality(Bool_t isXi = ChosenParticleXi,
                     Int_t Choice = 0,
                     Int_t part = ExtrParticle,
                     Int_t BkgType = ExtrBkgType,
                     Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn;
  TFile *fileIn[numCent + 1];

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = "Pzs2VsCentrality" + NameAnalysis[!isV2] + "_";
  stringout += SinputFileName;
  stringout += "_" + ParticleName[!isXi] + ChargeName[ExtrCharge + 1];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[BkgType];
  stringout += "_" + TypeHisto[Choice];
  if (isApplyWeights)
    stringout += "_Weighted";
  if (!useCommonBDTValue)
    stringout += "_BDTCentDep";
  if (isRun2Binning)
    stringout += "_Run2Binning";
  stringoutpdf = stringout;
  stringout += ".root";

  // canvases
  gStyle->SetOptStat(0);
  TCanvas *canvasPzs = new TCanvas("canvasPzs", "canvasPzs", 900, 700);
  StyleCanvas(canvasPzs, 0.05, 0.15, 0.15, 0.05);
  TH1F *fHistPzs = new TH1F("fHistPzs", "fHistPzs", numCent, fCentFT0C);

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *LegendTitle;
  if (Choice == 2)
    LegendTitle = new TLegend(0.54, 0.55, 0.95, 0.72);
  else
    LegendTitle = new TLegend(0.54, 0.75, 0.95, 0.92);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(33);
  LegendTitle->SetTextSize(0.04);
  LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  LegendTitle->AddEntry("", "PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
  LegendTitle->AddEntry("", ParticleName[part] + ChargeName[ExtrCharge + 1] + " |#it{y}| < 0.5", "");

  TLine *lineat0 = new TLine(0, 0, 100, 0);
  lineat0->SetLineColor(1);
  lineat0->SetLineStyle(2);

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  TH1F *fHistSpectrum[numCent + 1];
  TString Smolt[numCent + 1];
  TString SmoltBis[numCent + 1];
  // get spectra in multiplicity classes
  for (Int_t m = 0; m < numCent; m++)
  {
    if (m == numCent)
    { // 0-80%
      CentFT0CMin = 0;
      CentFT0CMax = 80;
    }
    else
    {
      CentFT0CMin = CentFT0C[m];
      CentFT0CMax = CentFT0C[m + 1];
    }

    PathIn = "OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_";
    PathIn += SinputFileName;
    PathIn += "_" + ParticleName[!isXi] + ChargeName[ExtrCharge + 1];
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[BkgType];
    Smolt[m] += Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);
    SmoltBis[m] += Form("%i#minus%i", CentFT0CMin, CentFT0CMax);
    PathIn += Smolt[m];
    if (isApplyWeights)
      PathIn += "_Weighted";
    if (!useCommonBDTValue)
      PathIn += "_BDTCentDep";
    if (isRun2Binning)
      PathIn += "_Run2Binning";
    PathIn += ".root";
    cout << "Path in : " << PathIn << endl;
    fileIn[m] = TFile::Open(PathIn);
    fHistSpectrum[m] = (TH1F *)fileIn[m]->Get("histoPzs2PtIntMixed");
    if (!fHistSpectrum[m])
    {
      cout << " no hist " << endl;
      return;
    }
    fHistSpectrum[m]->SetName("histoPzs2_" + Smolt[m]);
    fHistPzs->SetBinContent(m + 1, fHistSpectrum[m]->GetBinContent(1));
    fHistPzs->SetBinError(m + 1, fHistSpectrum[m]->GetBinError(1));
    cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << " Pzs2: " << fHistSpectrum[m]->GetBinContent(1) << endl;
  } // end loop on mult

  Float_t xTitle = 30;
  Float_t xOffset = 1.5;
  Float_t yTitle = 30;
  Float_t yOffset = 2;

  Float_t xLabel = 30;
  Float_t yLabel = 30;
  Float_t xLabelOffset = 0.03;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.03;
  Float_t tickY = 0.042;

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasPzs->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleXCent, TitleYPzs, "", 1, 1.15, 1.6);
  StyleHistoYield(fHistPzs, YLow[part], YUp[part], ColorPart[part], MarkerPart[part], TitleXCent, TitleYPzs, "", MarkerPartSize[part], 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->Draw("");
  fHistPzs->Draw("same");
  LegendTitle->Draw("");

  TFile *fileout = new TFile(stringout, "RECREATE");
  fHistPzs->Write();
  fileout->Close();

  canvasPzs->SaveAs(stringoutpdf + ".pdf");
  canvasPzs->SaveAs(stringoutpdf + ".png");
  cout << "\nStarting from the files (for the different mult): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
}
