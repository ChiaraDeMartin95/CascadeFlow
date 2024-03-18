// This macro was originally written by:
// chiara.de.martin@cern.ch

#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TPad.h"
#include "TSpline.h"

void StyleHisto(TH1D &histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRange,
                Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize)
{
  histo.GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)
    histo.GetXaxis()->SetRangeUser(XLow, XUp);
  histo.SetLineColor(color);
  histo.SetMarkerColor(color);
  histo.SetMarkerStyle(style);
  histo.SetMarkerSize(mSize);
  histo.GetXaxis()->SetTitle(titleX);
  histo.GetXaxis()->SetLabelSize(0.05);
  histo.GetXaxis()->SetTitleSize(0.05);
  histo.GetXaxis()->SetTitleOffset(xOffset);
  histo.GetYaxis()->SetTitle(titleY);
  histo.GetYaxis()->SetTitleSize(0.05);
  histo.GetYaxis()->SetLabelSize(0.05);
  histo.GetYaxis()->SetTitleOffset(yOffset);
  histo.SetTitle(title);
}

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

const Int_t numPart = 7;
const Int_t numChoice = 5; // mean, sigma, purity, yield, efficiency for MC
const Int_t numCent = 3;//10;
Int_t CentFT0C[numCent + 1] = {10, 30, 50, 90};
//Int_t CentFT0C[numCent + 1] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
Float_t ParticleMassPDG[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
using namespace ROOT;

void ProcessTree(Bool_t isXi = 1, TString inputFileName = "16March_New")
{
  TString TreeName = "O2cascanalysis";

  TString inputFile = "TreeForAnalysis/AnalysisResults_trees_" + inputFileName + ".root";
  RDataFrame d1(TreeName, inputFile);
  auto h = d1.Histo1D("fPt");
  h->Draw();

  // invariant mass histograms
  auto mass_Xi_Bef = d1.Histo1D({"mass_Xi_Bef", "Invariant mass of #Lambda#pi", 100, 1.28, 1.36}, "fMassXi");
  auto mass_Omega_Bef = d1.Histo1D({"mass_Omega_Bef", "Invariant mass of #LambdaK", 100, 1.6, 1.73}, "fMassOmega");
  // BDT response histogram
  auto BDT_response_Bef = d1.Histo1D({"BDT_response_Bef", "BDT response", 100, 0, 1}, "fBDTResponseXi");
  if (!isXi)
    BDT_response_Bef = d1.Histo1D({"BDT_response_Bef", "BDT response", 100, 0, 1}, "fBDTResponseOmega");
  // apply BDT selection
  auto d2 = d1.Filter("fBDTResponseXi < 0.98");
  if (!isXi)
    d2 = d1.Filter("fBDTResponseOmega > 0.9");

  auto BDT_response = d2.Histo1D({"BDT_response", "BDT response", 100, 0, 1}, "fBDTResponseXi");
  if (!isXi)
    BDT_response = d2.Histo1D({"BDT_response", "BDT response", 100, 0, 1}, "fBDTResponseOmega");

  // invariant mass histograms
  auto mass_Xi = d2.Histo1D({"mass_Xi", "Invariant mass of #Lambda#pi", 100, 1.28, 1.36}, "fMassXi");
  auto mass_Omega = d2.Histo1D({"mass_Omega", "Invariant mass of #LambdaK", 100, 1.6, 1.73}, "fMassOmega");

  // 3D histograms
  TFile *file = new TFile("OutputAnalysis/Output_" + inputFileName + ".root", "RECREATE");

  for (Int_t cent = 0; cent < numCent; cent++)
  {
    auto dcent = d2.Filter(Form("fCentFT0C>%i && fCentFT0C<%i", CentFT0C[cent], CentFT0C[cent + 1]));
    if (isXi)
    {
      auto massVsPtVsV2C = dcent.Histo3D({Form("massVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "Invariant mass vs Pt vs V2C", 100, 1.28, 1.36, 100, 0, 10, 100, 0, 1}, "fMassXi", "fPt", "fV2C");
      massVsPtVsV2C->Write();
    }
    else
    {
      auto massVsPtVsV2C = dcent.Histo3D({Form("massVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "Invariant mass vs Pt vs V2C", 100, 1.6, 1.73, 100, 0, 10, 100, 0, 1}, "fMassOmega", "fPt", "fV2C");
      massVsPtVsV2C->Write();
    }
  }
  // draw histograms
  TCanvas *cMassXi = new TCanvas("cMassXi", "cMassXi", 900, 600);
  StyleCanvas(cMassXi, 0.1, 0.1, 0.03, 0.1);
  StyleHisto(*mass_Xi_Bef, 0, 1.2 * mass_Xi_Bef->GetMaximum(), kRed, 20, "M_{#Lambda#pi}", "Counts", "", kTRUE, 1.2, 1.4, 1.2, 1.2, 0.7);
  StyleHisto(*mass_Xi, 0, 1.2 * mass_Xi_Bef->GetMaximum(), kBlue, 20, "M_{#Lambda#pi}", "Counts", "", kTRUE, 1.2, 1.4, 1.2, 1.2, 0.7);
  mass_Xi_Bef->Draw("E");
  mass_Xi->Draw("E SAME");

  TCanvas *cMassOmega = new TCanvas("cMassOmega", "cMassOmega", 900, 600);
  StyleCanvas(cMassOmega, 0.1, 0.1, 0.03, 0.1);
  StyleHisto(*mass_Omega_Bef, 0, 1.2 * mass_Omega_Bef->GetMaximum(), kRed, 20, "M_{#LambdaK}", "Counts", "", kTRUE, 1.5, 1.7, 1.2, 1.2, 0.7);
  StyleHisto(*mass_Omega, 0, 1.2 * mass_Omega_Bef->GetMaximum(), kBlue, 20, "M_{#LambdaK}", "Counts", "", kTRUE, 1.5, 1.7, 1.2, 1.2, 0.7);
  mass_Omega_Bef->Draw("E");
  mass_Omega->Draw("E SAME");

  // create THnSparse

  // create TH3 histograms

  cMassXi->Write();
  cMassOmega->Write();
  mass_Xi_Bef->Write();
  mass_Omega_Bef->Write();
  BDT_response_Bef->Write();
  mass_Xi->Write();
  mass_Omega->Write();
  BDT_response->Write();
  file->Close();

  cout << "I created the file " << file->GetName() << endl;
}
