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
const Int_t numPtBins = 4;
Float_t PtBins[numPtBins + 1] = {0, 1.2, 2, 3, 4};
const Int_t numCent = 3;//10;
Int_t CentFT0C[numCent + 1] = {10, 30, 50, 90};

Float_t ParticleMassPDG[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
using namespace ROOT;

void ComputeV2(Bool_t isXi = 1, TString inputFileName = "16March_New")
{

  TString SinputFile = "OutputAnalysis/Output_" + inputFileName + ".root";
  TFile *inputFile = new TFile(SinputFile);
  TH3D *hmassVsPtVsV2C[numCent];
  TH2F *hmassVsPt[numCent];
  TH2F *hmassVsV2C[numCent][numPtBins];
  TH1F *hmass[numCent][numPtBins];
  TH1F *hV2C[numCent][numPtBins];
  TString hName = "";
  TString hNameMass = "";
  TString hNameV2C = "";
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    hName = Form("massVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]);
    hmassVsPtVsV2C[cent] = (TH3D *)inputFile->Get(hName);
    hmassVsPt[cent] = (TH2F *)hmassVsPtVsV2C[cent]->Project3D("yx");
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      hNameMass = Form("mass_cent%i-%i_pt%i", CentFT0C[cent], CentFT0C[cent + 1], pt);
      hNameV2C = Form("V2C_cent%i-%i_pt%i", CentFT0C[cent], CentFT0C[cent + 1], pt);
      // hmassVsPtVsV2C[cent]->GetXaxis()->SetRangeUser(1.28, 1.36);
      hmassVsPtVsV2C[cent]->GetYaxis()->SetRangeUser(PtBins[pt], PtBins[pt + 1]);
      hmassVsV2C[cent][pt] = (TH2F *)hmassVsPtVsV2C[cent]->Project3D("xz");
      hmassVsV2C[cent][pt]->SetName(hNameV2C);
      hmass[cent][pt] = (TH1F *)hmassVsPtVsV2C[cent]->Project3D("x");
      hmass[cent][pt]->SetName(hNameMass);
      hV2C[cent][pt] = (TH1F *)hmass[cent][pt]->Clone(hNameV2C);
      for (Int_t bin = 0; bin < hmass[cent][pt]->GetNbinsX(); bin++)
      {
        // hV2C[cent]->SetBinContent(bin, hmassVsV2C[cent]->GetXaxis()->GetMean());
      }
    }
  }

  TFile *file = new TFile("OutputAnalysis/V2_" + inputFileName + ".root", "RECREATE");
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      hmass[cent][pt]->Write();
      hV2C[cent][pt]->Write();
      hmassVsV2C[cent][pt]->Write();
    }
    hmassVsPtVsV2C[cent]->Write();
    hmassVsPt[cent]->Write();
  }
  file->Close();
  cout << "I created the file " << file->GetName() << endl;
}
