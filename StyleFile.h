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

