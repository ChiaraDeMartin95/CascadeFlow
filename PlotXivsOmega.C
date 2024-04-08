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
#include "CommonVar.h"

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRange,
                Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)
    histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetLabelSize(0.07);
  histo->GetXaxis()->SetTitleSize(0.07);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.09);
  histo->GetYaxis()->SetLabelSize(0.07);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->SetTitle(title);
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

TString NameRun2[4] = {"10-20 %", "20-30 %", "30-40 %", "40-50 %"};
Int_t Color[12] = {kYellow + 2, kGreen + 2, kBlack, kBlue + 1, kYellow + 2, kGreen + 2, kBlack, kBlue + 1, kYellow + 2, kGreen + 2, kBlack, kBlue + 1};
Int_t Marker[] = {1, 25, 30, 27, 24};
Int_t MarkerRun3[12] = {20, 21, 29, 33, 20, 21, 29, 33, 20, 21, 29, 33};
Float_t MarkerSize[] = {1.5, 1.5, 1.5, 2., 1.5, 1.5, 1.5, 2., 1.5, 1.5, 1.5, 2.};

void PlotXivsOmega(Bool_t isPlotRatio = 1,
                   TString inputFileName = SinputFileName,
                   Bool_t UseTwoGauss = ExtrUseTwoGauss,
                   Int_t BkgType = ExtrBkgType)
{

  TString SinputFile = "";
  TFile *inputFile;
  TH1F *hV2C[numCent][numPart];
  TH1F *hV2CRatio[numCent];

  TH1F *histoDummy = new TH1F("histoDummy", "histoDummy", 100, 0, 6);
  if (isPlotRatio)
    StyleHisto(histoDummy, 0. + 1e-4, 2. - 1e-4, 1, 20, "#it{p}_{T} (GeV/#it{c})", "v_{2} #Xi / v_{2} #Omega", "", kFALSE, 0, 100, 1, 1, 0.05);
  else
    StyleHisto(histoDummy, -0.2 + 1e-4, 0.8 - 1e-4, 1, 20, "#it{p}_{T} (GeV/#it{c})", "v_{2}", "", kFALSE, 0, 100, 1, 1, 0.05);

  for (Int_t cent = 0; cent < numCent; cent++)
  {
    for (Int_t part = 0; part < numPart; part++)
    {
      SinputFile = "OutputAnalysis/FitV2_" + inputFileName + "_" + ParticleName[part];
      SinputFile += IsOneOrTwoGauss[UseTwoGauss];
      SinputFile += SIsBkgParab[BkgType];
      SinputFile += Form("_Cent%i-%i.root", CentFT0C[cent], CentFT0C[cent + 1]);
      inputFile = new TFile(SinputFile);
      cout << "Input file with Run 3 v2: " << SinputFile << endl;

      hV2C[cent][part] = (TH1F *)inputFile->Get("histoV2");
      hV2C[cent][part]->SetName(Form("hV2C_%i_part%i", cent, part));
      hV2C[cent][part]->SetMarkerStyle(MarkerPart[part]);
      hV2C[cent][part]->SetMarkerSize(MarkerPartSize[part]);
      hV2C[cent][part]->SetMarkerColor(ColorPart[part]);
      hV2C[cent][part]->SetLineColor(ColorPart[part]);
    }
  }

  TLegend *LegendTitle = new TLegend(0.15, 0.81, 0.42, 0.91);
  LegendTitle->SetBorderSize(0);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextSize(0.04);
  LegendTitle->SetMargin(0.);

  TLegend *legendPart = new TLegend(0.25, 0.7, 0.6, 0.9);
  legendPart->SetBorderSize(0);
  legendPart->SetFillStyle(0);
  legendPart->SetTextSize(0.08);

  TFile *file = new TFile("OutputAnalysis/PlotXivsOmega_" + inputFileName + ".root", "RECREATE");

  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 1200);
  Float_t PadWidth = 0.54;
  Float_t PadHeight = 2. / numCent;
  TPad *pad[numCent];
  StyleCanvas(canvas, 0.1, 0.05, 0.05, 0.15);
  gStyle->SetOptStat(0);

  Float_t plotHeight = 0.238;
  Float_t MinY = 0;
  Float_t MaxY = 0;
  TF1* fLine = new TF1("fLine", "1", 0, 6);
  fLine->SetLineColor(kBlack);
  fLine->SetLineStyle(10);
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    if (cent == 0 || cent == 1)
    {
      MinY = 0.28 + 2 * plotHeight;
      MaxY = 1;
    }
    else if (cent == 2 || cent == 3)
    {
      MinY = 0.28 + plotHeight;
      MaxY = 0.28 + 2 * plotHeight;
    }
    else if (cent == 4 || cent == 5)
    {
      MinY = 0.28;
      MaxY = 0.28 + plotHeight;
    }
    else
    {
      MinY = 0;
      MaxY = 0.28;
    }
    if (cent % 2 == 0)
    {
      pad[cent] = new TPad(Form("pad%i", cent), Form("pad%i", cent), 0, MinY, PadWidth, MaxY);
      if (cent == 0)
        StylePad(pad[cent], 0.18, 0.0, 0.03, 0.); // L, R, T, B
      else if (cent == 6)
        StylePad(pad[cent], 0.18, 0.0, 0, 0.15); // L, R, T, B
      else
        StylePad(pad[cent], 0.18, 0.0, 0, 0); // L, R, T, B
    }
    else
    {
      pad[cent] = new TPad(Form("pad%i", cent), Form("pad%i", cent), PadWidth, MinY, 1, MaxY);
      if (cent == 1)
        StylePad(pad[cent], 0., 0.02, 0.03, 0.); // L, R, T, B
      else if (cent == 7)
        StylePad(pad[cent], 0., 0.02, 0, 0.15); // L, R, T, B
      else
        StylePad(pad[cent], 0., 0.02, 0, 0.); // L, R, T, B
    }
    canvas->cd();
    pad[cent]->Draw();
    pad[cent]->cd();
    histoDummy->Draw("");

    TLegend *legendMult;
    if (cent % 2 == 0)
      legendMult = new TLegend(0.7, 0.7, 0.8, 0.9);
    else
      legendMult = new TLegend(0.65, 0.7, 0.8, 0.9);
    legendMult->SetBorderSize(0);
    legendMult->SetFillStyle(0);
    legendMult->SetTextSize(0.08);
    legendMult->SetHeader(Form("%i-%i %%", CentFT0C[cent], CentFT0C[cent + 1]));
    legendMult->Draw();

    for (Int_t part = 0; part < numPart; part++)
    {

      if (part == 0)
        hV2CRatio[cent] = (TH1F *)hV2C[cent][part]->Clone(Form("hV2CRatio_%i", cent));
      else
      {
        hV2CRatio[cent]->Divide(hV2C[cent][part]);
        for (Int_t bin = 1; bin <= hV2CRatio[cent]->GetNbinsX(); bin++)
        {
          //cout << hV2CRatio[cent]->GetBinCenter(bin) << " " << hV2C[cent][part]->GetBinCenter(bin) << endl;
          //cout << hV2CRatio[cent]->GetBinContent(bin) << " +- " << hV2CRatio[cent]->GetBinError(bin) << endl;
        }
      }
      if (isPlotRatio)
      {
        if (part == 1){
          hV2CRatio[cent]->DrawCopy("same");
          fLine->Draw("same");
        }
      }
      else
      {
        hV2C[cent][part]->Draw("same");
      }

      if (cent == 0)
        legendPart->AddEntry(hV2C[cent][part], ParticleName[part], "p");
    }
    if (cent == 0 && !isPlotRatio)
      legendPart->Draw();
  }
  LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}, Pb-Pb 5.36 TeV", "");

  canvas->SaveAs("OutputAnalysis/PlotXivsOmega_" + inputFileName + ".pdf");

  file->Close();
  cout << "I created the file " << file->GetName() << endl;
}
