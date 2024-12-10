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

// isSPReso:
// 1 --> resolution for SP method
// 0 --> resolution for EP method

void Resolution(Bool_t isSPReso = 1, Bool_t isLFReso = 1)
{

  TString ComputeResoFileName = "TreeForAnalysis/AnalysisResults_";
  if (isLFReso == 1)
    ComputeResoFileName += inputFileResoLF + ".root";
  else
    ComputeResoFileName += inputFileResoCFW + ".root";
  if (!isLFReso && (ComputeResoFileName.Index("CFW") == -1))
  {
    cout << "You are looking for file: " << ComputeResoFileName << endl;
    cout << "Error: CFW resolution not available" << endl;
    return;
  }
  TFile *inputfile = new TFile(ComputeResoFileName, "");
  TString nameQT0CTPCA = "QVectorsT0CTPCA";
  TString nameQT0CTPCC = "QVectorsT0CTPCC";
  TString nameQTPCAC = "QVectorsTPCAC";
  if (!isSPReso)
  {
    nameQT0CTPCA = "QVectorsNormT0CTPCA";
    nameQT0CTPCC = "QVectorsNormT0CTPCC";
    nameQTPCAC = "QVectorsNormTPCAC";
  }

  cout << "InputFile: " << ComputeResoFileName << endl;
  if (!inputfile)
  {
    cout << "Error: input file not found" << endl;
    return;
  }
  TDirectory *dir = (TDirectory *)inputfile->Get("lf-cascade-flow/resolution");
  if (!dir)
  {
    cout << "Error: input directory not found" << endl;
    return;
  }
  TH2F *hQT0CTPCA = (TH2F *)dir->Get(nameQT0CTPCA);
  TH2F *hQT0CTPCC = (TH2F *)dir->Get(nameQT0CTPCC);
  TH2F *hQTPCAC = (TH2F *)dir->Get(nameQTPCAC);
  if (!hQT0CTPCA || !hQT0CTPCC || !hQTPCAC)
  {
    cout << "Error: input histograms not found" << endl;
    return;
  }

  TH1D *hReso = new TH1D("hReso", "hReso", numCent, fCentFT0C);
  TH1D *hReso080 = new TH1D("hReso080", "hReso080", 1, 0, 1);
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
    cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << "%" << endl;
    hQT0CTPCA->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    hQT0CTPCC->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    hQTPCAC->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    TH1F *h1QT0CTPCA = (TH1F *)hQT0CTPCA->ProjectionX("h1QT0CTPCA");
    TH1F *h1QT0CTPCC = (TH1F *)hQT0CTPCC->ProjectionX("h1QT0CTPCC");
    TH1F *h1QTPCAC = (TH1F *)hQTPCAC->ProjectionX("h1QTPCAC");
    Float_t MeanT0CTPCA = h1QT0CTPCA->GetMean();
    Float_t MeanT0CTPCC = h1QT0CTPCC->GetMean();
    Float_t MeanTPCAC = h1QTPCAC->GetMean();
    Float_t ErrMeanT0CTPCA = h1QT0CTPCA->GetMeanError();
    Float_t ErrMeanT0CTPCC = h1QT0CTPCC->GetMeanError();
    Float_t ErrMeanTPCAC = h1QTPCAC->GetMeanError();
    cout << "MeanT0CTPCA: " << MeanT0CTPCA << " MeanT0CTPCC: " << MeanT0CTPCC << " MeanTPCAC: " << MeanTPCAC << " Reso: " << sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC) << endl;
    // cout << "Sourav resolution (EP): " << ftcReso[cent] << endl;
    Float_t RelErr2 = pow(ErrMeanT0CTPCA / MeanT0CTPCA, 2) + pow(ErrMeanT0CTPCC / MeanT0CTPCC, 2) + pow(ErrMeanTPCAC / MeanTPCAC, 2);
    Float_t ErrReso = sqrt(RelErr2 * (MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
    if (cent == numCent)
    {
      hReso080->SetBinContent(1, sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
      hReso080->SetBinError(1, ErrReso);
    }
    else
    {
      hReso->SetBinContent(cent + 1, sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
      hReso->SetBinError(cent + 1, ErrReso);
    }
  }

  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  StyleCanvas(canvas, 0.05, 0.1, 0.1, 0.05);
  hReso->GetYaxis()->SetRangeUser(0., 1.5);
  hReso->SetLineColor(kBlack);
  hReso->SetMarkerColor(kBlack);
  hReso->SetMarkerStyle(20);
  hReso->SetMarkerSize(1.5);
  hReso->GetXaxis()->SetTitle("Centrality (%)");
  hReso->GetYaxis()->SetTitle("Resolution");
  hReso->GetXaxis()->SetTitleSize(0.05);
  hReso->GetXaxis()->SetLabelSize(0.05);
  hReso->GetXaxis()->SetTitleOffset(1.);
  hReso->GetYaxis()->SetTitleSize(0.05);
  hReso->GetYaxis()->SetLabelSize(0.05);
  hReso->GetYaxis()->SetTitleOffset(1.);
  hReso->SetTitle("");
  hReso->Draw("");

  TH1F *hResoSourav = (TH1F *)hReso->Clone("hResoSourav");
  hResoSourav->SetLineColor(kRed);
  hResoSourav->SetMarkerColor(kRed);
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    hResoSourav->SetBinContent(cent + 1, ftcResoSourav[cent]);
    hResoSourav->SetBinError(cent + 1, 0);
  }
  // if (!isSPReso)
  // hResoSourav->Draw("same");

  TLegend *legend = new TLegend(0.15, 0.67, 0.75, 0.93);
  legend->SetFillStyle(0);
  legend->SetMargin(0);
  legend->SetTextSize(0.05);
  legend->SetTextAlign(12);
  legend->AddEntry("", "#bf{ALICE Performance}", "");
  legend->AddEntry("", "Run 3, Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
  legend->AddEntry("", "T0C (#minus3.3 < #it{#eta} < #minus2.1) and TPC (0.1 < |#it{#eta}| < 0.8)", "");
  legend->Draw();

  TString Soutputfile = "";
  if (!isSPReso)
  { // event plane method
    if (isLFReso)
      Soutputfile = ResoFileName_EPLF;
    else
      Soutputfile = ResoFileName_EPCFW;
  }
  else
  {
    if (isLFReso)
      Soutputfile = ResoFileName_SPLF;
    else
      Soutputfile = ResoFileName_SPCFW;
  }
  canvas->SaveAs(Soutputfile + ".pdf");
  canvas->SaveAs(Soutputfile + ".png");
  canvas->SaveAs(Soutputfile + ".eps");
  TFile *outputfile = new TFile(Soutputfile + ".root", "RECREATE");
  hReso->Write();
  hReso080->Write();
  outputfile->Close();
  cout << "\nOutputFile: " << Soutputfile << endl;
}
