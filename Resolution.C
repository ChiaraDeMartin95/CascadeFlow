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
    // return;
  }
  TFile *inputfile = new TFile(ComputeResoFileName, "");
  TString nameQT0CTPCA = "QVectorsT0CTPCA";
  TString nameQT0CTPCC = "QVectorsT0CTPCC";
  TString nameQTPCAC = "QVectorsTPCAC";
  TString nameQT0CV0A = "QVectorsT0CV0A";
  TString nameQV0ATPCA = "QVectorsV0ATPCA";
  TString nameQV0ATPCC = "QVectorsV0ATPCC";
  if (!isSPReso)
  {
    nameQT0CTPCA = "QVectorsNormT0CTPCA";
    nameQT0CTPCC = "QVectorsNormT0CTPCC";
    nameQTPCAC = "QVectorsNormTPCAC";
    nameQT0CV0A = "QVectorsNormT0CV0A";
    nameQV0ATPCA = "QVectorsNormV0ATPCA";
    nameQV0ATPCC = "QVectorsNormV0ATPCC";
  }

  cout << "InputFile: " << ComputeResoFileName << endl;
  if (!inputfile)
  {
    cout << "Error: input file not found" << endl;
    return;
  }
  TDirectory *dir = (TDirectory *)inputfile->Get("lf-cascade-flow/resolution");
  // TDirectory *dir = (TDirectory *)inputfile->Get("lf-cascade-flow");
  if (!dir)
  {
    cout << "Error: input directory not found" << endl;
    return;
  }
  TH2F *hQT0CTPCA = (TH2F *)dir->Get(nameQT0CTPCA);
  TH2F *hQT0CTPCC = (TH2F *)dir->Get(nameQT0CTPCC);
  TH2F *hQTPCAC = (TH2F *)dir->Get(nameQTPCAC);
  TH2F *hQT0CV0A = (TH2F *)dir->Get(nameQT0CV0A);
  TH2F *hQV0ATPCA = (TH2F *)dir->Get(nameQV0ATPCA);
  TH2F *hQV0ATPCC = (TH2F *)dir->Get(nameQV0ATPCC);
  if (!hQT0CTPCA || !hQT0CTPCC || !hQTPCAC)
  {
    cout << "Error: input histograms not found" << endl;
    return;
  }
  if (!hQT0CV0A)
  {
    hQT0CV0A = (TH2F *)dir->Get(nameQT0CTPCA);
    cout << "Error: input histograms T0CV0A not found, taking T0CTPCA instead" << endl;
  }
  if (!hQV0ATPCA)
  {
    hQV0ATPCA = (TH2F *)dir->Get(nameQT0CTPCC);
    cout << "Error: input histograms V0ATPCA not found, taking T0CTPCC instead" << endl;
  }
  if (!hQV0ATPCC)
  {
    hQV0ATPCC = (TH2F *)dir->Get(nameQTPCAC);
    cout << "Error: input histograms V0ATPCC not found, taking TPCAC instead" << endl;
  }

  TH1D *hReso = new TH1D("hReso", "hReso", numCent, fCentFT0C);
  if (isOOCentrality)
    hReso = new TH1D("hReso", "hReso", numCentLambdaOO, fCentFT0CLambdaOO);
  TH1D *hResoV0A = (TH1D*)hReso->Clone("hResoV0A");  
  TH1D *hResoPerCentBins = new TH1D("hResoPerCentBins", "hResoPerCentBins", 100, 0, 100);
  TH1D *hResoPerCentBinsV0A = new TH1D("hResoPerCentBinsV0A", "hResoPerCentBinsV0A", 100, 0, 100);
  TH1D *hReso080 = new TH1D("hReso080", "hReso080", 1, 0, 1);
  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  Int_t commonnumCent = numCent;
  if (isOOCentrality)
    commonnumCent = numCentLambdaOO;
  for (Int_t cent = 0; cent < commonnumCent + 1; cent++)
  {

    if (cent == commonnumCent)
    { // 0-80%
      CentFT0CMin = 0;
      CentFT0CMax = 80;
      if (isOOCentrality)
      {
        CentFT0CMin = 0;
        CentFT0CMax = 90;
      }
    }
    else
    {
      CentFT0CMin = CentFT0C[cent];
      CentFT0CMax = CentFT0C[cent + 1];
      if (isOOCentrality)
      {
        CentFT0CMin = CentFT0CLambdaOO[cent];
        CentFT0CMax = CentFT0CLambdaOO[cent + 1];
      }
    }
    cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << "%" << endl;
    hQT0CTPCA->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    hQT0CTPCC->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    hQTPCAC->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    hQT0CV0A->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    hQV0ATPCA->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    hQV0ATPCC->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    TH1F *h1QT0CTPCA = (TH1F *)hQT0CTPCA->ProjectionX("h1QT0CTPCA");
    TH1F *h1QT0CTPCC = (TH1F *)hQT0CTPCC->ProjectionX("h1QT0CTPCC");
    TH1F *h1QTPCAC = (TH1F *)hQTPCAC->ProjectionX("h1QTPCAC");
    TH1F *h1QT0CV0A = (TH1F *)hQT0CV0A->ProjectionX("h1QT0CV0A");
    TH1F *h1QV0ATPCA = (TH1F *)hQV0ATPCA->ProjectionX("h1QV0ATPCA");
    TH1F *h1QV0ATPCC = (TH1F *)hQV0ATPCC->ProjectionX("h1QV0ATPCC");
    Float_t MeanT0CTPCA = h1QT0CTPCA->GetMean();
    Float_t MeanT0CTPCC = h1QT0CTPCC->GetMean();
    Float_t MeanTPCAC = h1QTPCAC->GetMean();
    Float_t MeanT0CV0A = h1QT0CV0A->GetMean();
    Float_t MeanV0ATPCA = h1QV0ATPCA->GetMean();
    Float_t MeanV0ATPCC = h1QV0ATPCC->GetMean();
    Float_t ErrMeanT0CTPCA = h1QT0CTPCA->GetMeanError();
    Float_t ErrMeanT0CTPCC = h1QT0CTPCC->GetMeanError();
    Float_t ErrMeanTPCAC = h1QTPCAC->GetMeanError();
    Float_t ErrMeanT0CV0A = h1QT0CV0A->GetMeanError();
    Float_t ErrMeanV0ATPCA = h1QV0ATPCA->GetMeanError();
    Float_t ErrMeanV0ATPCC = h1QV0ATPCC->GetMeanError();
    cout << "MeanT0CTPCA: " << MeanT0CTPCA << " MeanT0CTPCC: " << MeanT0CTPCC << " MeanTPCAC: " << MeanTPCAC << " Reso: " << sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC) << endl;
    cout << "MeanT0CV0A: " << MeanT0CV0A << " MeanV0ATPCA: " << MeanV0ATPCA << " MeanV0ATPCC: " << MeanV0ATPCC << "Reso: " << sqrt(MeanT0CV0A * MeanV0ATPCA / MeanV0ATPCA) << endl;
    // cout << "Sourav resolution (EP): " << ftcReso[cent] << endl;
    Float_t RelErr2 = pow(ErrMeanT0CTPCA / MeanT0CTPCA, 2) + pow(ErrMeanT0CTPCC / MeanT0CTPCC, 2) + pow(ErrMeanTPCAC / MeanTPCAC, 2);
    Float_t ErrReso = sqrt(RelErr2 * (MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
    Float_t RelErr2V0A = pow(ErrMeanV0ATPCA / MeanV0ATPCA, 2) + pow(ErrMeanT0CTPCA / MeanT0CTPCA, 2) + pow(ErrMeanT0CV0A / MeanT0CV0A, 2);
    Float_t ErrResoV0A = sqrt(RelErr2V0A * (MeanT0CV0A * MeanT0CTPCA / MeanV0ATPCA));
    if (cent == commonnumCent)
    {
      hReso080->SetBinContent(1, sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
      hReso080->SetBinError(1, ErrReso);
    }
    else
    {
      hReso->SetBinContent(cent + 1, sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
      hReso->SetBinError(cent + 1, ErrReso);
      hResoV0A->SetBinContent(cent + 1, sqrt(MeanT0CV0A * MeanT0CTPCA / MeanTPCAC));
      hResoV0A->SetBinError(cent + 1, ErrResoV0A);
    }
  }
  // per cent bins
  for (Int_t cent = 0; cent < 100; cent++)
  {
    hQT0CTPCA->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    hQT0CTPCC->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    hQTPCAC->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    hQT0CV0A->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    TH1F *h1QT0CTPCA = (TH1F *)hQT0CTPCA->ProjectionX("h1QT0CTPCA");
    TH1F *h1QT0CTPCC = (TH1F *)hQT0CTPCC->ProjectionX("h1QT0CTPCC");
    TH1F *h1QTPCAC = (TH1F *)hQTPCAC->ProjectionX("h1QTPCAC");
    TH1F *h1QT0CV0A = (TH1F *)hQT0CV0A->ProjectionX("h1QT0CV0A");
    Float_t MeanT0CTPCA = h1QT0CTPCA->GetMean();
    Float_t MeanT0CTPCC = h1QT0CTPCC->GetMean();
    Float_t MeanTPCAC = h1QTPCAC->GetMean();
    Float_t MeanT0CV0A = h1QT0CV0A->GetMean();
    Float_t MeanV0ATPC = h1QT0CV0A->GetMean();
    Float_t ErrMeanT0CTPCA = h1QT0CTPCA->GetMeanError();
    Float_t ErrMeanT0CTPCC = h1QT0CTPCC->GetMeanError();
    Float_t ErrMeanTPCAC = h1QTPCAC->GetMeanError();
    Float_t ErrMeanT0CV0A = h1QT0CV0A->GetMeanError();
    Float_t ErrMeanV0ATPC = h1QT0CV0A->GetMeanError();
    Float_t RelErr2 = pow(ErrMeanT0CTPCA / MeanT0CTPCA, 2) + pow(ErrMeanT0CTPCC / MeanT0CTPCC, 2) + pow(ErrMeanTPCAC / MeanTPCAC, 2);
    Float_t ErrReso = sqrt(RelErr2 * (MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
    Float_t ErrResoV0A = sqrt(RelErr2 * (MeanT0CV0A * MeanT0CTPCC / MeanTPCAC));
    hResoPerCentBins->SetBinContent(cent + 1, sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
    hResoPerCentBinsV0A->SetBinContent(cent + 1, sqrt(MeanV0ATPC * MeanT0CTPCC / MeanTPCAC));
    cout << "Cent: " << cent << " Reso: " << sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC) << endl;
    cout << MeanT0CTPCA << " " << MeanT0CTPCC << " " << MeanTPCAC << " " << ErrMeanT0CTPCA << " " << ErrMeanT0CTPCC << " " << ErrMeanTPCAC << endl;
    hResoPerCentBins->SetBinError(cent + 1, ErrReso);
    hResoPerCentBinsV0A->SetBinError(cent + 1, ErrResoV0A);
  }

  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  StyleCanvas(canvas, 0.05, 0.1, 0.1, 0.05);
  hReso->GetYaxis()->SetRangeUser(0., 1.5);
  if (isOOCentrality)
    hReso->GetYaxis()->SetRangeUser(0., 0.65);
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
  hResoV0A->SetLineColor(kGreen + 2);
  hResoV0A->SetMarkerColor(kGreen + 2);
  hResoV0A->SetMarkerStyle(23);
  hResoV0A->SetMarkerSize(1.5);
  hResoV0A->Draw("same");
  hResoPerCentBins->SetLineColor(kBlue);
  hResoPerCentBins->SetMarkerColor(kBlue);
  hResoPerCentBins->SetMarkerStyle(21);
  hResoPerCentBins->SetMarkerSize(0.8);
  //hResoPerCentBins->Draw("same");
  hResoPerCentBinsV0A->SetLineColor(kRed);
  hResoPerCentBinsV0A->SetMarkerColor(kRed);
  hResoPerCentBinsV0A->SetMarkerStyle(22);
  hResoPerCentBinsV0A->SetMarkerSize(0.8);
  hResoPerCentBinsV0A->Draw("same");

  TH1F *hResoSourav = (TH1F *)hReso->Clone("hResoSourav");
  hResoSourav->SetLineColor(kRed);
  hResoSourav->SetMarkerColor(kRed);
  for (Int_t cent = 0; cent < commonnumCent; cent++)
  {
    // hResoSourav->SetBinContent(cent + 1, ftcResoSourav[cent]);
    // hResoSourav->SetBinError(cent + 1, 0);
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
  hResoPerCentBins->Write();
  outputfile->Close();
  cout << "\nOutputFile: " << Soutputfile << endl;
}
