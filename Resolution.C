#include "Riostream.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TH3F.h>
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TProfile.h"
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TFile.h>
// #include "CommonVar.h"
#include "CommonVarLambda.h"

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

// isSPReso:
// 1 --> resolution for SP method
// 0 --> resolution for EP method

void Resolution(Bool_t isSPReso = 1, Bool_t isLFReso = 1)
{

  TString ComputeResoFileName = "../TreeForAnalysis/AnalysisResults_";
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
  TString nameQT0CT0A = "QVectorsT0CT0A";
  TString nameQT0ATPCA = "QVectorsT0ATPCA";
  TString nameQT0ATPCC = "QVectorsT0ATPCC";
  if (!isSPReso)
  {
    // these are probably normalised in the wrong way if CFW vectors are used!
    // nameQT0CTPCA = "QVectorsNormT0CTPCA";
    // nameQT0CTPCC = "QVectorsNormT0CTPCC";
    // nameQTPCAC = "QVectorsNormTPCAC";
    // nameQT0CV0A = "QVectorsNormT0CV0A";
    // nameQV0ATPCA = "QVectorsNormV0ATPCA";
    // nameQV0ATPCC = "QVectorsNormV0ATPCC";
    // the following names are to compute resolution AFTER shift correction (the right thing to do!)
    nameQT0CTPCA = "QVectorsT0CTPCA_Shifted";
    nameQT0CTPCC = "QVectorsT0CTPCC_Shifted";
    nameQTPCAC = "QVectorsTPCAC_Shifted";
    nameQT0CV0A = "QVectorsT0CV0A_Shifted";
    nameQV0ATPCA = "QVectorsV0ATPCA_Shifted";
    nameQV0ATPCC = "QVectorsV0ATPCC_Shifted";
    nameQT0CT0A = "QVectorsT0CT0A_Shifted";
    nameQT0ATPCA = "QVectorsT0ATPCA_Shifted";
    nameQT0ATPCC = "QVectorsT0ATPCC_Shifted";
    // the following names are to compute resolution BEFORE shift correction (not ideal)
    // nameQT0CTPCA = "EP_T0CTPCA";
    // nameQT0CTPCC = "EP_T0CTPCC";
    // nameQTPCAC = "EP_TPCAC";
    // nameQT0CV0A = "EP_T0CV0A";
    // nameQV0ATPCA = "EP_V0ATPCA";
    // nameQV0ATPCC = "EP_V0ATPCC";
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
  TH2F *hQT0CT0A = (TH2F *)dir->Get(nameQT0CT0A);
  TH2F *hQT0ATPCA = (TH2F *)dir->Get(nameQT0ATPCA);
  TH2F *hQT0ATPCC = (TH2F *)dir->Get(nameQT0ATPCC);
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
  if (!hQT0CT0A)
  {
    hQT0CT0A = (TH2F *)dir->Get(nameQT0CTPCA);
    cout << "Error: input histograms T0CT0A not found, taking T0CTPCA instead" << endl;
  }
  if (!hQT0ATPCA)
  {
    hQT0ATPCA = (TH2F *)dir->Get(nameQT0CTPCC);
    cout << "Error: input histograms T0ATPCA not found, taking T0CTPCC instead" << endl;
  }
  if (!hQT0ATPCC)
  {
    hQT0ATPCC = (TH2F *)dir->Get(nameQTPCAC);
    cout << "Error: input histograms T0ATPCC not found, taking TPCAC instead" << endl;
  }

  TH1D *hReso = new TH1D("hReso", "hReso", numCent, fCentFT0C);
  if (isOOCentrality)
    hReso = new TH1D("hReso", "hReso", numCentLambdaOO, fCentFT0CLambdaOO);
  TH1D *hResoV0ATPCA = (TH1D *)hReso->Clone("hResoV0ATPCA");
  TH1D *hResoV0ATPCC = (TH1D *)hReso->Clone("hResoV0ATPCC");
  TH1D *hResoT0ATPCA = (TH1D *)hReso->Clone("hResoT0ATPCA");
  TH1D *hResoT0ATPCC = (TH1D *)hReso->Clone("hResoT0ATPCC");
  TH1D *hReso2SubEvents = (TH1D *)hReso->Clone("hReso2SubEvents");
  TH1D *hReso2SubEventsT0A = (TH1D *)hReso->Clone("hReso2SubEventsT0A");
  TH1D *hReso2T0CTPCA = (TH1D *)hReso->Clone("hReso2T0CTPCA");
  TH1D *hReso2T0CTPCC = (TH1D *)hReso->Clone("hReso2T0CTPCC");
  TH1D *hReso2TPCAC = (TH1D *)hReso->Clone("hReso2TPCAC");
  TH1D *hReso2V0ATPCA = (TH1D *)hReso->Clone("hReso2V0ATPCA");
  TH1D *hReso2V0ATPCC = (TH1D *)hReso->Clone("hReso2V0ATPCC");
  TH1D *hReso2T0ATPCA = (TH1D *)hReso->Clone("hReso2T0ATPCA");
  TH1D *hReso2T0ATPCC = (TH1D *)hReso->Clone("hReso2T0ATPCC");
  TH1D *hResoPerCentBins = new TH1D("hResoPerCentBins", "hResoPerCentBins", 100, 0, 100);
  TH1D *hResoPerCentBinsV0A = new TH1D("hResoPerCentBinsV0A", "hResoPerCentBinsV0A", 100, 0, 100);
  TH1D *hResoPerCentBinsT0ATPCC = new TH1D("hResoPerCentBinsT0ATPCC", "hResoPerCentBinsT0ATPCC", 100, 0, 100);
  TH1D *hReso080 = new TH1D("hReso080", "hReso080", 1, 0, 1);
  TH1D *hResoV0ATPCA080 = new TH1D("hResoV0ATPCA080", "hResoV0ATPCA080", 1, 0, 1);
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
        CentFT0CMax = CentFT0CMaxLambdaOO;
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
    hQT0CT0A->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    hQV0ATPCA->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    hQV0ATPCC->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    hQT0ATPCA->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    hQT0ATPCC->GetYaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    TH1F *h1QT0CTPCA = (TH1F *)hQT0CTPCA->ProjectionX("h1QT0CTPCA");
    TH1F *h1QT0CTPCC = (TH1F *)hQT0CTPCC->ProjectionX("h1QT0CTPCC");
    TH1F *h1QTPCAC = (TH1F *)hQTPCAC->ProjectionX("h1QTPCAC");
    TH1F *h1QT0CV0A = (TH1F *)hQT0CV0A->ProjectionX("h1QT0CV0A");
    TH1F *h1QT0CT0A = (TH1F *)hQT0CT0A->ProjectionX("h1QT0CT0A");
    TH1F *h1QV0ATPCA = (TH1F *)hQV0ATPCA->ProjectionX("h1QV0ATPCA");
    TH1F *h1QV0ATPCC = (TH1F *)hQV0ATPCC->ProjectionX("h1QV0ATPCC");
    TH1F *h1QT0ATPCA = (TH1F *)hQT0ATPCA->ProjectionX("h1QT0ATPCA");
    TH1F *h1QT0ATPCC = (TH1F *)hQT0ATPCC->ProjectionX("h1QT0ATPCC");

    Float_t MeanT0CTPCA = h1QT0CTPCA->GetMean();
    Float_t MeanT0CTPCC = h1QT0CTPCC->GetMean();
    Float_t MeanTPCAC = h1QTPCAC->GetMean();
    Float_t MeanT0CV0A = h1QT0CV0A->GetMean();
    Float_t MeanT0CT0A = h1QT0CT0A->GetMean();
    Float_t MeanV0ATPCA = h1QV0ATPCA->GetMean();
    Float_t MeanV0ATPCC = h1QV0ATPCC->GetMean();
    Float_t MeanT0ATPCA = h1QT0ATPCA->GetMean();
    Float_t MeanT0ATPCC = h1QT0ATPCC->GetMean();

    Float_t ErrMeanT0CTPCA = h1QT0CTPCA->GetMeanError();
    Float_t ErrMeanT0CTPCC = h1QT0CTPCC->GetMeanError();
    Float_t ErrMeanTPCAC = h1QTPCAC->GetMeanError();
    Float_t ErrMeanT0CV0A = h1QT0CV0A->GetMeanError();
    Float_t ErrMeanT0CT0A = h1QT0CT0A->GetMeanError();
    Float_t ErrMeanV0ATPCA = h1QV0ATPCA->GetMeanError();
    Float_t ErrMeanV0ATPCC = h1QV0ATPCC->GetMeanError();
    Float_t ErrMeanT0ATPCA = h1QT0ATPCA->GetMeanError();
    Float_t ErrMeanT0ATPCC = h1QT0ATPCC->GetMeanError();
    cout << "MeanT0CTPCA: " << MeanT0CTPCA << " MeanT0CTPCC: " << MeanT0CTPCC << " MeanTPCAC: " << MeanTPCAC << " Reso: " << sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC) << endl;
    cout << "MeanT0CV0A: " << MeanT0CV0A << " MeanV0ATPCA: " << MeanV0ATPCA << " MeanV0ATPCC: " << MeanV0ATPCC << "Reso: " << sqrt(MeanT0CV0A * MeanV0ATPCA / MeanV0ATPCA) << endl;
    // cout << "Sourav resolution (EP): " << ftcReso[cent] << endl;
    Float_t RelErr2 = pow(ErrMeanT0CTPCA / MeanT0CTPCA, 2) + pow(ErrMeanT0CTPCC / MeanT0CTPCC, 2) + pow(ErrMeanTPCAC / MeanTPCAC, 2);
    Float_t ErrReso = sqrt(RelErr2 * (MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
    Float_t RelErr2V0ATPCA = pow(ErrMeanV0ATPCA / MeanV0ATPCA, 2) + pow(ErrMeanT0CTPCA / MeanT0CTPCA, 2) + pow(ErrMeanT0CV0A / MeanT0CV0A, 2);
    Float_t RelErr2V0ATPCC = pow(ErrMeanV0ATPCC / MeanV0ATPCC, 2) + pow(ErrMeanT0CTPCC / MeanT0CTPCC, 2) + pow(ErrMeanT0CV0A / MeanT0CV0A, 2);
    Float_t RelErr2T0ATPCA = pow(ErrMeanT0ATPCA / MeanT0ATPCA, 2) + pow(ErrMeanT0CTPCA / MeanT0CTPCA, 2) + pow(ErrMeanT0CT0A / MeanT0CT0A, 2);
    Float_t RelErr2T0ATPCC = pow(ErrMeanT0ATPCC / MeanT0ATPCC, 2) + pow(ErrMeanT0CTPCC / MeanT0CTPCC, 2) + pow(ErrMeanT0CT0A / MeanT0CT0A, 2);
    Float_t ErrResoV0ATPCA = sqrt(RelErr2V0ATPCA * (MeanT0CV0A * MeanT0CTPCA / MeanV0ATPCA));
    Float_t ErrResoV0ATPCC = sqrt(RelErr2V0ATPCC * (MeanT0CV0A * MeanT0CTPCC / MeanV0ATPCC));
    Float_t ErrResoT0ATPCA = sqrt(RelErr2T0ATPCA * (MeanT0CT0A * MeanT0CTPCA / MeanT0ATPCA));
    Float_t ErrResoT0ATPCC = sqrt(RelErr2T0ATPCC * (MeanT0CT0A * MeanT0CTPCC / MeanT0ATPCC));

    if (cent == commonnumCent)
    {
      hReso080->SetBinContent(1, sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
      hReso080->SetBinError(1, ErrReso);
      hResoV0ATPCA080->SetBinContent(1, sqrt(MeanT0CV0A * MeanT0CTPCA / MeanV0ATPCA));
      hResoV0ATPCA080->SetBinError(1, ErrResoV0ATPCA);
    }
    else
    {
      hReso->SetBinContent(cent + 1, sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
      hReso->SetBinError(cent + 1, ErrReso);
      hReso2SubEvents->SetBinContent(cent + 1, sqrt(MeanT0CV0A));
      hReso2SubEvents->SetBinError(cent + 1, ErrMeanT0CV0A / 2 / sqrt(MeanT0CV0A));

      hReso2SubEventsT0A->SetBinContent(cent + 1, sqrt(MeanT0CT0A));
      hReso2SubEventsT0A->SetBinError(cent + 1, ErrMeanT0CT0A / 2 / sqrt(MeanT0CT0A));

      hResoV0ATPCA->SetBinContent(cent + 1, sqrt(MeanT0CV0A * MeanT0CTPCA / MeanV0ATPCA));
      hResoV0ATPCA->SetBinError(cent + 1, ErrResoV0ATPCA);
      hResoV0ATPCC->SetBinContent(cent + 1, sqrt(MeanT0CV0A * MeanT0CTPCC / MeanV0ATPCC));
      hResoV0ATPCC->SetBinError(cent + 1, ErrResoV0ATPCC);

      hResoT0ATPCA->SetBinContent(cent + 1, sqrt(MeanT0CT0A * MeanT0CTPCA / MeanT0ATPCA));
      hResoT0ATPCA->SetBinError(cent + 1, ErrResoT0ATPCA);
      hResoT0ATPCC->SetBinContent(cent + 1, sqrt(MeanT0CT0A * MeanT0CTPCC / MeanV0ATPCC));
      hResoT0ATPCC->SetBinError(cent + 1, ErrResoT0ATPCC);

      hReso2T0CTPCA->SetBinContent(cent + 1, sqrt(MeanT0CTPCA));
      hReso2T0CTPCA->SetBinError(cent + 1, ErrMeanT0CTPCA / 2 / sqrt(MeanT0CTPCA));
      hReso2T0CTPCC->SetBinContent(cent + 1, sqrt(MeanT0CTPCC));
      hReso2T0CTPCC->SetBinError(cent + 1, ErrMeanT0CTPCC / 2 / sqrt(MeanT0CTPCC));
      hReso2TPCAC->SetBinContent(cent + 1, sqrt(MeanTPCAC));
      hReso2TPCAC->SetBinError(cent + 1, ErrMeanTPCAC / 2 / sqrt(MeanTPCAC));
      hReso2V0ATPCA->SetBinContent(cent + 1, sqrt(MeanV0ATPCA));
      hReso2V0ATPCA->SetBinError(cent + 1, ErrMeanV0ATPCA / 2 / sqrt(MeanV0ATPCA));
      hReso2V0ATPCC->SetBinContent(cent + 1, sqrt(MeanV0ATPCC));
      hReso2V0ATPCC->SetBinError(cent + 1, ErrMeanV0ATPCC / 2 / sqrt(MeanV0ATPCC));
      hReso2T0ATPCA->SetBinContent(cent + 1, sqrt(MeanT0ATPCA));
      hReso2T0ATPCA->SetBinError(cent + 1, ErrMeanT0ATPCA / 2 / sqrt(MeanT0ATPCA));
      hReso2T0ATPCC->SetBinContent(cent + 1, sqrt(MeanT0ATPCC));
      hReso2T0ATPCC->SetBinError(cent + 1, ErrMeanT0ATPCC / 2 / sqrt(MeanT0ATPCC));
    }
  }
  // per cent bins
  for (Int_t cent = 0; cent < 100; cent++)
  {
    hQT0CTPCA->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    hQT0CTPCC->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    hQTPCAC->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    hQT0CV0A->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    hQV0ATPCA->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    hQV0ATPCC->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    hQT0ATPCA->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    hQT0ATPCC->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    hQT0CT0A->GetYaxis()->SetRangeUser(cent + 0.001, cent + 0.001);
    TH1F *h1QT0CTPCA = (TH1F *)hQT0CTPCA->ProjectionX("h1QT0CTPCA");
    TH1F *h1QT0CTPCC = (TH1F *)hQT0CTPCC->ProjectionX("h1QT0CTPCC");
    TH1F *h1QTPCAC = (TH1F *)hQTPCAC->ProjectionX("h1QTPCAC");
    TH1F *h1QT0CV0A = (TH1F *)hQT0CV0A->ProjectionX("h1QT0CV0A");
    TH1F *h1QV0ATPCA = (TH1F *)hQV0ATPCA->ProjectionX("h1QV0ATPCA");
    TH1F *h1QV0ATPCC = (TH1F *)hQV0ATPCC->ProjectionX("h1QV0ATPCC");
    TH1F *h1QT0ATPCA = (TH1F *)hQT0ATPCA->ProjectionX("h1QT0ATPCA");
    TH1F *h1QT0ATPCC = (TH1F *)hQT0ATPCC->ProjectionX("h1QT0ATPCC");
    TH1F *h1QT0CT0A = (TH1F *)hQT0CT0A->ProjectionX("h1QT0CT0A");
    Float_t MeanT0CTPCA = h1QT0CTPCA->GetMean();
    Float_t MeanT0CTPCC = h1QT0CTPCC->GetMean();
    Float_t MeanTPCAC = h1QTPCAC->GetMean();
    Float_t MeanT0CV0A = h1QT0CV0A->GetMean();
    Float_t MeanV0ATPCA = h1QV0ATPCA->GetMean();
    Float_t MeanV0ATPCC = h1QV0ATPCC->GetMean();
    Float_t MeanT0CT0A = h1QT0CT0A->GetMean();
    Float_t MeanT0ATPCA = h1QT0ATPCA->GetMean();
    Float_t MeanT0ATPCC = h1QT0ATPCC->GetMean();
    Float_t ErrMeanT0CTPCA = h1QT0CTPCA->GetMeanError();
    Float_t ErrMeanT0CTPCC = h1QT0CTPCC->GetMeanError();
    Float_t ErrMeanTPCAC = h1QTPCAC->GetMeanError();
    Float_t ErrMeanT0CV0A = h1QT0CV0A->GetMeanError();
    Float_t ErrMeanV0ATPCA = h1QV0ATPCA->GetMeanError();
    Float_t ErrMeanV0ATPCC = h1QV0ATPCC->GetMeanError();
    Float_t ErrMeanT0CT0A = h1QT0CT0A->GetMeanError();
    Float_t ErrMeanT0ATPCA = h1QT0ATPCA->GetMeanError();
    Float_t ErrMeanT0ATPCC = h1QT0ATPCC->GetMeanError();
    Float_t RelErr2 = pow(ErrMeanT0CTPCA / MeanT0CTPCA, 2) + pow(ErrMeanT0CTPCC / MeanT0CTPCC, 2) + pow(ErrMeanTPCAC / MeanTPCAC, 2);
    Float_t ErrReso = sqrt(RelErr2 * (MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
    Float_t RelErr2V0A = pow(ErrMeanV0ATPCA / MeanV0ATPCA, 2) + pow(ErrMeanT0CTPCA / MeanT0CTPCA, 2) + pow(ErrMeanT0CV0A / MeanT0CV0A, 2);
    Float_t ErrResoV0A = sqrt(RelErr2V0A * (MeanT0CV0A * MeanT0CTPCA / MeanV0ATPCA));
    Float_t RelErr2T0ATPCC = pow(ErrMeanT0ATPCC / MeanT0ATPCC, 2) + pow(ErrMeanT0CTPCC / MeanT0CTPCC, 2) + pow(ErrMeanT0CT0A / MeanT0CT0A, 2);
    Float_t ErrResoT0ATPCC = sqrt(RelErr2T0ATPCC * (MeanT0CT0A * MeanT0CTPCC / MeanT0ATPCC));

    hResoPerCentBins->SetBinContent(cent + 1, sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC));
    hResoPerCentBinsV0A->SetBinContent(cent + 1, sqrt(MeanT0CV0A * MeanT0CTPCA / MeanV0ATPCA));
    hResoPerCentBinsT0ATPCC->SetBinContent(cent + 1, sqrt(MeanT0CT0A * MeanT0CTPCC / MeanT0ATPCC));
    cout << "Cent: " << cent << " Reso: " << sqrt(MeanT0CTPCA * MeanT0CTPCC / MeanTPCAC) << endl;
    cout << MeanT0CTPCA << " " << MeanT0CTPCC << " " << MeanTPCAC << " " << ErrMeanT0CTPCA << " " << ErrMeanT0CTPCC << " " << ErrMeanTPCAC << endl;
    cout << "hResoPerCentBinsV0A: " << hResoPerCentBinsV0A->GetBinContent(cent + 1) << endl;
    cout << MeanV0ATPCA << " " << MeanT0CV0A << " " << MeanT0CTPCA << " " << ErrMeanV0ATPCA << " " << ErrMeanT0CV0A << " " << ErrMeanT0CTPCA << endl;
    cout << "hResoPerCentBins: " << hResoPerCentBins->GetBinContent(cent + 1) << endl;
    if (MeanT0CV0A < 0)
    {
      hResoPerCentBinsV0A->SetBinContent(cent + 1, 0);
      ErrResoV0A = 0;
    }
    hResoPerCentBins->SetBinError(cent + 1, ErrReso);
    hResoPerCentBinsV0A->SetBinError(cent + 1, ErrResoV0A);
    hResoPerCentBinsT0ATPCC->SetBinError(cent + 1, ErrResoT0ATPCC);
  }

  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  StyleCanvas(canvas, 0.05, 0.1, 0.1, 0.05);
  hReso->GetXaxis()->SetRangeUser(0., 80.);
  hReso->GetYaxis()->SetRangeUser(0., 1.5);
  if (isOOCentrality)
  {
    hReso->GetYaxis()->SetRangeUser(0., 0.65);
    hReso->GetXaxis()->SetRangeUser(0., 100.);
  }
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
  hReso2SubEvents->GetXaxis()->SetTitle("Centrality (%)");
  hReso2SubEvents->GetYaxis()->SetTitle("Resolution");
  hReso2SubEvents->GetXaxis()->SetTitleSize(0.05);
  hReso2SubEvents->GetXaxis()->SetLabelSize(0.05);
  hReso2SubEvents->GetXaxis()->SetTitleOffset(1.);
  hReso2SubEvents->GetYaxis()->SetTitleSize(0.05);
  hReso2SubEvents->GetYaxis()->SetLabelSize(0.05);
  hReso2SubEvents->GetYaxis()->SetTitleOffset(1.);
  hReso2SubEvents->SetTitle("");
  if (isOOCentrality)
  {
    hReso2SubEvents->GetYaxis()->SetRangeUser(0., 0.4);
    hReso2SubEvents->GetXaxis()->SetRangeUser(0., 100.);
  }
  else
  {
    hReso2SubEvents->GetYaxis()->SetRangeUser(0., 1.5);
    hReso2SubEvents->GetXaxis()->SetRangeUser(0., 80.);
  }
  hReso->Draw("");
  hReso2SubEvents->SetLineColor(kOrange + 7);
  hReso2SubEvents->SetMarkerColor(kOrange + 7);
  hReso2SubEvents->SetMarkerStyle(25);
  hReso2SubEvents->SetMarkerSize(1.5);
  // hReso2SubEvents->Draw("same");
  hReso2SubEventsT0A->SetLineColor(kOrange + 3);
  hReso2SubEventsT0A->SetMarkerColor(kOrange + 3);
  hReso2SubEventsT0A->SetMarkerStyle(28);
  hReso2SubEventsT0A->SetMarkerSize(1.5);
  // hReso2SubEventsT0A->Draw("same");
  hResoV0ATPCA->SetLineColor(kGreen + 2);
  hResoV0ATPCA->SetMarkerColor(kGreen + 2);
  hResoV0ATPCA->SetMarkerStyle(23);
  hResoV0ATPCA->SetMarkerSize(1.5);
  // hResoV0ATPCA->Draw("same");
  hResoV0ATPCC->SetLineColor(kGreen + 3);
  hResoV0ATPCC->SetMarkerColor(kGreen + 3);
  hResoV0ATPCC->SetMarkerStyle(24);
  hResoV0ATPCC->SetMarkerSize(1.5);
  // hResoV0ATPCC->Draw("same");
  hResoT0ATPCA->SetLineColor(kRed + 1);
  hResoT0ATPCA->SetMarkerColor(kRed + 1);
  hResoT0ATPCA->SetMarkerStyle(23);
  hResoT0ATPCA->SetMarkerSize(1.5);
  // hResoT0ATPCA->Draw("same");
  hResoT0ATPCC->SetLineColor(kRed + 2);
  hResoT0ATPCC->SetMarkerColor(kRed + 2);
  hResoT0ATPCC->SetMarkerStyle(24);
  hResoT0ATPCC->SetMarkerSize(1.5);
  // hResoT0ATPCC->Draw("same");
  hResoPerCentBins->SetLineColor(kBlue);
  hResoPerCentBins->SetMarkerColor(kBlue);
  hResoPerCentBins->SetMarkerStyle(21);
  hResoPerCentBins->SetMarkerSize(0.8);
  // hResoPerCentBins->Draw("same");
  hResoPerCentBinsV0A->SetLineColor(kRed);
  hResoPerCentBinsV0A->SetMarkerColor(kRed);
  hResoPerCentBinsV0A->SetMarkerStyle(22);
  hResoPerCentBinsV0A->SetMarkerSize(0.8);
  // hResoPerCentBinsV0A->Draw("same");
  hReso2T0CTPCA->SetLineColor(kCyan + 2);
  hReso2T0CTPCA->SetMarkerColor(kCyan + 2);
  hReso2T0CTPCA->SetMarkerStyle(26);
  hReso2T0CTPCA->SetMarkerSize(1.5);
  // hReso2T0CTPCA->Draw("same");
  hReso2T0CTPCC->SetLineColor(kCyan + 3);
  hReso2T0CTPCC->SetMarkerColor(kCyan + 3);
  hReso2T0CTPCC->SetMarkerStyle(27);
  hReso2T0CTPCC->SetMarkerSize(1.5);
  // hReso2T0CTPCC->Draw("same");
  hReso2TPCAC->SetLineColor(kCyan + 4);
  hReso2TPCAC->SetMarkerColor(kCyan + 4);
  hReso2TPCAC->SetMarkerStyle(28);
  hReso2TPCAC->SetMarkerSize(1.5);
  // hReso2TPCAC->Draw("same");
  hReso2V0ATPCA->SetLineColor(kBlue + 2);
  hReso2V0ATPCA->SetMarkerColor(kBlue + 2);
  hReso2V0ATPCA->SetMarkerStyle(29);
  hReso2V0ATPCA->SetMarkerSize(1.5);
  // hReso2V0ATPCA->Draw("same");
  hReso2V0ATPCC->SetLineColor(kBlue + 3);
  hReso2V0ATPCC->SetMarkerColor(kBlue + 3);
  hReso2V0ATPCC->SetMarkerStyle(30);
  hReso2V0ATPCC->SetMarkerSize(1.5);
  // hReso2V0ATPCC->Draw("same");
  hReso2T0ATPCA->SetLineColor(kMagenta + 2);
  hReso2T0ATPCA->SetMarkerColor(kMagenta + 2);
  hReso2T0ATPCA->SetMarkerStyle(29);
  hReso2T0ATPCA->SetMarkerSize(1.5);
  // hReso2T0ATPCA->Draw("same");
  hReso2T0ATPCC->SetLineColor(kMagenta + 3);
  hReso2T0ATPCC->SetMarkerColor(kMagenta + 3);
  hReso2T0ATPCC->SetMarkerStyle(30);
  hReso2T0ATPCC->SetMarkerSize(1.5);
  // hReso2T0ATPCC->Draw("same");

  TLegend *legendRes = new TLegend(0.6, 0.6, 0.9, 0.93);
  legendRes->SetFillStyle(0);
  legendRes->SetMargin(0);
  legendRes->SetTextSize(0.03);
  legendRes->SetTextAlign(32);
  // legendRes->AddEntry(hReso, "EP T0C-TPC", "lp");
  legendRes->AddEntry(hResoV0ATPCA, "EP 3 sub-events T0C-V0A-TPCA", "lp");
  // legendRes->AddEntry(hResoV0ATPCC, "EP 3 sub-events T0C-V0A-TPCC", "lp");
  // legendRes->AddEntry(hResoT0ATPCA, "EP 3 sub-events T0C-T0A-TPCA", "lp");
  // legendRes->AddEntry(hResoT0ATPCC, "EP 3 sub-events T0C-T0A-TPCC", "lp");
  legendRes->AddEntry(hReso2SubEvents, "2 sub-events T0C-V0A", "lp");
  legendRes->AddEntry(hReso2SubEventsT0A, "2 sub-events T0C-T0A", "lp");
  // legendRes->AddEntry(hReso2T0CTPCA, "2 sub-events T0C-TPC A", "lp");
  // legendRes->AddEntry(hReso2T0CTPCC, "2 sub-events T0C-TPC C", "lp");
  // legendRes->AddEntry(hReso2TPCAC, "2 sub-events TPC A-C", "lp");
  legendRes->AddEntry(hReso2V0ATPCA, "2 sub-events V0A-TPC A", "lp");
  legendRes->AddEntry(hReso2V0ATPCC, "2 sub-events V0A-TPC C", "lp");
  legendRes->AddEntry(hReso2T0ATPCA, "2 sub-events T0A-TPC A", "lp");
  legendRes->AddEntry(hReso2T0ATPCC, "2 sub-events T0A-TPC C", "lp");
  // legendRes->Draw();

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

  TString Soutputfile = "../";
  if (!isSPReso)
  { // event plane method
    if (isLFReso)
      Soutputfile += ResoFileName_EPLF;
    else
      Soutputfile += ResoFileName_EPCFW;
  }
  else
  {
    if (isLFReso)
      Soutputfile += ResoFileName_SPLF;
    else
      Soutputfile += ResoFileName_SPCFW;
  }
  canvas->SaveAs(Soutputfile + ".pdf");
  canvas->SaveAs(Soutputfile + ".png");
  canvas->SaveAs(Soutputfile + ".eps");
  TFile *outputfile = new TFile(Soutputfile + ".root", "RECREATE");
  hReso->Write();
  hReso080->Write();
  hResoV0ATPCA080->Write();
  hResoPerCentBins->Write();
  hResoV0ATPCA->Write();
  hResoV0ATPCC->Write();
  hResoT0ATPCA->Write();
  hResoT0ATPCC->Write();
  hResoPerCentBinsV0A->Write();
  hResoPerCentBinsT0ATPCC->Write();
  TList *listReso = new TList();
  listReso->Add(hResoPerCentBinsV0A);
  listReso->Add(hResoPerCentBinsT0ATPCC);
  listReso->Write("ccdb_object", TObject::kSingleKey);
  outputfile->Close();
  cout << "\nOutputFile: " << Soutputfile << endl;
}
