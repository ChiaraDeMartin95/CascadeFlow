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
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->SetTitle(title);
}

TString NameRun2[4] = {"10-20 %", "20-30 %", "30-40 %", "40-50 %"};
Int_t Color[12] = {kYellow + 2, kGreen + 2, kBlack, kBlue + 1, kYellow + 2, kGreen + 2, kBlack, kBlue + 1, kYellow + 2, kGreen + 2, kBlack, kBlue + 1};
Int_t Marker[4] = {24, 25, 30, 27};
Int_t MarkerRun3[12] = {20, 21, 29, 33, 20, 21, 29, 33, 20, 21, 29, 33};
Float_t MarkerSize[12] = {1.5, 1.5, 1.5, 2., 1.5, 1.5, 1.5, 2., 1.5, 1.5, 1.5, 2.};

void CompareWPublished(Bool_t isXi = ChosenParticleXi,
                       TString inputFileName = SinputFileName,
                       Bool_t UseTwoGauss = ExtrUseTwoGauss,
                       Int_t BkgType = ExtrBkgType)
{

  TString SinputFile = "";
  TFile *inputFile;
  TH1F *hV2C[numCent];

  for (Int_t cent = 0; cent < numCent; cent++)
  //for (Int_t cent = 0; cent < 4; cent++)
  {
    SinputFile = "OutputAnalysis/FitV2_" + inputFileName + "_" + ParticleName[!isXi];
    SinputFile += IsOneOrTwoGauss[UseTwoGauss];
    SinputFile += SIsBkgParab[BkgType];
    SinputFile += Form("_Cent%.i-%.i.root", CentFT0C[cent], CentFT0C[cent + 1]);
    inputFile = new TFile(SinputFile);
    cout << "Input file with Run 3 v2: " << SinputFile << endl;

    hV2C[cent] = (TH1F *)inputFile->Get("histoV2");
    hV2C[cent]->SetName(Form("hV2C_%i", cent));
    hV2C[cent]->SetMarkerStyle(MarkerRun3[cent]);
    hV2C[cent]->SetMarkerSize(MarkerSize[cent]);
    hV2C[cent]->SetMarkerColor(ColorMult[cent]);
    hV2C[cent]->SetLineColor(ColorMult[cent]);
    hV2C[cent]->SetLineWidth(2);
  }

  TString SPublishedFileRun2 = "Run2Results/HEPData-ins2093750-v1-root.root";
  TFile *PublishedFileRun2 = new TFile(SPublishedFileRun2);
  TString Tables[4] = {"Table 8", "Table 17", "Table 26", "Table 35"};
  if (!isXi)
  {
    Tables[0] = "Table 9";
    Tables[1] = "Table 18";
    Tables[2] = "Table 27";
    Tables[3] = "Table 36";
  }

  TDirectoryFile *dirPublished;
  TGraphErrors *gV2Run2[4];
  for (Int_t i = 0; i < 4; i++)
  {
    dirPublished = (TDirectoryFile *)PublishedFileRun2->Get(Tables[i]);
    gV2Run2[i] = (TGraphErrors *)dirPublished->Get("Graph1D_y1");
    gV2Run2[i]->SetName(NameRun2[i]);
    gV2Run2[i]->SetMarkerStyle(Marker[i]);
    gV2Run2[i]->SetMarkerSize(MarkerSize[i]);
    gV2Run2[i]->SetMarkerColor(Color[i]);
    gV2Run2[i]->SetLineColor(Color[i]);
    gV2Run2[i]->SetLineWidth(2);
  }

  TLegend *legendRun3 = new TLegend(0.13, 0.6, 0.4, 0.88);
  legendRun3->SetBorderSize(0);
  legendRun3->SetFillStyle(0);
  legendRun3->SetTextSize(0.03);
  legendRun3->SetHeader("Run 3");

  TLegend *legendRun2 = new TLegend(0.71, 0.12, 0.81, 0.30);
  legendRun2->SetBorderSize(0);
  legendRun2->SetFillStyle(0);
  legendRun2->SetTextSize(0.03);
  legendRun2->SetHeader("Run 2 published");

  TH1F *histoDummy = new TH1F("histoDummy", "histoDummy", 100, 0, 6);
  StyleHisto(histoDummy, -0.2, 0.8, 1, 20, "#it{p}_{T} (GeV/#it{c})", "v_{2}", "", kFALSE, 0, 100, 1, 1, 0.05);
  TFile *file = new TFile("OutputAnalysis/CompareWPublished_" + inputFileName + ".root", "RECREATE");
  TCanvas *canvasvsRun2 = new TCanvas("canvasvsRun2", "canvasvsRun2", 1200, 800);
  canvasvsRun2->cd();
  histoDummy->Draw("");
  
  for (Int_t i = 0; i < numCent; i++)
  {
    legendRun3->AddEntry(hV2C[i], Form("%i-%i", CentFT0C[i], CentFT0C[i+1]), "P");
    if (i<4) {
      gV2Run2[i]->DrawClone("same P");
      legendRun2->AddEntry(gV2Run2[i], NameRun2[i], "P");
    }
    hV2C[i]->Draw("same");
  }
  legendRun3->Draw();
  legendRun2->Draw();
  canvasvsRun2->Write();
  file->Close();
  cout << "I created the file " << file->GetName() << endl;
}
