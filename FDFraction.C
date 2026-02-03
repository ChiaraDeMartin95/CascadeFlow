#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "CommonVarLambda.h"

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

void FDFraction()
{

  TString PathInFD = "../TreeForAnalysis/AnalysisResults_" + SinputFileNameFDFraction + ".root";
  TFile *fileInFD = new TFile(PathInFD, "READ");
  TDirectoryFile *dirFD = (TDirectoryFile *)fileInFD->Get("lf-cascade-flow/histos");
  if (!dirFD)
  {
    cout << "Directory lf-cascade-flow/histos not found in file " << PathInFD << endl;
    return;
  }
  TH3F *hFDLambda = (TH3F *)dirFD->Get("hCentvsPtvsPrimaryFracLambda");
  if (!hFDLambda)
  {
    cout << "Histogram hCentvsPtvsPrimFracLambda not found in file " << PathInFD << endl;
    return;
  }

  TH1F *hFDFractionLambdavsCent = new TH1F("hFDFractionLambdavsCent", "hFDFractionLambdavsCent", numCentLambdaOO, 0, 100);
  TH1F *hFDFractionALambdavsCent = new TH1F("hFDFractionALambdavsCent", "hFDFractionALambdavsCent", numCentLambdaOO, 0, 100);
  TH1F *hFDFractionAllLambdavsCent = new TH1F("hFDFractionAllLambdavsCent", "hFDFractionAllLambdavsCent", numCentLambdaOO, 0, 100);
  TCanvas *canvasFDLambda = new TCanvas("canvasFDLambda", "canvasFDLambda", 900, 700);
  gStyle->SetOptStat(0);

  hFDLambda->GetYaxis()->SetRangeUser(MinPt[ChosenParticle], MaxPt[ChosenParticle]); // pt range
  TH2F *hFDLambdaProj2D = (TH2F *)hFDLambda->Project3D("zxe");                       // particle vs centrality
  hFDLambdaProj2D->SetName("FDFraction2D");

  for (Int_t mul = 0; mul < numCentLambdaOO; mul++)
  {
    Int_t CentFT0CMax = 0;
    Int_t CentFT0CMin = 0;
    if (mul == numCentLambdaOO)
    {
      CentFT0CMin = 0;
      CentFT0CMax = CentFT0CMaxLambdaOO;
    }
    else
    {
      CentFT0CMin = CentFT0CLambdaOO[mul];
      CentFT0CMax = CentFT0CLambdaOO[mul + 1];
    }

    hFDLambdaProj2D->GetXaxis()->SetRange(hFDLambdaProj2D->GetXaxis()->FindBin(CentFT0CMin + 0.1), hFDLambdaProj2D->GetXaxis()->FindBin(CentFT0CMax - 0.1));
    TH1F *hFDLambdaProj = (TH1F *)hFDLambdaProj2D->ProjectionY(Form("FDFraction_%d_%d", CentFT0CMin, CentFT0CMax));
    hFDFractionLambdavsCent->SetBinContent(mul + 1, hFDLambdaProj->GetBinContent(2) / hFDLambdaProj->GetBinContent(1));
    hFDFractionALambdavsCent->SetBinContent(mul + 1, hFDLambdaProj->GetBinContent(4) / hFDLambdaProj->GetBinContent(3));
    hFDFractionAllLambdavsCent->SetBinContent(mul + 1, (hFDLambdaProj->GetBinContent(2) + hFDLambdaProj->GetBinContent(4)) / (hFDLambdaProj->GetBinContent(1) + hFDLambdaProj->GetBinContent(3)));
    cout << "Centrality " << CentFT0CMin << "-" << CentFT0CMax << "%: FD Lambda Fraction = " << hFDFractionLambdavsCent->GetBinContent(mul + 1) << ", FD Anti-Lambda Fraction = " << hFDFractionALambdavsCent->GetBinContent(mul + 1) << endl;
    hFDFractionLambdavsCent->SetBinError(mul + 1, sqrt(1. / hFDLambdaProj->GetBinContent(2) + 1. / hFDLambdaProj->GetBinContent(1)) * hFDFractionLambdavsCent->GetBinContent(mul + 1));
    hFDFractionALambdavsCent->SetBinError(mul + 1, sqrt(1. / hFDLambdaProj->GetBinContent(4) + 1. / hFDLambdaProj->GetBinContent(3)) * hFDFractionALambdavsCent->GetBinContent(mul + 1));
    hFDFractionAllLambdavsCent->SetBinError(mul + 1, sqrt(1. / (hFDLambdaProj->GetBinContent(2) + hFDLambdaProj->GetBinContent(4)) + 1. / (hFDLambdaProj->GetBinContent(1) + hFDLambdaProj->GetBinContent(3))) * hFDFractionAllLambdavsCent->GetBinContent(mul + 1));

    // StyleHistoYield(hFDLambdaProj, 0, 1, kBlue + 1, 20, "p_{T} (GeV/c)", "Fraction of primary #Lambda", "", 1.2, 1.5, 1.7);
    canvasFDLambda->cd();
    // hFDLambdaProj2D->Draw("same");
    hFDLambdaProj2D->Draw("same");
    // canvasFDLambda->SaveAs(Form("../FDFraction/Canvas_FDFraction_Lambda_CentBin%d.png", mul));
    // canvasFDLambda->SaveAs(Form("../FDFraction/Canvas_FDFraction_Lambda_CentBin%d.pdf", mul));
  }

  Float_t xTitle = 30;
  Float_t xOffset = 1.3;
  Float_t yTitle = 38;
  Float_t yOffset = 1.6;

  Float_t xLabel = 30;
  Float_t yLabel = 30;
  Float_t xLabelOffset = 0.015;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.03;
  Float_t tickY = 0.025;

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, -1000);
  SetFont(hDummy);
  StyleHistoYield(hDummy, 0.08, 0.12, 1, 1, TitleXCent, "Fraction of secondary #Lambda", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  if (ChosenParticle == 6)
    hDummy->GetXaxis()->SetRangeUser(0, 90);

  TCanvas *canvasFDFractionVsCent = new TCanvas("canvasFDFractionVsCent", "canvasFDFractionVsCent", 900, 700);
  StyleCanvas(canvasFDFractionVsCent, 0.06, 0.12, 0.15, 0.03);
  canvasFDFractionVsCent->cd();
  SetFont(hFDFractionLambdavsCent);
  StyleHistoYield(hFDFractionLambdavsCent, 0, 1, kBlue + 1, 20, "Centrality (%)", "Fraction of primary #Lambda", "", 1.5, 1.5, 1.7);
  StyleHistoYield(hFDFractionALambdavsCent, 0, 1, kGreen + 1, 20, "Centrality (%)", "Fraction of primary #Lambda", "", 1.5, 1.5, 1.7);
  StyleHistoYield(hFDFractionAllLambdavsCent, 0, 1, kRed + 1, 20, "Centrality (%)", "Fraction of primary #Lambda", "", 1.5, 1.5, 1.7);
  hDummy->Draw("");
  hFDFractionLambdavsCent->Draw("same");
  hFDFractionALambdavsCent->Draw("same");
  hFDFractionAllLambdavsCent->Draw("same");

  TLegend *legendFDFraction = new TLegend(0.2, 0.7, 0.5, 0.9);
  legendFDFraction->SetFillStyle(0);
  legendFDFraction->SetTextAlign(12);
  legendFDFraction->SetTextSize(0.048);
  legendFDFraction->AddEntry(hFDFractionLambdavsCent, "#Lambda", "pl");
  legendFDFraction->AddEntry(hFDFractionALambdavsCent, "#bar{#Lambda}", "pl");
  legendFDFraction->AddEntry(hFDFractionAllLambdavsCent, "#Lambda + #bar{#Lambda}", "pl");
  legendFDFraction->Draw("");

  canvasFDFractionVsCent->SaveAs("../FDFraction/FDFraction_Lambda_vs_Cent.pdf");
  canvasFDFractionVsCent->SaveAs("../FDFraction/FDFraction_Lambda_vs_Cent.png");

  Float_t FDFraction[numCentLambdaOO] = {0};
  Float_t FDFractionLambda[numCentLambdaOO] = {0};
  Float_t FDFractionALambda[numCentLambdaOO] = {0};
  Float_t Correction[numCentLambdaOO] = {0};
  Float_t CorrectionLambda[numCentLambdaOO] = {0};
  Float_t CorrectionALambda[numCentLambdaOO] = {0};
  TH1F *hFDFractionRelSistVsCent = new TH1F("hFDFractionRelSistVsCent", "hFDFractionRelSistVsCent", numCentLambdaOO, 0, 100);
  for (Int_t m = 0; m < numCentLambdaOO; m++)
  {
    FDFraction[m] = hFDFractionAllLambdavsCent->GetBinContent(m + 1);
    FDFractionLambda[m] = hFDFractionLambdavsCent->GetBinContent(m + 1);
    FDFractionALambda[m] = hFDFractionALambdavsCent->GetBinContent(m + 1);
    cout << "Centrality bin " << m << ": FD Fraction Lambda = " << FDFractionLambda[m] << ", FD Fraction Anti-Lambda = " << FDFractionALambda[m] << ", FD Fraction All = " << FDFraction[m] << endl;
    cout << "Max difference / 2 = " << abs(FDFractionLambda[m] - FDFractionALambda[m]) / 2 << endl;
    Correction[m] = 1. / (CXiToLambda * FDFraction[m] + (1 - FDFraction[m]));
    CorrectionLambda[m] = 1. / (CXiToLambda * FDFractionLambda[m] + (1 - FDFractionLambda[m]));
    CorrectionALambda[m] = 1. / (CXiToLambda * FDFractionALambda[m] + (1 - FDFractionALambda[m]));
    cout << "  -> Correction Lambda = " << CorrectionLambda[m] << ", Correction Anti-Lambda = " << CorrectionALambda[m] << ", Correction All = " << Correction[m] << endl;
    cout << "  -> Max difference / 2 = " << abs(CorrectionLambda[m] - CorrectionALambda[m]) / 2 << endl;
    cout << "Relative difference: " << abs(CorrectionLambda[m] - CorrectionALambda[m]) / (2 * Correction[m]) << endl;
    hFDFractionRelSistVsCent->SetBinContent(m + 1, abs(CorrectionLambda[m] - CorrectionALambda[m]) / 2);
  }
  TString stringout = "../LambdaFDFraction" + SinputFileNameFDFraction + ".root";
  TFile *fileout = new TFile(stringout, "RECREATE");
  hFDFractionLambdavsCent->Write();
  hFDFractionALambdavsCent->Write();
  hFDFractionAllLambdavsCent->Write();
  hFDFractionRelSistVsCent->Write();
  fileout->Close();

  cout << "\nStarting from the file: " << PathInFD << endl;
  cout << "\nI have created the file:\n " << stringout << endl;
}