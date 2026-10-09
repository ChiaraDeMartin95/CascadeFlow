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
#include "CommonVarPub.h"
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

  TH1F *hFDFractionLambdavsPt = new TH1F("hFDFractionLambdavsPt", "hFDFractionLambdavsPt", numPtBins, PtBins);
  TH1F *hFDFractionALambdavsPt = new TH1F("hFDFractionALambdavsPt", "hFDFractionALambdavsPt", numPtBins, PtBins);
  TH1F *hFDFractionAllLambdavsPt = new TH1F("hFDFractionAllLambdavsPt", "hFDFractionAllLambdavsPt", numPtBins, PtBins);

  TCanvas *canvasFDLambda = new TCanvas("canvasFDLambda", "canvasFDLambda", 900, 700);
  gStyle->SetOptStat(0);

  TH3F *hFDLambdaClone = (TH3F *)hFDLambda->Clone("hFDLambdaClone");

  hFDLambda->GetYaxis()->SetRangeUser(MinPt[ChosenParticle], MaxPt[ChosenParticle]); // pt range
  TH2F *hFDLambdaProj2D = (TH2F *)hFDLambda->Project3D("zxe");                       // particle vs centrality
  hFDLambdaProj2D->SetName("FDFraction2D");

  hFDLambdaClone->GetXaxis()->SetRangeUser(0, CentFT0CMaxLambdaOO);     // centrality range
  TH2F *hFDLambdaProj2DvsPt = (TH2F *)hFDLambdaClone->Project3D("zye"); // particle vs pt
  hFDLambdaProj2DvsPt->SetName("FDFraction2DvsPt");

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
    hFDFractionLambdavsCent->SetBinContent(mul + 1, hFDLambdaProj->GetBinContent(2) / (hFDLambdaProj->GetBinContent(1) + hFDLambdaProj->GetBinContent(2)));
    hFDFractionALambdavsCent->SetBinContent(mul + 1, hFDLambdaProj->GetBinContent(4) / (hFDLambdaProj->GetBinContent(3) + hFDLambdaProj->GetBinContent(4)));
    hFDFractionAllLambdavsCent->SetBinContent(mul + 1, (hFDLambdaProj->GetBinContent(2) + hFDLambdaProj->GetBinContent(4)) / (hFDLambdaProj->GetBinContent(1) + hFDLambdaProj->GetBinContent(3) + hFDLambdaProj->GetBinContent(2) + hFDLambdaProj->GetBinContent(4)));
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

  for (Int_t pt = 0; pt < numPtBins; pt++)
  {
    Float_t PtMin = PtBins[pt];
    Float_t PtMax = PtBins[pt + 1];
    hFDLambdaProj2DvsPt->GetXaxis()->SetRange(hFDLambdaProj2DvsPt->GetXaxis()->FindBin(PtMin + 0.01), hFDLambdaProj2DvsPt->GetXaxis()->FindBin(PtMax - 0.01));
    TH1F *hFDLambdaProjPt = (TH1F *)hFDLambdaProj2DvsPt->ProjectionY(Form("FDFraction_Pt_%d_%d", Int_t(PtMin * 10), Int_t(PtMax * 10)));
    hFDFractionLambdavsPt->SetBinContent(pt + 1, hFDLambdaProjPt->GetBinContent(2) / hFDLambdaProjPt->GetBinContent(1));
    hFDFractionALambdavsPt->SetBinContent(pt + 1, hFDLambdaProjPt->GetBinContent(4) / hFDLambdaProjPt->GetBinContent(3));
    hFDFractionAllLambdavsPt->SetBinContent(pt + 1, (hFDLambdaProjPt->GetBinContent(2) + hFDLambdaProjPt->GetBinContent(4)) / (hFDLambdaProjPt->GetBinContent(1) + hFDLambdaProjPt->GetBinContent(3)));
    cout << "Pt " << PtMin << "-" << PtMax << " GeV/c: FD Lambda Fraction = " << hFDFractionLambdavsPt->GetBinContent(pt + 1) << ", FD Anti-Lambda Fraction = " << hFDFractionALambdavsPt->GetBinContent(pt + 1) << endl;
    hFDFractionLambdavsPt->SetBinError(pt + 1, sqrt(1. / hFDLambdaProjPt->GetBinContent(2) + 1. / hFDLambdaProjPt->GetBinContent(1)) * hFDFractionLambdavsPt->GetBinContent(pt + 1));
    hFDFractionALambdavsPt->SetBinError(pt + 1, sqrt(1. / hFDLambdaProjPt->GetBinContent(4) + 1. / hFDLambdaProjPt->GetBinContent(3)) * hFDFractionALambdavsPt->GetBinContent(pt + 1));
    hFDFractionAllLambdavsPt->SetBinError(pt + 1, sqrt(1. / (hFDLambdaProjPt->GetBinContent(2) + hFDLambdaProjPt->GetBinContent(4)) + 1. / (hFDLambdaProjPt->GetBinContent(1) + hFDLambdaProjPt->GetBinContent(3))) * hFDFractionAllLambdavsPt->GetBinContent(pt + 1));
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
  StyleHistoYield(hDummy, 0, 0.3, 1, 1, TitleXCent, "Fraction of secondary #Lambda", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  if (ChosenParticle == 6)
    hDummy->GetXaxis()->SetRangeUser(0, 90);

  TCanvas *canvasFDFractionVsCent = new TCanvas("canvasFDFractionVsCent", "canvasFDFractionVsCent", 900, 700);
  StyleCanvas(canvasFDFractionVsCent, 0.08, 0.12, 0.15, 0.03);
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

  TH1F *hDummyPt = new TH1F("hDummyPt", "hDummyPt", 10000, 0, 10);
  for (Int_t i = 1; i <= hDummyPt->GetNbinsX(); i++)
    hDummyPt->SetBinContent(i, -1000);
  SetFont(hDummyPt);
  StyleHistoYield(hDummyPt, 0, 0.3, 1, 1, TitleXPt, "Fraction of secondary #Lambda", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummyPt, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummyPt, tickX, tickY);
  StyleHistoYield(hFDFractionLambdavsPt, 0, 1, kBlue + 1, 20, TitleXPt, "Fraction of secondary #Lambda", "", 1.5, 1.5, 1.7);
  StyleHistoYield(hFDFractionALambdavsPt, 0, 1, kGreen + 1, 20, TitleXPt, "Fraction of secondary #Lambda", "", 1.5, 1.5, 1.7);
  StyleHistoYield(hFDFractionAllLambdavsPt, 0, 1, kRed + 1, 20, TitleXPt, "Fraction of secondary #Lambda", "", 1.5, 1.5, 1.7);
  TCanvas *canvasFDFractionVsPt = new TCanvas("canvasFDFractionVsPt", "canvasFDFractionVsPt", 900, 700);
  StyleCanvas(canvasFDFractionVsPt, 0.06, 0.12, 0.15, 0.03);
  canvasFDFractionVsPt->cd();
  hDummyPt->Draw("");
  hFDFractionLambdavsPt->Draw("same");
  hFDFractionALambdavsPt->Draw("same");
  hFDFractionAllLambdavsPt->Draw("same");
  legendFDFraction->Draw("");
  canvasFDFractionVsPt->SaveAs("../FDFraction/FDFraction_Lambda_vs_Pt.pdf");
  canvasFDFractionVsPt->SaveAs("../FDFraction/FDFraction_Lambda_vs_Pt.png");

  Float_t FDFraction[numCentLambdaOO] = {0};
  Float_t FDFractionLambda[numCentLambdaOO] = {0};
  Float_t FDFractionALambda[numCentLambdaOO] = {0};
  Float_t Correction[numCentLambdaOO] = {0};
  Float_t CorrectionLambda[numCentLambdaOO] = {0};
  Float_t CorrectionALambda[numCentLambdaOO] = {0};
  TH1F *hFDFractionRelSistVsCent = new TH1F("hFDFractionRelSistVsCent", "hFDFractionRelSistVsCent", numCentLambdaOO, 0, 100);
  TH1D *hRelSystV1 = new TH1D("hRelSystV1", "hRelSystV1", numCentLambdaOO, 0, 100);
  TH1D *hRelSystV2 = new TH1D("hRelSystV2", "hRelSystV2", numCentLambdaOO, 0, 100);
  TH1F *hRelSystV1Plus = new TH1F("hRelSystV1Plus", "hRelSystV1Plus", numCentLambdaOO, 0, 100);
  TH1F *hRelSystV1Minus = new TH1F("hRelSystV1Minus", "hRelSystV1Minus", numCentLambdaOO, 0, 100);
  TH1D *hRelSystSouravV1 = new TH1D("hRelSystSouravV1", "hRelSystSouravV1", numCentLambdaOO, 0, 100);
  TH1D *hRelSystSouravV1Plus = new TH1D("hRelSystSouravV1Plus", "hRelSystSouravV1Plus", numCentLambdaOO, 0, 100);
  TH1D *hRelSystSouravV1Minus = new TH1D("hRelSystSouravV1Minus", "hRelSystSouravV1Minus", numCentLambdaOO, 0, 100);
  // rel syst version 1: fFD |C * R -1|
  // rel syst version 2: fFD |C * R -1| / (1-fFD +C *  fFD* R)
  Float_t XiLambdaRatio = 0.98; // R
  Float_t XiLambdaRatioErr = 0.25;
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
    hRelSystV1->SetBinContent(m + 1, FDFraction[m] * std::abs(CXiToLambda * XiLambdaRatio - 1));
    hRelSystV1Plus->SetBinContent(m + 1, FDFraction[m] * std::abs(CXiToLambda * (XiLambdaRatio + XiLambdaRatioErr) - 1));
    hRelSystV1Minus->SetBinContent(m + 1, FDFraction[m] * std::abs(CXiToLambda * (XiLambdaRatio - XiLambdaRatioErr) - 1));
    hRelSystV2->SetBinContent(m + 1, FDFraction[m] * std::abs(CXiToLambda * XiLambdaRatio - 1) / (1 - FDFraction[m] + CXiToLambda * FDFraction[m] * XiLambdaRatio));
    hRelSystSouravV1->SetBinContent(m + 1, FDFraction[m] * std::abs(CXiToLambda * XiLambdaRatio - 1) / (1 - FDFraction[m] * CXiToLambda * XiLambdaRatio));
    hRelSystSouravV1Plus->SetBinContent(m + 1, FDFraction[m] * std::abs(CXiToLambda * (XiLambdaRatio + XiLambdaRatioErr) - 1) / (1 - FDFraction[m] * CXiToLambda * (XiLambdaRatio + XiLambdaRatioErr)));
    hRelSystSouravV1Minus->SetBinContent(m + 1, FDFraction[m] * std::abs(CXiToLambda * (XiLambdaRatio - XiLambdaRatioErr) - 1) / (1 - FDFraction[m] * CXiToLambda * (XiLambdaRatio - XiLambdaRatioErr)));
  }
  TCanvas *canvasSystFDFractionvsCent = new TCanvas("canvasSystFDFractionvsCent", "canvasSystFDFractionvsCent", 900, 700);
  StyleCanvas(canvasSystFDFractionvsCent, 0.06, 0.15, 0.13, 0.03);
  canvasSystFDFractionvsCent->cd();
  TH1F *hDummyRelError = new TH1F("hDummyRelError", "hDummyRelError", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummyRelError->GetNbinsX(); i++)
    hDummyRelError->SetBinContent(i, -1000);
  SetFont(hDummyRelError);
  StyleHistoYield(hDummyRelError, 0., 0.08, 1, 1, TitleXCent, "Rel. error", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummyRelError, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummyRelError, tickX, tickY);
  hDummyRelError->Draw();
  // hFDFractionRelSistVsCent->Draw("same");
  hRelSystV1->SetLineColor(kRed);
  // hRelSystV1->Draw("same");
  //  hRelSystV2->SetLineColor(kBlue);
  //  hRelSystV2->Draw("same");
  hRelSystV1Plus->SetLineColor(kGreen + 2);
  // hRelSystV1Plus->Draw("same");
  hRelSystV1Minus->SetLineColor(kMagenta);
  // hRelSystV1Minus->Draw("same");
  hRelSystSouravV1->SetLineColor(kRed);
  hRelSystSouravV1->SetLineStyle(3);
  hRelSystSouravV1->SetLineWidth(3);
  hRelSystSouravV1->Draw("same");
  hRelSystSouravV1Plus->SetLineColor(kGreen + 2);
  hRelSystSouravV1Plus->SetLineStyle(3);
  hRelSystSouravV1Plus->SetLineWidth(3);
  hRelSystSouravV1Plus->Draw("same");
  hRelSystSouravV1Minus->SetLineColor(kMagenta);
  hRelSystSouravV1Minus->SetLineStyle(3);
  hRelSystSouravV1Minus->SetLineWidth(3);
  hRelSystSouravV1Minus->Draw("same");
  canvasSystFDFractionvsCent->SaveAs("../FDFraction/SystFDFraction_vs_Cent.pdf");
  canvasSystFDFractionvsCent->SaveAs("../FDFraction/SystFDFraction_vs_Cent.png");

  TCanvas *canvasSystFDFractionVsPt = new TCanvas("canvasSystFDFractionVsPt", "canvasSystFDFractionVsPt", 900, 700);
  StyleCanvas(canvasSystFDFractionVsPt, 0.06, 0.15, 0.13, 0.03);
  canvasSystFDFractionVsPt->cd();
  TH1F *hDummyPtRelError = new TH1F("hDummyPtRelError", "hDummyPtRelError", 10000, 0, 10);
  for (Int_t i = 1; i <= hDummyPtRelError->GetNbinsX(); i++)
    hDummyPtRelError->SetBinContent(i, -1000);
  SetFont(hDummyPtRelError);
  StyleHistoYield(hDummyPtRelError, 0, 0.08, 1, 1, TitleXPt, "Rel. error", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummyPtRelError, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummyPtRelError, tickX, tickY);
  TH1F *hRelSystV1Pt = (TH1F *)hFDFractionAllLambdavsPt->Clone("hRelSystV1Pt");
  TH1F *hRelSystV1PlusPt = (TH1F *)hFDFractionAllLambdavsPt->Clone("hRelSystV1PlusPt");
  TH1F *hRelSystV1MinusPt = (TH1F *)hFDFractionAllLambdavsPt->Clone("hRelSystV1MinusPt");
  TH1F *hRelSystSouravV1Pt = (TH1F *)hFDFractionAllLambdavsPt->Clone("hRelSystSouravV1Pt");
  TH1F *hRelSystSouravV1PlusPt = (TH1F *)hFDFractionAllLambdavsPt->Clone("hRelSystSouravV1PlusPt");
  TH1F *hRelSystSouravV1MinusPt = (TH1F *)hFDFractionAllLambdavsPt->Clone("hRelSystSouravV1MinusPt");
  for (Int_t i = 1; i <= hRelSystV1Pt->GetNbinsX(); i++)
  {
    hRelSystV1Pt->SetBinContent(i, hFDFractionAllLambdavsPt->GetBinContent(i) * std::abs(CXiToLambda * XiLambdaRatio - 1));
    hRelSystV1PlusPt->SetBinContent(i, hFDFractionAllLambdavsPt->GetBinContent(i) * std::abs(CXiToLambda * (XiLambdaRatio + XiLambdaRatioErr) - 1));
    hRelSystV1MinusPt->SetBinContent(i, hFDFractionAllLambdavsPt->GetBinContent(i) * std::abs(CXiToLambda * (XiLambdaRatio - XiLambdaRatioErr) - 1));
    hRelSystSouravV1Pt->SetBinContent(i, hFDFractionAllLambdavsPt->GetBinContent(i) * std::abs(CXiToLambda * XiLambdaRatio - 1) / (1 - hFDFractionAllLambdavsPt->GetBinContent(i) * CXiToLambda * XiLambdaRatio));
    hRelSystSouravV1PlusPt->SetBinContent(i, hFDFractionAllLambdavsPt->GetBinContent(i) * std::abs(CXiToLambda * (XiLambdaRatio + XiLambdaRatioErr) - 1) / (1 - hFDFractionAllLambdavsPt->GetBinContent(i) * CXiToLambda * (XiLambdaRatio + XiLambdaRatioErr)));
    hRelSystSouravV1MinusPt->SetBinContent(i, hFDFractionAllLambdavsPt->GetBinContent(i) * std::abs(CXiToLambda * (XiLambdaRatio - XiLambdaRatioErr) - 1) / (1 - hFDFractionAllLambdavsPt->GetBinContent(i) * CXiToLambda * (XiLambdaRatio - XiLambdaRatioErr)));
  }
  hDummyPtRelError->Draw();
  hRelSystV1Pt->SetLineColor(kRed);
  hRelSystV1Pt->SetMarkerColor(kRed);
  //hRelSystV1Pt->Draw("same");
  hRelSystV1PlusPt->SetLineColor(kGreen + 2);
  hRelSystV1PlusPt->SetMarkerColor(kGreen + 2);
  //hRelSystV1PlusPt->Draw("same");
  hRelSystV1MinusPt->SetLineColor(kMagenta);
  hRelSystV1MinusPt->SetMarkerColor(kMagenta);
  //hRelSystV1MinusPt->Draw("same");

  hRelSystSouravV1Pt->SetLineColor(kRed);
  hRelSystSouravV1Pt->SetMarkerColor(kRed);
  hRelSystSouravV1Pt->SetLineStyle(3);
  hRelSystSouravV1Pt->SetLineWidth(3);
  hRelSystSouravV1Pt->Draw("same");
  hRelSystSouravV1PlusPt->SetLineColor(kGreen + 2);
  hRelSystSouravV1PlusPt->SetMarkerColor(kGreen + 2);
  hRelSystSouravV1PlusPt->SetLineStyle(3);
  hRelSystSouravV1PlusPt->SetLineWidth(3);
  hRelSystSouravV1PlusPt->Draw("same");
  hRelSystSouravV1MinusPt->SetLineColor(kMagenta);
  hRelSystSouravV1MinusPt->SetMarkerColor(kMagenta);
  hRelSystSouravV1MinusPt->SetLineStyle(3);
  hRelSystSouravV1MinusPt->SetLineWidth(3);
  hRelSystSouravV1MinusPt->Draw("same");
  canvasSystFDFractionVsPt->SaveAs("../FDFraction/SystFDFraction_vs_Pt.pdf");
  canvasSystFDFractionVsPt->SaveAs("../FDFraction/SystFDFraction_vs_Pt.png");

  TString stringout = "../LambdaFDFraction" + SinputFileNameFDFraction + ".root";
  TFile *fileout = new TFile(stringout, "RECREATE");
  hFDFractionLambdavsCent->Write();
  hFDFractionALambdavsCent->Write();
  hFDFractionAllLambdavsCent->Write();
  hFDFractionRelSistVsCent->Write();
  hRelSystV1MinusPt->Write();
  hRelSystSouravV1MinusPt->Write();

  fileout->Close();

  cout << "\nStarting from the file: " << PathInFD << endl;
  cout << "\nI have created the file:\n " << stringout << endl;
}