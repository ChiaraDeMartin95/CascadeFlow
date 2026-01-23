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
#include "TFitResult.h"

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

void MultVsCent()
{

  TH1F *hMultVsCent = new TH1F("hMultVsCent", "hMultVsCent", numCentLambdaOO, 0, 100);

  for (Int_t mul = 0; mul < numCentLambdaOO - 1; mul++)
  {
    Int_t CentFT0CMax = 0;
    Int_t CentFT0CMin = 0;
    CentFT0CMin = CentFT0CLambdaOO[mul];
    CentFT0CMax = CentFT0CLambdaOO[mul + 1];

    hMultVsCent->SetBinContent(mul + 1, dNdEtaOOPrel[mul]);
    hMultVsCent->SetBinError(mul + 1, dNdEtaOOErrPrel[mul]);

    if (CentFT0CMax > 60)
    {
      hMultVsCent->SetBinContent(mul + 1, -10);
      hMultVsCent->SetBinError(mul + 1, 0);
    }
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

  gStyle->SetOptStat(0);
  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, -1000);
  SetFont(hDummy);
  StyleHistoYield(hDummy, 0, 150, 1, 1, TitleXCent, "<dN/d#eta>", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  StyleHistoYield(hMultVsCent, 0, 150, kBlue + 1, 20, TitleXCent, "<dN/d#eta>", "", 1.5, 1.15, 1.6);

  TF1 *fitExpo = new TF1("fitExpo", "expo", 0, 60);
  TFitResultPtr fFitResult = hMultVsCent->Fit(fitExpo, "SR+");
  TMatrixDSym cov = fFitResult->GetCovarianceMatrix();
  Double_t covFit = cov[0][1];

  TF1 *fitExpoFull = new TF1("fitExpoFull", "expo", 0, 100);
  fitExpoFull->SetLineColor(kBlack);
  //hMultVsCent->Fit(fitExpoFull, "R+");

  TF1 *fitExpoFinal = new TF1("fitExpoFinal", "expo", 0, 100);
  fitExpoFinal->SetParameters(fitExpo->GetParameter(0), fitExpo->GetParameter(1));
  fitExpoFinal->SetLineColor(kRed);

  TCanvas *canvasMultVsCent = new TCanvas("canvasMultVsCent", "canvasMultVsCent", 900, 700);
  StyleCanvas(canvasMultVsCent, 0.08, 0.12, 0.15, 0.05);
  canvasMultVsCent->cd();
  hDummy->Draw("");
  hMultVsCent->Draw("same");
  fitExpoFinal->Draw("same");
  canvasMultVsCent->SaveAs("../MultVsCent/Canvas_MultVsCent_LambdaOO.png");
  canvasMultVsCent->SaveAs("../MultVsCent/Canvas_MultVsCent_LambdaOO.pdf");

  cout << "\nhola " << sqrt(cov[0][0]) << " " << fitExpo->GetParError(0) << endl;

  cout << "Fit results for <dN/d#eta> vs Centrality:" << endl;
  for (Int_t mul = 0; mul < numCentLambdaOO - 1; mul++)
  {
    Float_t centValue = (CentFT0CLambdaOO[mul] + CentFT0CLambdaOO[mul + 1]) / 2.;
    Float_t fitValue = fitExpoFinal->Eval(centValue);
    Float_t fitError = fitValue * TMath::Sqrt(pow(fitExpo->GetParError(0), 2) + pow(centValue * fitExpo->GetParError(1), 2) + 2 * centValue * covFit);
    cout << "Centrality " << CentFT0CLambdaOO[mul] << "-" << CentFT0CLambdaOO[mul + 1]
         << "%: Fit Value = " << fitValue << " +- " << fitError
         << ", Original Value = " << hMultVsCent->GetBinContent(mul + 1) << "+- " << hMultVsCent->GetBinError(mul + 1) << endl;
  }
}