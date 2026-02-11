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
// #include "CommonVar.h"
#include "CommonVarLambda.h"
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
void QCTopoVariables()
{

  gStyle->SetOptStat(0);

  TString SinputFile = "../OutputAnalysis/Output_LHC25_OO_pass2_Train589711_Lambda_CentWeighted_Eta08_isLoosest_isOOCentrality_ResoOnTheFly_Nvar2.root";
  cout << "Input file: " << SinputFile << endl;

  TFile *inputFile = new TFile(SinputFile);
  TH1F *hV0Radius = (TH1F *)inputFile->Get("histoBefV0Radius");
  if (!hV0Radius)
  {
    cout << "histoBefV0Radius not found in " << SinputFile << endl;
    return;
  }
  TH1F *hCosPA = (TH1F *)inputFile->Get("histoBefV0CosPA");
  if (!hCosPA)
  {
    cout << "histoBefV0CosPA not found in " << SinputFile << endl;
    return;
  }
  TH1F *hDCA = (TH1F *)inputFile->Get("histoBefDcaV0Daughters");
  if (!hDCA)
  {
    cout << "histoBefDcaV0Daughters not found in " << SinputFile << endl;
    return;
  }
  TH1F *hDCAPosToPV = (TH1F *)inputFile->Get("histoBefDcaPosToPV");
  if (!hDCAPosToPV)
  {
    cout << "histoBefDcaPosToPV not found in " << SinputFile << endl;
    return;
  }
  TH1F *hDCANegToPV = (TH1F *)inputFile->Get("histoBefDcaNegToPV");
  if (!hDCANegToPV)
  {
    cout << "histoBefDcaNegToPV not found in " << SinputFile << endl;
    return;
  }

  TCanvas *cV0Radius = new TCanvas("cV0Radius", "cV0Radius", 800, 600);
  StyleCanvas(cV0Radius, 0.12, 0.1, 0.12, 0.1);
  hV0Radius->GetYaxis()->SetRangeUser(1, 1.2 * hV0Radius->GetMaximum());
  TLine *lineV0Radius = new TLine(DefaultV0RadiusCut, 0, DefaultV0RadiusCut, hV0Radius->GetMaximum());
  lineV0Radius->SetLineColor(kBlack);
  TLine *lineV0RadiusTight = new TLine(LowerlimitV0RadiusCut, 0, LowerlimitV0RadiusCut, hV0Radius->GetMaximum());
  lineV0RadiusTight->SetLineColor(kRed);
  TLine *lineV0RadiusLoose = new TLine(UpperlimitV0RadiusCut, 0, UpperlimitV0RadiusCut, hV0Radius->GetMaximum());
  lineV0RadiusLoose->SetLineColor(kGreen + 2);
  hV0Radius->Draw();
  lineV0Radius->Draw("SAME");
  lineV0RadiusTight->Draw("SAME");
  lineV0RadiusLoose->Draw("SAME");
  cV0Radius->SaveAs("../Systematics/QCTopoVariableV0Radius.pdf");

  TCanvas *cCosPA = new TCanvas("cCosPA", "cCosPA", 800, 600);
  StyleCanvas(cCosPA, 0.12, 0.1, 0.12, 0.1);
  hCosPA->GetYaxis()->SetRangeUser(1, 2 * hCosPA->GetMaximum());
  TLine *lineCosPA = new TLine(DefaultV0CosPA, 1, DefaultV0CosPA, hCosPA->GetMaximum());
  lineCosPA->SetLineColor(kBlack);
  TLine *lineCosPATight = new TLine(UpperlimitV0CosPA, 1, UpperlimitV0CosPA, hCosPA->GetMaximum());
  lineCosPATight->SetLineColor(kRed);
  TLine *lineCosPALoose = new TLine(LowerlimitV0CosPA, 1, LowerlimitV0CosPA, hCosPA->GetMaximum());
  lineCosPALoose->SetLineColor(kGreen + 2);
  gPad->SetLogy();
  hCosPA->Draw();
  lineCosPA->Draw("SAME");
  lineCosPATight->Draw("SAME");
  lineCosPALoose->Draw("SAME");
  cCosPA->SaveAs("../Systematics/QCTopoVariableCosPA.pdf");

  TCanvas *cDCA = new TCanvas("cDCA", "cDCA", 800, 600);
  StyleCanvas(cDCA, 0.12, 0.1, 0.12, 0.1);
  hDCA->GetYaxis()->SetRangeUser(1, 2 * hDCA->GetMaximum());
  TLine *lineDCA = new TLine(DefaultDcaV0DauCut, 1, DefaultDcaV0DauCut, hDCA->GetMaximum());
  lineDCA->SetLineColor(kBlack);
  TLine *lineDCATight = new TLine(LowerlimitDcaV0DauCut, 1, LowerlimitDcaV0DauCut, hDCA->GetMaximum());
  lineDCATight->SetLineColor(kRed);
  TLine *lineDCALoose = new TLine(UpperlimitDcaV0DauCut, 1, UpperlimitDcaV0DauCut, hDCA->GetMaximum());
  lineDCALoose->SetLineColor(kGreen + 2);
  gPad->SetLogy();
  hDCA->Draw();
  lineDCA->Draw("SAME");
  lineDCATight->Draw("SAME");
  lineDCALoose->Draw("SAME");
  cDCA->SaveAs("../Systematics/QCTopoVariableDCA.pdf");

  TCanvas *cDCAPosToPV = new TCanvas("cDCAPosToPV", "cDCAPosToPV", 800, 600);
  StyleCanvas(cDCAPosToPV, 0.12, 0.1, 0.12, 0.1);
  hDCAPosToPV->GetYaxis()->SetRangeUser(0, 1.2 * hDCAPosToPV->GetMaximum());
  hDCAPosToPV->GetXaxis()->SetRangeUser(0, 1);
  TLine *lineDCAPosToPV = new TLine(DefaultDcaPosToPV, 0, DefaultDcaPosToPV, hDCAPosToPV->GetMaximum());
  lineDCAPosToPV->SetLineColor(kBlack);
  TLine *lineDCAPosToPVTight = new TLine(UpperlimitDcaPosToPV, 0, UpperlimitDcaPosToPV, hDCAPosToPV->GetMaximum());
  lineDCAPosToPVTight->SetLineColor(kRed);
  TLine *lineDCAPosToPVLoose = new TLine(LowerlimitDcaPosToPV, 0, LowerlimitDcaPosToPV, hDCAPosToPV->GetMaximum());
  lineDCAPosToPVLoose->SetLineColor(kGreen + 2);
  hDCAPosToPV->Draw();
  lineDCAPosToPV->Draw("SAME");
  lineDCAPosToPVTight->Draw("SAME");
  lineDCAPosToPVLoose->Draw("SAME");
  cDCAPosToPV->SaveAs("../Systematics/QCTopoVariableDCAPosToPV.pdf");

  TCanvas *cDCANegToPV = new TCanvas("cDCANegToPV", "cDCANegToPV", 800, 600);
  StyleCanvas(cDCANegToPV, 0.12, 0.1, 0.12, 0.1);
  hDCANegToPV->GetYaxis()->SetRangeUser(0, 1.2 * hDCANegToPV->GetMaximum());
  hDCANegToPV->GetXaxis()->SetRangeUser(0, 1);
  TLine *lineDCANegToPV = new TLine(DefaultDcaNegToPV, 0, DefaultDcaNegToPV, hDCANegToPV->GetMaximum());
  lineDCANegToPV->SetLineColor(kBlack);
  TLine *lineDCANegToPVTight = new TLine(UpperlimitDcaNegToPV, 0, UpperlimitDcaNegToPV, hDCANegToPV->GetMaximum());
  lineDCANegToPVTight->SetLineColor(kRed);
  TLine *lineDCANegToPVLoose = new TLine(LowerlimitDcaNegToPV, 0, LowerlimitDcaNegToPV, hDCANegToPV->GetMaximum());
  lineDCANegToPVLoose->SetLineColor(kGreen + 2);
  hDCANegToPV->Draw();
  lineDCANegToPV->Draw("SAME");
  lineDCANegToPVTight->Draw("SAME");
  lineDCANegToPVLoose->Draw("SAME");
  cDCANegToPV->SaveAs("../Systematics/QCTopoVariableDCANegToPV.pdf");
}
