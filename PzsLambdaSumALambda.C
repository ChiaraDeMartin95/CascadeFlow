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
#include "TGraphErrors.h"
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

void PzsLambdaSumALambda()
{

  TString SfileLambda = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2_Train589711_LambdaPart_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly.root";
  TFile *fLambda = TFile::Open(SfileLambda);
  TH1F *fHistPzsLambda = (TH1F *)fLambda->Get("fHistPzs");
  fHistPzsLambda->SetName("fHistPzsLambda");

  TString SfileAntiLambda = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2_Train589711_AntiLambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly.root";
  TFile *fAntiLambda = TFile::Open(SfileAntiLambda);
  TH1F *fHistPzsAntiLambda = (TH1F *)fAntiLambda->Get("fHistPzs");
  fHistPzsAntiLambda->SetName("fHistPzsAntiLambda");

  TString SfileAllLambda = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2_Train589711_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly.root";
  TFile *fAllLambda = TFile::Open(SfileAllLambda);
  TH1F *fHistPzsAllLambda = (TH1F *)fAllLambda->Get("fHistPzs");
  fHistPzsAllLambda->SetName("fHistPzsAllLambda");

  TH1F *hSumLambdaAntiLambda = (TH1F *)fHistPzsLambda->Clone("hSumLambdaAntiLambda");
  for (Int_t i = 1; i <= hSumLambdaAntiLambda->GetNbinsX(); i++)
  {
    float x1 = fHistPzsLambda->GetBinContent(i);
    float x2 = fHistPzsAntiLambda->GetBinContent(i);
    float err1 = fHistPzsLambda->GetBinError(i);
    float err2 = fHistPzsAntiLambda->GetBinError(i);
    hSumLambdaAntiLambda->SetBinContent(i, (x1 / pow(err1, 2) + x2 / pow(err2, 2)) / (1. / pow(err1, 2) + 1. / pow(err2, 2)));
    hSumLambdaAntiLambda->SetBinError(i, TMath::Sqrt(1. / (1. / pow(err1, 2) + 1. / pow(err2, 2))));
  }
  // hSumLambdaAntiLambda->Add(fHistPzsAntiLambda);

  TCanvas *cPzsSumALambda = new TCanvas("cPzsSumALambda", "cPzsSumALambda", 800, 600);
  StyleCanvas(cPzsSumALambda, 0.05, 0.12, 0.12, 0.05);
  StylePad(cPzsSumALambda, 0.12, 0.05, 0.05, 0.12);
  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10, 0, 90);
  StyleHistoYield(hDummy, 0, 0.03, kWhite, 20, "FT0C centrality (%)", "Pzs2 (GeV/c)", "", 1.5, 1.2, 1.3);
  hDummy->Draw();
  StyleHisto(hSumLambdaAntiLambda, 0, 0.03, kRed + 1, 20, "FT0C centrality (%)", "Pzs2 (GeV/c)", "Sum #Lambda + #bar{#Lambda}");
  StyleHisto(fHistPzsAllLambda, 0, 0.03, kBlue + 1, 21, "FT0C centrality (%)", "Pzs2 (GeV/c)", "#Lambda + #bar{#Lambda} from single file");
  hSumLambdaAntiLambda->Draw("E1");
  fHistPzsAllLambda->Draw("E1 SAME");
  TLegend *legPzsSumALambda = new TLegend(0.5, 0.7, 0.9, 0.9);
  legPzsSumALambda->SetBorderSize(0);
  legPzsSumALambda->SetTextSize(0.04);
  legPzsSumALambda->SetFillStyle(0);
  legPzsSumALambda->AddEntry(hSumLambdaAntiLambda, "Sum #Lambda + #bar{#Lambda}", "lep");
  legPzsSumALambda->AddEntry(fHistPzsAllLambda, "#Lambda + #bar{#Lambda} from single file", "lep");
  legPzsSumALambda->Draw();

  TString SfileOutput = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2_Train589711_LambdaPlusAntiLambda.root";
  TFile *fOutput = TFile::Open(SfileOutput, "RECREATE");
  fOutput->cd();
  hSumLambdaAntiLambda->Write();
  fOutput->Close();
}
