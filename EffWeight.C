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
//#include "CommonVar.h"
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
void EffWeight()
{

  TFile *fEff = new TFile("../" + SinputFileNameEfficiency, "READ");
  gStyle->SetOptStat(0);
  TDirectoryFile *dirEff = (TDirectoryFile *)fEff->Get("efficiencyHistograms_CENT_10_20"); //no centrality dependence expected
  if (!dirEff)
  {
    cout << "Error: Directory not found in the file!" << endl;
    return;
  }
  TH1F *hEfficiency = (TH1F *)dirEff->FindObject("hLambdaEff_CENT_10.00_20.00");
  if (!hEfficiency)
  {
    cout << "Error: Efficiency histogram not found in the file!" << endl;
    return;
  }
  
  TH1F *hEffWeight = (TH1F*)hEfficiency->Clone("hEffWeight");
  for (Int_t b = 1; b <= hEfficiency->GetNbinsX(); b++)
  {
    if (hEfficiency->GetBinContent(b) > 0)
      hEffWeight->SetBinContent(b, 1. / (hEfficiency->GetBinContent(b)));
  }

  TFile *fout = new TFile("../EfficiencyWeight.root", "RECREATE");
  hEffWeight->Write();
  TList *list = new TList();
  list->Add(hEffWeight);
  list->Write("ccdb_object", TObject::kSingleKey);
  fout->Close();

  TCanvas *cEff = new TCanvas("cEff", "cEff", 800, 600);
  StyleCanvas(cEff, 0.02, 0.13, 0.1, 0.03);
  cEff->cd();
  StyleHistoYield(hEfficiency, 0, 1, 1, 20, TitleXCent, "Efficiency", "", 1, 1.15, 1);
  hEfficiency->SetTitle("");
  hEfficiency->Draw("");
  cEff->SaveAs("../QCPlots/hEff.png");
  cEff->SaveAs("../QCPlots/hEff.pdf");

  TCanvas *cEffWeight = new TCanvas("cEffWeight", "cEffWeight", 800, 600);
  StyleCanvas(cEffWeight, 0.02, 0.13, 0.1, 0.03);
  cEffWeight->cd();
  StyleHistoYield(hEffWeight, 1, 15, 1, 20, TitleXCent, "Efficiency weight", "", 1, 1.15, 1);
  hEffWeight->SetTitle("");
  hEffWeight->Draw("");
  cEffWeight->SaveAs("../QCPlots/hEffWeight.png");
  cEffWeight->SaveAs("../QCPlots/hEffWeight.pdf");
}
