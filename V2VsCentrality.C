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

Float_t YLow = {-0.001};
Float_t YUp = {0.011};
Float_t xTitle = 30;
Float_t xOffset = 1.3;
Float_t yTitle = 38;   // 30
Float_t yOffset = 1.1; // 2.2 if setmaxdigits not set

Float_t xLabel = 30;
Float_t yLabel = 30;
Float_t xLabelOffset = 0.015;
Float_t yLabelOffset = 0.01;

Float_t tickX = 0.03;
Float_t tickY = 0.025;

void V2VsCentrality()
{

  TString SfileinPbPb = "../HEPData-ins1723697-v1-Table_31.root";
  TFile *fileinPbPb = new TFile(SfileinPbPb);
  if (!fileinPbPb || fileinPbPb->IsZombie())
  {
    std::cout << "Error opening file: " << SfileinPbPb << std::endl;
    return;
  }
  TDirectoryFile *dir = (TDirectoryFile *)fileinPbPb->Get("Table 31");
  if (!dir)
  {
    std::cout << "Error: Directory 'Table 31' not found in file: " << SfileinPbPb << std::endl;
    return;
  }
  TH1F *histoV2PbPb = (TH1F *)dir->Get("Hist1D_y1");
  if (!histoV2PbPb)
  {
    std::cout << "Error: Histogram 'Hist1D_y1' not found in directory 'Table 31'." << std::endl;
    return;
  }
  TH1F *histoV2PbPbStatErr = (TH1F *)dir->Get("Hist1D_y1_e1");
  if (!histoV2PbPbStatErr)
  {
    std::cout << "Error: Histogram 'Hist1D_y1_e1' not found in directory 'Table 31'." << std::endl;
    return;
  }
  TH1F *histoV2PbPbSystErr = (TH1F *)dir->Get("Hist1D_y1_e2");
  if (!histoV2PbPbSystErr)
  {
    std::cout << "Error: Histogram 'Hist1D_y1_e2' not found in directory 'Table 31'." << std::endl;
    return;
  }

  TString SfileinpPb = "../HEPData-ins1723697-v1-Table_10.root";
  TFile *fileinpPb = new TFile(SfileinpPb);
  if (!fileinpPb || fileinpPb->IsZombie())
  {
    std::cout << "Error opening file: " << SfileinpPb << std::endl;
    return;
  }
  TDirectoryFile *dirpPb = (TDirectoryFile *)fileinpPb->Get("Table 10");
  if (!dirpPb)
  {
    std::cout << "Error: Directory 'Table 10' not found in file: " << SfileinpPb << std::endl;
    return;
  } 
  TH1F *histoV2pPb = (TH1F *)dirpPb->Get("Hist1D_y1");
  if (!histoV2pPb)
  {
    std::cout << "Error: Histogram 'Hist1D_y1' not found in directory 'Table 10'." << std::endl;
    return;
  }
  TH1F *histoV2pPbStatErr = (TH1F *)dirpPb->Get("Hist1D_y1_e1");
  if (!histoV2pPbStatErr)
  {
    std::cout << "Error: Histogram 'Hist1D_y1_e1' not found in directory 'Table 10'." << std::endl;
    return;
  }
  TH1F *histoV2pPbSystErr = (TH1F *)dirpPb->Get("Hist1D_y1_e2");
  if (!histoV2pPbSystErr)
  {
    std::cout << "Error: Histogram 'Hist1D_y1_e2' not found in directory 'Table 10'." << std::endl;
    return;
  }

  TGraphAsymmErrors *gV2PbPb = new TGraphAsymmErrors(histoV2PbPb->GetNbinsX());
  TGraphAsymmErrors *gV2PbPbSyst = new TGraphAsymmErrors(histoV2PbPb->GetNbinsX());
  for (Int_t i = 1; i <= histoV2PbPb->GetNbinsX(); i++)
  {
    if (histoV2PbPb->GetBinContent(i) <= 0) continue;
    gV2PbPb->SetPoint(i - 1, histoV2PbPb->GetBinCenter(i) / 1.6, histoV2PbPb->GetBinContent(i));
    cout << "histoV2PbPb->GetBinCenter(i) " << histoV2PbPb->GetBinCenter(i) / 1.6 << " histoV2PbPb->GetBinContent(i) " << histoV2PbPb->GetBinContent(i) << endl;
    gV2PbPb->SetPointError(i - 1, abs(histoV2PbPb->GetBinCenter(i) - histoV2PbPb->GetBinLowEdge(i)) / 1.6, abs(histoV2PbPb->GetBinCenter(i) - histoV2PbPb->GetBinLowEdge(i + 1)) / 1.6, histoV2PbPbStatErr->GetBinContent(i), histoV2PbPbStatErr->GetBinContent(i));
    gV2PbPbSyst->SetPoint(i - 1, histoV2PbPb->GetBinCenter(i) / 1.6, histoV2PbPb->GetBinContent(i));
    gV2PbPbSyst->SetPointError(i - 1, abs(histoV2PbPb->GetBinCenter(i) - histoV2PbPb->GetBinLowEdge(i)) / 1.6, abs(histoV2PbPb->GetBinCenter(i) - histoV2PbPb->GetBinLowEdge(i + 1)) / 1.6, histoV2PbPbSystErr->GetBinContent(i), histoV2PbPbSystErr->GetBinContent(i));
  }

  TGraphAsymmErrors *gV2pPb = new TGraphAsymmErrors(histoV2pPb->GetNbinsX());
  TGraphAsymmErrors *gV2pPbSyst = new TGraphAsymmErrors(histoV2pPb->GetNbinsX());
  for (Int_t i = 1; i <= histoV2pPb->GetNbinsX(); i++) 
  {
    if (histoV2pPb->GetBinContent(i) <= 0) continue;
    gV2pPb->SetPoint(i - 1, histoV2pPb->GetBinCenter(i) / 1.6, histoV2pPb->GetBinContent(i));
    cout << "histoV2pPb->GetBinCenter(i) " << histoV2pPb->GetBinCenter(i) / 1.6 << " histoV2pPb->GetBinContent(i) " << histoV2pPb->GetBinContent(i) << endl;
    gV2pPb->SetPointError(i - 1, abs(histoV2pPb->GetBinCenter(i) - histoV2pPb->GetBinLowEdge(i)) / 1.6, abs(histoV2pPb->GetBinCenter(i) - histoV2pPb->GetBinLowEdge(i + 1)) / 1.6, histoV2pPbStatErr->GetBinContent(i), histoV2pPbStatErr->GetBinContent(i));
    gV2pPbSyst->SetPoint(i - 1, histoV2pPb->GetBinCenter(i) / 1.6, histoV2pPb->GetBinContent(i));
    gV2pPbSyst->SetPointError(i - 1, abs(histoV2pPb->GetBinCenter(i) - histoV2pPb->GetBinLowEdge(i)) / 1.6, abs(histoV2pPb->GetBinCenter(i) - histoV2pPb->GetBinLowEdge(i + 1)) / 1.6, histoV2pPbSystErr->GetBinContent(i), histoV2pPbSystErr->GetBinContent(i));
  }

  TGraphAsymmErrors *gV2OO = new TGraphAsymmErrors(numV2OOPubCent);
  TGraphAsymmErrors *gV2OOSyst = new TGraphAsymmErrors(numV2OOPubCent);
  for (Int_t i = 0; i < numV2OOPubCent; i++)
  {
    gV2OO->SetPoint(i, V2OOPubCentMid[i], V2OOPub[i]);
    cout << "V2OOPubCentMid[i] " << V2OOPubCentMid[i] << " V2OOPub[i] " << V2OOPub[i] << endl;
    gV2OO->SetPointError(i, V2OOPubCentMidErr[i], V2OOPubCentMidErr[i], V2OOPubErrStat[i], V2OOPubErrStat[i]);
    gV2OOSyst->SetPoint(i, V2OOPubCentMid[i], V2OOPub[i]);
    gV2OOSyst->SetPointError(i, V2OOPubCentMidErr[i], V2OOPubCentMidErr[i], V2OOPubErrSyst[i], V2OOPubErrSyst[i]);
  }

  TCanvas *canvasV2VsMultiplicity = new TCanvas("canvasV2VsMultiplicity", "canvasV2VsMultiplicity", 900, 700);
  StyleCanvas(canvasV2VsMultiplicity, 0.06, 0.15, 0.15, 0.03);
  canvasV2VsMultiplicity->cd();
  gPad->SetLogx();
  gStyle->SetOptStat(0);

  TH1F *hDummyVsMultiplicity = new TH1F("hDummyVsMultiplicity", "hDummyVsMultiplicity", 1500, 0, 1500);
  for (Int_t i = 1; i <= hDummyVsMultiplicity->GetNbinsX(); i++)
    hDummyVsMultiplicity->SetBinContent(i, -1000);
  canvasV2VsMultiplicity->cd();
  SetFont(hDummyVsMultiplicity);
  TString titledNdeta = "#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";
  StyleHistoYield(hDummyVsMultiplicity, 0, 0.15, 1, 1, titledNdeta, "v_{2}", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummyVsMultiplicity, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummyVsMultiplicity, tickX, tickY);
  hDummyVsMultiplicity->GetXaxis()->SetRangeUser(7, 4000);
  hDummyVsMultiplicity->GetYaxis()->SetTitleOffset(1.7);
  hDummyVsMultiplicity->Draw("");

  gV2PbPb->SetLineColor(kAzure - 3);
  gV2PbPb->SetMarkerColor(kAzure - 3);
  gV2PbPb->SetMarkerStyle(20);
  gV2PbPb->SetMarkerSize(1.5);
  gV2PbPb->Draw("same p");
  gV2PbPbSyst->SetMarkerStyle(20);
  gV2PbPbSyst->SetFillStyle(0);
  gV2PbPbSyst->SetMarkerColor(kAzure - 3);
  gV2PbPbSyst->SetFillColor(kAzure - 3);
  gV2PbPbSyst->SetLineColor(kAzure - 3);
  gV2PbPbSyst->Draw("same p2");
  gV2pPb->SetLineColor(kBlack);
  gV2pPb->SetMarkerColor(kBlack);
  gV2pPb->SetMarkerStyle(21);
  gV2pPb->SetMarkerSize(1.5);
  //gV2pPb->Draw("same p");
  gV2pPbSyst->SetMarkerStyle(21);
  gV2pPbSyst->SetFillStyle(0);
  gV2pPbSyst->SetMarkerColor(kBlack);
  gV2pPbSyst->SetFillColor(kBlack);
  gV2pPbSyst->SetLineColor(kBlack);
  //gV2pPbSyst->Draw("same p2");
  gV2OO->SetLineColor(kMagenta +1);
  gV2OO->SetMarkerColor(kMagenta + 1);
  gV2OO->SetMarkerStyle(20);
  gV2OO->SetMarkerSize(1.5);
  gV2OO->Draw("same p");
  gV2OOSyst->SetMarkerStyle(20);
  gV2OOSyst->SetFillStyle(0);
  gV2OOSyst->SetMarkerColor(kMagenta + 1);
  gV2OOSyst->SetFillColor(kMagenta + 1);
  gV2OOSyst->SetLineColor(kMagenta + 1);
  gV2OOSyst->Draw("same p2");

  TLegend *legendSystem = new TLegend(0.2, 0.7, 0.60, 0.9);
  legendSystem->SetFillStyle(0);
  legendSystem->SetTextAlign(12);
  legendSystem->SetTextSize(0.035);
  legendSystem->AddEntry("", "V_{2}{2, |#Delta#it{#eta} | > 1.4}, |#it{#eta} | < 0.8, 0.2 < #it{p}_{T} < 3 GeV/#it{c}", "");
  legendSystem->AddEntry(gV2PbPb, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  //legendSystem->AddEntry(gV2pPb, "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  legendSystem->AddEntry(gV2OO, "O-O #sqrt{#it{s}_{NN}} = 5.36 TeV");

  legendSystem->Draw("");
  canvasV2VsMultiplicity->SaveAs("../V2VsMultiplicity.pdf");
  canvasV2VsMultiplicity->SaveAs("../V2VsMultiplicity.png");
  canvasV2VsMultiplicity->SaveAs("../V2VsMultiplicity.eps");

}
