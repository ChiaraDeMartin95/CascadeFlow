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
void QCPlots(TString inputFileName = SinputFileName, Bool_t isAfterEPSel = 1)
{

  gStyle->SetOptStat(0);

  TString SAfterEventsel = "";
  if (isAfterEPSel)
    SAfterEventsel = "AfterEP";

  // TString SinputFile = "TreeForAnalysis/AnalysisResults_" + SinputFileName;
  TString SinputFile = "TreeForAnalysis/AnalysisResults_3004.root";
  TFile *file = new TFile(SinputFile, "READ");
  TDirectoryFile *dir = (TDirectoryFile *)file->Get("lf-cascade-flow");

  // event selection
  TH1F *hNEvents = (TH1F *)dir->Get("hNEvents");
  TH2F *hPVContribvsFT0CBefSel = (TH2F *)dir->Get("hEventPVcontributorsVsCentralityBefCuts");
  TH2F *hPVContribvsFT0C = (TH2F *)dir->Get("hEventPVcontributorsVsCentrality" + SAfterEventsel);
  TH2F *hGlobalTrkvsFT0CBefSel = (TH2F *)dir->Get("hEventGlobalTracksVsCentralityBefCuts");
  TH2F *hGlobalTrkvsFT0C = (TH2F *)dir->Get("hEventGlobalTracksVsCentrality" + SAfterEventsel);
  TH2F *hGloablTrkvsPVContribBefSel = (TH2F *)dir->Get("hEventNchCorrelationBefCuts");
  TH2F *hGloablTrkvsPVContrib = (TH2F *)dir->Get("hEventNchCorrelation" + SAfterEventsel);

  // event plane vs FT0C
  TH2F *hPsiT0CvsFT0C = (TH2F *)dir->Get("hPsiT0CvsCentFT0C");

  // v2 vs FT0C
  TH2F *hv2CEPvsFT0C = (TH2F *)dir->Get("hv2CEPvsFT0C");

  TCanvas *c[11];
  for (int i = 0; i < 11; i++)
  {
    c[i] = new TCanvas(Form("c%d", i), Form("c%d", i), 800, 600);
    if (i == 0)
      StyleCanvas(c[i], 0.05, 0.1, 0.15, 0.05);
    else
      StyleCanvas(c[i], 0.05, 0.1, 0.10, 0.1);
  }
  c[0]->cd();
  hNEvents->Scale(1. / hNEvents->GetBinContent(1));
  StyleHisto(hNEvents, 0, 1.2, kBlack, 1, "", "Fraction of selected events", "");
  hNEvents->Draw("l");
  c[0]->SaveAs("QCPlots/hNEvents.pdf");
  c[0]->SaveAs("QCPlots/hNEvents.png");

  c[1]->cd();
  hPVContribvsFT0CBefSel->GetXaxis()->SetTitle("FT0C(%)");
  hPVContribvsFT0CBefSel->GetYaxis()->SetTitle("PV contributors");
  hPVContribvsFT0CBefSel->SetTitle("");
  hPVContribvsFT0CBefSel->Draw("colz");
  c[1]->SaveAs("QCPlots/hPVContribvsFT0C_BefSel.pdf");
  c[1]->SaveAs("QCPlots/hPVContribvsFT0C_BefSel.png");

  c[2]->cd();
  hPVContribvsFT0C->GetXaxis()->SetTitle("FT0C(%)");
  hPVContribvsFT0C->GetYaxis()->SetTitle("PV contributors");
  hPVContribvsFT0C->SetTitle("");
  hPVContribvsFT0C->Draw("colz");
  c[2]->SaveAs("QCPlots/hPVContribvsFT0C" + SAfterEventsel + ".pdf");
  c[2]->SaveAs("QCPlots/hPVContribvsFT0C" + SAfterEventsel + ".png");

  c[3]->cd();
  hGlobalTrkvsFT0CBefSel->GetXaxis()->SetTitle("FT0C(%)");
  hGlobalTrkvsFT0CBefSel->GetYaxis()->SetTitle("Global tracks");
  hGlobalTrkvsFT0CBefSel->SetTitle("");
  hGlobalTrkvsFT0CBefSel->Draw("colz");
  c[3]->SaveAs("QCPlots/hGlobalTrkvsFT0C_BefSel.pdf");
  c[3]->SaveAs("QCPlots/hGlobalTrkvsFT0C_BefSel.png");

  c[4]->cd();
  hGlobalTrkvsFT0C->GetXaxis()->SetTitle("FT0C(%)");
  hGlobalTrkvsFT0C->GetYaxis()->SetTitle("Global tracks");
  hGlobalTrkvsFT0C->SetTitle("");
  hGlobalTrkvsFT0C->Draw("colz");
  c[4]->SaveAs("QCPlots/hGlobalTrkvsFT0C" + SAfterEventsel + ".pdf");
  c[4]->SaveAs("QCPlots/hGlobalTrkvsFT0C" + SAfterEventsel + ".png");

  c[5]->cd();
  hGloablTrkvsPVContribBefSel->GetXaxis()->SetTitle("PV contributors");
  hGloablTrkvsPVContribBefSel->GetYaxis()->SetTitle("Global tracks");
  hGloablTrkvsPVContribBefSel->SetTitle("");
  hGloablTrkvsPVContribBefSel->Draw("colz");
  c[5]->SaveAs("QCPlots/hGloablTrkvsPVContrib_BefSel.pdf");
  c[5]->SaveAs("QCPlots/hGloablTrkvsPVContrib_BefSel.png");

  c[6]->cd();
  hGloablTrkvsPVContrib->GetXaxis()->SetTitle("PV contributors");
  hGloablTrkvsPVContrib->GetYaxis()->SetTitle("Global tracks");
  hGloablTrkvsPVContrib->SetTitle("");
  hGloablTrkvsPVContrib->Draw("colz");
  c[6]->SaveAs("QCPlots/hGloablTrkvsPVContrib" + SAfterEventsel + ".pdf");
  c[6]->SaveAs("QCPlots/hGloablTrkvsPVContrib" + SAfterEventsel + ".png");

  c[7]->cd();
  hPsiT0CvsFT0C->GetXaxis()->SetTitle("FT0C(%)");
  hPsiT0CvsFT0C->GetYaxis()->SetTitle("#Psi_{T0C}");
  hPsiT0CvsFT0C->SetTitle("");
  hPsiT0CvsFT0C->Draw("colz");
  c[7]->SaveAs("QCPlots/hPsiT0CvsFT0C.pdf");
  c[7]->SaveAs("QCPlots/hPsiT0CvsFT0C.png");

  TH1D *hPsiT0C[numCent];
  c[8]->cd();
  for (Int_t mult = 0; mult < numCent; mult++)
  {
    if (mult == 0 || mult > (numCent - 3))
      continue;
    hPsiT0C[mult] = hPsiT0CvsFT0C->ProjectionY(Form("hPsiT0C%d", mult), hPsiT0CvsFT0C->GetXaxis()->FindBin(CentFT0C[mult] + 0.001), hPsiT0CvsFT0C->GetXaxis()->FindBin(CentFT0C[mult + 1] - 0.001));
    hPsiT0C[mult]->Scale(1. / hPsiT0C[mult]->Integral());
    hPsiT0C[mult]->SetLineColor(ColorMult[mult]);
    hPsiT0C[mult]->SetMarkerColor(ColorMult[mult]);
    hPsiT0C[mult]->SetMarkerStyle(MarkerMult[mult]);
    hPsiT0C[mult]->SetMarkerSize(SizeMult[mult]);
    hPsiT0C[mult]->GetXaxis()->SetTitle("#Psi_{T0C}");
    hPsiT0C[mult]->Draw("same hist");
  }

  TLegend *leg = new TLegend(0.35, 0.2, 0.65, 0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.03);
  for (Int_t mult = 0; mult < numCent; mult++)
  {
    if (mult == 0 || mult > (numCent - 3))
      continue;
    leg->AddEntry(hPsiT0C[mult], Form("%d-%d%%", CentFT0C[mult], CentFT0C[mult + 1]), "lp");
  }
  leg->Draw();
  c[8]->SaveAs("QCPlots/hPsiT0C.pdf");
  c[8]->SaveAs("QCPlots/hPsiT0C.png");

  c[9]->cd();
  hv2CEPvsFT0C->GetXaxis()->SetTitle("FT0C(%)");
  hv2CEPvsFT0C->GetYaxis()->SetTitle("v_{2}");
  hv2CEPvsFT0C->SetTitle("");
  hv2CEPvsFT0C->Draw("colz");
  c[9]->SaveAs("QCPlots/hv2CEPvsFT0C.pdf");
  c[9]->SaveAs("QCPlots/hv2CEPvsFT0C.pdg");

  TH1D *hv2CEP[numCent];
  c[10]->cd();
  for (Int_t mult = 0; mult < numCent; mult++)
  {
    if (mult == 0 || mult > (numCent - 3))
      continue;
    hv2CEP[mult] = hv2CEPvsFT0C->ProjectionY(Form("hv2CEP%d", mult), hv2CEPvsFT0C->GetXaxis()->FindBin(CentFT0C[mult] + 0.001), hv2CEPvsFT0C->GetXaxis()->FindBin(CentFT0C[mult + 1] - 0.001));
    hv2CEP[mult]->Scale(1. / hv2CEP[mult]->Integral());
    hv2CEP[mult]->SetLineColor(ColorMult[mult]);
    hv2CEP[mult]->SetMarkerColor(ColorMult[mult]);
    hv2CEP[mult]->SetMarkerStyle(MarkerMult[mult]);
    hv2CEP[mult]->SetMarkerSize(SizeMult[mult]);
    hv2CEP[mult]->GetXaxis()->SetTitle("v_{2}");
    hv2CEP[mult]->Draw("same hist");
  }
  TLegend *leg2 = new TLegend(0.6, 0.6, 0.85, 0.9);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.03);
  for (Int_t mult = 0; mult < numCent; mult++)
  {
    if (mult == 0 || mult > (numCent - 3))
      continue;
    leg2->AddEntry(hPsiT0C[mult], Form("%d-%d%%", CentFT0C[mult], CentFT0C[mult + 1]), "lp");
  }
  leg2->Draw();

  TString OutputFile = "QCPlots/QCPlots_" + inputFileName;
  for (Int_t i = 0; i < 11; i++)
  {
    if (i == 0)
      c[i]->SaveAs(OutputFile + ".pdf(");
    else if (i == 10)
      c[i]->SaveAs(OutputFile + ".pdf)");
    else
      c[i]->SaveAs(OutputFile + ".pdf");
  }

  cout << "\nStarting from the files (for the different mult): " << SinputFile << endl;
}
