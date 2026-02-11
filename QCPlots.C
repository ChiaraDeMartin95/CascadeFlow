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
void QCPlots(Bool_t isEff = 0, Bool_t isAfterEPSel = 0)
{

  TString inputFileName = SinputFileName;
  //if (isEff)
  //  inputFileName = SinputFileNameEff;
  gStyle->SetOptStat(0);
  Int_t nrebinx = 4;
  Int_t nrebiny = 4;
  const Int_t nCanvas = 17;

  TString SAfterEventsel = "";
  if (isAfterEPSel)
    SAfterEventsel = "AfterEP";

  TString SinputFile = "../TreeForAnalysis";
  if (isEff)
    SinputFile = "FileForEfficiency";
  SinputFile += "/AnalysisResults_" + inputFileName + ".root";
  cout << "Input file: " << SinputFile << endl;
  // TString SinputFile = "TreeForTrainingSignal/AnalysisResults_" + SinputFileName + ".root";
  // TString SinputFile = "TreeForAnalysis/AnalysisResults_3004.root";
  //  TString SinputFile = "TreeForAnalysis/AnalysisResults_LHC23zzh_pass3_DerivedStrangeness_Train205658_Test.root";
  TFile *file = new TFile(SinputFile, "READ");
  if (!file)
  {
    std::cout << "Error opening file: " << SinputFile << std::endl;
    return;
  }
  TDirectoryFile *dir = (TDirectoryFile *)file->Get("lf-cascade-flow");
  if (!dir)
  {
    std::cout << "Directory lf-cascade-flow not found!" << std::endl;
    return;
  }

  TDirectoryFile *dirHistos;
  dirHistos = (TDirectoryFile *)dir->Get("histos");
  if (!dirHistos)
  {
    std::cout << "Directory histos not found!" << std::endl;
    return;
  }
  TH1F *hNEvents = (TH1F *)dirHistos->Get("hNEvents");
  TH2F *hPVContribvsFT0CBefSel = (TH2F *)dirHistos->Get("hEventPVcontributorsVsCentralityBefCuts");
  TH2F *hPVContribvsFT0C = (TH2F *)dirHistos->Get("hEventPVcontributorsVsCentrality" + SAfterEventsel);
  TH2F *hGlobalTrkvsFT0CBefSel = (TH2F *)dirHistos->Get("hEventGlobalTracksVsCentralityBefCuts");
  TH2F *hGlobalTrkvsFT0C = (TH2F *)dirHistos->Get("hEventGlobalTracksVsCentrality" + SAfterEventsel);
  TH2F *hGlobalTrkvsPVContribBefSel = (TH2F *)dirHistos->Get("hEventNchCorrelationBefCuts");
  TH2F *hGlobalTrkvsPVContrib = (TH2F *)dirHistos->Get("hEventNchCorrelation" + SAfterEventsel);
  TH1F *hCentrality = (TH1F *)dirHistos->Get("hEventCentrality");
  TH1F *hVertexZ = (TH1F *)dirHistos->Get("hEventVertexZ");

  // event plane vs FT0C
  TH2F *hPsiT0CvsFT0C = (TH2F *)dirHistos->Get("hPsiT0CvsCentFT0C");

  // event planes after shift correction (only available in recent files)
  TH2F *hPsiV0AvsFT0C_Shifted = (TH2F *)dirHistos->Get("Psi_EP_FV0A_shifted");
  if (!hPsiV0AvsFT0C_Shifted)
  {
    cout << "No shifted EP histos found, taking in input a random histo" << endl;
    hPsiV0AvsFT0C_Shifted = (TH2F *)hPsiT0CvsFT0C->Clone("hPsiT0CvsCentFT0C_dummy1");
    // return;
  }
  TH2F *hPsiT0AvsFT0C_Shifted = (TH2F *)dirHistos->Get("Psi_EP_FT0A_shifted");
  if (!hPsiT0AvsFT0C_Shifted)
  {
    cout << "No shifted EP histos found, taking in input a random histo" << endl;
    hPsiT0AvsFT0C_Shifted = (TH2F *)hPsiT0CvsFT0C->Clone("hPsiT0CvsCentFT0C_dummy1");
    // return;
  }
  TH2F *hPsiTPCLvsFT0C_Shifted = (TH2F *)dirHistos->Get("Psi_EP_TPCA_shifted");
  if (!hPsiTPCLvsFT0C_Shifted)
  {
    cout << "No shifted EP histos found, taking in input a random histo" << endl;
    hPsiTPCLvsFT0C_Shifted = (TH2F *)hPsiT0CvsFT0C->Clone("hPsiT0CvsCentFT0C_dummy2");
    // return;
  }
  TH2F *hPsiTPCRvsFT0C_Shifted = (TH2F *)dirHistos->Get("Psi_EP_TPCC_shifted");
  if (!hPsiTPCRvsFT0C_Shifted)
  {
    cout << "No shifted EP histos found, taking in input a random histo" << endl;
    hPsiTPCRvsFT0C_Shifted = (TH2F *)hPsiT0CvsFT0C->Clone("hPsiT0CvsCentFT0C_dummy3");
    // return;
  }

  // v2 vs FT0C
  TH2F *hv2CEPvsFT0C = (TH2F *)dirHistos->Get("hv2CEPvsFT0C");

  TCanvas *c[nCanvas];
  for (int i = 0; i < nCanvas; i++)
  {
    c[i] = new TCanvas(Form("c%d", i), Form("c%d", i), 800, 600);
    if (i == 0)
      StyleCanvas(c[i], 0.05, 0.1, 0.15, 0.05);
    else
      StyleCanvas(c[i], 0.05, 0.1, 0.10, 0.1);
  }
  c[0]->cd();
  for (Int_t b = 1; b <= hNEvents->GetNbinsX(); b++)
  {
    cout << "Number of events: " << hNEvents->GetXaxis()->GetBinLabel(b) << " : " << hNEvents->GetBinContent(b) << endl;
  }
  hNEvents->Scale(1. / hNEvents->GetBinContent(1));
  StyleHisto(hNEvents, 0, 1.2, kBlack, 1, "", "Fraction of selected events", "");
  hNEvents->Draw("l");
  c[0]->SaveAs("../QCPlots/hNEvents" + inputFileName + ".pdf");
  c[0]->SaveAs("../QCPlots/hNEvents" + inputFileName + ".png");

  cout << "Hello 4" << endl;
  c[1]->cd();
  gPad->SetLogz();
  hPVContribvsFT0CBefSel->GetXaxis()->SetTitle("FT0C(%)");
  hPVContribvsFT0CBefSel->GetYaxis()->SetTitle("PV contributors");
  hPVContribvsFT0CBefSel->SetTitle("");
  hPVContribvsFT0CBefSel->Rebin2D(nrebinx, nrebiny);
  if (isOOCentrality)
  {
    hPVContribvsFT0CBefSel->GetYaxis()->SetRangeUser(0, 700);
  }
  hPVContribvsFT0CBefSel->Draw("colz");
  c[1]->SaveAs("../QCPlots/hPVContribvsFT0C_BefSel" + inputFileName + ".pdf");
  c[1]->SaveAs("../QCPlots/hPVContribvsFT0C_BefSel" + inputFileName + ".png");

  c[2]->cd();
  gPad->SetLogz();
  hPVContribvsFT0C->GetXaxis()->SetTitle("FT0C(%)");
  hPVContribvsFT0C->GetYaxis()->SetTitle("PV contributors");
  hPVContribvsFT0C->SetTitle("");
  hPVContribvsFT0C->Rebin2D(nrebinx, nrebiny);
  if (isOOCentrality)
  {
    hPVContribvsFT0C->GetYaxis()->SetRangeUser(0, 700);
  }
  hPVContribvsFT0C->Draw("colz");
  c[2]->SaveAs("../QCPlots/hPVContribvsFT0C" + SAfterEventsel + inputFileName + ".pdf");
  c[2]->SaveAs("../QCPlots/hPVContribvsFT0C" + SAfterEventsel + inputFileName + ".png");

  c[3]->cd();
  gPad->SetLogz();
  hGlobalTrkvsFT0CBefSel->GetXaxis()->SetTitle("FT0C(%)");
  hGlobalTrkvsFT0CBefSel->GetYaxis()->SetTitle("Global tracks");
  hGlobalTrkvsFT0CBefSel->SetTitle("");
  hGlobalTrkvsFT0CBefSel->Rebin2D(nrebinx, nrebiny);
  if (isOOCentrality)
  {
    hGlobalTrkvsFT0CBefSel->GetYaxis()->SetRangeUser(0, 500);
  }
  hGlobalTrkvsFT0CBefSel->Draw("colz");
  c[3]->SaveAs("../QCPlots/hGlobalTrkvsFT0C_BefSel" + inputFileName + ".pdf");
  c[3]->SaveAs("../QCPlots/hGlobalTrkvsFT0C_BefSel" + inputFileName + ".png");

  c[4]->cd();
  gPad->SetLogz();
  hGlobalTrkvsFT0C->GetXaxis()->SetTitle("FT0C(%)");
  hGlobalTrkvsFT0C->GetYaxis()->SetTitle("Global tracks");
  hGlobalTrkvsFT0C->SetTitle("");
  hGlobalTrkvsFT0C->Rebin2D(nrebinx, nrebiny);
  if (isOOCentrality)
  {
    hGlobalTrkvsFT0C->GetYaxis()->SetRangeUser(0, 500);
  }
  hGlobalTrkvsFT0C->Draw("colz");
  c[4]->SaveAs("../QCPlots/hGlobalTrkvsFT0C" + SAfterEventsel + inputFileName + ".pdf");
  c[4]->SaveAs("../QCPlots/hGlobalTrkvsFT0C" + SAfterEventsel + inputFileName + ".png");

  c[5]->cd();
  gPad->SetLogz();
  hGlobalTrkvsPVContribBefSel->GetXaxis()->SetTitle("PV contributors");
  hGlobalTrkvsPVContribBefSel->GetYaxis()->SetTitle("Global tracks");
  hGlobalTrkvsPVContribBefSel->SetTitle("");
  hGlobalTrkvsPVContribBefSel->Rebin2D(4 * nrebinx, 4 * nrebiny);
  if (isOOCentrality)
  {
    hGlobalTrkvsPVContribBefSel->GetYaxis()->SetRangeUser(0, 500);
    hGlobalTrkvsPVContribBefSel->GetXaxis()->SetRangeUser(0, 700);
  }
  hGlobalTrkvsPVContribBefSel->Draw("colz");
  c[5]->SaveAs("../QCPlots/hGlobalTrkvsPVContrib_BefSel" + SinputFileName + ".pdf");
  c[5]->SaveAs("../QCPlots/hGlobalTrkvsPVContrib_BefSel" + SinputFileName + ".png");

  c[6]->cd();
  gPad->SetLogz();
  hGlobalTrkvsPVContrib->GetXaxis()->SetTitle("PV contributors");
  hGlobalTrkvsPVContrib->GetYaxis()->SetTitle("Global tracks");
  hGlobalTrkvsPVContrib->SetTitle("");
  hGlobalTrkvsPVContrib->Rebin2D(4 * nrebinx, 4 * nrebiny);
  if (isOOCentrality)
  {
    hGlobalTrkvsPVContrib->GetYaxis()->SetRangeUser(0, 500);
    hGlobalTrkvsPVContrib->GetXaxis()->SetRangeUser(0, 700);
  }
  hGlobalTrkvsPVContrib->Draw("colz");
  c[6]->SaveAs("../QCPlots/hGlobalTrkvsPVContrib" + SAfterEventsel + inputFileName + ".pdf");
  c[6]->SaveAs("../QCPlots/hGlobalTrkvsPVContrib" + SAfterEventsel + inputFileName + ".png");

  c[7]->cd();
  hPsiT0CvsFT0C->GetXaxis()->SetTitle("FT0C(%)");
  hPsiT0CvsFT0C->GetYaxis()->SetTitle("#Psi_{T0C}");
  hPsiT0CvsFT0C->SetTitle("");
  hPsiT0CvsFT0C->Draw("colz");
  c[7]->SaveAs("../QCPlots/hPsiT0CvsFT0C" + inputFileName + ".pdf");
  c[7]->SaveAs("../QCPlots/hPsiT0CvsFT0C" + inputFileName + ".png");

  TH1D *hPsiT0C[commonNumCent];
  c[8]->cd();
  for (Int_t mult = 0; mult < commonNumCent; mult++)
  {
    if (isOOCentrality)
    {
      hPsiT0C[mult] = hPsiT0CvsFT0C->ProjectionY(Form("hPsiT0C%d", mult), hPsiT0CvsFT0C->GetXaxis()->FindBin(CentFT0CLambdaOO[mult] + 0.001), hPsiT0CvsFT0C->GetXaxis()->FindBin(CentFT0CLambdaOO[mult + 1] - 0.001));
    }
    else
      hPsiT0C[mult] = hPsiT0CvsFT0C->ProjectionY(Form("hPsiT0C%d", mult), hPsiT0CvsFT0C->GetXaxis()->FindBin(CentFT0C[mult] + 0.001), hPsiT0CvsFT0C->GetXaxis()->FindBin(CentFT0C[mult + 1] - 0.001));
    hPsiT0C[mult]->Scale(1. / hPsiT0C[mult]->Integral());
    hPsiT0C[mult]->SetLineColor(ColorMult[mult]);
    hPsiT0C[mult]->SetMarkerColor(ColorMult[mult]);
    hPsiT0C[mult]->SetMarkerStyle(MarkerMult[mult]);
    hPsiT0C[mult]->SetMarkerSize(SizeMult[mult]);
    hPsiT0C[mult]->GetXaxis()->SetTitle("#Psi_{T0C}");
    if (mult < (commonNumCent - 2))
      hPsiT0C[mult]->Draw("same hist");
  }

  TLegend *leg = new TLegend(0.35, 0.2, 0.65, 0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.03);
  for (Int_t mult = 0; mult < commonNumCent; mult++)
  {
    // if (mult == 0 || mult > (commonNumCent - 1))
    //   continue;
    if (isOOCentrality)
      leg->AddEntry(hPsiT0C[mult], Form("%d-%d%%", CentFT0CLambdaOO[mult], CentFT0CLambdaOO[mult + 1]), "lp");
    else
      leg->AddEntry(hPsiT0C[mult], Form("%d-%d%%", CentFT0C[mult], CentFT0C[mult + 1]), "lp");
  }
  leg->Draw();
  c[8]->SaveAs("../QCPlots/hPsiT0C" + inputFileName + ".pdf");
  c[8]->SaveAs("../QCPlots/hPsiT0C" + inputFileName + ".png");

  c[9]->cd();
  hv2CEPvsFT0C->GetXaxis()->SetTitle("FT0C(%)");
  hv2CEPvsFT0C->GetYaxis()->SetTitle("v_{2}");
  hv2CEPvsFT0C->SetTitle("");
  hv2CEPvsFT0C->Draw("colz");
  c[9]->SaveAs("../QCPlots/hv2CEPvsFT0C" + inputFileName + ".pdf");
  c[9]->SaveAs("../QCPlots/hv2CEPvsFT0C" + inputFileName + ".pdg");

  TH1D *hv2CEP[commonNumCent];
  c[10]->cd();
  for (Int_t mult = 0; mult < commonNumCent; mult++)
  {
    if (mult == 0 || mult > (commonNumCent - 3))
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
  for (Int_t mult = 0; mult < commonNumCent; mult++)
  {
    if (mult == 0 || mult > (commonNumCent - 3))
      continue;
    leg2->AddEntry(hPsiT0C[mult], Form("%d-%d%%", CentFT0C[mult], CentFT0C[mult + 1]), "lp");
  }
  leg2->Draw();

  c[11]->cd();
  TF1 *pol0 = new TF1("pol0", "pol0", 0, 60);
  hCentrality->Fit("pol0", "R");
  hCentrality->GetXaxis()->SetTitle("FT0C(%)");
  hCentrality->SetTitle("");
  hCentrality->Draw("");
  c[11]->SaveAs("../QCPlots/hCentrality" + inputFileName + ".pdf");
  c[11]->SaveAs("../QCPlots/hCentality" + inputFileName + ".png");

  c[12]->cd();
  hVertexZ->GetXaxis()->SetTitle("z (cm)");
  hVertexZ->SetTitle("");
  hVertexZ->Draw("");
  c[12]->SaveAs("../QCPlots/hVertexZ" + inputFileName + ".pdf");
  c[12]->SaveAs("../QCPlots/hVertexZ" + inputFileName + ".png");

  TH1D *hPsiV0A[commonNumCent];
  c[13]->cd();
  for (Int_t mult = 0; mult < commonNumCent; mult++)
  {
    if (isOOCentrality)
    {
      hPsiV0A[mult] = hPsiV0AvsFT0C_Shifted->ProjectionY(Form("hPsiV0A%d", mult), hPsiV0AvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0CLambdaOO[mult] + 0.001), hPsiV0AvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0CLambdaOO[mult + 1] - 0.001));
    }
    else
      hPsiV0A[mult] = hPsiV0AvsFT0C_Shifted->ProjectionY(Form("hPsiV0A%d", mult), hPsiV0AvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0C[mult] + 0.001), hPsiV0AvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0C[mult + 1] - 0.001));
    hPsiV0A[mult]->Scale(1. / hPsiV0A[mult]->Integral());
    hPsiV0A[mult]->SetLineColor(ColorMult[mult]);
    hPsiV0A[mult]->SetMarkerColor(ColorMult[mult]);
    hPsiV0A[mult]->SetMarkerStyle(MarkerMult[mult]);
    hPsiV0A[mult]->SetMarkerSize(SizeMult[mult]);
    hPsiV0A[mult]->GetXaxis()->SetTitle("#Psi_{V0A}");
    hPsiV0A[mult]->SetTitle("EP V0A");
    hPsiV0A[mult]->SetTitle("");
    if (mult < (commonNumCent - 2))
      hPsiV0A[mult]->Draw("same hist");
  }
  leg->Draw();
  c[13]->SaveAs("../QCPlots/hPsiV0A_Shifted" + inputFileName + ".pdf");
  c[13]->SaveAs("../QCPlots/hPsiV0A_Shifted" + inputFileName + ".png");

  TH1D *hPsiTPCL[commonNumCent];
  c[14]->cd();
  for (Int_t mult = 0; mult < commonNumCent; mult++)
  {
    if (isOOCentrality)
    {
      hPsiTPCL[mult] = hPsiTPCLvsFT0C_Shifted->ProjectionY(Form("hPsiTPCL%d", mult), hPsiTPCLvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0CLambdaOO[mult] + 0.001), hPsiTPCLvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0CLambdaOO[mult + 1] - 0.001));
    }
    else
      hPsiTPCL[mult] = hPsiTPCLvsFT0C_Shifted->ProjectionY(Form("hPsiTPCL%d", mult), hPsiTPCLvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0C[mult] + 0.001), hPsiTPCLvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0C[mult + 1] - 0.001));
    hPsiTPCL[mult]->Scale(1. / hPsiTPCL[mult]->Integral());
    hPsiTPCL[mult]->SetLineColor(ColorMult[mult]);
    hPsiTPCL[mult]->SetMarkerColor(ColorMult[mult]);
    hPsiTPCL[mult]->SetMarkerStyle(MarkerMult[mult]);
    hPsiTPCL[mult]->SetMarkerSize(SizeMult[mult]);
    hPsiTPCL[mult]->GetXaxis()->SetTitle("#Psi_{TPC-L}");
    hPsiTPCL[mult]->SetTitle("EP TPC-L");
    hPsiTPCL[mult]->SetTitle("");
    if (mult < (commonNumCent - 2))
      hPsiTPCL[mult]->Draw("same hist");
  }
  leg->Draw();
  c[14]->SaveAs("../QCPlots/hPsiTPCL_Shifted" + inputFileName + ".pdf");
  c[14]->SaveAs("../QCPlots/hPsiTPCL_Shifted" + inputFileName + ".png");

  TH1D *hPsiTPCR[commonNumCent];
  c[15]->cd();
  for (Int_t mult = 0; mult < commonNumCent; mult++)
  {
    if (isOOCentrality)
    {
      hPsiTPCR[mult] = hPsiTPCRvsFT0C_Shifted->ProjectionY(Form("hPsiTPCR%d", mult), hPsiTPCRvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0CLambdaOO[mult] + 0.001), hPsiTPCRvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0CLambdaOO[mult + 1] - 0.001));
    }
    else
      hPsiTPCR[mult] = hPsiTPCRvsFT0C_Shifted->ProjectionY(Form("hPsiTPCR%d", mult), hPsiTPCRvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0C[mult] + 0.001), hPsiTPCRvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0C[mult + 1] - 0.001));
    hPsiTPCR[mult]->Scale(1. / hPsiTPCR[mult]->Integral());
    hPsiTPCR[mult]->SetLineColor(ColorMult[mult]);
    hPsiTPCR[mult]->SetMarkerColor(ColorMult[mult]);
    hPsiTPCR[mult]->SetMarkerStyle(MarkerMult[mult]);
    hPsiTPCR[mult]->SetMarkerSize(SizeMult[mult]);
    hPsiTPCR[mult]->GetXaxis()->SetTitle("#Psi_{TPC-R}");
    hPsiTPCR[mult]->SetTitle("EP TPC-R");
    hPsiTPCR[mult]->SetTitle("");
    if (mult < (commonNumCent - 2))
      hPsiTPCR[mult]->Draw("same hist");
  }
  leg->Draw();
  c[15]->SaveAs("../QCPlots/hPsiTPCR_Shifted" + inputFileName + ".pdf");
  c[15]->SaveAs("../QCPlots/hPsiTPCR_Shifted" + inputFileName + ".png");

  TH1D *hPsiT0A[commonNumCent];
  c[16]->cd();
  for (Int_t mult = 0; mult < commonNumCent; mult++)
  {
    if (isOOCentrality)
    {
      hPsiT0A[mult] = hPsiT0AvsFT0C_Shifted->ProjectionY(Form("hPsiT0A%d", mult), hPsiT0AvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0CLambdaOO[mult] + 0.001), hPsiT0AvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0CLambdaOO[mult + 1] - 0.001));
    }
    else
      hPsiT0A[mult] = hPsiT0AvsFT0C_Shifted->ProjectionY(Form("hPsiT0A%d", mult), hPsiT0AvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0C[mult] + 0.001), hPsiT0AvsFT0C_Shifted->GetXaxis()->FindBin(CentFT0C[mult + 1] - 0.001));
    hPsiT0A[mult]->Scale(1. / hPsiT0A[mult]->Integral());
    hPsiT0A[mult]->SetLineColor(ColorMult[mult]);
    hPsiT0A[mult]->SetMarkerColor(ColorMult[mult]);
    hPsiT0A[mult]->SetMarkerStyle(MarkerMult[mult]);
    hPsiT0A[mult]->SetMarkerSize(SizeMult[mult]);
    hPsiT0A[mult]->GetXaxis()->SetTitle("#Psi_{T0A}");
    hPsiT0A[mult]->SetTitle("EP T0A");
    hPsiT0A[mult]->SetTitle("");
    if (mult < (commonNumCent - 2))
      hPsiT0A[mult]->Draw("same hist");
  }
  leg->Draw();
  c[16]->SaveAs("../QCPlots/hPsiT0A_Shifted" + inputFileName + ".pdf");
  c[16]->SaveAs("../QCPlots/hPsiT0A_Shifted" + inputFileName + ".png");


  TString OutputFile = "../QCPlots/QCPlots_" + inputFileName;
  for (Int_t i = 0; i < nCanvas; i++)
  {
    if (i == 0)
      c[i]->SaveAs(OutputFile + ".pdf(");
    else if (i == nCanvas - 1)
      c[i]->SaveAs(OutputFile + ".pdf)");
    else
      c[i]->SaveAs(OutputFile + ".pdf");
  }

  TH1F *hCentWeight = new TH1F("hCentWeight", "hCentWeight", 100, 0, 100);
  for (Int_t b = 1; b <= hCentrality->GetNbinsX(); b++)
  {
    if (hCentrality->GetBinContent(b) > 0)
      hCentWeight->SetBinContent(b, 1. / (hCentrality->GetBinContent(b) / pol0->GetParameter(0)));
  }

  TFile *fout = new TFile("../CentralityWeight_" + inputFileName + ".root", "RECREATE");
  hCentWeight->Write();
  TList *list = new TList();
  list->Add(hCentWeight);
  list->Write("ccdb_object", TObject::kSingleKey);
  fout->Close();

  TCanvas *cCentWeight = new TCanvas("cCentWeight", "cCentWeight", 800, 600);
  StyleCanvas(cCentWeight, 0.05, 0.1, 0.15, 0.05);
  cCentWeight->cd();
  hCentWeight->GetXaxis()->SetTitle("FT0C(%)");
  hCentWeight->GetYaxis()->SetTitle("Weight");
  hCentWeight->SetTitle("");
  hCentWeight->Draw("");
  cCentWeight->SaveAs("../QCPlots/hCentWeight_" + inputFileName + ".png");
  cCentWeight->SaveAs("../QCPlots/hCentWeight_" + inputFileName + ".pdf");

  cout << "\nStarting from the files (for the different mult): " << SinputFile << endl;
}
