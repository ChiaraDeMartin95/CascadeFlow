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
#include "TGraphErrors.h"

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
void StyleCanvas(TCanvas *canvas, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(LMargin);
  canvas->SetRightMargin(RMargin);
  canvas->SetTopMargin(TMargin);
  canvas->SetBottomMargin(BMargin);
  gStyle->SetOptStat(0);
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

TString NameRun2[5] = {"10-20 %", "20-30 %", "30-40 %", "40-50 %", ""};
TString NameRun1[5] = {"10-20 %", "20-30 %", "30-40 %", "40-50 %", "50-60 %"};
Int_t Color[12] = {kYellow + 2, kGreen + 2, kBlack, kBlue + 1, kYellow + 2, kGreen + 2, kBlack, kBlue + 1, kYellow + 2, kGreen + 2, kBlack, kBlue + 1};
Int_t Marker[] = {1, 25, 30, 27, 24, 25};
Int_t MarkerRun3[12] = {20, 21, 29, 33, 20, 21, 29, 33, 20, 21, 29, 33};
Float_t MarkerSize[] = {1.5, 1.5, 1.5, 2., 1.5, 1.5, 1.5, 2., 1.5, 1.5, 1.5, 2.};

void CompareWPublished(Bool_t isRun2Comparison = 1, // 0 for Run1 comparison, 1 for Run2 comparison
                       Bool_t isXi = ChosenParticleXi,
                       TString inputFileName = SinputFileName,
                       Bool_t UseTwoGauss = ExtrUseTwoGauss,
                       Int_t BkgType = ExtrBkgType)
{

  if (isRun2Comparison && numPtBins != 6)
  {
    cout << "The pt binning you chose is not the one used in Run 2, please change in CommonVar.h file" << endl;
    return;
  }
  TString SinputFile = "";
  TFile *inputFile;
  TH1F *hV2C[numCent];
  TH1F *hV2CRatio[numCent];

  TH1F *histoDummy = new TH1F("histoDummy", "histoDummy", 100, 0, 6);
  StyleHisto(histoDummy, -0.2 + 1e-4, 0.8 - 1e-4, 1, 20, "#it{p}_{T} (GeV/#it{c})", "v_{2}", "", kFALSE, 0, 100, 1, 1, 0.05);
  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 100, 0, 6);
  StyleHisto(hDummyRatio, 0 + 1e-4, 2 - 1e-4, 1, 20, "#it{p}_{T} (GeV/#it{c})", "", "", kFALSE, 0, 100, 1, 1, 0.05);
  hDummyRatio->GetXaxis()->SetLabelSize(0.1);
  hDummyRatio->GetYaxis()->SetLabelSize(0.1);
  hDummyRatio->GetXaxis()->SetTitleSize(0.1);
  hDummyRatio->GetYaxis()->SetTitleSize(0.1);
  hDummyRatio->GetYaxis()->SetTitleOffset(0.5);

  for (Int_t cent = 0; cent < numCent; cent++)
  // for (Int_t cent = 0; cent < 4; cent++)
  {
    SinputFile = "OutputAnalysis/FitV2_" + inputFileName + "_" + ParticleName[!isXi];
    SinputFile += IsOneOrTwoGauss[UseTwoGauss];
    SinputFile += SIsBkgParab[BkgType];
    SinputFile += Form("_Cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]);
    if (isApplyWeights)
      SinputFile += "_Weighted";
    if (v2type == 1)
      SinputFile += "_SP";
    if (!useCommonBDTValue)
      SinputFile += "_BDTCentDep";
    if (isRun2Binning)
      SinputFile += "_Run2Binning";
    SinputFile += ".root";
    inputFile = new TFile(SinputFile);
    cout << "Input file with Run 3 v2: " << SinputFile << endl;

    hV2C[cent] = (TH1F *)inputFile->Get("histoV2Mixed");
    hV2C[cent]->SetName(Form("hV2C_%i", cent));
    hV2C[cent]->SetMarkerStyle(MarkerRun3[cent]);
    hV2C[cent]->SetMarkerSize(MarkerSize[cent]);
    hV2C[cent]->SetMarkerColor(ColorMult[cent]);
    hV2C[cent]->SetLineColor(ColorMult[cent]);
    hV2C[cent]->SetLineWidth(2);
    hV2CRatio[cent] = (TH1F *)hV2C[cent]->Clone(Form("hV2CRatio_%i", cent));
  }

  TString SPublishedFileRun2 = "Run2Results/HEPData-ins2093750-v1-root.root";
  TString SPublishedFileRun1 = "Run1Results/HEPData-ins1297103-v1-root.root";
  // in classes 10-20%, 20-30%, 30-40%, 40-50%
  //TString TablesRun2[4] = {"Table 8", "Table 17", "Table 26", "Table 35"}; // v2 measured with 2-particle correlations
  TString TablesRun2[4] = {"Table 79", "Table 87", "Table 95", "Table 103"}; // mean value of v2
  if (!isXi)
  {
    TablesRun2[0] = "Table 80";
    TablesRun2[1] = "Table 88";
    TablesRun2[2] = "Table 96";
    TablesRun2[3] = "Table 104";
  }
  // in classes 10-20%, 20-30%, 30-40%, 40-50%, 50-60%
  TString TablesRun1[5] = {"Table 50", "Table 51", "Table 52", "Table 53", "Table 54"}; // v2 from SP, !y! < 0.5
  if (!isXi)
  {
    TablesRun1[0] = "Table 56";
    TablesRun1[1] = "Table 57";
    TablesRun1[2] = "Table 58";
    TablesRun1[3] = "Table 59";
    TablesRun1[4] = "Table 60";
  }

  TString SPublishedFile = "";
  if (isRun2Comparison)
    SPublishedFile = SPublishedFileRun2;
  else
    SPublishedFile = SPublishedFileRun1;
  TFile *PublishedFile = new TFile(SPublishedFile);

  TString NamePublished[5];
  for (Int_t i = 0; i < 5; i++)
  {
    if (isRun2Comparison)
      NamePublished[i] = NameRun2[i];
    else
      NamePublished[i] = NameRun1[i];
  }

  TDirectoryFile *dirPublished;
  TGraphErrors *gV2Run2[5];
  for (Int_t i = 0; i < 5; i++)
  {
    if (isRun2Comparison && i == 4)
      continue;
    if (isRun2Comparison)
      dirPublished = (TDirectoryFile *)PublishedFile->Get(TablesRun2[i]);
    else
      dirPublished = (TDirectoryFile *)PublishedFile->Get(TablesRun1[i]);
    gV2Run2[i] = (TGraphErrors *)dirPublished->Get("Graph1D_y1");
    gV2Run2[i]->SetName(NameRun2[i]);
    gV2Run2[i]->SetMarkerStyle(Marker[i + 1]);
    gV2Run2[i]->SetMarkerSize(MarkerSize[i + 1]);
    gV2Run2[i]->SetMarkerColor(ColorMult[i + 1]);
    gV2Run2[i]->SetLineColor(ColorMult[i + 1]);
    gV2Run2[i]->SetLineWidth(2);
  }

  TLegend *LegendTitle = new TLegend(0.18, 0.81, 0.42, 0.91); // 0.15
  LegendTitle->SetBorderSize(0);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextSize(0.04);
  LegendTitle->SetMargin(0.);

  TLegend *legendRun3 = new TLegend(0.18, 0.57, 0.42, 0.77); // 0.15
  legendRun3->SetBorderSize(0);
  legendRun3->SetFillStyle(0);
  legendRun3->SetTextSize(0.03);
  legendRun3->SetHeader("PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV");

  TLegend *legendRun2 = new TLegend(0.43, 0.57, 0.6, 0.77); // 0.4
  legendRun2->SetBorderSize(0);
  legendRun2->SetFillStyle(0);
  legendRun2->SetTextSize(0.03);
  legendRun2->SetHeader("PbPb, #sqrt{#it{s}_{NN}} = 5.02 TeV, JHEP 05 (2023) 243");

  TFile *file = new TFile("OutputAnalysis/CompareWPublished_" + inputFileName + "_" + ParticleName[!isXi] + ".root", "RECREATE");
  TCanvas *canvasvsRun2 = new TCanvas("canvasvsRun2", "canvasvsRun2", 1200, 800);
  StyleCanvas(canvasvsRun2, 0.1, 0.05, 0.05, 0.15);
  gStyle->SetOptStat(0);
  canvasvsRun2->cd();
  histoDummy->Draw("");

  for (Int_t i = 0; i < numCent; i++)
  {
    if (i < 4)
    {
      gV2Run2[i]->DrawClone("same P");
      legendRun2->AddEntry(gV2Run2[i], NameRun2[i], "P");
    }
    if (CentFT0C[i] == 10 || CentFT0C[i] == 20 || CentFT0C[i] == 30 || CentFT0C[i] == 40)
    {
      legendRun3->AddEntry(hV2C[i], Form("%i-%i %%", CentFT0C[i], CentFT0C[i + 1]), "P");
      hV2C[i]->Draw("same");
    }
    if (!isRun2Comparison)
    {
      if (CentFT0C[i] == 50)
      {
        legendRun3->AddEntry(hV2C[i], Form("%i-%i %%", CentFT0C[i], CentFT0C[i + 1]), "P");
        hV2C[i]->Draw("same");
      }
      if (i == 4)
      {
        legendRun2->AddEntry(gV2Run2[4], NameRun1[4], "P");
        gV2Run2[4]->DrawClone("same P");
      }
    }
  }
  LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  LegendTitle->AddEntry("", ParticleNameLegend[!isXi] + ", |#it{#eta}| < 0.8", "");

  LegendTitle->Draw();
  legendRun3->Draw();
  legendRun2->Draw();

  canvasvsRun2->SaveAs("OutputAnalysis/CompareWPublished_" + inputFileName + "_" + ParticleName[!isXi] + ".pdf");
  canvasvsRun2->SaveAs("OutputAnalysis/CompareWPublished_" + inputFileName + "_" + ParticleName[!isXi] + ".png");
  //--------------------------------------
  if (!isXi)
  {
    canvasvsRun2->Write();
    file->Close();
    cout << "I created the file " << file->GetName() << endl;
    return;
  }
  // Ratio Run3 / Run2 for Xi only (for Omega, the pt binning is different)
  Int_t index = 0;
  for (Int_t i = 0; i < numCent; i++)
  {
    if (CentFT0C[i] == 10 || CentFT0C[i] == 20 || CentFT0C[i] == 30 || CentFT0C[i] == 40 || CentFT0C[i] == 50)
    {
      cout << "CentFT0C[i]: " << CentFT0C[i] << endl;
      if (CentFT0C[i] == 10)
        index = 0;
      if (CentFT0C[i] == 20)
        index = 1;
      if (CentFT0C[i] == 30)
        index = 2;
      if (CentFT0C[i] == 40)
        index = 3;
      if (CentFT0C[i] == 50)
        index = 4;
      if (isRun2Comparison && index == 4)
        continue;
      for (Int_t b = 1; b <= hV2CRatio[i]->GetNbinsX(); b++)
      {
        hV2CRatio[i]->SetBinContent(b, hV2CRatio[i]->GetBinContent(b) / gV2Run2[index]->GetPointY(b - 1));

        // cout << "V2 Run3: " << hV2C[i]->GetBinCenter(b) << endl;
        // cout << "V2 Run2: " << gV2Run2[index]->GetPointX(b - 1) << endl;
        // cout << "V2 Run3/Run2: " << hV2CRatio[i]->GetBinContent(b) << " +/- " << hV2CRatio[i]->GetBinError(b) << endl;
        hV2CRatio[i]->SetBinError(b, sqrt(pow(hV2C[i]->GetBinError(b) / hV2C[i]->GetBinContent(b), 2) + pow(gV2Run2[index]->GetErrorY(b - 1) / gV2Run2[index]->GetPointY(b - 1), 2)) * hV2CRatio[i]->GetBinContent(b));
      }
    }
  }

  // Plot + ratio Run3/Run2
  TCanvas *canvaswRatio = new TCanvas("canvaswRatio", "canvaswRatio", 700, 900);
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);
  StylePad(pad1, 0.13, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.13, 0.01, 0.02, 0.3); // L, R, T, B

  canvaswRatio->cd();
  pad1->Draw();
  pad1->cd();
  histoDummy->Draw("same");

  for (Int_t i = 0; i < numCent; i++)
  {
    if (i < 4)
    {
      gV2Run2[i]->Draw("same P");
    }
    if (CentFT0C[i] == 10 || CentFT0C[i] == 20 || CentFT0C[i] == 30 || CentFT0C[i] == 40)
    {
      hV2C[i]->Draw("same");
    }
    if (!isRun2Comparison)
    {
      if (CentFT0C[i] == 50)
      {
        hV2C[i]->Draw("same");
      }
      if (i == 4)
      {
        gV2Run2[i]->Draw("same P");
      }
    }
  }
  LegendTitle->Draw();
  legendRun3->Draw();
  legendRun2->Draw();

  canvaswRatio->cd();
  padL1->Draw();
  padL1->cd();
  // SetFont(hDummyRatio);
  // SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  // SetTickLength(hDummyRatio, tickX, tickY);
  hDummyRatio->Draw("same");
  for (Int_t i = 0; i < numCent; i++)
  {
    if (CentFT0C[i] == 10 || CentFT0C[i] == 20 || CentFT0C[i] == 30 || CentFT0C[i] == 40)
    {
      hV2CRatio[i]->Draw("same");
    }
  }
  TF1 *f1 = new TF1("f1", "1", 0, 6);
  f1->SetLineColor(kBlack);
  f1->SetLineStyle(2);
  f1->Draw("same");
  canvaswRatio->SaveAs("OutputAnalysis/CompareWPublished_" + inputFileName + "_" + ParticleName[!isXi] + "_Ratio.pdf");
  canvaswRatio->SaveAs("OutputAnalysis/CompareWPublished_" + inputFileName + "_" + ParticleName[!isXi] + "_Ratio.png");
  canvasvsRun2->Write();
  file->Close();
  cout << "I created the file " << file->GetName() << endl;
}
