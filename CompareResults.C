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

Float_t YLow = 0;
Float_t YUp = 0;
Float_t YLowRatio = 0;
Float_t YUpRatio = 0;

TString hTitleX = "";
TString hTitleY = "";

Float_t xTitle = 15;
Float_t xOffset = 4;
Float_t yTitle = 30;
Float_t yOffset = 2;

Float_t xLabel = 30;
Float_t yLabel = 30;
Float_t xLabelOffset = 0.05;
Float_t yLabelOffset = 0.01;

Float_t tickX = 0.03;
Float_t tickY = 0.042;

Float_t LimSupMultRatio = 5.1;
Float_t LimInfMultRatio = 1e-2;
Float_t YoffsetSpectraRatio = 1.1;
Float_t xTitleR = 35;
Float_t xOffsetR = 1;
Float_t yTitleR = 30;
Float_t yOffsetR = 2;

Float_t xLabelR = 25;
Float_t yLabelR = 25;
Float_t xLabelOffsetR = 0.02;
Float_t yLabelOffsetR = 0.04;

TString Sinputfile = "";
TString namehisto[2] = {""};
TString CommonFileName = "";

Int_t numOptions = 0;
TH1F *hDef;
TString fileName[10] = {""};
TH1F *h[10];
TH1F *hRatio[10];
TString sleg[10] = {""};

void CompareResults(Int_t TypeComp = 0, Int_t mult = 0)
{

  gStyle->SetOptStat(0);

  if (mult > numCent)
  {
    cout << "Multiplciity out of range" << endl;
    return;
  }

  // TypeComp = 0 --> weighted vs unweighted v2
  // TypeComp = 1 --> LHC23 vs LHC23zzh
  // TypeComp = 2 --> SP vs EP
  // TypeComp = 3 --> BDT cut 0.7 vs default
  // TypeComp = 4 --> v2 from fit vs v2 from histo
  // TypeComp = 5 --> resolution LF vs CFW
  // TypeComp = 6 --> v2 LF vs CFW
  // TypeComp = 7 --> v2 with and without occupancy event selection
  // TypeComp = 8 --> v2 pass4 vs pass3 (LHCzzh)

  // TypeComp = 0 --> weighted vs unweighted v2
  if (TypeComp == 0)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_" + SinputFileName + "_" + ParticleName[!ChosenParticleXi] + ChargeName[ExtrCharge + 1] + SEtaSysChoice[ExtrEtaSysChoice];
    CommonFileName += IsOneOrTwoGauss[ExtrUseTwoGauss];
    CommonFileName += SIsBkgParab[ExtrBkgType];
    CommonFileName += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] = "";
    fileName[1] = "_Weighted";
    namehisto[0] = "histoV2";
    namehisto[1] = "histoV2";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "Default";
    sleg[1] = "Weighted";
  }
  // TypeComp = 1 --> LHC23 vs LHC23zzh
  else if (TypeComp == 1)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23_PbPb_pass3_Train218607_";
    fileName[1] = "LHC23zzh_pass3_Train224930_";
    fileName[0] += ParticleName[!ChosenParticleXi] + ChargeName[ExtrCharge + 1] + SEtaSysChoice[ExtrEtaSysChoice];
    fileName[0] += IsOneOrTwoGauss[ExtrUseTwoGauss];
    fileName[0] += SIsBkgParab[ExtrBkgType];
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] += ParticleName[!ChosenParticleXi] + ChargeName[ExtrCharge + 1] + SEtaSysChoice[ExtrEtaSysChoice];
    fileName[1] += IsOneOrTwoGauss[ExtrUseTwoGauss];
    fileName[1] += SIsBkgParab[ExtrBkgType];
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";

    namehisto[0] = "histoV2NoFit";
    namehisto[1] = "histoV2NoFit";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "LHC23_PbPb_pass3";
    sleg[1] = "LHC23zzh_pass3";
  }
  // TypeComp = 2 --> SP vs EP
  else if (TypeComp == 2)
  {
    numOptions = 2;
    // CommonFileName = "OutputAnalysis/FitV2_" + SinputFileName + "_" + ParticleName[!ChosenParticleXi] + ChargeName[ExtrCharge + 1] + SEtaSysChoice[ExtrEtaSysChoice];
    // CommonFileName += IsOneOrTwoGauss[ExtrUseTwoGauss];
    // CommonFileName += SIsBkgParab[ExtrBkgType];
    // CommonFileName += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    // CommonFileName += "_Weighted";
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23zzh_pass3_Train226234_CFW_Xi_BkgParab";
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] = fileName[0];
    fileName[0] += "";
    fileName[1] += "_SP";
    fileName[0] += "_Run2Binning";
    fileName[1] += "_Run2Binning";
    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "EP";
    sleg[1] = "SP";
  }
  // TypeComp = 3 --> BDT cut 0.7 vs default
  else if (TypeComp == 3)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23_PbPb_pass3_Train218607_";
    fileName[1] = "LHC23zzh_pass3_Train225737_BDT0.7_";
    fileName[0] += ParticleName[!ChosenParticleXi] + ChargeName[ExtrCharge + 1] + SEtaSysChoice[ExtrEtaSysChoice];
    fileName[0] += IsOneOrTwoGauss[ExtrUseTwoGauss];
    fileName[0] += SIsBkgParab[ExtrBkgType];
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted_BDTCentDep";
    fileName[1] += ParticleName[!ChosenParticleXi] + ChargeName[ExtrCharge + 1] + SEtaSysChoice[ExtrEtaSysChoice];
    fileName[1] += IsOneOrTwoGauss[ExtrUseTwoGauss];
    fileName[1] += SIsBkgParab[ExtrBkgType];
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";

    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "";
    sleg[1] = "BDTcut > 0.97";
  }
  // v2 from fit vs v2 from histo
  else if (TypeComp == 4)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_" + SinputFileName + "_" + ParticleName[!ChosenParticleXi] + ChargeName[ExtrCharge + 1] + SEtaSysChoice[ExtrEtaSysChoice];
    CommonFileName += IsOneOrTwoGauss[ExtrUseTwoGauss];
    CommonFileName += SIsBkgParab[ExtrBkgType];
    CommonFileName += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    CommonFileName += "_Weighted";
    fileName[0] = "";
    fileName[0] = "";
    namehisto[0] = "histoV2";
    namehisto[1] = "histoV2NoFit";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "v2 from fit";
    sleg[1] = "v2 no fit";
  }
  // resolution LF vs CFW
  else if (TypeComp == 5)
  {
    numOptions = 2;
    fileName[0] = "Resolution/Resolution_EP_LF";
    fileName[1] = "Resolution/Resolution_EP_CFW";
    namehisto[0] = "hReso";
    namehisto[1] = "hReso";
    hTitleY = "Resolution";
    hTitleX = "Centrality (%)";
    YLow = 0;
    YUp = 1.8;
    YLowRatio = 0.6;
    YUpRatio = 1.2;
    sleg[0] = "LF";
    sleg[1] = "central FW";
  }
  // v2 LF vs CFW
  else if (TypeComp == 6)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23zzh_pass3_Train226234_CFW_Xi_BkgParab";
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] = "LHC23zzh_pass3_Train224930_Xi_BkgParab";
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";
    // fileName[0] += "_SP";
    // fileName[1] += "_SP";
    fileName[0] += "_Run2Binning";
    fileName[1] += "_Run2Binning";
    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "CFW";
    sleg[1] = "LF";
  }
  // v2 with and without occupancy event selection
  else if (TypeComp == 7)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23_PbPb_pass3_Train218607_Xi_BkgParab";
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] = "LHC23_PbPb_pass3_Train231308_Xi_BkgParab";
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";

    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "w/o occupancy sel.";
    sleg[1] = "w/ occupancy sel.";
  }
  // v2 with pass4 vs pass3
  else if (TypeComp == 8)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23zzh_pass3_Train224930_Xi_BkgParab";
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] = "LHC23zzh_pass4_test3_Train232412_Xi_BkgParab";
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";
    fileName[0] += "_Run2Binning";
    fileName[1] += "_Run2Binning";
    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.;
    YUpRatio = 2.;
    sleg[0] = "pass 3";
    sleg[1] = "pass 4";
  }
  // v2 with pass4; test 5 vs test 3
  else if (TypeComp == 9)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23zzh_pass4_test3_Train232412_Xi_BkgParab";
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] = "LHC23zzh_pass4_test5_Train235645_Xi_BkgParab";
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";
    fileName[0] += "_Run2Binning";
    fileName[1] += "_Run2Binning";
    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.;
    YUpRatio = 2.;
    sleg[0] = "test 3";
    sleg[1] = "test 5";
  }

  TString hTitleYRatio = "Ratio to " + sleg[0];
  for (Int_t i = 0; i < numOptions; i++)
  {
    Sinputfile = CommonFileName + fileName[i] + ".root";
    cout << "Input file: " << Sinputfile << endl;
    TFile *inputFile = new TFile(Sinputfile);
    if (i == 0)
    {
      hDef = (TH1F *)inputFile->Get(namehisto[i]);
      if (!hDef)
      {
        cout << "Histogram not found" << endl;
        return;
      }
      hDef->SetName("hDefault");
    }
    else
    {
      h[i] = (TH1F *)inputFile->Get(namehisto[i]);
      if (!h[i])
      {
        cout << "Histogram not found" << endl;
        return;
      }
      h[i]->SetName(Form("hVar_%i", i));
      hRatio[i] = (TH1F *)h[i]->Clone(Form("hRatio_%i", i));
      hRatio[i]->Divide(hDef);
      ErrRatioCorr(h[i], hDef, hRatio[i], 0);
    }
  }

  TCanvas *canvas = new TCanvas("canvas", "canvas", 700, 900);
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B
  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvas->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow, YUp, 1, 1, hTitleX, hTitleY, "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(MinPt[!ChosenParticleXi], MaxPt[!ChosenParticleXi]);
  if (TypeComp == 5)
    hDummy->GetXaxis()->SetRangeUser(0, 80);
  pad1->Draw();
  pad1->cd();
  hDummy->Draw("same");

  hDef->SetMarkerColor(ColorMult[5]);
  hDef->SetLineColor(ColorMult[5]);
  hDef->SetMarkerStyle(MarkerMult[0]);
  hDef->SetMarkerSize(0.6 * SizeMult[0]);
  hDef->Draw("same");
  for (Int_t i = 1; i < numOptions; i++)
  {
    h[i]->SetMarkerColor(ColorMult[1]);
    h[i]->SetLineColor(ColorMult[1]);
    h[i]->SetMarkerStyle(MarkerMult[1]);
    h[i]->SetMarkerSize(0.6 * SizeMult[1]);
    h[i]->Draw("same");
  }

  TLegend *leg = new TLegend(0.3, 0.7, 0.9, 0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  if (TypeComp == 5)
    leg->AddEntry("", ParticleName[!ChosenParticleXi] + " Pb-Pb 5.36 TeV", "");
  else
    leg->AddEntry("", ParticleName[!ChosenParticleXi] + Form("  Pb-Pb 5.36 TeV, FT0C %i-%i", CentFT0C[mult], CentFT0C[mult + 1]) + "%", "");
  leg->AddEntry(hDef, sleg[0], "lp");
  for (Int_t i = 1; i < numOptions; i++)
  {
    leg->AddEntry(h[i], sleg[1], "lp");
  }
  leg->Draw("same");

  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummyRatio->GetNbinsX(); i++)
    hDummyRatio->SetBinContent(i, 1e-12);
  SetFont(hDummyRatio);
  StyleHistoYield(hDummyRatio, YLowRatio, YUpRatio, 1, 1, hTitleX, hTitleYRatio, "", 1, 1.15, YoffsetSpectraRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetTickLength(hDummyRatio, tickX, tickY);
  hDummyRatio->GetXaxis()->SetRangeUser(MinPt[!ChosenParticleXi], MaxPt[!ChosenParticleXi]);
  if (TypeComp == 5)
    hDummyRatio->GetXaxis()->SetRangeUser(0, 80);
  canvas->cd();
  padL1->Draw();
  padL1->cd();
  hDummyRatio->Draw("same");
  for (Int_t i = 1; i < numOptions; i++)
  {
    hRatio[i]->SetMarkerColor(ColorMult[1]);
    hRatio[i]->SetLineColor(ColorMult[1]);
    hRatio[i]->SetMarkerStyle(MarkerMult[1]);
    hRatio[i]->SetMarkerSize(0.6 * SizeMult[1]);
    hRatio[i]->Draw("same");
  }
  TF1 *line = new TF1("line", "1", MinPt[!ChosenParticleXi], MaxPt[!ChosenParticleXi]);
  if (TypeComp == 5)
    line = new TF1("line", "1", 0, 80);
  line->SetLineColor(kBlack);
  line->SetLineStyle(9);
  line->Draw("same");

  canvas->SaveAs(Form("Canvas_CompareResults%i.png", TypeComp));
}
