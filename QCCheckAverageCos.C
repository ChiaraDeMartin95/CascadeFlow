// This macro was originally written by:
// chiara.de.martin@cern.ch
#include "TStyle.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPad.h"
// #include "CommonVar_v2.h"
#include "CommonVarPub.h"
// #include "CommonVarOmega.h"
// #include "CommonVarXi.h"
#include "CommonVarLambda.h"

void StyleHisto(TH1D *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title)
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
void StyleHistoYield(TH1D *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset)
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

void SetFont(TH1D *histo)
{
  histo->GetXaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetYaxis()->SetLabelFont(43);
}
void SetTickLength(TH1D *histo, Float_t TickLengthX, Float_t TickLengthY)
{
  histo->GetXaxis()->SetTickLength(TickLengthX);
  histo->GetYaxis()->SetTickLength(TickLengthY);
}
void SetHistoTextSize(TH1D *histo, Float_t XSize, Float_t XLabelSize, Float_t XOffset, Float_t XLabelOffset, Float_t YSize, Float_t YLabelSize, Float_t YOffset, Float_t YLabelOffset)
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

void QCCheckAverageCos(Int_t indexMultTrial = 0,
                       Int_t ChosenPart = ChosenParticle,
                       Bool_t isMC = 0,
                       Bool_t isMassCutForAcceptance = 1, // default is 1, with 0 we have no mass cut so to be able to do the fit of acceptance vs mass; only for LAMBDA
                       Bool_t isRapiditySel = ExtrisRapiditySel,
                       TString inputFileName = SinputFileName,
                       Int_t RebinFactor = 1,
                       Int_t EtaSysChoice = ExtrEtaSysChoice,
                       Bool_t isSysMultTrial = ExtrisSysMultTrial)
{

  if (isMC)
  {
    inputFileName = SinputFileNameMC;
  }

  if (isOOCentrality && commonNumCent != numCentLambdaOO)
  {
    cout << "O-O centrality is selected, but commonNumCent is not set to numCentLambdaOO. Please check the settings." << endl;
    return;
  }
  if (!isOOCentrality && commonNumCent != numCent && commonNumCent != numCentOmega && commonNumCent != numCentOmegaRed && commonNumCent != numCentXiRed)
  {
    cout << "Pb-Pb centrality is selected, but commonNumCent is not set to numCent or numCentOmega. Please check the settings." << endl;
    return;
  }
  if (isReducedPtBins && numPtBins != numPtBinsReduced)
  {
    cout << "Reduced pt bins are selected, but numPtBins is not set to numPtBinsReduced. Please check the settings." << endl;
    return;
  }

  if (isSysMultTrial && !isMC)
    inputFileName = SinputFileNameSyst;
  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trialsBDT)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;
  TString SBDT = "";
  if ((BDTscoreCut != DefaultBDTscoreCut || isSysMultTrial) && ChosenPart != 6 && ChosenPart != 7 && ChosenPart != 8)
    SBDT = Form("_BDT%.3f", BDTscoreCut);

  Int_t Part = 0;
  if ((ChosenPart == 1) || (ChosenPart == 4) || (ChosenPart == 5))
  {
    Part = 1; // Omega
  }

  TString SinputFile = "../OutputAnalysis/Output" + STHN[ExtrisFromTHN] + "_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice]; // + SBDT;
  if (isApplyWeights)
    SinputFile += "_Weighted";
  if (isApplyCentWeight)
    SinputFile += "_CentWeighted";
  if (v2type == 1)
    SinputFile += "_SP";
  if (!useCommonBDTValue)
    SinputFile += "_BDTCentDep";
  if (isRun2Binning)
    SinputFile += "_Run2Binning";
  if (!isRapiditySel)
    SinputFile += "_Eta08";
  SinputFile += SBDT;
  if (ChosenPart >= 6 && ExtrisSysLambdaMultTrial)
  {
    if (isLoosest)
      SinputFile += "_isLoosest";
    else if (isTightest)
      SinputFile += "_isTightest";
    // SinputFile += Form("_SysMultTrial_%i", indexMultTrial);
  }
  if (isOOCentrality)
    SinputFile += "_isOOCentrality";
  if (ExtrisApplyResoOnTheFly)
    SinputFile += "_ResoOnTheFly";
  if (ExtrisApplyEffWeights && ChosenPart >= 6)
    SinputFile += "_EffWeighted";
  if (ChosenPart >= 6 && !ExtrisFromTHN)
  {
    SinputFile += "_Nvar1";
    if (SinputFileName == "LHC25_OO_pass2_Train598890_MyEff")
      SinputFile += "_050PtCut";
  }
  // SinputFile += "_Nvar1_TestMoreBins";
  if (ChosenPart >= 6 && !isMassCutForAcceptance)
    SinputFile += "_NoMassCutForAcceptance";
  if (ExtrisCentOmegaRed && Part == 1)
    SinputFile += "_OmegaRedCent";
  if (ExtrisCentXiRed && Part == 0)
    SinputFile += "_XiRedCent";
  SinputFile += "_NewHistos.root";
  cout << "Input file: " << SinputFile << endl;

  TFile *inputFile = new TFile(SinputFile);

  TH1D *hCos2Phi[commonNumCent + 1];
  TH1D *hCos2Psi[commonNumCent + 1];
  TH1D *hSin2Phi[commonNumCent + 1];
  TH1D *hSin2Psi[commonNumCent + 1];
  TH1D *hCosTheta[commonNumCent + 1];

  int nBinsCent;
  double *binsCent;
  if (isOOCentrality)
  {
    nBinsCent = numCentLambdaOO;
    binsCent = fCentFT0CLambdaOO;
  }
  else if (Part == 1)
  { // Omega in Pb-Pb
    nBinsCent = numCentOmega;
    binsCent = fCentFT0COmega;
    if (ExtrisCentOmegaRed)
    {
      nBinsCent = numCentOmegaRed;
      binsCent = fCentFT0COmegaRed;
    }
  }
  else
  {
    nBinsCent = numCent;
    binsCent = fCentFT0C;
  }

  TH1D *hCos2PhiVsCent = new TH1D("hCos2PhiVsCent", "Average cos(2#phi) vs Centrality", nBinsCent, binsCent);
  TH1D *hCos2PsiVsCent = new TH1D("hCos2PsiVsCent", "Average cos(2#psi) vs Centrality", nBinsCent, binsCent);
  TH1D *hSin2PhiVsCent = new TH1D("hSin2PhiVsCent", "Average sin(2#phi) vs Centrality", nBinsCent, binsCent);
  TH1D *hSin2PsiVsCent = new TH1D("hSin2PsiVsCent", "Average sin(2#psi) vs Centrality", nBinsCent, binsCent);
  TH1D *hCosThetaVsCent = new TH1D("hCosThetaVsCent", "Average cos(#theta) vs Centrality", nBinsCent, binsCent);

  TProfile *pCos2PsiVsCent = (TProfile*)inputFile->Get("profileCos2Psi");
  if (!pCos2PsiVsCent)
    return;
  TProfile *pSin2PsiVsCent = (TProfile*)inputFile->Get("profileSin2Psi");
  if (!pSin2PsiVsCent)
    return;
  TProfile *pSin2PsiDiff = (TProfile*)inputFile->Get("profileSin2PsiDiff");
  if (!pSin2PsiDiff)
    return;
  TProfile *pCosTheta = (TProfile*)inputFile->Get("profileCosTheta");
  if (!pCosTheta)
    return;

  if (!pCos2PsiVsCent)
    return;

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;

  for (Int_t cent = 0; cent < commonNumCent + 1; cent++)
  {
    if (cent == numCent)
    { // 0-80%
      CentFT0CMin = 0;
      CentFT0CMax = CentFT0CMaxPbPb;
    }
    else
    {
      CentFT0CMin = CentFT0C[cent];
      CentFT0CMax = CentFT0C[cent + 1];
    }
    if (ExtrisCentXiRed)
    {
      if (cent == numCentXiRed)
      {
        CentFT0CMin = 0;
        CentFT0CMax = CentFT0CMaxXiRed;
      }
      else
      {
        CentFT0CMin = CentFT0CXiRed[cent];
        CentFT0CMax = CentFT0CXiRed[cent + 1];
      }
    }
    if (Part == 1) // Omega in Pb-Pb
    {
      if (cent == numCentOmega)
      {
        CentFT0CMin = 0;
        CentFT0CMax = CentFT0CMaxOmega;
      }
      else
      {
        CentFT0CMin = CentFT0COmega[cent];
        CentFT0CMax = CentFT0COmega[cent + 1];
      }
      if (ExtrisCentOmegaRed)
      {
        if (cent == numCentOmegaRed)
        {
          CentFT0CMin = 0;
          CentFT0CMax = CentFT0CMaxOmegaRed;
        }
        else
        {
          CentFT0CMin = CentFT0COmegaRed[cent];
          CentFT0CMax = CentFT0COmegaRed[cent + 1];
        }
      }
    }
    if (isOOCentrality)
    {
      if (cent == numCentLambdaOO)
      {
        CentFT0CMin = 0;
        CentFT0CMax = CentFT0CMaxLambdaOO;
      }
      else
      {
        CentFT0CMin = CentFT0CLambdaOO[cent];
        CentFT0CMax = CentFT0CLambdaOO[cent + 1];
      }
    }
    hCos2Phi[cent] = (TH1D *)inputFile->Get(Form("cos2Phi_cent%i-%i", CentFT0CMin, CentFT0CMax));
    hCos2Psi[cent] = (TH1D *)inputFile->Get(Form("cos2Psi_cent%i-%i", CentFT0CMin, CentFT0CMax));
    hSin2Phi[cent] = (TH1D *)inputFile->Get(Form("sin2Phi_cent%i-%i", CentFT0CMin, CentFT0CMax));
    hSin2Psi[cent] = (TH1D *)inputFile->Get(Form("sin2Psi_cent%i-%i", CentFT0CMin, CentFT0CMax));
    hCosTheta[cent] = (TH1D *)inputFile->Get(Form("cosTheta_cent%i-%i", CentFT0CMin, CentFT0CMax));

    //cout << "Mean cos(2#phi) in cent: " << CentFT0CMin << "-" << CentFT0CMax << " = " << hCos2Psi[cent]->GetMean() << "+-" << hCos2Psi[cent]->GetMeanError() << endl;
    //hSin2Psi[cent]->GetXaxis()->SetRangeUser(-0.8, 0.8);
    //hCos2Psi[cent]->GetXaxis()->SetRangeUser(-0.8, 0.8);
    //cout << "Mean cos(2#phi) in cent: " << CentFT0CMin << "-" << CentFT0CMax << " = " << hCos2Psi[cent]->GetMean() << "+-" << hCos2Psi[cent]->GetMeanError() << endl;
    hCos2PhiVsCent->SetBinContent(cent + 1, hCos2Phi[cent]->GetMean());
    hCos2PsiVsCent->SetBinContent(cent + 1, hCos2Psi[cent]->GetMean());
    hSin2PhiVsCent->SetBinContent(cent + 1, hSin2Phi[cent]->GetMean());
    hSin2PsiVsCent->SetBinContent(cent + 1, hSin2Psi[cent]->GetMean());
    hCosThetaVsCent->SetBinContent(cent + 1, hCosTheta[cent]->GetMean());

    hCos2PhiVsCent->SetBinError(cent + 1, hCos2Phi[cent]->GetMeanError());
    hCos2PsiVsCent->SetBinError(cent + 1, hCos2Psi[cent]->GetMeanError());
    hSin2PhiVsCent->SetBinError(cent + 1, hSin2Phi[cent]->GetMeanError());
    hSin2PsiVsCent->SetBinError(cent + 1, hSin2Psi[cent]->GetMeanError());
    hCosThetaVsCent->SetBinError(cent + 1, hCosTheta[cent]->GetMeanError());
  }
  gStyle->SetOptStat(0);
  TH1D *hDummy = new TH1D("hDummy", "hDummy", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, -1000);
  SetFont(hDummy);
  StyleHistoYield(hDummy, -0.01, 0.01, 1, 1, TitleXCent, "", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);

  TF1 *lineAtzero = new TF1("lineAtzero", "0", 0, 100);
  lineAtzero->SetLineColor(kBlack);
  lineAtzero->SetLineWidth(2);
  lineAtzero->SetLineStyle(2);

  TCanvas *canvasCos2Phi = new TCanvas("canvasCos2Phi", "Cos2Phi vs Centrality", 800, 600);
  StyleCanvas(canvasCos2Phi, 0.06, 0.15, 0.15, 0.03);
  hDummy->Draw("");
  StyleHistoYield(hCos2PhiVsCent, 0, 1, kBlue + 1, 20, "Centrality (%)", "Cos(2#varphi)", "", 1.5, 1.5, 1.7);
  hCos2PhiVsCent->Draw("same");
  lineAtzero->Draw("same");

  TCanvas *canvasCos2Psi = new TCanvas("canvasCos2Psi", "Cos2Psi vs Centrality", 800, 600);
  StyleCanvas(canvasCos2Psi, 0.06, 0.15, 0.15, 0.03);
  hDummy->Draw("");
  StyleHistoYield(hCos2PsiVsCent, 0, 1, kRed + 1, 20, "Centrality (%)", "Cos(2#Psi)", "", 1.5, 1.5, 1.7);
  pCos2PsiVsCent->SetLineColor(kRed + 1);
  pCos2PsiVsCent->SetLineWidth(2);
  pCos2PsiVsCent->SetLineStyle(2);
  hCos2PsiVsCent->Draw("same");
  pCos2PsiVsCent->Draw("same");
  lineAtzero->Draw("same");

  TCanvas *canvasSin2Phi = new TCanvas("canvasSin2Phi", "Sin2Phi vs Centrality", 800, 600);
  StyleCanvas(canvasSin2Phi, 0.06, 0.15, 0.15, 0.03);
  hDummy->Draw("");
  StyleHistoYield(hSin2PhiVsCent, 0, 1, kGreen + 1, 20, "Centrality (%)", "Sin(2#varphi)", "", 1.5, 1.5, 1.7);
  hSin2PhiVsCent->Draw("same");
  lineAtzero->Draw("same");

  TCanvas *canvasSin2Psi = new TCanvas("canvasSin2Psi", "Sin2Psi vs Centrality", 800, 600);
  StyleCanvas(canvasSin2Psi, 0.06, 0.15, 0.15, 0.03);
  hDummy->GetYaxis()->SetRangeUser(-0.001, 0.001);
  hDummy->Draw("");
  StyleHistoYield(hSin2PsiVsCent, 0, 1, kMagenta + 1, 20, "Centrality (%)", "Sin(2#Psi)", "", 1.5, 1.5, 1.7);
  hSin2PsiVsCent->Draw("same");
  pSin2PsiVsCent->SetLineColor(kMagenta + 1);
  pSin2PsiVsCent->SetLineWidth(2);
  pSin2PsiVsCent->SetLineStyle(2);
  pSin2PsiVsCent->Draw("same");
  lineAtzero->Draw("same");

  TCanvas *canvasSin2PsiDiff = new TCanvas("canvasSin2PsiDiff", "Sin2Psi Difference vs Centrality", 800, 600);
  StyleCanvas(canvasSin2PsiDiff, 0.06, 0.15, 0.15, 0.03);
  hDummy->Draw("");
  StyleHistoYield(pSin2PsiDiff, 0, 1, kCyan + 1, 20, "Centrality (%)", "Sin(2#phi-2#psi)", "", 1.5, 1.5, 1.7);
  pSin2PsiDiff->Draw("same");
  lineAtzero->Draw("same");

  TCanvas *canvasCosTheta = new TCanvas("canvasCosTheta", "CosTheta vs Centrality", 800, 600);
  StyleCanvas(canvasCosTheta, 0.06, 0.15, 0.15, 0.03);
  hDummy->Draw("");
  StyleHistoYield(hCosThetaVsCent, 0, 1, kOrange + 1, 20, "Centrality (%)", "Cos(#theta)", "", 1.5, 1.5, 1.7);
  hCosThetaVsCent->Draw("same");
  StyleHistoYield(pCosTheta, 0, 1, kOrange + 1, 20, "Centrality (%)", "Cos(#theta)", "", 1.5, 1.5, 1.7);
  pCosTheta->SetLineColor(kOrange + 1);
  pCosTheta->SetLineWidth(2);
  pCosTheta->SetLineStyle(2);
  pCosTheta->Draw("same");
  lineAtzero->Draw("same");

  TH1D *hCosThetaTimesSinPsiDiff = new TH1D("hCosThetaTimesSinPsiDiff", "CosTheta * Sin2Psi Difference vs Centrality", 10, 0, 100);
  for (int cent = 0; cent < commonNumCent; cent++) {
    hCosThetaTimesSinPsiDiff->SetBinContent(cent + 1, hCosThetaVsCent->GetBinContent(cent + 1) * pSin2PsiDiff->GetBinContent(cent + 1));
    hCosThetaTimesSinPsiDiff->SetBinError(cent + 1, sqrt(pow(hCosThetaVsCent->GetBinError(cent + 1)/hCosThetaVsCent->GetBinContent(cent + 1), 2) + pow(pSin2PsiDiff->GetBinError(cent + 1)/pSin2PsiDiff->GetBinContent(cent + 1), 2)) * hCosThetaTimesSinPsiDiff->GetBinContent(cent + 1));
  }
  StyleHistoYield(hCosThetaTimesSinPsiDiff, -0.0000003, 0.0000003, kBlue + 1, 20, "Centrality (%)", "#LT cos(#theta*_{p}) sin(2#varphi - 2#Psi) #GT", "", 1.5, 1.5, 1.7);
  hCosThetaTimesSinPsiDiff->GetYaxis()->SetTitleOffset(1.2);
  TCanvas *canvasCosThetaTimesSinPsiDiff = new TCanvas("canvasCosThetaTimesSinPsiDiff", "CosTheta * Sin2Psi Difference vs Centrality", 800, 600);
  StyleCanvas(canvasCosThetaTimesSinPsiDiff, 0.06, 0.15, 0.15, 0.03);
  hCosThetaTimesSinPsiDiff->Draw("same");
  lineAtzero->Draw("same");
  canvasCosThetaTimesSinPsiDiff->SaveAs("../LowerOrderCumulants.png");
}
