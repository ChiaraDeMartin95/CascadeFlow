// This macro was originally written by:
// chiara.de.martin@cern.ch
#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "THn.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TString.h"
#include "TPad.h"
#include "StyleFile.h"
#include "CommonVar.h"
#include "TRandom3.h"
#include <ROOT/RDataFrame.hxx>

using namespace ROOT;
using namespace std;

void ProcessTHN(Int_t indexMultTrial = 0,
                Int_t ChosenPart = ChosenParticle,
                TString inputFileName = SinputFileName,
                Int_t EtaSysChoice = ExtrEtaSysChoice,
                Bool_t isSysMultTrial = ExtrisSysMultTrial)
{

  // phi and efficiency weights should be applied in task
  Int_t Part = 0;
  if ((ChosenPart == 1) || (ChosenPart == 4) || (ChosenPart == 5))
  {
    Part = 1; // Omega
  }
  else if (ChosenPart == 6)
  {
    Part = 6; // Lambda
  }

  if (isSysMultTrial)
  {
    inputFileName = SinputFileNameSyst;
  }

  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trialsBDT)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;

  // inputFileName = "TestTHN";
  TString SinputFile = "TreeForAnalysis";
  SinputFile += "/AnalysisResults_" + inputFileName + ".root";

  cout << "BDT score cut: " << BDTscoreCut << endl;
  cout << "ChosenPart: " << ParticleName[ChosenPart] << endl;
  cout << "EtaSysChoice: " << EtaSysChoice << endl;
  cout << "Use common BDT value " << useCommonBDTValue << endl;

  // get THNs
  TFile *inputFile = new TFile(SinputFile, "READ");
  if (!inputFile)
  {
    cout << "File not found" << endl;
    return;
  }
  TDirectoryFile *dir = (TDirectoryFile *)inputFile->Get("lf-cascade-flow");
  TDirectoryFile *dir1 = (TDirectoryFile *)dir->Get("histos");
  THn *hV2 = (THnF *)dir1->Get("h" + ParticleName[Part] + "V2");
  // axes: thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassCasc, thnAxisBDTScore, thnAxisV2
  if (!hV2)
  {
    cout << "THn hV2 not found" << endl;
    // return;
  }
  THn *hXiPzs2 = (THnF *)dir1->Get("h" + ParticleName[Part] + "Pzs2");
  // axes: thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassCasc, thnAxisBDTScore, thnAxisPzs2Xi
  if (!hXiPzs2)
  {
    cout << "THn hXiPzs2 not found" << endl;
    // return;
  }
  THn *hXiPzs2FromLambda = (THnF *)dir1->Get("h" + ParticleName[Part] + "Pzs2FromLambda");
  // axes: thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassCasc, thnAxisBDTScore, thnAxisPzs2Lambda
  if (!hXiPzs2FromLambda)
  {
    cout << "THn hXiPzs2FromLambda not found" << endl;
    // return;
  }
  THn *hXiCos2Theta = (THnF *)dir1->Get("h" + ParticleName[Part] + "Cos2Theta");
  // axes: thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassCasc, thnAxisBDTScore, thnAxisCos2Theta
  if (!hXiCos2Theta)
  {
    cout << "THn hXiCos2Theta not found" << endl;
    // return;
  }

  THn *hXiCos2ThetaFromLambda = (THnF *)dir1->Get("h" + ParticleName[Part] + "Cos2ThetaFromLambda");
  // axes: thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassCasc, thnAxisBDTScore, thnAxisCos2ThetaLambda
  // new axes: thnAxisFT0C, thnAxisCharge, Pt Xi, Pt Xi, Mass Xi, thnAxisBDTScore, thnAxisCos2ThetaLambda
  if (!hXiCos2ThetaFromLambda)
  {
    cout << "THn hXiCos2ThetaFromLambda not found" << endl;
    // return;
  }
  THn *hXiCos2ThetaFromLambdaL = (THnF *)dir1->Get("h" + ParticleName[Part] + "Cos2ThetaFromLambdaL");
  // axes: thnAxisFT0C, thnAxisCharge, eta Lambda, pt Lambda, Mass Lambda, thnAxisBDTScore, thnAxisCos2ThetaLambda
  if (!hXiCos2ThetaFromLambdaL)
  {
    cout << "THn hXiCos2ThetaFromLambdaL not found" << endl;
    // return;
  }

  if (!hV2)
  {
    // hV2 = (THnF *)hXiPzs2->Clone("hV2");
    if (isProducedAcceptancePlots)
    {
      if (ChosenPart == 6)
        hV2 = (THnF *)hXiCos2Theta->Clone("hV2");
      else
        hV2 = (THnF *)hXiCos2ThetaFromLambdaL->Clone("hV2");
    }
    else
    {
      if (ChosenPart == 6)
        hV2 = (THnF *)hXiPzs2->Clone("hV2");
      else
        hV2 = (THnF *)hXiPzs2FromLambda->Clone("hV2");
    }
  }
  if (!hXiPzs2)
  {
    hXiPzs2 = (THnF *)hV2->Clone("hXiPzs2");
    if (!hXiPzs2)
    {
      cout << "!!!Not found" << endl;
      return;
    }
  }
  if (!hXiPzs2FromLambda)
  {
    hXiPzs2FromLambda = (THnF *)hV2->Clone("hXiPzs2FromLambda");
    if (!hXiPzs2FromLambda)
      return;
  }
  if (!hXiCos2Theta)
  {
    hXiCos2Theta = (THnF *)hV2->Clone("hXiCos2Theta");
    if (!hXiCos2Theta)
      return;
  }
  if (!hXiCos2ThetaFromLambda)
  {
    hXiCos2ThetaFromLambda = (THnF *)hV2->Clone("hXiCos2ThetaFromLambda");
    if (!hXiCos2ThetaFromLambda)
      return;
  }
  if (!hXiCos2ThetaFromLambdaL)
  {
    hXiCos2ThetaFromLambdaL = (THnF *)hXiCos2ThetaFromLambda->Clone("hXiCos2ThetaFromLambdaL");
    if (!hXiCos2ThetaFromLambdaL)
      return;
  }

  // V2, Pzs2, Cos2Theta ***
  TH3D *hmassVsPtVsV2C[numCent + 1];
  TH3D *hmassVsPtVsPzs2[numCent + 1];
  TH3D *hmassVsPtVsPzs2LambdaFromC[numCent + 1];
  TH3D *hmassVsPsiVsPz[numCent + 1];
  TH3D *hmassVsPsiVsPzLambdaFromC[numCent + 1];
  TString hName[numCent + 1] = {""};
  TString hNamePzs2_3D[numCent + 1] = {""};
  TString hNamePzs2LambdaFromC_3D[numCent + 1] = {""};
  TString hNamePzVsPsi_3D[numCent + 1] = {""};
  TString hNamePzVsPsiLambdaFromC_3D[numCent + 1] = {""};

  // acceptance correction ***
  TH3D *hmassVsPtVsCos2Theta[numCent + 1];
  TH3D *hmassVsPtVsCos2ThetaLambdaFromC[numCent + 1];
  TH3D *hmassVsPsiVsCos2Theta[numCent + 1];
  TH3D *hmassVsPsiVsCos2ThetaLambdaFromC[numCent + 1];
  TH3D *hEtaVsPtVsCos2ThetaLambdaFromC[numCent + 1];
  TString hNameCos2Theta_3D[numCent + 1] = {""};
  TString hNameCos2ThetaLambdaFromC_3D[numCent + 1] = {""};
  TString hNameCos2ThetaLambdaFromC_Eta3D[numCent + 1] = {""};
  TString hNameCos2ThetaVsPsi_3D[numCent + 1] = {""};
  TString hNameCos2ThetaVsPsiLambdaFromC_3D[numCent + 1] = {""};

  TH2F *hPhiCentHisto[numCent];

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;

  // QCPlots
  TCanvas *canvasQC = new TCanvas("canvasQC", "canvasQC", 1400, 800);
  canvasQC->Divide(3, 2);

  TH1F *hCentrality[numCent + 1];
  TH1F *hCharge[numCent + 1];
  TH1F *hPt[numCent + 1];
  TH1F *hMass[numCent + 1];
  TH1F *hBDT[numCent + 1];
  TH1F *hMassLambda[numCent + 1];

  TH1F *hDummyCentrality = (TH1F *)hV2->Projection(0);
  TH1F *hDummyCharge = (TH1F *)hV2->Projection(1);
  TH1F *hDummyPt = (TH1F *)hV2->Projection(3);
  TH1F *hDummyMass = (TH1F *)hV2->Projection(2);
  TH1F *hDummyBDT = (TH1F *)hV2->Projection(4);
  TH1F *hDummyMassLambda = (TH1F *)hXiCos2ThetaFromLambdaL->Projection(4);
  // TH1F *hDummyPt = (TH1F *)hXiCos2ThetaFromLambdaL->Projection(3);
  // TH1F *hDummyMass = (TH1F *)hXiCos2ThetaFromLambda->Projection(4);
  // TH1F *hDummyBDT = (TH1F *)hV2->Projection(5);
  // TH1F *hDummyMassLambda = (TH1F *)hXiCos2ThetaFromLambdaL->Projection(4);
  hDummyCentrality->Reset();
  hDummyCharge->Reset();
  hDummyPt->Reset();
  hDummyMass->Reset();
  hDummyBDT->Reset();
  hDummyMassLambda->Reset();
  hDummyCentrality->GetYaxis()->SetRangeUser(0, 1.2);
  hDummyCharge->GetYaxis()->SetRangeUser(0, 1);
  hDummyPt->GetYaxis()->SetRangeUser(0, 0.25);
  hDummyMass->GetYaxis()->SetRangeUser(0, 0.15);
  hDummyBDT->GetYaxis()->SetRangeUser(0, 1.2);
  hDummyMassLambda->GetYaxis()->SetRangeUser(0, 0.1);

  TH1F *hBDTSelection = new TH1F("hBDTSelection", "BDT selection", 8, 0, 80);

  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    if (cent == numCent)
    { // 0-80%
      CentFT0CMin = 0;
      CentFT0CMax = 80;
    }
    else
    {
      CentFT0CMin = CentFT0C[cent];
      CentFT0CMax = CentFT0C[cent + 1];
    }
    cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << endl;
    // Selection of centrality range
    hV2->GetAxis(0)->SetRange(hV2->GetAxis(0)->FindBin(CentFT0CMin + 0.001), hV2->GetAxis(0)->FindBin(CentFT0CMax - 0.001));
    hXiPzs2->GetAxis(0)->SetRange(hXiPzs2->GetAxis(0)->FindBin(CentFT0CMin + 0.001), hXiPzs2->GetAxis(0)->FindBin(CentFT0CMax - 0.001));
    hXiCos2Theta->GetAxis(0)->SetRange(hXiCos2Theta->GetAxis(0)->FindBin(CentFT0CMin + 0.001), hXiCos2Theta->GetAxis(0)->FindBin(CentFT0CMax - 0.001));
    hXiPzs2FromLambda->GetAxis(0)->SetRange(hXiPzs2FromLambda->GetAxis(0)->FindBin(CentFT0CMin + 0.001), hXiPzs2FromLambda->GetAxis(0)->FindBin(CentFT0CMax - 0.001));
    hXiCos2ThetaFromLambda->GetAxis(0)->SetRange(hXiCos2ThetaFromLambda->GetAxis(0)->FindBin(CentFT0CMin + 0.001), hXiCos2ThetaFromLambda->GetAxis(0)->FindBin(CentFT0CMax - 0.001));
    hXiCos2ThetaFromLambdaL->GetAxis(0)->SetRange(hXiCos2ThetaFromLambda->GetAxis(0)->FindBin(CentFT0CMin + 0.001), hXiCos2ThetaFromLambda->GetAxis(0)->FindBin(CentFT0CMax - 0.001));

    // Selection of charge
    if (ChosenParticle == 2 || ChosenParticle == 4)
    {
      hV2->GetAxis(1)->SetRange(1, 1);          // Charge < 0
      hXiPzs2->GetAxis(1)->SetRange(1, 1);      // Charge < 0
      hXiCos2Theta->GetAxis(1)->SetRange(1, 1); // Charge < 0
      hXiPzs2FromLambda->GetAxis(1)->SetRange(1, 1);
      hXiCos2ThetaFromLambda->GetAxis(1)->SetRange(1, 1);
      hXiCos2ThetaFromLambdaL->GetAxis(1)->SetRange(1, 1);
    }
    else if (ChosenParticle == 3 || ChosenParticle == 5)
    {
      hV2->GetAxis(1)->SetRange(2, 2);          // Charge > 0
      hXiPzs2->GetAxis(1)->SetRange(2, 2);      // Charge > 0
      hXiCos2Theta->GetAxis(1)->SetRange(2, 2); // Charge > 0
      hXiPzs2FromLambda->GetAxis(1)->SetRange(2, 2);
      hXiCos2ThetaFromLambda->GetAxis(1)->SetRange(2, 2);
      hXiCos2ThetaFromLambdaL->GetAxis(1)->SetRange(2, 2);
    }
    else
    {
      hV2->GetAxis(1)->SetRange(1, 2);          // All charges
      hXiPzs2->GetAxis(1)->SetRange(1, 2);      // All charges
      hXiCos2Theta->GetAxis(1)->SetRange(1, 2); // All charges
      hXiPzs2FromLambda->GetAxis(1)->SetRange(1, 2);
      hXiCos2ThetaFromLambda->GetAxis(1)->SetRange(1, 2);
      hXiCos2ThetaFromLambdaL->GetAxis(1)->SetRange(1, 2);
    }

    // Selection on BDT score -- be careful, in the 0-80% case, the BDT score is not defined and I apply the same cut as in the 0-10% case
    if (!useCommonBDTValue && !isSysMultTrial)
      BDTscoreCut = bdtCut[cent];
    if (ChosenPart != 6)
      cout << "BDTscoreCut: " << BDTscoreCut << endl;
    if (ChosenPart != 6)
    { // No BDT axis for Lambda
      hV2->GetAxis(4)->SetRange(hV2->GetAxis(4)->FindBin(BDTscoreCut + 0.001), hV2->GetAxis(4)->FindBin(1 - 0.001));
      hXiPzs2->GetAxis(4)->SetRange(hXiPzs2->GetAxis(4)->FindBin(BDTscoreCut + 0.001), hXiPzs2->GetAxis(4)->FindBin(1 - 0.001));
      hXiCos2Theta->GetAxis(4)->SetRange(hXiCos2Theta->GetAxis(4)->FindBin(BDTscoreCut + 0.001), hXiCos2Theta->GetAxis(4)->FindBin(1 - 0.001));
      hXiPzs2FromLambda->GetAxis(4)->SetRange(hXiPzs2FromLambda->GetAxis(4)->FindBin(BDTscoreCut + 0.001), hXiPzs2FromLambda->GetAxis(4)->FindBin(1 - 0.001));

      hXiCos2ThetaFromLambda->GetAxis(5)->SetRange(hXiCos2ThetaFromLambda->GetAxis(5)->FindBin(BDTscoreCut + 0.001), hXiCos2ThetaFromLambda->GetAxis(5)->FindBin(1 - 0.001));
      hXiCos2ThetaFromLambdaL->GetAxis(5)->SetRange(hXiCos2ThetaFromLambdaL->GetAxis(5)->FindBin(BDTscoreCut + 0.001), hXiCos2ThetaFromLambdaL->GetAxis(5)->FindBin(1 - 0.001));
    }

    // LambdaMass selection
    hXiCos2ThetaFromLambdaL->GetAxis(4)->SetRange(hXiCos2ThetaFromLambdaL->GetAxis(4)->FindBin(1.112 + 0.00001), hXiCos2ThetaFromLambdaL->GetAxis(4)->FindBin(1.119 - 0.00001));

    // Lambda Eta selection
    hXiCos2ThetaFromLambdaL->GetAxis(2)->SetRange(hXiCos2ThetaFromLambdaL->GetAxis(2)->FindBin(-0.8 + 0.00001), hXiCos2ThetaFromLambdaL->GetAxis(2)->FindBin(0.8 - 0.00001));

    if (ChosenPart == 6)
    {
      hXiCos2Theta->GetAxis(4)->SetRange(hXiCos2Theta->GetAxis(4)->FindBin(1.112 + 0.00001), hXiCos2Theta->GetAxis(4)->FindBin(1.119 - 0.00001));
      hXiCos2Theta->GetAxis(2)->SetRange(hXiCos2Theta->GetAxis(2)->FindBin(-0.8 + 0.00001), hXiCos2Theta->GetAxis(2)->FindBin(0.8 - 0.00001));
    }

    // QCPlots
    hBDTSelection->SetBinContent(cent + 1, BDTscoreCut);
    hBDTSelection->SetBinError(cent + 1, 0);

    cout << "QC plots" << endl;
    gStyle->SetOptStat(0);
    canvasQC->cd(1);
    hCentrality[cent] = (TH1F *)hV2->Projection(0);
    hCentrality[cent]->Scale(1. / hCentrality[cent]->Integral());
    hCentrality[cent]->SetLineColor(ColorMult[cent]);
    hCentrality[cent]->SetMarkerColor(ColorMult[cent]);
    hCentrality[cent]->GetXaxis()->SetRangeUser(0, 80);
    if (cent == 0)
      hDummyCentrality->Draw();
    hCentrality[cent]->Draw("same");

    canvasQC->cd(2);
    hCharge[cent] = (TH1F *)hV2->Projection(1);
    hCharge[cent]->Scale(1. / hCharge[cent]->Integral());
    hCharge[cent]->SetLineColor(ColorMult[cent]);
    hCharge[cent]->SetMarkerColor(ColorMult[cent]);
    hCharge[cent]->GetXaxis()->SetRangeUser(0, 2);
    if (cent == 0)
      hDummyCharge->Draw();
    hCharge[cent]->Draw("same");

    canvasQC->cd(3);
    // hPt[cent] = (TH1F *)hXiPzs2FromLambda->Projection(2);
    hPt[cent] = (TH1F *)hXiCos2Theta->Projection(3); // pt
    hPt[cent]->Scale(1. / hPt[cent]->Integral());
    hPt[cent]->SetLineColor(ColorMult[cent]);
    hPt[cent]->SetMarkerColor(ColorMult[cent]);
    hPt[cent]->GetXaxis()->SetRangeUser(0, 10);
    if (cent == 0)
      hDummyPt->Draw();
    hPt[cent]->Draw("same");

    canvasQC->cd(4);
    // hMass[cent] = (TH1F *)hXiPzs2FromLambda->Projection(3);
    hMass[cent] = (TH1F *)hXiCos2Theta->Projection(2);
    hMass[cent]->Scale(1. / hMass[cent]->Integral());
    hMass[cent]->SetLineColor(ColorMult[cent]);
    hMass[cent]->SetMarkerColor(ColorMult[cent]);
    hMass[cent]->GetXaxis()->SetRangeUser(1.3, 1.345);
    if (cent == 0)
      hDummyMass->Draw();
    hMass[cent]->Draw("same");

    canvasQC->cd(5);
    hBDT[cent] = (TH1F *)hXiPzs2FromLambda->Projection(4);
    hBDT[cent]->Scale(1. / hBDT[cent]->Integral());
    hBDT[cent]->SetLineColor(ColorMult[cent]);
    hBDT[cent]->SetMarkerColor(ColorMult[cent]);
    hBDT[cent]->GetXaxis()->SetRangeUser(0, 1);
    hBDT[cent]->GetYaxis()->SetRangeUser(0, 0.4);
    if (cent == 0)
      hDummyBDT->Draw();
    hBDT[cent]->Draw("same");

    canvasQC->cd(6);
    // hMassLambda[cent] = (TH1F *)hXiCos2ThetaFromLambdaL->Projection(4);
    hMassLambda[cent] = (TH1F *)hXiPzs2FromLambda->Projection(3);
    hMassLambda[cent]->Scale(1. / hMassLambda[cent]->Integral());
    hMassLambda[cent]->SetLineColor(ColorMult[cent]);
    hMassLambda[cent]->SetMarkerColor(ColorMult[cent]);
    hMassLambda[cent]->GetXaxis()->SetRangeUser(1.1, 1.13);
    if (cent == 0)
      hDummyMassLambda->Draw();
    hMassLambda[cent]->Draw("same");

    // Histo name definition
    hName[cent] = Form("massVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax);
    if (ExtrisApplyEffWeights)
      hName[cent] = Form("massVsPtVsV2CWeighted_cent%i-%i", CentFT0CMin, CentFT0CMax);

    hNamePzs2_3D[cent] = Form("massVsPtVsPzs2_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNamePzs2LambdaFromC_3D[cent] = Form("massVsPtVsPzs2LambdaFromC_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNamePzVsPsi_3D[cent] = Form("massVsPsiVsPz_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNamePzVsPsiLambdaFromC_3D[cent] = Form("massVsPsiVsPzLambdaFromC_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);

    hNameCos2Theta_3D[cent] = Form("massVsPtVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNameCos2ThetaLambdaFromC_3D[cent] = Form("massVsPtVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNameCos2ThetaLambdaFromC_Eta3D[cent] = Form("etaVsPtVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNameCos2ThetaVsPsi_3D[cent] = Form("massVsPsiVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNameCos2ThetaVsPsiLambdaFromC_3D[cent] = Form("massVsPsiVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax);

    // Projections
    Int_t AxisNumber = 5;
    if (ChosenPart == 6) // Lambda
      AxisNumber = 4;    // no BDT axis for Lambda

    hmassVsPtVsV2C[cent] = (TH3D *)hV2->Projection(3, 2, AxisNumber); // mass, pt, v2
    hmassVsPtVsV2C[cent]->SetName(hName[cent]);
    hmassVsPtVsPzs2[cent] = (TH3D *)hXiPzs2->Projection(3, 2, AxisNumber); // mass, pt, pzs2
    hmassVsPtVsPzs2[cent]->SetName(hNamePzs2_3D[cent]);
    hmassVsPtVsPzs2LambdaFromC[cent] = (TH3D *)hXiPzs2FromLambda->Projection(3, 2, AxisNumber); // mass, pt, pzs2LambdaFromC
    hmassVsPtVsPzs2LambdaFromC[cent]->SetName(hNamePzs2LambdaFromC_3D[cent]);

    // BE CAREFUL: AXES TO BE UPDATED
    hmassVsPsiVsPz[cent] = (TH3D *)hXiPzs2->Projection(3, 2, AxisNumber); // mass, psi, cosThetaLambda
    hmassVsPsiVsPz[cent]->SetName(hNamePzVsPsi_3D[cent]);
    hmassVsPsiVsPzLambdaFromC[cent] = (TH3D *)hXiPzs2->Projection(3, 2, AxisNumber); // mass, psi, cosThetaProton
    hmassVsPsiVsPzLambdaFromC[cent]->SetName(hNamePzVsPsiLambdaFromC_3D[cent]);

    hmassVsPtVsCos2Theta[cent] = (TH3D *)hXiCos2Theta->Projection(3, 2, AxisNumber); // mass, pt, cos2ThetaLambda
    hmassVsPtVsCos2Theta[cent]->SetName(hNameCos2Theta_3D[cent]);
    hmassVsPtVsCos2ThetaLambdaFromC[cent] = (TH3D *)hXiCos2ThetaFromLambda->Projection(3, 2, AxisNumber); // mass, pt, cos2ThetaProton
    hmassVsPtVsCos2ThetaLambdaFromC[cent]->SetName(hNameCos2ThetaLambdaFromC_3D[cent]);

    if (isProducedAcceptancePlots)
      hEtaVsPtVsCos2ThetaLambdaFromC[cent] = (TH3D *)hXiCos2ThetaFromLambdaL->Projection(2, 3, AxisNumber + 1); // eta, pt, cos2ThetaProton
    else
      hEtaVsPtVsCos2ThetaLambdaFromC[cent] = (TH3D *)hXiCos2ThetaFromLambdaL->Projection(2, 3, AxisNumber); // WRONG!
    hEtaVsPtVsCos2ThetaLambdaFromC[cent]->SetName(hNameCos2ThetaLambdaFromC_Eta3D[cent]);
    // BE CAREFUL: AXES TO BE UPDATED
    // hmassVsPsiVsCos2Theta[cent] = (TH3D *)hXiCos2Theta->Projection(3, 2, AxisNumber + 1); // mass, psi, cos2ThetaLambda
    hmassVsPsiVsCos2Theta[cent] = (TH3D *)hXiCos2Theta->Projection(3, 2, AxisNumber); // mass, psi, cos2ThetaLambda
    hmassVsPsiVsCos2Theta[cent]->SetName(hNameCos2ThetaVsPsi_3D[cent]);
    // hmassVsPsiVsCos2ThetaLambdaFromC[cent] = (TH3D *)hXiCos2Theta->Projection(3, 2, AxisNumber + 1); // mass, psi, cos2ThetaProton
    hmassVsPsiVsCos2ThetaLambdaFromC[cent] = (TH3D *)hXiCos2Theta->Projection(3, 2, AxisNumber); // mass, psi, cos2ThetaProton
    hmassVsPsiVsCos2ThetaLambdaFromC[cent]->SetName(hNameCos2ThetaVsPsiLambdaFromC_3D[cent]);
  }

  // canvasBDT
  TCanvas *canvasBDT = new TCanvas("canvasBDT", "canvasBDT", 800, 800);
  StyleCanvas(canvasBDT, 0.1, 0.05, 0.05, 0.15);
  hBDTSelection->Draw("HIST");

  // create output file
  TString SBDT = "";
  if (BDTscoreCut != DefaultBDTscoreCut || isSysMultTrial)
    SBDT = Form("_BDT%.3f", BDTscoreCut);
  TString OutputFileName = "OutputAnalysis/Output_FromTHN_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice];
  if (isApplyWeights)
    OutputFileName += "_Weighted";
  if (v2type == 1)
    OutputFileName += "_SP";
  if (!useCommonBDTValue)
    OutputFileName += "_BDTCentDep";
  if (isRun2Binning)
    OutputFileName += "_Run2Binning";
  if (!ExtrisRapiditySel)
    OutputFileName += "_Eta08";
  if (ChosenPart != 6)
    OutputFileName += SBDT;
  // if (isRedCentrality) OutputFileName += "_isRedCent";
  OutputFileName += ".root";
  TFile *file = new TFile(OutputFileName, "RECREATE");

  hBDTSelection->Write();
  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    hmassVsPtVsV2C[cent]->Write();
    hmassVsPtVsPzs2[cent]->Write();
    hmassVsPtVsPzs2LambdaFromC[cent]->Write();
    hmassVsPsiVsPz[cent]->Write();
    hmassVsPsiVsPzLambdaFromC[cent]->Write();
    hmassVsPtVsCos2Theta[cent]->Write();
    hmassVsPtVsCos2ThetaLambdaFromC[cent]->Write();
    hmassVsPsiVsCos2Theta[cent]->Write();
    hmassVsPsiVsCos2ThetaLambdaFromC[cent]->Write();
    hEtaVsPtVsCos2ThetaLambdaFromC[cent]->Write();
  }

  file->Close();

  cout << "I created the file " << file->GetName() << endl;
}
