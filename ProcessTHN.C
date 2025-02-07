// This macro was originally written by:
// chiara.de.martin@cern.ch
#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "THnSparse.h"
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

void ProcessTHN(Bool_t isEff = 0,
                Int_t indexMultTrial = 0,
                Int_t ChosenPart = ChosenParticle,
                TString inputFileName = SinputFileName,
                Int_t EtaSysChoice = ExtrEtaSysChoice,
                Bool_t isSysMultTrial = ExtrisSysMultTrial)
{

  // phi and efficiency weights should be applied in task

  Bool_t isXi = 1;
  Int_t Part = 0;
  if ((ChosenPart == 1) || (ChosenPart == 4) || (ChosenPart == 5))
  {
    isXi = 0; // Omega
    Part = 1;
  }
  if (isEff)
    inputFileName = SinputFileNameEff;

  if (isSysMultTrial)
  {
    inputFileName = SinputFileNameSyst;
    if (isEff)
      inputFileName = SinputFileNameEff;
  }

  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trialsBDT)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;

  inputFileName = "TestTHN";
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
  THnSparseF *hV2 = (THnSparseF *)dir1->Get("h" + ParticleName[Part] + "V2");
  // axes: thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassCasc, thnAxisBDTScore, thnAxisV2
  THnSparseF *hXiPzs2 = (THnSparseF *)dir1->Get("h" + ParticleName[Part] + "Pzs2");
  // axes: thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassCasc, thnAxisMassLambda, thnAxisBDTScore, thnAxisPzs2Xi, thnAxisPzs2Lambda
  THnSparseF *hXiCos2Theta = (THnSparseF *)dir1->Get("h" + ParticleName[Part] + "Cos2Theta");
  // axes: thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassCasc, thnAxisMassLambda, thnAxisBDTScore, thnAxisCos2Theta, thnAxisCos2Theta

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
  TString hNameCos2Theta_3D[numCent + 1] = {""};
  TString hNameCos2ThetaLambdaFromC_3D[numCent + 1] = {""};
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
  TH1F *hDummyPt = (TH1F *)hV2->Projection(2);
  TH1F *hDummyMass = (TH1F *)hV2->Projection(3);
  TH1F *hDummyBDT = (TH1F *)hV2->Projection(4);
  TH1F *hDummyMassLambda = (TH1F *)hXiPzs2->Projection(4);
  hDummyCentrality->Reset();
  hDummyCharge->Reset();
  hDummyPt->Reset();
  hDummyMass->Reset();
  hDummyBDT->Reset();
  hDummyMassLambda->Reset();
  hDummyCentrality->GetYaxis()->SetRangeUser(0, 0.2);
  hDummyCharge->GetYaxis()->SetRangeUser(0, 1);
  hDummyPt->GetYaxis()->SetRangeUser(0, 0.25);
  hDummyMass->GetYaxis()->SetRangeUser(0, 0.15);
  hDummyBDT->GetYaxis()->SetRangeUser(0, 0.4);
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
    // Selection of centrality range
    hV2->GetAxis(0)->SetRange(hV2->GetAxis(0)->FindBin(CentFT0CMin + 0.001), hV2->GetAxis(0)->FindBin(CentFT0CMax - 0.001));
    hXiPzs2->GetAxis(0)->SetRange(hXiPzs2->GetAxis(0)->FindBin(CentFT0CMin + 0.001), hXiPzs2->GetAxis(0)->FindBin(CentFT0CMax - 0.001));
    hXiCos2Theta->GetAxis(0)->SetRange(hXiCos2Theta->GetAxis(0)->FindBin(CentFT0CMin + 0.001), hXiCos2Theta->GetAxis(0)->FindBin(CentFT0CMax - 0.001));

    // Selection of charge
    if (ChosenParticle == 2 || ChosenParticle == 4)
    {
      hV2->GetAxis(1)->SetRange(1, 1);          // Charge < 0
      hXiPzs2->GetAxis(1)->SetRange(1, 1);      // Charge < 0
      hXiCos2Theta->GetAxis(1)->SetRange(1, 1); // Charge < 0
    }
    else if (ChosenParticle == 3 || ChosenParticle == 5)
    {
      hV2->GetAxis(1)->SetRange(2, 2);          // Charge > 0
      hXiPzs2->GetAxis(1)->SetRange(2, 2);      // Charge > 0
      hXiCos2Theta->GetAxis(1)->SetRange(2, 2); // Charge > 0
    }
    else
    {
      hV2->GetAxis(1)->SetRange(1, 2);          // All charges
      hXiPzs2->GetAxis(1)->SetRange(1, 2);      // All charges
      hXiCos2Theta->GetAxis(1)->SetRange(1, 2); // All charges
    }

    // Selection on BDT score -- be careful, in the 0-80% case, the BDT score is not defined and I apply the same cut as in the 0-10% case
    Float_t BDTscoreCut = DefaultBDTscoreCut;
    if (!useCommonBDTValue)
      BDTscoreCut = bdtCut[cent];
    hV2->GetAxis(4)->SetRange(hV2->GetAxis(4)->FindBin(BDTscoreCut + 0.001), hV2->GetAxis(4)->FindBin(1 - 0.001));

    // QCPlots
    hBDTSelection->SetBinContent(cent + 1, BDTscoreCut);
    hBDTSelection->SetBinError(cent + 1, 0);

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
    hPt[cent] = (TH1F *)hV2->Projection(2);
    hPt[cent]->Scale(1. / hPt[cent]->Integral());
    hPt[cent]->SetLineColor(ColorMult[cent]);
    hPt[cent]->SetMarkerColor(ColorMult[cent]);
    hPt[cent]->GetXaxis()->SetRangeUser(0, 10);
    if (cent == 0)
      hDummyPt->Draw();
    hPt[cent]->Draw("same");

    canvasQC->cd(4);
    hMass[cent] = (TH1F *)hV2->Projection(3);
    hMass[cent]->Scale(1. / hMass[cent]->Integral());
    hMass[cent]->SetLineColor(ColorMult[cent]);
    hMass[cent]->SetMarkerColor(ColorMult[cent]);
    hMass[cent]->GetXaxis()->SetRangeUser(1.3, 1.345);
    if (cent == 0)
      hDummyMass->Draw();
    hMass[cent]->Draw("same");

    canvasQC->cd(5);
    hBDT[cent] = (TH1F *)hV2->Projection(4);
    hBDT[cent]->Scale(1. / hBDT[cent]->Integral());
    hBDT[cent]->SetLineColor(ColorMult[cent]);
    hBDT[cent]->SetMarkerColor(ColorMult[cent]);
    hBDT[cent]->GetXaxis()->SetRangeUser(0.4, 1);
    hBDT[cent]->GetYaxis()->SetRangeUser(0, 0.4);
    if (cent == 0)
      hDummyBDT->Draw();
    hBDT[cent]->Draw("same");

    canvasQC->cd(6);
    hMassLambda[cent] = (TH1F *)hXiPzs2->Projection(4);
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
    hNameCos2ThetaVsPsi_3D[cent] = Form("massVsPsiVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNameCos2ThetaVsPsiLambdaFromC_3D[cent] = Form("massVsPsiVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax);

    // Projections
    hmassVsPtVsV2C[cent] = (TH3D *)hV2->Projection(3, 2, 5); // mass, pt, v2
    hmassVsPtVsV2C[cent]->SetName(hName[cent]);

    hmassVsPtVsPzs2[cent] = (TH3D *)hXiPzs2->Projection(3, 2, 6); // mass, pt, pzs2
    hmassVsPtVsPzs2[cent]->SetName(hNamePzs2_3D[cent]);
    hmassVsPtVsPzs2LambdaFromC[cent] = (TH3D *)hXiPzs2->Projection(3, 2, 7); // mass, pt, pzs2LambdaFromC
    hmassVsPtVsPzs2LambdaFromC[cent]->SetName(hNamePzs2LambdaFromC_3D[cent]);

    // BE CAREFUL: AXES TO BE UPDATED
    hmassVsPsiVsPz[cent] = (TH3D *)hXiPzs2->Projection(3, 2, 6); // mass, psi, cosThetaLambda
    hmassVsPsiVsPz[cent]->SetName(hNamePzVsPsi_3D[cent]);
    hmassVsPsiVsPzLambdaFromC[cent] = (TH3D *)hXiPzs2->Projection(3, 2, 6); // mass, psi, cosThetaProton
    hmassVsPsiVsPzLambdaFromC[cent]->SetName(hNamePzVsPsiLambdaFromC_3D[cent]);

    hmassVsPtVsCos2Theta[cent] = (TH3D *)hXiCos2Theta->Projection(3, 2, 6); // mass, pt, cos2ThetaLambda
    hmassVsPtVsCos2Theta[cent]->SetName(hNameCos2Theta_3D[cent]);
    hmassVsPtVsCos2ThetaLambdaFromC[cent] = (TH3D *)hXiCos2Theta->Projection(3, 2, 7); // mass, pt, cos2ThetaProton
    hmassVsPtVsCos2ThetaLambdaFromC[cent]->SetName(hNameCos2ThetaLambdaFromC_3D[cent]);

    // BE CAREFUL: AXES TO BE UPDATED
    hmassVsPsiVsCos2Theta[cent] = (TH3D *)hXiCos2Theta->Projection(3, 2, 6); // mass, psi, cos2ThetaLambda
    hmassVsPsiVsCos2Theta[cent]->SetName(hNameCos2ThetaVsPsi_3D[cent]);
    hmassVsPsiVsCos2ThetaLambdaFromC[cent] = (TH3D *)hXiCos2Theta->Projection(3, 2, 6); // mass, psi, cos2ThetaProton
    hmassVsPsiVsCos2ThetaLambdaFromC[cent]->SetName(hNameCos2ThetaVsPsiLambdaFromC_3D[cent]);
  }

  // canvasBDT
  TCanvas *canvasBDT = new TCanvas("canvasBDT", "canvasBDT", 800, 800);
  StyleCanvas(canvasBDT, 0.1, 0.05, 0.05, 0.15);
  hBDTSelection->Draw("HIST");

  // create output file
  TString SBDT = "";
  if (BDTscoreCut != DefaultBDTscoreCut)
    SBDT = Form("_BDT%.3f", BDTscoreCut);
  TString OutputFileName = "OutputAnalysis/Output_FromTHN_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (isApplyWeights)
    OutputFileName += "_Weighted";
  if (v2type == 1)
    OutputFileName += "_SP";
  if (!useCommonBDTValue)
    OutputFileName += "_BDTCentDep";
  if (isRun2Binning)
    OutputFileName += "_Run2Binning";
  OutputFileName += "_Eta08";
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
  }

  file->Close();

  cout << "I created the file " << file->GetName() << endl;
}
