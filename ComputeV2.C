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
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TPad.h"
#include "CommonVar.h"
#include "StyleFile.h"

void ComputeV2(Int_t indexMultTrial = 0, Bool_t isXi = ChosenParticleXi, TString inputFileName = SinputFileName, Int_t RebinFactor = 1, Int_t EtaSysChoice = ExtrEtaSysChoice, Bool_t isSysMultTrial = ExtrisSysMultTrial)
{

  if (isSysMultTrial)
    inputFileName = SinputFileNameSyst;
  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trialsBDT)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;
  TString SBDT = "";
  if (BDTscoreCut != DefaultBDTscoreCut)
    SBDT = Form("_BDT%.3f", BDTscoreCut);

  TString SinputFile = "OutputAnalysis/Output_" + inputFileName + "_" + ParticleName[!isXi] + ChargeName[ExtrCharge + 1] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (isApplyWeights)
    SinputFile += "_Weighted";
  if (v2type == 1)
    SinputFile += "_SP";
  if (!useCommonBDTValue)
    SinputFile += "_BDTCentDep";
  if (isRun2Binning)
    SinputFile += "_Run2Binning";
  SinputFile += ".root";
  cout << "Input file: " << SinputFile << endl;
  TFile *inputFile = new TFile(SinputFile);
  TH3D *hmassVsPtVsV2C[numCent + 1];
  TH3D *hmassVsPtVsPzs2[numCent + 1];
  TProfile2D *profmassVsPt[numCent + 1];
  TH2F *hmassVsPt[numCent + 1];
  TH2F *hmassVsV2C[numCent + 1][numPtBins];
  TH2F *hmassVsPzs2[numCent + 1][numPtBins];
  TH1F *hmass[numCent + 1][numPtBins];
  TH1F *hV2C[numCent + 1][numPtBins];
  TH1F *hPzs2[numCent + 1][numPtBins];
  TH1F *hmassVsV2Cx[numCent + 1][numPtBins];
  TH1F *hV2CFromProfile[numCent + 1][numPtBins];
  TProfile *pV2C[numCent + 1][numPtBins];
  TProfile *pPzs2[numCent + 1][numPtBins];
  TH2F *hPhiCentHisto[numCent];
  TH1F *hPhiCentHisto1D[numCent][numPtBins];
  TString hName = "";
  TString hNamePzs2_3D = "";
  TString profName = "";
  TString hNameMass = "";
  TString hNameV2C = "";
  TString hNamePzs2 = "";
  TString hNameMassV2C = "";
  TString hNameMassPzs2 = "";
  TString hNameV2CFromProfile2D = "";
  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;

  // QC: phi distribution of selected candidates in centrality classes
  gStyle->SetOptStat(0);
  TCanvas *QCPhi = new TCanvas("QCPhi", "QCPhi", 1400, 1200);
  QCPhi->Divide(4, 4);
  std::vector<double> centBins;
  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    centBins.push_back(static_cast<double>(CentFT0C[cent]));
  }

  TH3D *weights{nullptr};
  for (Int_t pt = 0; pt < numPtBins; pt++)
  {
    QCPhi->cd(pt + 1);
    for (Int_t cent = 0; cent < numCent; cent++)
    {
      // QCPlot
      hPhiCentHisto[cent] = (TH2F *)inputFile->Get(Form("PhiHist_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
      if (!hPhiCentHisto[cent])
      {
        cout << "Histogram hPhiCentHisto not available" << endl;
        return;
      }
      hPhiCentHisto[cent]->SetName(Form("hPhiCentHisto_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
      if (!weights)
      {
        std::vector<double> phiBins(hPhiCentHisto[cent]->GetYaxis()->GetNbins() + 1, 0.);
        for (uint32_t bin = 1; bin < phiBins.size() - 1; bin++)
        {
          phiBins[bin] = bin * 2 * TMath::Pi() / (phiBins.size() - 1);
        }
        phiBins.back() = 2 * TMath::Pi();
        weights = new TH3D("weights", "weights", numCent, centBins.data(), numPtBins, PtBins, hPhiCentHisto[cent]->GetYaxis()->GetNbins(), phiBins.data());
      }

      hPhiCentHisto[cent]->GetXaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
      hPhiCentHisto1D[cent][pt] = (TH1F *)hPhiCentHisto[cent]->ProjectionY(Form("Weights_cent%i-%i_pt%.3f-%.3f", CentFT0C[cent], CentFT0C[cent + 1], PtBins[pt], PtBins[pt + 1]));
      hPhiCentHisto1D[cent][pt]->SetMarkerStyle(20);
      hPhiCentHisto1D[cent][pt]->SetMarkerSize(0.5);
      hPhiCentHisto1D[cent][pt]->SetMarkerColor(ColorMult[cent]);
      hPhiCentHisto1D[cent][pt]->SetLineColor(ColorMult[cent]);
      hPhiCentHisto1D[cent][pt]->SetTitle(Form("%.2f < p_{T} < %.2f GeV/c", PtBins[pt], PtBins[pt + 1]));
      // hPhiCentHisto1D[cent][pt]->Scale(1. / hPhiCentHisto1D[cent]->Integral());
      hPhiCentHisto1D[cent][pt]->Scale(1. / hPhiCentHisto1D[cent][pt]->GetBinContent(hPhiCentHisto1D[cent][pt]->GetMinimumBin()));
      hPhiCentHisto1D[cent][pt]->GetYaxis()->SetRangeUser(0, 1.2 * hPhiCentHisto1D[cent][pt]->GetBinContent(hPhiCentHisto1D[cent][pt]->GetMaximumBin()));
      hPhiCentHisto1D[cent][pt]->GetXaxis()->SetRangeUser(0, 2 * TMath::Pi());
      hPhiCentHisto1D[cent][pt]->Draw("same hist");

      for (Int_t bin = 0; bin < hPhiCentHisto1D[cent][pt]->GetNbinsX(); bin++)
      {
        weights->Fill(CentFT0C[cent] + 0.0001, PtBins[pt] + 0.0001, hPhiCentHisto1D[cent][pt]->GetBinCenter(bin + 1), hPhiCentHisto1D[cent][pt]->GetBinContent(bin + 1));
      }
    }
  }

  // average pt for each pt interval: an estimate based on all selected candidates
  TCanvas *cAvgPt = new TCanvas("cAvgPt", "cAvgPt", 1400, 1200);
  TCanvas *cPtDeviation = new TCanvas("cPtDeviation", "cPtDeviation", 1400, 1200);
  TH1F *hHistoPt[numCent + 1];
  TH1F *hAvgPt[numCent + 1];
  TH1F *hPtDeviation[numCent + 1];
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
    hPhiCentHisto[cent] = (TH2F *)inputFile->Get(Form("PhiHist_cent%i-%i", CentFT0CMin, CentFT0CMax));
    if (!hPhiCentHisto[cent])
    {
      cout << "Histogram hPhiCentHisto not available" << endl;
      return;
    }
    hPhiCentHisto[cent]->GetYaxis()->SetRangeUser(0, 2 * TMath::Pi());
    hHistoPt[cent] = (TH1F *)hPhiCentHisto[cent]->ProjectionX(Form("PtHist_cent%i-%i", 0, -1));
    hAvgPt[cent] = new TH1F(Form("AvgPt_cent%i-%i", CentFT0CMin, CentFT0CMax), Form("AvgPt_cent%i-%i", CentFT0CMin, CentFT0CMax), numPtBins, PtBins);
    hPtDeviation[cent] = new TH1F(Form("PtDeviation_cent%i-%i", CentFT0CMin, CentFT0CMax), Form("PtDeviation_cent%i-%i", CentFT0CMin, CentFT0CMax), numPtBins, PtBins);
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      hHistoPt[cent]->GetXaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
      hAvgPt[cent]->SetBinContent(pt + 1, hHistoPt[cent]->GetMean());
      hAvgPt[cent]->SetBinError(pt + 1, hHistoPt[cent]->GetMeanError());
      hPtDeviation[cent]->SetBinContent(pt + 1, (hHistoPt[cent]->GetMean() - (PtBins[pt] + PtBins[pt + 1]) / 2) / ((PtBins[pt] + PtBins[pt + 1]) / 2));
      hPtDeviation[cent]->SetBinError(pt + 1, 0);
    }
    cAvgPt->cd();
    hAvgPt[cent]->SetMarkerStyle(20);
    hAvgPt[cent]->SetMarkerSize(0.5);
    hAvgPt[cent]->SetMarkerColor(ColorMult[cent]);
    hAvgPt[cent]->SetLineColor(ColorMult[cent]);
    hAvgPt[cent]->Draw("same");

    cPtDeviation->cd();
    hPtDeviation[cent]->SetMarkerStyle(20);
    hPtDeviation[cent]->SetMarkerSize(0.5);
    hPtDeviation[cent]->SetMarkerColor(ColorMult[cent]);
    hPtDeviation[cent]->SetLineColor(ColorMult[cent]);
    hPtDeviation[cent]->Draw("same");
  }

  // v2 and polarization computation
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
    hName = Form("massVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax);
    profName = Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNamePzs2_3D = Form("massVsPtVsPzs2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hmassVsPtVsV2C[cent] = (TH3D *)inputFile->Get(hName);
    hmassVsPt[cent] = (TH2F *)hmassVsPtVsV2C[cent]->Project3D("yx");
    profmassVsPt[cent] = (TProfile2D *)inputFile->Get(profName);
    hmassVsPtVsPzs2[cent] = (TH3D *)inputFile->Get(hNamePzs2_3D);
    if (!hmassVsPtVsPzs2[cent])
    {
      cout << "Histogram hmassVsPtVsPzs2 not available" << endl;
      return;
    }
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      hNameMass = Form("mass_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameV2C = Form("V2C_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNamePzs2 = Form("Pzs2_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameV2CFromProfile2D = Form("V2CFromProfile_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameMassV2C = Form("MassvsV2C_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameMassPzs2 = Form("MassvsPzs2_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);

      hmassVsPtVsV2C[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);

      profmassVsPt[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
      hV2CFromProfile[cent][pt] = (TH1F *)profmassVsPt[cent]->ProjectionX(hNameV2CFromProfile2D, 0, -1, "e"); // v2C from TProfile2D

      hmassVsV2C[cent][pt] = (TH2F *)hmassVsPtVsV2C[cent]->Project3D("xze"); // mass vs V2C //"e" option does not change the results
      hmassVsV2C[cent][pt]->SetName(hNameMassV2C);
      hmassVsV2C[cent][pt]->RebinY(RebinFactor);

      hmassVsPtVsPzs2[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
      hmassVsPzs2[cent][pt] = (TH2F *)hmassVsPtVsPzs2[cent]->Project3D("xze"); // mass vs Pzs2
      hmassVsPzs2[cent][pt]->SetName(hNameMassPzs2);
      hmassVsPzs2[cent][pt]->RebinY(RebinFactor);

      hmass[cent][pt] = (TH1F *)hmassVsPtVsV2C[cent]->Project3D("xe"); // mass
      hmass[cent][pt]->SetName(hNameMass);
      hmass[cent][pt]->Rebin(RebinFactor);

      pV2C[cent][pt] = hmassVsV2C[cent][pt]->ProfileY(); // v2C //the error is the standard error of the mean
      pV2C[cent][pt]->SetName(hNameV2C + "_Profile");

      pPzs2[cent][pt] = hmassVsPzs2[cent][pt]->ProfileY(); // Pzs2 //the error is the standard error of the mean
      pPzs2[cent][pt]->SetName(hNamePzs2 + "_Profile");

      hV2C[cent][pt] = (TH1F *)hmass[cent][pt]->Clone(hNameV2C);
      hPzs2[cent][pt] = (TH1F *)hmass[cent][pt]->Clone(hNamePzs2);
      for (Int_t bin = 0; bin < hmass[cent][pt]->GetNbinsX(); bin++)
      {
        // hmassVsV2Cx[cent][pt] = (TH1F *)hmassVsV2C[cent][pt]->ProjectionX("", bin + 1, bin + 1);
        // hmassVsV2Cx[cent][pt]->SetName(Form("MassvsV2_%i_%i_%i", cent, pt, bin));
        // hmassVsV2Cx[cent][pt]->ResetStats();
        // hV2C[cent][pt]->SetBinContent(bin + 1, hmassVsV2Cx[cent][pt]->GetMean());
        // hV2C[cent][pt]->SetBinError(bin + 1, hmassVsV2Cx[cent][pt]->GetMeanError());
        hV2C[cent][pt]->SetBinContent(bin + 1, hmassVsV2C[cent][pt]->ProjectionX("", bin + 1, bin + 1)->GetMean());
        hV2C[cent][pt]->SetBinError(bin + 1, hmassVsV2C[cent][pt]->ProjectionX("", bin + 1, bin + 1)->GetMeanError());
        hPzs2[cent][pt]->SetBinContent(bin + 1, hmassVsPzs2[cent][pt]->ProjectionX("", bin + 1, bin + 1)->GetMean());
        hPzs2[cent][pt]->SetBinError(bin + 1, hmassVsPzs2[cent][pt]->ProjectionX("", bin + 1, bin + 1)->GetMeanError());
      }
    }
  }

  QCPhi->SaveAs("QCPlots/QCPhiCasc.png");

  TString SOutputFile = "OutputAnalysis/V2_" + inputFileName + "_" + ParticleName[!isXi] + ChargeName[ExtrCharge + 1] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (isApplyWeights)
    SOutputFile += "_Weighted";
  if (v2type == 1)
    SOutputFile += "_SP";
  if (!useCommonBDTValue)
    SOutputFile += "_BDTCentDep";
  if (isRun2Binning)
    SOutputFile += "_Run2Binning";
  // SOutputFile += "_TestV2InRestrictedRange.root";
  SOutputFile += ".root";
  TFile *file = new TFile(SOutputFile, "RECREATE");
  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      hmass[cent][pt]->Write();
      hV2C[cent][pt]->Write();
      hPzs2[cent][pt]->Write();
      hmassVsV2C[cent][pt]->Write();
      hmassVsPzs2[cent][pt]->Write();
      pV2C[cent][pt]->Write();
      pPzs2[cent][pt]->Write();
      hV2CFromProfile[cent][pt]->Write();
      if (cent != numCent)
        hPhiCentHisto1D[cent][pt]->Write();
    }
    hmassVsPtVsV2C[cent]->Write();
    hmassVsPtVsPzs2[cent]->Write();
    hmassVsPt[cent]->Write();
    hAvgPt[cent]->Write();
  }
  file->Close();
  TString SweightsFile = "PhiWeights/Weights_" + inputFileName + "_" + ParticleName[!isXi] + SEtaSysChoice[EtaSysChoice] + SBDT + ".root";
  TFile *weightsFile;
  if (!isApplyWeights) // weights computed only once (when we apply weights, they have already been created!)
  {
    weightsFile = new TFile(SweightsFile, "RECREATE");
    weightsFile->cd();
    weights->Write();
    weightsFile->Close();
  }
  cout << "I created the file " << file->GetName() << endl;
  if (!isApplyWeights)
    cout << " and the file with weights: " << SweightsFile << endl;
}
