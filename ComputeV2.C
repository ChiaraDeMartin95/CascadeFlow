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

void ComputeV2(Int_t indexMultTrial = 0,
               Int_t ChosenPart = ChosenParticle,
               Bool_t isRapiditySel = ExtrisRapiditySel,
               TString inputFileName = SinputFileName,
               Int_t RebinFactor = 1,
               Int_t EtaSysChoice = ExtrEtaSysChoice,
               Bool_t isSysMultTrial = ExtrisSysMultTrial)
{
  if (isReducedPtBins && numPtBins != numPtBinsReduced)
  {
    cout << "Reduced pt bins are selected, but numPtBins is not set to numPtBinsReduced. Please check the settings." << endl;
    return;
  }
  if (isSysMultTrial)
    inputFileName = SinputFileNameSyst;
  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trialsBDT)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;
  TString SBDT = "";
  if (BDTscoreCut != DefaultBDTscoreCut || isSysMultTrial)
    SBDT = Form("_BDT%.3f", BDTscoreCut);

  TString SinputFile = "OutputAnalysis/Output" + STHN[ExtrisFromTHN] + "_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice]; // + SBDT;
  if (isApplyWeights)
    SinputFile += "_Weighted";
  if (v2type == 1)
    SinputFile += "_SP";
  if (!useCommonBDTValue)
    SinputFile += "_BDTCentDep";
  if (isRun2Binning)
    SinputFile += "_Run2Binning";
  if (!isRapiditySel)
    SinputFile += "_Eta08";
  SinputFile += SBDT;
  if (ChosenPart == 6 && isSysMultTrial)
    SinputFile += Form("_SysMultTrial_%i", indexMultTrial);
  SinputFile += ".root";
  cout << "Input file: " << SinputFile << endl;
  TFile *inputFile = new TFile(SinputFile);
  TH3D *hmassVsPtVsV2C[numCent + 1];

  TH3D *hmassVsPtVsPzs2[numCent + 1];
  TH3D *hmassVsPtVsPzs2LambdaFromC[numCent + 1];
  TH3D *hmassVsPsiVsPz[numCent + 1];
  TH3D *hmassVsPsiVsPzLambdaFromC[numCent + 1];

  TProfile2D *profmassVsPt[numCent + 1];
  TH2F *hmassVsPt[numCent + 1];
  TH2F *hmassVsV2C[numCent + 1][numPtBins + 1];

  TH2F *hmassVsPzs2[numCent + 1][numPtBins + 1];
  TH2F *hmassVsPzs2LambdaFromC[numCent + 1][numPtBins + 1];
  TH2F *hmassVsPz[numCent + 1][numPsiBins + 1];
  TH2F *hmassVsPzLambdaFromC[numCent + 1][numPsiBins + 1];

  TH1F *hmass[numCent + 1][numPtBins + 1];
  TH1F *hmassPsi[numCent + 1][numPsiBins + 1];
  TH1F *hV2C[numCent + 1][numPtBins + 1];

  TH1F *hPzs2[numCent + 1][numPtBins + 1];
  TH1F *hPzs2LambdaFromC[numCent + 1][numPtBins + 1];
  TH1F *hPz[numCent + 1][numPsiBins + 1];
  TH1F *hPzLambdaFromC[numCent + 1][numPsiBins + 1];

  TH1F *hmassVsV2Cx[numCent + 1][numPtBins + 1];
  TH1F *hV2CFromProfile[numCent + 1][numPtBins + 1];
  TProfile *pV2C[numCent + 1][numPtBins + 1];
  TProfile *pPzs2[numCent + 1][numPtBins + 1];
  TProfile *pPzs2LambdaFromC[numCent + 1][numPtBins + 1];
  TProfile *pPz[numCent + 1][numPsiBins + 1];
  TProfile *pPzLambdaFromC[numCent + 1][numPsiBins + 1];
  TH2F *hPhiCentHisto[numCent];
  TH1F *hPhiCentHisto1D[numCent][numPtBins + 1];

  TString hName[numCent + 1] = {""};
  TString hNamePzs2_3D[numCent + 1] = {""};
  TString hNamePzs2LambdaFromC_3D[numCent + 1] = {""};
  TString hNamePzVsPsi_3D[numCent + 1] = {""};
  TString hNamePzVsPsiLambdaFromC_3D[numCent + 1] = {""};
  TString profName[numCent + 1] = {""};

  TString hNameMass[numCent + 1][numPtBins + 1] = {""};
  TString hNameV2C[numCent + 1][numPtBins + 1] = {""};
  TString hNamePzs2[numCent + 1][numPtBins + 1] = {""};
  TString hNamePzs2LambdaFromC[numCent + 1][numPtBins + 1] = {""};
  TString hNamePz[numCent + 1][numPsiBins + 1] = {""};
  TString hNamePzLambdaFromC[numCent + 1][numPsiBins + 1] = {""};

  TString hNameMassV2C[numCent + 1][numPtBins + 1] = {""};
  TString hNameMassPzs2[numCent + 1][numPtBins + 1] = {""};
  TString hNameMassPzs2LambdaFromC[numCent + 1][numPtBins + 1] = {""};
  TString hNameMassPz[numCent + 1][numPsiBins + 1] = {""};
  TString hNameMassPzLambdaFromC[numCent + 1][numPsiBins + 1] = {""};

  TString hNameMassPsi[numCent + 1][numPsiBins + 1] = {""};
  TString hNameV2CFromProfile2D[numCent + 1][numPtBins + 1] = {""};

  // acceptance correction ***
  TH3D *hmassVsPtVsCos2Theta[numCent + 1];
  TH3D *hmassVsPtVsCos2ThetaLambdaFromC[numCent + 1];
  TH3D *hmassVsPsiVsCos2Theta[numCent + 1];
  TH3D *hmassVsPsiVsCos2ThetaLambdaFromC[numCent + 1];
  TString hNameCos2Theta_3D[numCent + 1] = {""};
  TString hNameCos2ThetaLambdaFromC_3D[numCent + 1] = {""};
  TString hNameCos2ThetaVsPsi_3D[numCent + 1] = {""};
  TString hNameCos2ThetaVsPsiLambdaFromC_3D[numCent + 1] = {""};

  TH2F *hmassVsCos2Theta[numCent + 1][numPtBins + 1];
  TH2F *hmassVsCos2ThetaLambdaFromC[numCent + 1][numPtBins + 1];
  TH2F *hmassVsCos2ThetaPsi[numCent + 1][numPsiBins + 1];
  TH2F *hmassVsCos2ThetaPsiLambdaFromC[numCent + 1][numPsiBins + 1];
  TString hNameMassCos2Theta[numCent + 1][numPtBins + 1] = {""};
  TString hNameMassCos2ThetaLambdaFromC[numCent + 1][numPtBins + 1] = {""};
  TString hNameMassCos2ThetaPsi[numCent + 1][numPsiBins + 1] = {""};
  TString hNameMassCos2ThetaPsiLambdaFromC[numCent + 1][numPsiBins + 1] = {""};

  TH1F *hCos2Theta[numCent + 1][numPtBins + 1];
  TH1F *hCos2ThetaLambdaFromC[numCent + 1][numPtBins + 1];
  TH1F *hCos2ThetaPsi[numCent + 1][numPsiBins + 1];
  TH1F *hCos2ThetaPsiLambdaFromC[numCent + 1][numPsiBins + 1];
  TString hNameCos2Theta[numCent + 1][numPtBins + 1] = {""};
  TString hNameCos2ThetaLambdaFromC[numCent + 1][numPtBins + 1] = {""};
  TString hNameCos2ThetaPsi[numCent + 1][numPsiBins + 1] = {""};
  TString hNameCos2ThetaPsiLambdaFromC[numCent + 1][numPsiBins + 1] = {""};
  // end of acceptance correction ***

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;

  // QC: phi distribution of selected candidates in centrality classes
  gStyle->SetOptStat(0);
  TCanvas *QCPhi = new TCanvas("QCPhi", "QCPhi", 1400, 1200);
  QCPhi->Divide(4, 4);
  std::vector<double> centBins;
  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    if (isOOCentrality && cent > numCentLambdaOO)
    {
      continue;
    }
    centBins.push_back(static_cast<double>(CentFT0C[cent]));
  }

  TH3D *weights{nullptr};
  if (!ExtrisFromTHN)
  {
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      QCPhi->cd(pt + 1);
      for (Int_t cent = 0; cent < numCent; cent++)
      {
        if (isOOCentrality && cent > (numCentLambdaOO - 1))
        {
          continue;
        }
        // QCPlot
        if (ChosenPart == 6)
        {
          hPhiCentHisto[cent] = (TH2F *)inputFile->Get(Form("PhiHist_cent%i-%i", CentFT0CLambdaOO[cent], CentFT0CLambdaOO[cent + 1]));
          cout << "Using Lambda centrality bins: " << CentFT0CLambdaOO[cent] << " - " << CentFT0CLambdaOO[cent + 1] << endl;
        }
        else
          hPhiCentHisto[cent] = (TH2F *)inputFile->Get(Form("PhiHist_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
        if (!hPhiCentHisto[cent])
        {
          cout << "Histogram hPhiCentHisto not available" << endl;
          return;
        }
        hPhiCentHisto[cent]->SetName(Form("hPhiCentHisto_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
        if (ChosenPart == 6)
          hPhiCentHisto[cent]->SetName(Form("hPhiCentHisto_cent%i-%i", CentFT0CLambdaOO[cent], CentFT0CLambdaOO[cent + 1]));
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
  }

  // average pt for each pt interval: an estimate based on all selected candidates
  TCanvas *cAvgPt = new TCanvas("cAvgPt", "cAvgPt", 1400, 1200);
  TCanvas *cPtDeviation = new TCanvas("cPtDeviation", "cPtDeviation", 1400, 1200);
  TH1F *hHistoPt[numCent + 1];
  TH1F *hAvgPt[numCent + 1];
  TH1F *hPtDeviation[numCent + 1];
  if (!ExtrisFromTHN)
  {
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
      if (isOOCentrality)
      {
        if (cent > numCentLambdaOO)
          continue;
        if (cent == numCentLambdaOO)
        {
          CentFT0CMin = 0;
          CentFT0CMax = 90;
        }
        else
        {
          CentFT0CMin = CentFT0CLambdaOO[cent];
          CentFT0CMax = CentFT0CLambdaOO[cent + 1];
        }
      }
      hPhiCentHisto[cent] = (TH2F *)inputFile->Get(Form("PhiHist_cent%i-%i", CentFT0CMin, CentFT0CMax));
      cout << "CentMin " << CentFT0CMin << " CentMax " << CentFT0CMax << endl;
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
    if (isOOCentrality)
    {
      if (cent > numCentLambdaOO)
        continue;
      if (cent == numCentLambdaOO)
      {
        CentFT0CMin = 0;
        CentFT0CMax = 90;
      }
      else
      {
        CentFT0CMin = CentFT0CLambdaOO[cent];
        CentFT0CMax = CentFT0CLambdaOO[cent + 1];
      }
    }
    hName[cent] = Form("massVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax);
    if (ExtrisApplyEffWeights)
      hName[cent] = Form("massVsPtVsV2CWeighted_cent%i-%i", CentFT0CMin, CentFT0CMax);
    profName[cent] = Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax);

    hNamePzs2_3D[cent] = Form("massVsPtVsPzs2_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);
    if (ChosenPart == 6)
      hNamePzs2_3D[cent] = Form("massVsPtVsPzs2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNamePzs2LambdaFromC_3D[cent] = Form("massVsPtVsPzs2LambdaFromC_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);
    if (ChosenPart == 6)
      hNamePzs2LambdaFromC_3D[cent] = Form("massVsPtVsPzs2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNamePzVsPsi_3D[cent] = Form("massVsPsiVsPz_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);
    if (ChosenPart == 6)
      hNamePzVsPsi_3D[cent] = Form("massVsPsiVsPz_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNamePzVsPsiLambdaFromC_3D[cent] = Form("massVsPsiVsPzLambdaFromC_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);
    if (ChosenPart == 6)
      hNamePzVsPsiLambdaFromC_3D[cent] = Form("massVsPsiVsPz_cent%i-%i", CentFT0CMin, CentFT0CMax);

    hNameCos2Theta_3D[cent] = Form("massVsPtVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNameCos2ThetaLambdaFromC_3D[cent] = Form("massVsPtVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax);
    if (ChosenPart == 6)
      hNameCos2ThetaLambdaFromC_3D[cent] = Form("massVsPtVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNameCos2ThetaVsPsi_3D[cent] = Form("massVsPsiVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNameCos2ThetaVsPsiLambdaFromC_3D[cent] = Form("massVsPsiVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax);
    if (ChosenPart == 6)
      hNameCos2ThetaVsPsiLambdaFromC_3D[cent] = Form("massVsPsiVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hmassVsPtVsV2C[cent] = (TH3D *)inputFile->Get(hName[cent]);
    hmassVsPt[cent] = (TH2F *)hmassVsPtVsV2C[cent]->Project3D("yx");
    hmassVsPt[cent]->SetName(Form("massVsPt_cent%i-%i", CentFT0CMin, CentFT0CMax));
    if (!ExtrisFromTHN && ChosenPart != 6)
      profmassVsPt[cent] = (TProfile2D *)inputFile->Get(profName[cent]);
    if (!profmassVsPt[cent] && !ExtrisFromTHN && ChosenPart != 6)
    {
      cout << "Profile profmassVsPt not available" << endl;
      return;
    }
    hmassVsPtVsPzs2[cent] = (TH3D *)inputFile->Get(hNamePzs2_3D[cent]);
    if (!hmassVsPtVsPzs2[cent])
    {
      cout << "Histogram hmassVsPtVsPzs2 not available" << endl;
      return;
    }
    hmassVsPtVsPzs2LambdaFromC[cent] = (TH3D *)inputFile->Get(hNamePzs2LambdaFromC_3D[cent]);
    if (!hmassVsPtVsPzs2LambdaFromC[cent])
    {
      cout << "Histogram hmassVsPtVsPzs2LambdaFromC not available" << endl;
      return;
    }
    hmassVsPtVsCos2Theta[cent] = (TH3D *)inputFile->Get(hNameCos2Theta_3D[cent]);
    if (!hmassVsPtVsCos2Theta[cent])
    {
      cout << "Histogram hmassVsPtVsCos2Theta not available" << endl;
      return;
    }
    hmassVsPtVsCos2ThetaLambdaFromC[cent] = (TH3D *)inputFile->Get(hNameCos2ThetaLambdaFromC_3D[cent]);
    if (!hmassVsPtVsCos2ThetaLambdaFromC[cent])
    {
      cout << "Histogram hmassVsPtVsCos2ThetaLambdaFromC not available" << endl;
      return;
    }
    hmassVsPsiVsPz[cent] = (TH3D *)inputFile->Get(hNamePzVsPsi_3D[cent]);
    if (!hmassVsPsiVsPz[cent])
    {
      cout << "Histogram hmassVsPsiVsPz not available" << endl;
      return;
    }
    hmassVsPsiVsPzLambdaFromC[cent] = (TH3D *)inputFile->Get(hNamePzVsPsiLambdaFromC_3D[cent]);
    if (!hmassVsPsiVsPzLambdaFromC[cent])
    {
      cout << "Histogram hmassVsPsiVsPzLambdaFromC not available" << endl;
      return;
    }
    hmassVsPsiVsCos2Theta[cent] = (TH3D *)inputFile->Get(hNameCos2ThetaVsPsi_3D[cent]);
    if (!hmassVsPsiVsCos2Theta[cent])
    {
      cout << "Histogram hmassVsPsiVsCos2Theta not available" << endl;
      return;
    }
    hmassVsPsiVsCos2ThetaLambdaFromC[cent] = (TH3D *)inputFile->Get(hNameCos2ThetaVsPsiLambdaFromC_3D[cent]);
    if (!hmassVsPsiVsCos2ThetaLambdaFromC[cent])
    {
      cout << "Histogram hmassVsPsiVsCos2ThetaLambdaFromC not available" << endl;
      return;
    }

    // psi binning for polarization
    Double_t PhiBins[numPsiBins + 1];
    for (Int_t psi = 0; psi < numPsiBins; psi++)
    {
      PhiBins[psi] = psi * 2 * TMath::Pi() / numPsiBins;
      hNameMassPsi[cent][psi] = Form("mass_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, psi);
      hNameMassPz[cent][psi] = Form("MassvsPz_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, psi);
      hNameMassPzLambdaFromC[cent][psi] = Form("MassvsPzLambdaFromC_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, psi);
      hNameMassCos2ThetaPsi[cent][psi] = Form("MassvsCos2Theta_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, psi);
      hNameMassCos2ThetaPsiLambdaFromC[cent][psi] = Form("MassvsCos2ThetaLambdaFromC_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, psi);

      hmassVsPsiVsPz[cent]->GetYaxis()->SetRangeUser(PhiBins[psi] + 0.0001, PhiBins[psi] + 2 * TMath::Pi() / numPsiBins - 0.0001);
      hmassPsi[cent][psi] = (TH1F *)hmassVsPsiVsPz[cent]->Project3D("xe"); // mass
      hmassPsi[cent][psi]->SetName(hNameMassPsi[cent][psi]);
      hmassPsi[cent][psi]->Rebin(RebinFactor);

      hmassVsPz[cent][psi] = (TH2F *)hmassVsPsiVsPz[cent]->Project3D("xze"); // mass & Pz 2D histo
      hmassVsPz[cent][psi]->SetName(hNameMassPz[cent][psi]);

      hmassVsPsiVsPzLambdaFromC[cent]->GetYaxis()->SetRangeUser(PhiBins[psi] + 0.0001, PhiBins[psi] + 2 * TMath::Pi() / numPsiBins - 0.0001);
      hmassVsPzLambdaFromC[cent][psi] = (TH2F *)hmassVsPsiVsPzLambdaFromC[cent]->Project3D("xze"); // mass & Pz 2D histo
      hmassVsPzLambdaFromC[cent][psi]->SetName(hNameMassPzLambdaFromC[cent][psi]);

      hmassVsPsiVsCos2Theta[cent]->GetYaxis()->SetRangeUser(PhiBins[psi] + 0.0001, PhiBins[psi] + 2 * TMath::Pi() / numPsiBins - 0.0001);
      hmassVsCos2ThetaPsi[cent][psi] = (TH2F *)hmassVsPsiVsCos2Theta[cent]->Project3D("xze"); // mass & cos2theta 2D histo
      hmassVsCos2ThetaPsi[cent][psi]->SetName(hNameMassCos2ThetaPsi[cent][psi]);

      hmassVsPsiVsCos2ThetaLambdaFromC[cent]->GetYaxis()->SetRangeUser(PhiBins[psi] + 0.0001, PhiBins[psi] + 2 * TMath::Pi() / numPsiBins - 0.0001);
      hmassVsCos2ThetaPsiLambdaFromC[cent][psi] = (TH2F *)hmassVsPsiVsCos2ThetaLambdaFromC[cent]->Project3D("xze"); // mass & cos2theta 2D histo
      hmassVsCos2ThetaPsiLambdaFromC[cent][psi]->SetName(hNameMassCos2ThetaPsiLambdaFromC[cent][psi]);

      hPz[cent][psi] = (TH1F *)hmassVsPsiVsPz[cent]->Project3D("xe"); // Pz vs mass
      hNamePz[cent][psi] = Form("Pz_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, psi);
      hPz[cent][psi]->SetName(hNamePz[cent][psi]);
      hPz[cent][psi]->Reset();

      hCos2ThetaPsi[cent][psi] = (TH1F *)hmassVsPsiVsCos2Theta[cent]->Project3D("xe"); // cos2theta vs mass
      hNameCos2ThetaPsi[cent][psi] = Form("Cos2Theta_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, psi);
      hCos2ThetaPsi[cent][psi]->SetName(hNameCos2ThetaPsi[cent][psi]);
      hCos2ThetaPsi[cent][psi]->Reset();

      for (Int_t bin = 0; bin < hmassVsPz[cent][psi]->GetNbinsY(); bin++)
      {
        TH1D *htemp = (TH1D *)hmassVsPz[cent][psi]->ProjectionX(Form("hhtemp_%i", bin), bin + 1, bin + 1);
        hPz[cent][psi]->SetBinContent(bin + 1, htemp->GetMean());
        hPz[cent][psi]->SetBinError(bin + 1, htemp->GetMeanError());
        TH1D *htemp2 = (TH1D *)hmassVsCos2ThetaPsi[cent][psi]->ProjectionX(Form("hhtemp2_%i", bin), bin + 1, bin + 1);
        hCos2ThetaPsi[cent][psi]->SetBinContent(bin + 1, htemp2->GetMean());
        hCos2ThetaPsi[cent][psi]->SetBinError(bin + 1, htemp2->GetMeanError());
      }
      pPz[cent][psi] = hmassVsPz[cent][psi]->ProfileY();
      pPz[cent][psi]->SetName(hNamePz[cent][psi] + "_Profile");

      hPzLambdaFromC[cent][psi] = (TH1F *)hmassVsPsiVsPzLambdaFromC[cent]->Project3D("xe"); // Pz vs mass
      hNamePzLambdaFromC[cent][psi] = Form("PzLambdaFromC_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, psi);
      hPzLambdaFromC[cent][psi]->SetName(hNamePzLambdaFromC[cent][psi]);
      hPzLambdaFromC[cent][psi]->Reset();

      hCos2ThetaPsiLambdaFromC[cent][psi] = (TH1F *)hmassVsPsiVsCos2ThetaLambdaFromC[cent]->Project3D("xe"); // cos2theta vs mass
      hNameCos2ThetaPsiLambdaFromC[cent][psi] = Form("Cos2ThetaLambdaFromC_cent%i-%i_psi%i", CentFT0CMin, CentFT0CMax, psi);
      hCos2ThetaPsiLambdaFromC[cent][psi]->SetName(hNameCos2ThetaPsiLambdaFromC[cent][psi]);
      hCos2ThetaPsiLambdaFromC[cent][psi]->Reset();

      for (Int_t bin = 0; bin < hmassVsPzLambdaFromC[cent][psi]->GetNbinsY(); bin++)
      {
        TH1D *htemp = (TH1D *)hmassVsPzLambdaFromC[cent][psi]->ProjectionX(Form("%i", bin), bin + 1, bin + 1);
        hPzLambdaFromC[cent][psi]->SetBinContent(bin + 1, htemp->GetMean());
        hPzLambdaFromC[cent][psi]->SetBinError(bin + 1, htemp->GetMeanError());
        TH1D *htemp2 = (TH1D *)hmassVsCos2ThetaPsiLambdaFromC[cent][psi]->ProjectionX(Form("%i", bin), bin + 1, bin + 1);
        hCos2ThetaPsiLambdaFromC[cent][psi]->SetBinContent(bin + 1, htemp2->GetMean());
        hCos2ThetaPsiLambdaFromC[cent][psi]->SetBinError(bin + 1, htemp2->GetMeanError());
      }
      pPzLambdaFromC[cent][psi] = hmassVsPzLambdaFromC[cent][psi]->ProfileY();
      pPzLambdaFromC[cent][psi]->SetName(hNamePzLambdaFromC[cent][psi] + "_Profile");
    }

    for (Int_t pt = 0; pt < numPtBins + 1; pt++)
    {
      hNameMass[cent][pt] = Form("mass_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameV2C[cent][pt] = Form("V2C_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNamePzs2[cent][pt] = Form("Pzs2_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNamePzs2LambdaFromC[cent][pt] = Form("Pzs2LambdaFromC_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameCos2Theta[cent][pt] = Form("Cos2Theta_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameCos2ThetaLambdaFromC[cent][pt] = Form("Cos2ThetaLambdaFromC_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameV2CFromProfile2D[cent][pt] = Form("V2CFromProfile_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameMassV2C[cent][pt] = Form("MassvsV2C_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameMassPzs2[cent][pt] = Form("MassvsPzs2_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameMassPzs2LambdaFromC[cent][pt] = Form("MassvsPzs2LambdaFromC_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameMassCos2Theta[cent][pt] = Form("MassvsCos2Theta_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameMassCos2ThetaLambdaFromC[cent][pt] = Form("MassvsCos2ThetaLambdaFromC_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);

      if (pt == numPtBins)
      {
        hmassVsPtVsV2C[cent]->GetYaxis()->SetRangeUser(PtBins[0] + 0.0001, PtBins[numPtBins] - 0.0001);
        if (!ExtrisFromTHN && ChosenPart != 6)
          profmassVsPt[cent]->GetYaxis()->SetRangeUser(PtBins[0] + 0.0001, PtBins[numPtBins] - 0.0001);
        hmassVsPtVsPzs2[cent]->GetYaxis()->SetRangeUser(PtBins[0] + 0.0001, PtBins[numPtBins] - 0.0001);
        hmassVsPtVsPzs2LambdaFromC[cent]->GetYaxis()->SetRangeUser(PtBins[0] + 0.0001, PtBins[numPtBins] - 0.0001);
        hmassVsPtVsCos2Theta[cent]->GetYaxis()->SetRangeUser(PtBins[0] + 0.0001, PtBins[numPtBins] - 0.0001);
        hmassVsPtVsCos2ThetaLambdaFromC[cent]->GetYaxis()->SetRangeUser(PtBins[0] + 0.0001, PtBins[numPtBins] - 0.0001);
      }
      else
      {
        hmassVsPtVsV2C[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
        if (!ExtrisFromTHN && ChosenPart != 6)
          profmassVsPt[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
        hmassVsPtVsPzs2[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
        hmassVsPtVsPzs2LambdaFromC[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
        hmassVsPtVsCos2Theta[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
        hmassVsPtVsCos2ThetaLambdaFromC[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
      }
      if (!ExtrisFromTHN && ChosenPart != 6)
        hV2CFromProfile[cent][pt] = (TH1F *)profmassVsPt[cent]->ProjectionX(hNameV2CFromProfile2D[cent][pt], 0, -1, "e"); // v2C from TProfile2D
      hmassVsV2C[cent][pt] = (TH2F *)hmassVsPtVsV2C[cent]->Project3D("xze");                                              // mass vs V2C //"e" option does not change the results
      hmassVsV2C[cent][pt]->SetName(hNameMassV2C[cent][pt]);
      hmassVsV2C[cent][pt]->RebinY(RebinFactor);

      hmassVsPzs2[cent][pt] = (TH2F *)hmassVsPtVsPzs2[cent]->Project3D("xze"); // mass vs Pzs2
      hmassVsPzs2[cent][pt]->SetName(hNameMassPzs2[cent][pt]);
      hmassVsPzs2[cent][pt]->RebinY(RebinFactor);

      hmassVsPzs2LambdaFromC[cent][pt] = (TH2F *)hmassVsPtVsPzs2LambdaFromC[cent]->Project3D("xze"); // mass vs Pzs2LambdaFromC
      hmassVsPzs2LambdaFromC[cent][pt]->SetName(hNameMassPzs2LambdaFromC[cent][pt]);
      hmassVsPzs2LambdaFromC[cent][pt]->RebinY(RebinFactor);

      hmassVsCos2Theta[cent][pt] = (TH2F *)hmassVsPtVsCos2Theta[cent]->Project3D("xze"); // mass vs cos2theta
      hmassVsCos2Theta[cent][pt]->SetName(hNameMassCos2Theta[cent][pt]);
      hmassVsCos2Theta[cent][pt]->RebinY(RebinFactor);

      hmassVsCos2ThetaLambdaFromC[cent][pt] = (TH2F *)hmassVsPtVsCos2ThetaLambdaFromC[cent]->Project3D("xze"); // mass vs cos2thetaLambdaFromC
      hmassVsCos2ThetaLambdaFromC[cent][pt]->SetName(hNameMassCos2ThetaLambdaFromC[cent][pt]);
      hmassVsCos2ThetaLambdaFromC[cent][pt]->RebinY(RebinFactor);

      hmass[cent][pt] = (TH1F *)hmassVsPtVsV2C[cent]->Project3D("xe"); // mass
      hmass[cent][pt]->SetName(hNameMass[cent][pt]);
      hmass[cent][pt]->Rebin(RebinFactor);

      pV2C[cent][pt] = hmassVsV2C[cent][pt]->ProfileY(); // v2C //the error is the standard error of the mean
      pV2C[cent][pt]->SetName(hNameV2C[cent][pt] + "_Profile");

      pPzs2[cent][pt] = hmassVsPzs2[cent][pt]->ProfileY(); // Pzs2 //the error is the standard error of the mean
      pPzs2[cent][pt]->SetName(hNamePzs2[cent][pt] + "_Profile");

      pPzs2LambdaFromC[cent][pt] = hmassVsPzs2LambdaFromC[cent][pt]->ProfileY(); // Pzs2LambdaFromC //the error is the standard error of the mean
      pPzs2LambdaFromC[cent][pt]->SetName(hNamePzs2LambdaFromC[cent][pt] + "_Profile");

      hV2C[cent][pt] = (TH1F *)hmass[cent][pt]->Clone(hNameV2C[cent][pt]);

      hPzs2[cent][pt] = (TH1F *)hmass[cent][pt]->Clone(hNamePzs2[cent][pt]);
      hPzs2LambdaFromC[cent][pt] = (TH1F *)hmass[cent][pt]->Clone(hNamePzs2LambdaFromC[cent][pt]);
      hCos2Theta[cent][pt] = (TH1F *)hmass[cent][pt]->Clone(hNameCos2Theta[cent][pt]);
      hCos2ThetaLambdaFromC[cent][pt] = (TH1F *)hmass[cent][pt]->Clone(hNameCos2ThetaLambdaFromC[cent][pt]);

      for (Int_t bin = 0; bin < hmass[cent][pt]->GetNbinsX(); bin++)
      {
        TH1D *htemp = (TH1D *)hmassVsV2C[cent][pt]->ProjectionX(Form("_htemp_%i", bin), bin + 1, bin + 1);
        hV2C[cent][pt]->SetBinContent(bin + 1, htemp->GetMean());
        hV2C[cent][pt]->SetBinError(bin + 1, htemp->GetMeanError());
        TH1D *htemp2 = (TH1D *)hmassVsPzs2[cent][pt]->ProjectionX(Form("_htemp2_%i", bin), bin + 1, bin + 1);
        hPzs2[cent][pt]->SetBinContent(bin + 1, htemp2->GetMean());
        hPzs2[cent][pt]->SetBinError(bin + 1, htemp2->GetMeanError());
        TH1D *htemp3 = (TH1D *)hmassVsPzs2LambdaFromC[cent][pt]->ProjectionX(Form("_htemp3_%i", bin), bin + 1, bin + 1);
        hPzs2LambdaFromC[cent][pt]->SetBinContent(bin + 1, htemp3->GetMean());
        hPzs2LambdaFromC[cent][pt]->SetBinError(bin + 1, htemp3->GetMeanError());
        TH1D *htemp4 = (TH1D *)hmassVsCos2Theta[cent][pt]->ProjectionX(Form("_htemp4_%i", bin), bin + 1, bin + 1);
        hCos2Theta[cent][pt]->SetBinContent(bin + 1, htemp4->GetMean());
        hCos2Theta[cent][pt]->SetBinError(bin + 1, htemp4->GetMeanError());
        TH1D *htemp5 = (TH1D *)hmassVsCos2ThetaLambdaFromC[cent][pt]->ProjectionX(Form("_htemp5_%i", bin), bin + 1, bin + 1);
        hCos2ThetaLambdaFromC[cent][pt]->SetBinContent(bin + 1, htemp5->GetMean());
        hCos2ThetaLambdaFromC[cent][pt]->SetBinError(bin + 1, htemp5->GetMeanError());
      }
    }
  }
  QCPhi->SaveAs("QCPlots/QCPhiCasc.png");

  TString SOutputFile = "OutputAnalysis/V2_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (isApplyWeights)
    SOutputFile += "_Weighted";
  if (v2type == 1)
    SOutputFile += "_SP";
  if (!useCommonBDTValue)
    SOutputFile += "_BDTCentDep";
  if (isRun2Binning)
    SOutputFile += "_Run2Binning";
  if (ExtrisApplyEffWeights)
  {
    SOutputFile += "_EffW";
  }
  SOutputFile += "_WithAlpha";
  if (!isRapiditySel)
    SOutputFile += "_Eta08";
  if (isReducedPtBins)
    SOutputFile += "_ReducedPtBins";
  SOutputFile += STHN[ExtrisFromTHN];
  if (ChosenPart == 6 && isSysMultTrial)
    SOutputFile += Form("_SysMultTrial_%i", indexMultTrial);
  SOutputFile += ".root";
  TFile *file = new TFile(SOutputFile, "RECREATE");
  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    cout << "hello" << cent << endl;
    if (isOOCentrality)
    {
      if (cent > numCentLambdaOO)
        continue;
    }
    hmassVsPtVsV2C[cent]->Write();
    hmassVsPtVsPzs2[cent]->Write();
    hmassVsPtVsPzs2LambdaFromC[cent]->Write();
    hmassVsPtVsCos2Theta[cent]->Write();
    hmassVsPtVsCos2ThetaLambdaFromC[cent]->Write();
    hmassVsPt[cent]->Write();
    if (!ExtrisFromTHN && ChosenPart != 6)
      hAvgPt[cent]->Write();
    for (Int_t pt = 0; pt < numPtBins + 1; pt++)
    {
      hmassVsV2C[cent][pt]->Write();
      hV2C[cent][pt]->Write();
      if (!ExtrisFromTHN && ChosenPart != 6)
        hV2CFromProfile[cent][pt]->Write();
      pV2C[cent][pt]->Write();
      hmassVsPzs2[cent][pt]->Write();
      hmassVsPzs2LambdaFromC[cent][pt]->Write();
      hPzs2[cent][pt]->Write();
      hPzs2LambdaFromC[cent][pt]->Write();
      pPzs2[cent][pt]->Write();
      pPzs2LambdaFromC[cent][pt]->Write();
      hmassVsCos2Theta[cent][pt]->Write();
      hmassVsCos2ThetaLambdaFromC[cent][pt]->Write();
      hCos2Theta[cent][pt]->Write();
      hCos2ThetaLambdaFromC[cent][pt]->Write();
      hmass[cent][pt]->Write();
      if (isOOCentrality)
      {
        if (cent != numCentLambdaOO && pt != numPtBins && !ExtrisFromTHN)
          hPhiCentHisto1D[cent][pt]->Write();
      }
      else if (cent != numCent && pt != numPtBins && !ExtrisFromTHN)
        hPhiCentHisto1D[cent][pt]->Write();
    }
    cout << "hello" << cent << endl;
    for (Int_t psi = 0; psi < numPsiBins; psi++)
    {
      //   hmassVsPz[cent][psi]->Write();
      // hmassVsPzLambdaFromC[cent][psi]->Write();
      // hPz[cent][psi]->Write();
      // hPzLambdaFromC[cent][psi]->Write();
      // hmassVsCos2ThetaPsi[cent][psi]->Write();
      // hmassVsCos2ThetaPsiLambdaFromC[cent][psi]->Write();
      // hCos2ThetaPsi[cent][psi]->Write();
      // hCos2ThetaPsiLambdaFromC[cent][psi]->Write();
      // pPz[cent][psi]->Write();
      // pPzLambdaFromC[cent][psi]->Write();
      // hmassPsi[cent][psi]->Write();
    }
  }
  file->Close();
  TString SweightsFile = "PhiWeights/Weights_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (!isRapiditySel || ExtrisFromTHN)
    SweightsFile += "_Eta08";
  SweightsFile += ".root";
  TFile *weightsFile;
  if (!isApplyWeights && !ExtrisFromTHN) // weights computed only once (when we apply weights, they have already been created!)
  {
    // weightsFile = new TFile(SweightsFile, "RECREATE");
    // weightsFile->cd();
    // weights->Write();
    // weightsFile->Close();
  }
  cout << "I created the file " << file->GetName() << endl;
  if (!isApplyWeights && !ExtrisFromTHN)
    cout << " and the file with weights: " << SweightsFile << endl;
}
