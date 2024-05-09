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
#include "StyleFile.h"

void ComputeV2(Int_t indexMultTrial = 0, Bool_t isXi = ChosenParticleXi, TString inputFileName = SinputFileName, Int_t RebinFactor = 2, Int_t EtaSysChoice = ExtrEtaSysChoice, Bool_t isSysMultTrial = ExtrisSysMultTrial)
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

  TString SinputFile = "OutputAnalysis/Output_" + inputFileName + "_" + ParticleName[!isXi] + SEtaSysChoice[EtaSysChoice] + SBDT + ".root";
  cout << "Input file: " << SinputFile << endl;
  TFile *inputFile = new TFile(SinputFile);
  TH3D *hmassVsPtVsV2C[numCent];
  TProfile2D *profmassVsPt[numCent];
  TH2F *hmassVsPt[numCent];
  TH2F *hmassVsV2C[numCent][numPtBins];
  TH1F *hmass[numCent][numPtBins];
  TH1F *hV2C[numCent][numPtBins];
  TH1F *hV2CFromProfile[numCent][numPtBins];
  TProfile *pV2C[numCent][numPtBins];
  TString hName = "";
  TString profName = "";
  TString hNameMass = "";
  TString hNameV2C = "";
  TString hNameMassV2C = "";
  TString hNameV2CFromProfile2D = "";
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    hName = Form("massVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]);
    profName = Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]);
    hmassVsPtVsV2C[cent] = (TH3D *)inputFile->Get(hName);
    profmassVsPt[cent] = (TProfile2D *)inputFile->Get(profName);
    if (!profmassVsPt[cent]){
      cout << "TProfile2D not found" << endl;
      return;
    }
    hmassVsPt[cent] = (TH2F *)hmassVsPtVsV2C[cent]->Project3D("yx");
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      hNameMass = Form("mass_cent%i-%i_pt%i", CentFT0C[cent], CentFT0C[cent + 1], pt);
      hNameV2C = Form("V2C_cent%i-%i_pt%i", CentFT0C[cent], CentFT0C[cent + 1], pt);
      hNameV2CFromProfile2D = Form("V2CFromProfile_cent%i-%i_pt%i", CentFT0C[cent], CentFT0C[cent + 1], pt);
      hNameMassV2C = Form("MassvsV2C_cent%i-%i_pt%i", CentFT0C[cent], CentFT0C[cent + 1], pt);

      hmassVsPtVsV2C[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);

      profmassVsPt[cent]->GetYaxis()->SetRangeUser(PtBins[pt] + 0.0001, PtBins[pt + 1] - 0.0001);
      hV2CFromProfile[cent][pt] = (TH1F *)profmassVsPt[cent]->ProjectionX(hNameV2CFromProfile2D, 0, -1, "e"); // v2C from TProfile2D

      hmassVsV2C[cent][pt] = (TH2F *)hmassVsPtVsV2C[cent]->Project3D("xze"); // mass vs V2C //"e" option does not change the results
      hmassVsV2C[cent][pt]->SetName(hNameMassV2C);
      hmassVsV2C[cent][pt]->RebinY(RebinFactor);

      hmass[cent][pt] = (TH1F *)hmassVsPtVsV2C[cent]->Project3D("xe"); // mass
      hmass[cent][pt]->SetName(hNameMass);
      hmass[cent][pt]->Rebin(RebinFactor);

      pV2C[cent][pt] = hmassVsV2C[cent][pt]->ProfileY(); // v2C //the error is the standard error of the mean
      pV2C[cent][pt]->SetName(hNameV2C + "_Profile");

      hV2C[cent][pt] = (TH1F *)hmass[cent][pt]->Clone(hNameV2C);
      for (Int_t bin = 0; bin < hmass[cent][pt]->GetNbinsX(); bin++)
      {
        hV2C[cent][pt]->SetBinContent(bin + 1, hmassVsV2C[cent][pt]->ProjectionX("", bin + 1, bin + 1)->GetMean());
        hV2C[cent][pt]->SetBinError(bin + 1, hmassVsV2C[cent][pt]->ProjectionX("", bin + 1, bin + 1)->GetMeanError());
      }
    }
  }

  TFile *file = new TFile("OutputAnalysis/V2_" + inputFileName + "_" + ParticleName[!isXi] + SEtaSysChoice[EtaSysChoice] + SBDT + ".root", "RECREATE");
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      hmass[cent][pt]->Write();
      hV2C[cent][pt]->Write();
      hmassVsV2C[cent][pt]->Write();
      pV2C[cent][pt]->Write();
      hV2CFromProfile[cent][pt]->Write();
    }
    hmassVsPtVsV2C[cent]->Write();
    hmassVsPt[cent]->Write();
  }
  file->Close();
  cout << "I created the file " << file->GetName() << endl;
}
