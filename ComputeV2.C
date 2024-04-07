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

void ComputeV2(Bool_t isXi = ChosenParticleXi, TString inputFileName = SinputFileName, Int_t RebinFactor = 2)
{

  TString SinputFile = "OutputAnalysis/Output_" + inputFileName + "_" + ParticleName[!isXi] + ".root";
  cout << "Input file: " << SinputFile << endl;
  TFile *inputFile = new TFile(SinputFile);
  TH3D *hmassVsPtVsV2C[numCent];
  TH2F *hmassVsPt[numCent];
  TH2F *hmassVsV2C[numCent][numPtBins];
  TH1F *hmass[numCent][numPtBins];
  TH1F *hV2C[numCent][numPtBins];
  TProfile *pV2C[numCent][numPtBins];
  TString hName = "";
  TString hNameMass = "";
  TString hNameV2C = "";
  TString hNameMassV2C = "";
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    hName = Form("massVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]);
    hmassVsPtVsV2C[cent] = (TH3D *)inputFile->Get(hName);
    hmassVsPt[cent] = (TH2F *)hmassVsPtVsV2C[cent]->Project3D("yx");
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      hNameMass = Form("mass_cent%i-%i_pt%i", CentFT0C[cent], CentFT0C[cent + 1], pt);
      hNameV2C = Form("V2C_cent%i-%i_pt%i", CentFT0C[cent], CentFT0C[cent + 1], pt);
      hNameMassV2C = Form("MassvsV2C_cent%i-%i_pt%i", CentFT0C[cent], CentFT0C[cent + 1], pt);

      hmassVsPtVsV2C[cent]->GetYaxis()->SetRangeUser(PtBins[pt], PtBins[pt + 1]);

      hmassVsV2C[cent][pt] = (TH2F *)hmassVsPtVsV2C[cent]->Project3D("xz");
      hmassVsV2C[cent][pt]->SetName(hNameMassV2C);
      hmassVsV2C[cent][pt]->RebinY(RebinFactor);

      hmass[cent][pt] = (TH1F *)hmassVsPtVsV2C[cent]->Project3D("x");
      hmass[cent][pt]->SetName(hNameMass);
      hmass[cent][pt]->Rebin(RebinFactor);

      hV2C[cent][pt] = (TH1F *)hmass[cent][pt]->Clone(hNameV2C);
      
      pV2C[cent][pt] = hmassVsV2C[cent][pt]->ProfileY();
      pV2C[cent][pt]->SetName(hNameV2C + "_Profile");

      for (Int_t bin = 0; bin < hmass[cent][pt]->GetNbinsX(); bin++)
      {
        hV2C[cent][pt]->SetBinContent(bin+1, hmassVsV2C[cent][pt]->ProjectionX("", bin+1, bin+2)->GetMean());
        hV2C[cent][pt]->SetBinError(bin+1, hmassVsV2C[cent][pt]->ProjectionX("", bin+1, bin+2)->GetMeanError());
      }
    }
  }

  TFile *file = new TFile("OutputAnalysis/V2_" + inputFileName + "_" + ParticleName[!isXi] + ".root", "RECREATE");
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      hmass[cent][pt]->Write();
      hV2C[cent][pt]->Write();
      hmassVsV2C[cent][pt]->Write();
      pV2C[cent][pt]->Write();
    }
    hmassVsPtVsV2C[cent]->Write();
    hmassVsPt[cent]->Write();
  }
  file->Close();
  cout << "I created the file " << file->GetName() << endl;
}
