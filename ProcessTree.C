// This macro was originally written by:
// chiara.de.martin@cern.ch
#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
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

using namespace ROOT;

void ProcessTree(Int_t indexMultTrial = 0, Bool_t isXi = ChosenParticleXi, TString inputFileName = SinputFileName, Int_t EtaSysChoice = ExtrEtaSysChoice, Bool_t isSysMultTrial = 1)
{

  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trials)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trials * indexMultTrial;

  cout << "Input file: " << inputFileName << endl;
  cout << "isXi: " << isXi << endl;
  cout << "EtaSysChoice: " << EtaSysChoice << endl;
  cout << "BDTscoreCut: " << BDTscoreCut << endl;

  TString TreeName = "O2cascanalysis";

  TString inputFile = "TreeForAnalysis/AnalysisResults_trees_" + inputFileName + "_New.root";
  RDataFrame d1(TreeName, inputFile);
  auto h = d1.Histo1D("fPt");
  h->Draw();

  // invariant mass histograms
  auto hmass_Bef = d1.Histo1D({"mass_Xi_Bef", "Invariant mass of #Lambda#pi", 100, 1.29, 1.35}, "fMassXi");
  if (!isXi)
    auto hmass_Bef = d1.Histo1D({"mass_Omega_Bef", "Invariant mass of #LambdaK", 100, 1.6, 1.73}, "fMassOmega");
  // BDT response histogram
  auto BDT_response_Bef = d1.Histo1D({"BDT_response_Bef", "BDT response", 100, 0, 1}, "fBDTResponseXi");
  if (!isXi)
    BDT_response_Bef = d1.Histo1D({"BDT_response_Bef", "BDT response", 100, 0, 1}, "fBDTResponseOmega");
  auto mass_vs_BDTResponse = d1.Histo2D({"mass_vs_BDTResponse", "Invariant mass vs BDT response", 100, 0, 1, 100, 1.28, 1.36}, "fBDTResponseXi", "fMassXi");
  if (!isXi)
    mass_vs_BDTResponse = d1.Histo2D({"mass_vs_BDTResponse", "Invariant mass vs BDT response", 100, 0, 1, 100, 1.6, 1.73}, "fBDTResponseOmega", "fMassOmega");

  // apply BDT selection
  string expression = Form("fBDTResponseXi > %.3f", BDTscoreCut);
  if (!isXi)
    expression = Form("fBDTResponseOmega > %.3f", BDTscoreCut);
  cout << "expression: " << expression << endl;
  auto d2 = d1.Filter(expression);
  
  // apply eta selection for systematic studies
  auto d3 = d2;
  if (EtaSysChoice == 0) // -0.8 < eta < 0.8
    d3 = d2.Filter("fEta > -0.8 && fEta < 0.8");
  else if (EtaSysChoice == 1) // eta > 0 && eta < 0.8
    d3 = d2.Filter("fEta > 0 && fEta < 0.8");
  else if (EtaSysChoice == 2) // eta < 0 && eta > -0.8
    d3 = d2.Filter("fEta < 0 && fEta > -0.8");

  auto BDT_response = d3.Histo1D({"BDT_response", "BDT response", 100, 0, 1}, "fBDTResponseXi");
  if (!isXi)
    BDT_response = d3.Histo1D({"BDT_response", "BDT response", 100, 0, 1}, "fBDTResponseOmega");

  // invariant mass histograms
  auto hmass = d3.Histo1D({"mass_Xi", "Invariant mass of #Lambda#pi", 100, 1.29, 1.35}, "fMassXi");
  if (!isXi)
    hmass = d3.Histo1D({"mass_Omega", "Invariant mass of #LambdaK", 100, 1.6, 1.73}, "fMassOmega");

  // eta distributions
  auto heta = d3.Histo1D({"eta", "Eta distribution of selected candidates", 200, -2, 2}, "fEta");

  // phi distributions
  // auto hphi = d3.Histo1D({"phi", "Phi distribution of selected candidates", 200, -3.2, 3.2}, "fPhi");

  // create output file
  Int_t ParticleIndex = 0; // 0 for Xi, 1 for Omega
  if (isXi)
    ParticleIndex = 0;
  else
    ParticleIndex = 1;

  TString SBDT = "";
  if (BDTscoreCut != DefaultBDTscoreCut)
    SBDT = Form("_BDT%.3f", BDTscoreCut);
  TString OutputFileName = "OutputAnalysis/Output_" + inputFileName + "_" + ParticleName[ParticleIndex] + SEtaSysChoice[EtaSysChoice] + SBDT + ".root";
  TFile *file = new TFile(OutputFileName, "RECREATE");

  // 3D histograms

  for (Int_t cent = 0; cent < numCent; cent++)
  {
    auto dcent = d3.Filter(Form("fCentFT0C>=%i && fCentFT0C<%i", CentFT0C[cent], CentFT0C[cent + 1]));
    auto v2C = dcent.Histo1D({Form("v2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "v2C", 240, -1.2, 1.2}, "fV2C");
    v2C->Write();
    if (isXi)
    {
      auto massVsPtVsV2C = dcent.Histo3D({Form("massVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "Invariant mass vs Pt vs V2C", 100, 1.28, 1.36, 100, 0, 10, 200, -1., 1.}, "fMassXi", "fPt", "fV2C");
      massVsPtVsV2C->Write();
      // profile: mean value of v2 vs mass and pt in centrality classes
      auto profile = dcent.Profile2D({Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "Mean invariant mass vs Pt vs V2C", 50, 1.28, 1.36, numPtBins, PtBins}, "fMassXi", "fPt", "fV2C");
      profile->Write();
    }
    else
    {
      auto massVsPtVsV2C = dcent.Histo3D({Form("massVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "Invariant mass vs Pt vs V2C", 100, 1.6, 1.73, 100, 0, 10, 200, -1., 1.}, "fMassOmega", "fPt", "fV2C");
      massVsPtVsV2C->Write();
      // profile: mean value of v2 vs mass and pt in centrality classes
      auto profile = dcent.Profile2D({Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "Mean invariant mass vs Pt vs V2C", 50, 1.6, 1.73, 100, 0, 10, -2, 2}, "fMassOmega", "fPt", "fV2C");
      profile->Write();
    }
  }

  // draw histograms
  TCanvas *cMass = new TCanvas("cMass", "cMass", 900, 600);
  TString TitleX = "M_{#Lambda#pi}";
  if (!isXi)
    TitleX = "M_{#LambdaK}";
  StyleCanvas(cMass, 0.1, 0.1, 0.03, 0.1);
  StyleHisto(*hmass_Bef, 0, 1.2 * hmass_Bef->GetMaximum(), kRed, 20, "M_{#Lambda#pi}", "Counts", "", kTRUE, 1.2, 1.4, 1.2, 1.2, 0.7);
  StyleHisto(*hmass, 0, 1.2 * hmass_Bef->GetMaximum(), kBlue, 20, "M_{#Lambda#pi}", "Counts", "", kTRUE, 1.2, 1.4, 1.2, 1.2, 0.7);
  hmass_Bef->Draw("E");
  hmass->Draw("E SAME");

  cMass->Write();
  hmass_Bef->Write();
  BDT_response_Bef->Write();
  mass_vs_BDTResponse->Write();
  hmass->Write();
  heta->Write();
  BDT_response->Write();
  file->Close();

  cout << "I created the file " << file->GetName() << endl;
}
