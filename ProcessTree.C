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
#include "TRandom3.h"
#include <ROOT/RDataFrame.hxx>


using namespace ROOT;
using namespace std;

constexpr double massSigmaParameters[4][2]{
    {4.9736e-3, 0.006815},
    {-2.39594, -2.257},
    {1.8064e-3, 0.00138},
    {1.03468e-1, 0.1898}};

void ProcessTree(Int_t indexMultTrial = 0, Bool_t isXi = ChosenParticleXi, TString inputFileName = SinputFileName, Int_t EtaSysChoice = ExtrEtaSysChoice, Bool_t isSysMultTrial = ExtrisSysMultTrial, Bool_t isApplyWeights = 1)
{

  ROOT::EnableImplicitMT();

  if (isSysMultTrial)
    inputFileName = SinputFileNameSyst;
  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trialsBDT)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;

  // weights to flatten phi distribution of cascades
  //TFile *fileWeights = new TFile(SfileWeights);

  cout << "Input file: " << inputFileName << endl;
  //if (isApplyWeights)
    //cout << "Weights applied from file " << SfileWeights << endl;
  cout << "isXi: " << isXi << endl;
  cout << "EtaSysChoice: " << EtaSysChoice << endl;
  cout << "BDTscoreCut: " << BDTscoreCut << endl;

  TString TreeName = "O2cascanalysis";

  TString inputFile = "TreeForAnalysis/AnalysisResults_trees_" + inputFileName + "_New.root";

  RDataFrame originalDF(TreeName, inputFile);

  std::vector<TRandom3> randoms;
  for (Int_t i = 0; i < originalDF.GetNSlots(); i++)
  {
    randoms.push_back(TRandom3(i + 1));
  }

  TFile *weightFile = new TFile(weightFileName, "READ");
  TH3D *weights{weightFileName ? (TH3D *)weightFile->Get("weights") : nullptr};

  auto d1 = originalDF.DefineSlot("random", [&randoms](unsigned int slot) -> float { return randoms[slot].Rndm(); }).Filter([&weights](float random, float cent, float pt, float phi) {
    if (!weights) return true;
    int centBin = weights->GetXaxis()->FindBin(cent);
    int ptBin = weights->GetYaxis()->FindBin(pt);
    int phiBin = weights->GetZaxis()->FindBin(phi);
    return random < weights->GetBinContent(centBin, ptBin, phiBin);
   }, {"random", "fCentFT0C", "fPt", "fPhi"});

  auto h = d1.Histo1D("fPt");

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

  // apply competing mass rejection for Omega
  // if (!isXi) d3 = d3.Filter("fMassXi > 1.34 || fMassXi < 1.3"); //rough mass cut to reject Xi
  string rejectMassXi = Form("abs(fMassXi - %.3f) > 5* (%.3f * exp(%.3f * fPt) + %.3f * exp(%.3f * fPt))", ParticleMassPDG[!isXi], massSigmaParameters[0][0], massSigmaParameters[1][0], massSigmaParameters[2][0], massSigmaParameters[3][0]);
  if (!isXi)
    d3 = d3.Filter(rejectMassXi);
  auto MassXi = d3.Histo2D({"mass_XivsPt", "Invariant mass of #Lambda#pi", 100, 1.29, 1.35, 100, 0, 10}, "fMassXi", "fPt");

  auto BDT_response = d3.Histo1D({"BDT_response", "BDT response", 100, 0, 1}, "fBDTResponseXi");
  if (!isXi)
    BDT_response = d3.Histo1D({"BDT_response", "BDT response", 100, 0, 1}, "fBDTResponseOmega");

  // invariant mass histograms
  auto hmass = d3.Histo1D({"mass_Xi", "Invariant mass of #Lambda#pi", 100, 1.29, 1.35}, "fMassXi");
  if (!isXi)
    hmass = d3.Histo1D({"mass_Omega", "Invariant mass of #LambdaK", 100, 1.6, 1.73}, "fMassOmega");

  // invariant mass histograms vs pt
  auto hmassvsPt = d3.Histo2D({"mass_XivsPt", "Invariant mass of #Lambda#pi", 100, 1.29, 1.35, 100, 0, 10}, "fMassXi", "fPt");
  if (!isXi)
    hmassvsPt = d3.Histo2D({"mass_OmegavsPt", "Invariant mass of #LambdaK", 100, 1.6, 1.73, 100, 0, 10}, "fMassOmega", "fPt");

  // eta distributions
  auto heta = d3.Histo1D({"eta", "Eta distribution of selected candidates", 200, -2, 2}, "fEta");

  // phi distributions
  auto hphi = d3.Histo1D({"phi", "Phi distribution of selected candidates", 200, 0, 2 * TMath::Pi()}, "fPhi");

  // eta - phi distributions
  auto hEtaPhi = d3.Histo2D({"PhivsEta", "Phi vs Eta distribution of selected candidates", 100, -1, 1, 200, 0, 2 * TMath::Pi()}, "fEta", "fPhi");

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
  TH1D *MassCutHisto[numCent];
  TH1D *v2CHisto[numCent];
  TH2D *hPhiCentHisto[numCent];
  TH3D *massVsPtVsV2CHisto[numCent];
  TProfile *profileHisto[numCent];
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    auto dcent = d3.Filter(Form("fCentFT0C>=%i && fCentFT0C<%i", CentFT0C[cent], CentFT0C[cent + 1]));
    string MassCut = "";
    if (isXi)
      MassCut = Form("abs(fMassXi - %.3f) < 3* (%.3f * exp(%.3f * fPt) + %.3f * exp(%.3f * fPt))", ParticleMassPDG[!isXi], massSigmaParameters[0][0], massSigmaParameters[1][0], massSigmaParameters[2][0], massSigmaParameters[3][0]);
    else
      MassCut = Form("abs(fMassOmega - %.3f) < 3* (%.3f * exp(%.3f * fPt) + %.3f * exp(%.3f * fPt))", ParticleMassPDG[!isXi], massSigmaParameters[0][1], massSigmaParameters[1][1], massSigmaParameters[2][1], massSigmaParameters[3][1]);
    auto dmasscut = dcent.Filter(MassCut);
    auto v2C = dcent.Histo1D({Form("v2CHist_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "v2C", 240, -1.2, 1.2}, "fV2C");
    auto hPhiCent = dmasscut.Histo2D({Form("PhiHist_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "Phi vs pt", 100, 0, 10, 100, 0, 2 * TMath::Pi()}, "fPt", "fPhi");
    if (isXi)
    {
      auto hMassCut = dmasscut.Histo1D({"massCut", "Invariant mass of #Lambda#pi", 100, 1.28, 1.36}, "fMassXi");
      auto massVsPtVsV2C = dcent.Histo3D({Form("massVsPtVsV2CHist_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "Invariant mass vs Pt vs V2C", 100, 1.28, 1.36, 100, 0, 10, 200, -1., 1.}, "fMassXi", "fPt", "fV2C");
      // profile: mean value of v2 vs mass and pt in centrality classes
      auto profile = dcent.Profile2D({Form("ProfilemassVsPtVsV2CHist_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "Mean invariant mass vs Pt vs V2C", 50, 1.28, 1.36, numPtBins, PtBins}, "fMassXi", "fPt", "fV2C");
      MassCutHisto[cent] = (TH1D *)hMassCut->Clone(Form("massCut_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
      massVsPtVsV2CHisto[cent] = (TH3D *)massVsPtVsV2C->Clone(Form("massVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
      profileHisto[cent] = (TProfile *)profile->Clone(Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
    }
    else
    {
      auto hMassCut = dmasscut.Histo1D({"massCut", "Invariant mass of #LambdaK", 100, 1.6, 1.73}, "fMassOmega");
      auto massVsPtVsV2C = dcent.Histo3D({Form("massVsPtVsV2CHist_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "Invariant mass vs Pt vs V2C", 100, 1.6, 1.73, 100, 0, 10, 200, -1., 1.}, "fMassOmega", "fPt", "fV2C");
      // profile: mean value of v2 vs mass and pt in centrality classes
      auto profile = dcent.Profile2D({Form("ProfilemassVsPtVsV2CHist_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]), "Mean invariant mass vs Pt vs V2C", 50, 1.6, 1.73, 100, 0, 10, -2, 2}, "fMassOmega", "fPt", "fV2C");
      MassCutHisto[cent] = (TH1D *)hMassCut->Clone(Form("massCut_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
      massVsPtVsV2CHisto[cent] = (TH3D *)massVsPtVsV2C->Clone(Form("massVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
      profileHisto[cent] = (TProfile *)profile->Clone(Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
    }
    v2CHisto[cent] = (TH1D *)v2C->Clone(Form("v2C_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
    hPhiCentHisto[cent] = (TH2D *)hPhiCent->Clone(Form("PhiHist_cent%i-%i", CentFT0C[cent], CentFT0C[cent + 1]));
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

  h->Write();
  cMass->Write();
  MassXi->Write();
  hmass_Bef->Write();
  BDT_response_Bef->Write();
  mass_vs_BDTResponse->Write();
  hmass->Write();
  hmassvsPt->Write();
  heta->Write();
  hphi->Write();
  hEtaPhi->Write();
  BDT_response->Write();
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    MassCutHisto[cent]->Write();
    v2CHisto[cent]->Write();
    hPhiCentHisto[cent]->Write();
    massVsPtVsV2CHisto[cent]->Write();
    profileHisto[cent]->Write();
  }
  file->Close();

  cout << "I created the file " << file->GetName() << endl;
}
