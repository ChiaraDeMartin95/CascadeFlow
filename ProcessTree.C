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

Float_t Minv2 = -1;
Float_t Maxv2 = 1;
Int_t Nv2 = 200;

Float_t MinPzs2 = -1;
Float_t MaxPzs2 = 1;
Int_t NPzs2 = 200;

Float_t MinPz = -1;
Float_t MaxPz = 1;
Int_t NPz = 200;

bool passbdtCut(float bdtscore, float cent)
{
  float sbdtCut = 0;
  if (useCommonBDTValue)
    sbdtCut = DefaultBDTscoreCut;
  else
    sbdtCut = bdtCut[int(cent / 10)];
  return (bdtscore > sbdtCut);
}

void ProcessTree(Bool_t isEff = 0,
                 Int_t indexMultTrial = 0,
                 Int_t ChosenPart = ChosenParticle,
                 TString inputFileName = SinputFileName,
                 Int_t EtaSysChoice = ExtrEtaSysChoice,
                 Bool_t isSysMultTrial = ExtrisSysMultTrial)
{
  
  Bool_t isApplyEffWeights = 1;
  if (isEff)
    isApplyEffWeights = 0; // just compute histos for efficiency, no use in applying weigths
  
  Bool_t isXi = 1;
  if ((ChosenPart == 1) || (ChosenPart == 4) || (ChosenPart == 5))
    isXi = 0; // Omega
  if (isEff)
    inputFileName = SinputFileNameEff;
  string v2Chosen = "fV2C";
  if (v2type == 1)
  {
    v2Chosen = "fV2CSP";
    Minv2 = -5;
    Maxv2 = 5;
  }
  else if (v2type == 2)
    v2Chosen = "fV2CEP";

  ROOT::EnableImplicitMT();

  string cosThetaStar = "fCosThetaStarLambdaFromXi";
  if (!isXi)
    cosThetaStar = "fCosThetaStarLambdaFromOmega";

  if (isSysMultTrial)
    inputFileName = SinputFileNameSyst;
  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trialsBDT)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;

  if (isApplyWeights)
    cout << "Weights applied from file " << weightFileName << endl;
  cout << "ChosenPart: " << ParticleName[ChosenPart] << endl;
  cout << "EtaSysChoice: " << EtaSysChoice << endl;
  cout << "Use common BDT value " << useCommonBDTValue << endl;

  TString TreeName = "O2cascanalysis";

  TString inputFile = "TreeForAnalysis";
  if (isEff)
    inputFile = "FileForEfficiency";
  inputFile += "/AnalysisResults_trees_" + inputFileName + "_New.root";

  RDataFrame originalDF(TreeName, inputFile);

  std::vector<TRandom3> randoms;
  for (int i = 0; i < (int)originalDF.GetNSlots(); i++)
  {
    randoms.push_back(TRandom3(i + 1));
  }

  TFile *weightFile = new TFile(weightFileName, "READ");
  TH3D *weights{weightFileName ? (TH3D *)weightFile->Get("weights") : nullptr};

  auto d1 = originalDF.DefineSlot("random", [&randoms](unsigned int slot) -> float
                                  { return randoms[slot].Rndm(); })
                .Filter([&weights](float random, float cent, float pt, float phi)
                        {
                          if (!weights) return true;
                          int centBin = weights->GetXaxis()->FindBin(cent);
                          int ptBin = weights->GetYaxis()->FindBin(pt);
                          int phiBin = weights->GetZaxis()->FindBin(phi);
                          if (isApplyWeights) return random < 1./weights->GetBinContent(centBin, ptBin, phiBin); 
                          else return random < 999; }, {"random", "fCentFT0C", "fPt", "fPhi"});

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

  // apply charge selection
  string chargecut = "abs(fSign) > 0";
  if (ChosenParticle == 3 || ChosenParticle == 5)
    chargecut = "fSign > 0";
  else if (ChosenParticle == 2 || ChosenParticle == 4)
    chargecut = "fSign < 0";
  auto d2 = d1.Filter(chargecut);

  // apply eta selection for systematic studies
  // auto d = d2;
  if (EtaSysChoice == 0) // -0.8 < eta < 0.8
    d2 = d2.Filter("fEta > -0.8 && fEta < 0.8");
  else if (EtaSysChoice == 1) // eta > 0 && eta < 0.8
    d2 = d2.Filter("fEta > 0 && fEta < 0.8");
  else if (EtaSysChoice == 2) // eta < 0 && eta > -0.8
    d2 = d2.Filter("fEta < 0 && fEta > -0.8");

  // apply competing mass rejection for Omega
  //   if (!isXi) d3 = d3.Filter("fMassXi > 1.34 || fMassXi < 1.3"); //rough mass cut to reject Xi
  string rejectMassXi = Form("abs(fMassXi - %.3f) > 5* (%.3f * exp(%.3f * fPt) + %.3f * exp(%.3f * fPt))", ParticleMassPDG[ChosenPart], massSigmaParameters[0][0], massSigmaParameters[1][0], massSigmaParameters[2][0], massSigmaParameters[3][0]);
  if (!isXi)
    d2 = d2.Filter(rejectMassXi);

  // apply PDG code selections in case it's MC
  if (isEff)
  {
    if (isXi)
      d2 = d2.Filter("abs(fMcPdgCode) == 3312");
    else
      d2 = d2.Filter("abs(fMcPdgCode) == 3334");
  }

  // define the rapidity
  string Common1 = Form("std::sqrt(fPt*fPt*std::cosh(fEta)*std::cosh(fEta) + %.3f*%.3f)", ParticleMassPDG[ChosenPart], ParticleMassPDG[ChosenPart]);
  string Common2 = "fPt*std::sinh(fEta)";
  string Num = Form("%s + %s", Common1.c_str(), Common2.c_str());
  string Denom = Form("%s - %s", Common1.c_str(), Common2.c_str());
  d2 = d2.Define("fRapidity", Form("0.5 * std::log((%s) / (%s))", Num.c_str(), Denom.c_str()));
  auto d2RapCut = d2.Filter("abs(fRapidity) < 0.5");
  cout << Form("0.5 * std::log((%s) / (%s))", Num.c_str(), Denom.c_str()) << endl;

  // pt vs centrality before BDT cut
  auto hPtvsCent_Bef = d2.Histo2D({"PtvsCent_BefBDT", "PtvsCent_BefBDT", 100, 0, 100, 200, 0, 20}, "fCentFT0C", "fPt");
  auto hPtvsCent_RapCut_Bef = d2RapCut.Histo2D({"PtvsCent_Y05_BefBDT", "PtvsCent_Y05_BefBDT", 100, 0, 100, 200, 0, 20}, "fCentFT0C", "fPt");

  //  apply BDT selection
  string cutvariable = "fBDTResponseXi";
  if (!isXi)
    cutvariable = "fBDTResponseOmega";
  auto d3 = d2.Filter(passbdtCut, {cutvariable, "fCentFT0C"});
  auto d3RapCut = d3.Filter(passbdtCut, {cutvariable, "fCentFT0C"});

  // pt vs centrality after BDT cut
  auto hPtvsCent_Aft = d3.Histo2D({"PtvsCent_AftBDT", "PtvsCent_AftBDT", 100, 0, 100, 200, 0, 20}, "fCentFT0C", "fPt");
  auto hPtvsCent_RapCut_Aft = d3RapCut.Histo2D({"PtvsCent_Y05_AftBDT", "PtvsCent_Y05_AftBDT", 100, 0, 100, 200, 0, 20}, "fCentFT0C", "fPt");

  auto MassXi = d3.Histo2D({"mass_XivsPt", "Invariant mass of #Lambda#pi", 100, 1.29, 1.35, 100, 0, 10}, "fMassXi", "fPt");

  auto BDT_response = d3.Histo2D({"BDT_response", "BDT response", 100, 0, 1, 100, 0, 100}, "fBDTResponseXi", "fCentFT0C");
  if (!isXi)
    BDT_response = d3.Histo2D({"BDT_response", "BDT response", 100, 0, 1, 100, 0, 100}, "fBDTResponseOmega", "fCentFT0C");

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

  TString SBDT = "";
  if (BDTscoreCut != DefaultBDTscoreCut)
    SBDT = Form("_BDT%.3f", BDTscoreCut);
  TString SEffWeights = "Efficiency/Efficiency_" + SinputFileNameEff + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (!useCommonBDTValue)
    SEffWeights += "_BDTCentDep";
  if (isRun2Binning)
    SEffWeights += "_Run2Binning";
  SEffWeights += ".root";
  TFile *fileEffWeights = new TFile(SEffWeights, "");
  if (isApplyEffWeights)
    cout << fileEffWeights->GetName() << endl;
  if (!fileEffWeights && isApplyEffWeights)
  {
    cout << "Efficiency file not found" << endl;
    return;
  }
  TF1 *lEffWeights[numPtBinsEff + 1];
  if (isApplyEffWeights)
  {
    for (Int_t i = 0; i < numPtBinsEff; i++)
    {
      lEffWeights[i] = (TF1 *)fileEffWeights->Get(Form("fitEffVsNch_%i", i));
      if (!lEffWeights[i])
      {
        cout << "Efficiency function not found" << endl;
        return;
      }
    }
  }

  // create output file
  TString OutputFileName = "OutputAnalysis/Output_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (isApplyWeights)
    OutputFileName += "_Weighted";
  if (v2type == 1)
    OutputFileName += "_SP";
  if (!useCommonBDTValue)
    OutputFileName += "_BDTCentDep";
  if (isRun2Binning)
    OutputFileName += "_Run2Binning";
  // OutputFileName += "_Effweight.root";
  OutputFileName += ".root";
  TFile *file = new TFile(OutputFileName, "RECREATE");
  cout << file->GetName() << endl;

  // 3D histograms
  TH1D *MassCutHisto[numCent + 1];
  TH1D *v2CHisto[numCent + 1];
  TH2D *hPhiCentHisto[numCent + 1];
  TH1D *hPsiCentHisto[numCent + 1];
  TH3D *massVsPtVsV2CHisto[numCent + 1];
  TH3D *massVsPtVsV2CWeightedHisto[numCent + 1];
  TH3D *massVsPtVsPzs2Histo[numCent + 1];
  TH3D *massVsPsiVsPzHisto[numCent + 1];
  TProfile *profileHisto[numCent + 1];
  TH1F *rapidityHisto[numCent + 1];
  TH1F *hEtaHisto[numCent + 1];
  TH1D *hPsiDiffHisto[numCent + 1];
  TH1D *hPsiDiff2Histo[numCent + 1];
  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
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
    auto dcent = d3.Filter(Form("fCentFT0C>=%.1f && fCentFT0C<%.1f", CentFT0CMin + 0.1, CentFT0CMax - 0.1));
    string MassCut = "";
    if (isXi)
      MassCut = Form("abs(fMassXi - %.3f) < 3* (%.3f * exp(%.3f * fPt) + %.3f * exp(%.3f * fPt))", ParticleMassPDG[ChosenPart], massSigmaParameters[0][0], massSigmaParameters[1][0], massSigmaParameters[2][0], massSigmaParameters[3][0]);
    else
      MassCut = Form("abs(fMassOmega - %.3f) < 3* (%.3f * exp(%.3f * fPt) + %.3f * exp(%.3f * fPt))", ParticleMassPDG[ChosenPart], massSigmaParameters[0][1], massSigmaParameters[1][1], massSigmaParameters[2][1], massSigmaParameters[3][1]);
    auto dmasscut = dcent.Filter(MassCut);

    auto v2C = dcent.Histo1D({Form("v2CHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "v2C", Nv2, Minv2, Maxv2}, v2Chosen);
    auto hPhiCent = dmasscut.Histo2D({Form("PhiHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Phi vs pt", 100, 0, 10, 100, 0, 2 * TMath::Pi()}, "fPt", "fPhi");
    auto hPsiCent = dmasscut.Histo1D({Form("PsiHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Psi", 100, -2 * TMath::Pi(), 2 * TMath::Pi()}, "fPsiT0C");
    dcent = dcent.Define("fPsiDiff", "if ((fPhi-fPsiT0C) < 0) return (fPhi-fPsiT0C+(float)TMath::Pi()); else if ((fPhi-fPsiT0C) > 2* TMath::Pi()) return (fPhi-fPsiT0C-2*(float)TMath::Pi()); else if ((fPhi-fPsiT0C) > TMath::Pi()) return (fPhi-fPsiT0C-(float)TMath::Pi()); else return (fPhi-fPsiT0C);");
    dcent = dcent.Define("f2PsiDiffCorr", "2*fPsiDiff");
    dcent = dcent.Define("f2PsiDiff", "2*fPhi-2*fPsiT0C");
    dcent = dcent.Define("Nch", Form("%.2f", dNdEtaAbhi[cent]));
    dcent = dcent.Define("v2Pub", Form("%.3f", v2PubRun2[cent]));
    auto dcentPzs2 = dcent.Filter("fRapidity > -0.5 && fRapidity < 0.5");
    // double numW = lEffWeights[0]->Eval(dNdEtaAbhi[cent]*(1+2*v2[cent]*cos(f2PsiDiffCorr)));
    // double denW = lEffWeights[0]->Eval(dNdEtaAbhi[cent]);
    float par1 = 2;
    float par2 = 1;
    // dcent = dcent.Define("numW", [&](double Nch, double v2, float y)
    //{ return par2 * std::exp(par1 * Nch * (1 + 2 * v2 * cos(y))); }, {"Nch", "v2Pub", "f2PsiDiffCorr"});
    dcent = dcent.Define("numW", [&weights](double Nch, double v2, float y, float cent, float pt)
                         { 
                                  int centBin = weights->GetXaxis()->FindBin(cent);
                                  int ptBin = weights->GetYaxis()->FindBin(pt);
                                  return std::exp(weights->GetBinContent(centBin, ptBin) * Nch * (1 + 2 * v2 * cos(y))); }, {"Nch", "v2Pub", "f2PsiDiffCorr", "fCentFT0C", "fPt"});

    dcent = dcent.Define("denW", [&](double Nch)
                         { return par2 * std::exp(par1 * Nch); }, {"Nch"});
    dcent = dcent.Define("fEffWeight", "numW/denW");
    if (isXi)
    {
      dcentPzs2 = dcentPzs2.Define("Pzs2Xi", "fCosThetaStarLambdaFromXi * sin(2*(fPhi-fPsiT0C))");
      auto hMassCut = dmasscut.Histo1D({"massCut", "Invariant mass of #Lambda#pi", 100, 1.28, 1.36}, "fMassXi");
      auto massVsPtVsV2C = dcent.Histo3D({Form("massVsPtVsV2CHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs V2C", 80, 1.28, 1.36, 100, 0, 10, Nv2, Minv2, Maxv2}, "fMassXi", "fPt", v2Chosen);
      auto massVsPtVsV2CWeighted = dcent.Histo3D({Form("massVsPtVsV2CHistWeighted_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs V2C", 80, 1.28, 1.36, 100, 0, 10, Nv2, Minv2, Maxv2}, "fMassXi", "fPt", v2Chosen, "fEffWeight");
      auto massVsPtVsPzs2 = dcentPzs2.Histo3D({Form("massVsPtVsPzs2Hist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Pzs2", 80, 1.28, 1.36, 100, 0, 10, NPzs2, MinPzs2, MaxPzs2}, "fMassXi", "fPt", "Pzs2Xi");
      auto massVsPsiVsPz = dcentPzs2.Histo3D({Form("massVsPsiVsPzHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Pz", 80, 1.28, 1.36, 20, 0, 2 * TMath::Pi(), NPz, MinPz, MaxPz}, "fMassXi", "f2PsiDiffCorr", "fCosThetaStarLambdaFromXi");
      auto PsiDiff = dcentPzs2.Histo1D({"2PsiDiffCorr", "2PsiDiffCorr", 100, 0, 2 * TMath::Pi()}, "f2PsiDiffCorr");
      auto PsiDiff2 = dcentPzs2.Histo1D({"2PsiDiff", "2PsiDiff", 100, -2 * TMath::Pi(), 2 * TMath::Pi()}, "f2PsiDiff");
      // profile: mean value of v2 vs mass and pt in centrality classes
      auto profile = dcent.Profile2D({Form("ProfilemassVsPtVsV2CHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Mean invariant mass vs Pt vs V2C", 80, 1.28, 1.36, numPtBins, PtBins}, "fMassXi", "fPt", v2Chosen);
      auto hrapidity = dcentPzs2.Histo1D({"rapidity", "Rapidity distribution of selected candidates", 200, -2, 2}, "fRapidity");
      auto hEta = dcent.Histo1D({"Eta", "Eta distribution of selected candidates", 200, -2, 2}, "fEta");

      MassCutHisto[cent] = (TH1D *)hMassCut->Clone(Form("massCut_cent%i-%i", CentFT0CMin, CentFT0CMax));
      massVsPtVsV2CHisto[cent] = (TH3D *)massVsPtVsV2C->Clone(Form("massVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax));
      massVsPtVsV2CWeightedHisto[cent] = (TH3D *)massVsPtVsV2CWeighted->Clone(Form("massVsPtVsV2CWeighted_cent%i-%i", CentFT0CMin, CentFT0CMax));
      massVsPtVsPzs2Histo[cent] = (TH3D *)massVsPtVsPzs2->Clone(Form("massVsPtVsPzs2_cent%i-%i", CentFT0CMin, CentFT0CMax));
      massVsPsiVsPzHisto[cent] = (TH3D *)massVsPsiVsPz->Clone(Form("massVsPsiVsPz_cent%i-%i", CentFT0CMin, CentFT0CMax));
      profileHisto[cent] = (TProfile *)profile->Clone(Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax));
      rapidityHisto[cent] = (TH1F *)hrapidity->Clone(Form("rapidity_cent%i-%i", CentFT0CMin, CentFT0CMax));
      hEtaHisto[cent] = (TH1F *)hEta->Clone(Form("Eta_cent%i-%i", CentFT0CMin, CentFT0CMax));
      hPsiDiffHisto[cent] = (TH1D *)PsiDiff->Clone(Form("2PsiDiffCorr_cent%i-%i", CentFT0CMin, CentFT0CMax));
      hPsiDiff2Histo[cent] = (TH1D *)PsiDiff2->Clone(Form("2PsiDiff_cent%i-%i", CentFT0CMin, CentFT0CMax));
    }
    else
    {
      dcentPzs2 = dcentPzs2.Define("Pzs2Omega", "fCosThetaStarLambdaFromOmega * sin(2*(fPhi-fPsiT0C))");
      auto hMassCut = dmasscut.Histo1D({"massCut", "Invariant mass of #LambdaK", 100, 1.6, 1.73}, "fMassOmega");
      auto massVsPtVsV2C = dcent.Histo3D({Form("massVsPtVsV2CHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs V2C", 80, 1.63, 1.726, 100, 0, 10, Nv2, Minv2, Maxv2}, "fMassOmega", "fPt", v2Chosen);
      auto massVsPtVsV2CWeighted = dcent.Histo3D({Form("massVsPtVsV2CHistWeighted_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs V2C", 80, 1.63, 1.726, 100, 0, 10, Nv2, Minv2, Maxv2}, "fMassOmega", "fPt", v2Chosen, "fEffWeight");
      auto massVsPtVsPzs2 = dcentPzs2.Histo3D({Form("massVsPtVsPzs2Hist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Pzs2", 80, 1.28, 1.36, 100, 0, 10, NPzs2, MinPzs2, MaxPzs2}, "fMassOmega", "fPt", "Pzs2Omega");
      auto massVsPsiVsPz = dcentPzs2.Histo3D({Form("massVsPsiVsPzHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Pz", 80, 1.63, 1.726, 20, 0, 2 * TMath::Pi(), NPz, MinPz, MaxPz}, "fMassOmega", "f2PsiDiffCorr", "fCosThetaStarLambdaFromOmega");
      // profile: mean value of v2 vs mass and pt in centrality classes
      auto profile = dcent.Profile2D({Form("ProfilemassVsPtVsV2CHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Mean invariant mass vs Pt vs V2C", 80, 1.63, 1.726, 100, 0, 10, -2, 2}, "fMassOmega", "fPt", v2Chosen);
      auto hrapidity = dcentPzs2.Histo1D({"rapidity", "Rapidity distribution of selected candidates", 200, -2, 2}, "fRapidity");
      auto hEta = dcent.Histo1D({"Eta", "Eta distribution of selected candidates", 200, -2, 2}, "fEta");
      auto PsiDiff = dcentPzs2.Histo1D({"2PsiDiffCorr", "2PsiDiffCorr", 100, 0, 2 * TMath::Pi()}, "f2PsiDiffCorr");
      auto PsiDiff2 = dcentPzs2.Histo1D({"2PsiDiff", "2PsiDiff", 100, -2 * TMath::Pi(), 2 * TMath::Pi()}, "f2PsiDiff");
      MassCutHisto[cent] = (TH1D *)hMassCut->Clone(Form("massCut_cent%i-%i", CentFT0CMin, CentFT0CMax));
      massVsPtVsV2CHisto[cent] = (TH3D *)massVsPtVsV2C->Clone(Form("massVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax));
      massVsPtVsV2CWeightedHisto[cent] = (TH3D *)massVsPtVsV2CWeighted->Clone(Form("massVsPtVsV2CWeighted_cent%i-%i", CentFT0CMin, CentFT0CMax));
      massVsPtVsPzs2Histo[cent] = (TH3D *)massVsPtVsPzs2->Clone(Form("massVsPtVsPzs2_cent%i-%i", CentFT0CMin, CentFT0CMax));
      massVsPsiVsPzHisto[cent] = (TH3D *)massVsPsiVsPz->Clone(Form("massVsPsiVsPz_cent%i-%i", CentFT0CMin, CentFT0CMax));
      profileHisto[cent] = (TProfile *)profile->Clone(Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax));
      rapidityHisto[cent] = (TH1F *)hrapidity->Clone(Form("rapidity_cent%i-%i", CentFT0CMin, CentFT0CMax));
      hEtaHisto[cent] = (TH1F *)hEta->Clone(Form("Eta_cent%i-%i", CentFT0CMin, CentFT0CMax));
      hPsiDiffHisto[cent] = (TH1D *)PsiDiff->Clone(Form("2PsiDiffCorr_cent%i-%i", CentFT0CMin, CentFT0CMax));
      hPsiDiff2Histo[cent] = (TH1D *)PsiDiff2->Clone(Form("2PsiDiff_cent%i-%i", CentFT0CMin, CentFT0CMax));
    }
    v2CHisto[cent] = (TH1D *)v2C->Clone(Form("v2C_cent%i-%i", CentFT0CMin, CentFT0CMax));
    hPhiCentHisto[cent] = (TH2D *)hPhiCent->Clone(Form("PhiHist_cent%i-%i", CentFT0CMin, CentFT0CMax));
    hPsiCentHisto[cent] = (TH1D *)hPsiCent->Clone(Form("PsiHist_cent%i-%i", CentFT0CMin, CentFT0CMax));
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
  hPtvsCent_Bef->Write();
  hPtvsCent_RapCut_Bef->Write();
  hPtvsCent_Aft->Write();
  hPtvsCent_RapCut_Aft->Write();
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
  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    MassCutHisto[cent]->Write();
    v2CHisto[cent]->Write();
    hPhiCentHisto[cent]->Write();
    hPsiCentHisto[cent]->Write();
    massVsPtVsV2CHisto[cent]->Write();
    massVsPtVsV2CWeightedHisto[cent]->Write();
    massVsPtVsPzs2Histo[cent]->Write();
    massVsPsiVsPzHisto[cent]->Write();
    profileHisto[cent]->Write();
    rapidityHisto[cent]->Write();
    hEtaHisto[cent]->Write();
    hPsiDiffHisto[cent]->Write();
    hPsiDiff2Histo[cent]->Write();
  }
  file->Close();

  cout << "I created the file " << file->GetName() << endl;
}
