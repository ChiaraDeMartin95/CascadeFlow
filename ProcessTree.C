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
Float_t MinPzs2WithAlphaXi = -2.8;
Float_t MinPzs2WithAlphaOmega = -65;
Float_t MaxPzs2 = 1;
Float_t MaxPzs2WithAlphaXi = 2.8;
Float_t MaxPzs2WithAlphaOmega = 65;
Int_t NPzs2 = 200;

Float_t MinPz = -1;
Float_t MinPzWithAlphaXi = -2.8;
Float_t MinPzWithAlphaOmega = -65;
Float_t MaxPz = 1;
Float_t MaxPzWithAlphaXi = 2.8;
Float_t MaxPzWithAlphaOmega = 65;
Int_t NPz = 200;

bool passbdtCut(float bdtscore, float cent, int indexMultTrial)
{
  float sbdtCut = 0;
  if (ExtrisSysMultTrial)
  {
    sbdtCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;
  }
  else
  {
    if (useCommonBDTValue)
      sbdtCut = DefaultBDTscoreCut;
    else
      sbdtCut = bdtCut[int(cent / 10)];
  }
  return (bdtscore > sbdtCut);
}

void ProcessTree(Bool_t isEff = 0,
                 Int_t indexMultTrial = 0,
                 Bool_t isRapiditySel = ExtrisRapiditySel,
                 Int_t ChosenPart = ChosenParticle,
                 Bool_t isApplyEffWeights = 0,
                 TString inputFileName = SinputFileName,
                 Int_t EtaSysChoice = ExtrEtaSysChoice,
                 Bool_t isSysMultTrial = ExtrisSysMultTrial)
{

  TString weightFileName = "PhiWeights/Weights_" + SinputFileName + "_" + ParticleName[ChosenParticle] + ".root";

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

  cout << "BDT score cut: " << BDTscoreCut << endl;
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
  if (ChosenPart == 3 || ChosenPart == 5)
    chargecut = "fSign > 0";
  else if (ChosenPart == 2 || ChosenPart == 4)
    chargecut = "fSign < 0";
  auto d2 = d1.Filter(chargecut);

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

  // define variable for syst studies
  d2 = d2.Define("BDTVarindex", Form("%i", indexMultTrial));

  // define the rapidity
  string Common1 = Form("std::sqrt(fPt*fPt*std::cosh(fEta)*std::cosh(fEta) + %.3f*%.3f)", ParticleMassPDG[ChosenPart], ParticleMassPDG[ChosenPart]);
  string Common2 = "fPt*std::sinh(fEta)";
  string Num = Form("%s + %s", Common1.c_str(), Common2.c_str());
  string Denom = Form("%s - %s", Common1.c_str(), Common2.c_str());
  d2 = d2.Define("fRapidity", Form("0.5 * std::log((%s) / (%s))", Num.c_str(), Denom.c_str()));

  // rapidity and eta selections
  auto d2RapCut = d2.Filter("abs(fRapidity) < 0.5"); // only |y| < 0.5 selection (no pseudorapidity selection)

  if (EtaSysChoice == 0) // -0.8 < eta < 0.8
    d2 = d2.Filter("fEta > -0.8 && fEta < 0.8");
  else if (EtaSysChoice == 1) // eta > 0 && eta < 0.8
    d2 = d2.Filter("fEta > 0 && fEta < 0.8");
  else if (EtaSysChoice == 2) // eta < 0 && eta > -0.8
    d2 = d2.Filter("fEta < 0 && fEta > -0.8");

  // pt vs centrality before BDT cut
  auto hPtvsCent_Bef = d2.Histo2D({"PtvsCent_BefBDT", "PtvsCent_BefBDT", 100, 0, 100, 400, 0, 20}, "fCentFT0C", "fPt");
  auto hPtvsCent_RapCut_Bef = d2RapCut.Histo2D({"PtvsCent_Y05_BefBDT", "PtvsCent_Y05_BefBDT", 100, 0, 100, 400, 0, 20}, "fCentFT0C", "fPt");
  auto BDT_response_Bef2 = d1.Histo1D({"BDT_response_Bef2", "BDT response", 100, 0, 1}, "fBDTResponseXi");
  if (!isXi)
    BDT_response_Bef2 = d1.Histo1D({"BDT_response_Bef2", "BDT response", 100, 0, 1}, "fBDTResponseOmega");

  //  apply BDT selection
  string cutvariable = "fBDTResponseXi";
  if (!isXi)
    cutvariable = "fBDTResponseOmega";

  auto d3 = d2.Filter(passbdtCut, {cutvariable, "fCentFT0C", "BDTVarindex"});
  auto d3RapCut = d2RapCut.Filter(passbdtCut, {cutvariable, "fCentFT0C", "BDTVarindex"});

  // pt vs centrality after BDT cut
  auto hPtvsCent_Aft = d3.Histo2D({"PtvsCent_AftBDT", "PtvsCent_AftBDT", 100, 0, 100, 400, 0, 20}, "fCentFT0C", "fPt");
  auto hPtvsCent_RapCut_Aft = d3RapCut.Histo2D({"PtvsCent_Y05_AftBDT", "PtvsCent_Y05_AftBDT", 100, 0, 100, 400, 0, 20}, "fCentFT0C", "fPt");

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
  // TString SEffWeights = "Efficiency/Efficiency_" + SinputFileNameEff + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice] + SBDT;
  TString SEffWeights = "Efficiency/Efficiency_" + SinputFileNameEff + "_" + ParticleName[!isXi] + SEtaSysChoice[EtaSysChoice] + SBDT;
  // if (!useCommonBDTValue)
  // SEffWeights += "_BDTCentDep";
  if (isRun2Binning)
    SEffWeights += "_Run2Binning";
  SEffWeights += "_" + RapidityCoverage[!isV2];
  SEffWeights += ".root";
  TFile *fileEffWeights = new TFile(SEffWeights, "");
  if (isApplyEffWeights)
    cout << fileEffWeights->GetName() << endl;
  if (!fileEffWeights && isApplyEffWeights)
  {
    cout << "Efficiency file not found" << endl;
    return;
  }
  TH1D *par0_Eff;
  TH1D *par1_Eff;
  if (isApplyEffWeights)
  {
    par0_Eff = (TH1D *)fileEffWeights->Get("hpar0_Eff");
    par1_Eff = (TH1D *)fileEffWeights->Get("hpar1_Eff");
    if (!par0_Eff || !par1_Eff)
    {
      cout << "Efficiency histograms not found" << endl;
      return;
    }
  }

  // Published v2 values vs cent and pt
  // TFile *v2PubFile = new TFile("v2Histo.root", "");
  // TH2F *v2Histo = (TH2F *)v2PubFile->Get("v2Histo");

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
  // OutputFileName += "_EffWBis.root";
  // OutputFileName += "_Prova1234.root";
  // OutputFileName += "_TestWithAlpha.root";
  if (!isRapiditySel)
    OutputFileName += "_Eta08";
  OutputFileName += ".root";
  TFile *file = new TFile(OutputFileName, "RECREATE");
  cout << file->GetName() << endl;

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;

  std::vector<ROOT::RDF::RResultPtr<TH1D>> hrapidityVector;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hEtaVector;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hEtaPzs2Vector;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> PsiDiffVector;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> PsiDiff2Vector;

  std::vector<ROOT::RDF::RResultPtr<TH1D>> hMassCutVector;
  std::vector<ROOT::RDF::RResultPtr<TH2D>> hPhiCentVector;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hPsiCentVector;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> v2CVector;
  std::vector<ROOT::RDF::RResultPtr<TProfile2D>> profileVector;

  std::vector<ROOT::RDF::RResultPtr<TH2D>> NchVarVector;
  std::vector<ROOT::RDF::RResultPtr<TH2D>> NchTimesV2Vector;
  std::vector<ROOT::RDF::RResultPtr<TH2D>> effWeightVector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> effWeight3DVector;

  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsV2CVector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsV2CWeightedVector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsPzs2Vector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsPzs2LambdaFromCVector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPsiVsPzVector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPsiVsPzLambdaFromCVector;

  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsPzs2VectorWithAlpha;            // decay parameter included in the calculation
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsPzs2LambdaFromCVectorWithAlpha; // decay parameter included in the calculation
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPsiVsPzVectorWithAlpha;             // decay parameter included in the calculation
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPsiVsPzLambdaFromCVectorWithAlpha;  // decay parameter included in the calculation

  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsCos2Vector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsCos2LambdaFromCVector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPsiVsCos2Vector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPsiVsCos2LambdaFromCVector;

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
    auto dcent = d3.Filter(Form("fCentFT0C>=%.1f && fCentFT0C<%.1f", CentFT0CMin + 0.001, CentFT0CMax - 0.001));
    string MassCut = "";
    if (isXi)
      MassCut = Form("abs(fMassXi - %.3f) < 3* (%.3f * exp(%.3f * fPt) + %.3f * exp(%.3f * fPt))", ParticleMassPDG[ChosenPart], massSigmaParameters[0][0], massSigmaParameters[1][0], massSigmaParameters[2][0], massSigmaParameters[3][0]);
    else
      MassCut = Form("abs(fMassOmega - %.3f) < 3* (%.3f * exp(%.3f * fPt) + %.3f * exp(%.3f * fPt))", ParticleMassPDG[ChosenPart], massSigmaParameters[0][1], massSigmaParameters[1][1], massSigmaParameters[2][1], massSigmaParameters[3][1]);
    auto dmasscut = dcent.Filter(MassCut);

    auto v2C = dcent.Histo1D({Form("v2C_cent%i-%i", CentFT0CMin, CentFT0CMax), "v2C", Nv2, Minv2, Maxv2}, v2Chosen);
    v2CVector.push_back(v2C);
    auto hPhiCent = dmasscut.Histo2D({Form("PhiHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Phi vs pt", 100, 0, 10, 100, 0, 2 * TMath::Pi()}, "fPt", "fPhi");
    hPhiCentVector.push_back(hPhiCent);
    auto hPsiCent = dmasscut.Histo1D({Form("PsiHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Psi", 100, -2 * TMath::Pi(), 2 * TMath::Pi()}, "fPsiT0C");
    hPsiCentVector.push_back(hPsiCent);

    dcent = dcent.Define("fPsiDiff", "if ((fPhi-fPsiT0C) < 0) return (fPhi-fPsiT0C+(float)TMath::Pi()); else if ((fPhi-fPsiT0C) > 2* TMath::Pi()) return (fPhi-fPsiT0C-2*(float)TMath::Pi()); else if ((fPhi-fPsiT0C) > TMath::Pi()) return (fPhi-fPsiT0C-(float)TMath::Pi()); else return (fPhi-fPsiT0C);");
    dcent = dcent.Define("f2PsiDiffCorr", "2*fPsiDiff");
    dcent = dcent.Define("f2PsiDiff", "2*fPhi-2*fPsiT0C");
    dcent = dcent.Define("Nch", Form("%.2f", dNdEtaAbhi[cent]));
    dcent = dcent.Define("v2Pub", Form("%.3f", v2PubRun2[cent]));

    // rapidity or eta selection
    auto dcentPzs2 = dcent;
    if (isRapiditySel)
    {
      dcentPzs2 = dcent.Filter("fRapidity > -0.5 && fRapidity < 0.5");
    }
    else
    {
      dcentPzs2 = dcent.Filter("fEta > -0.8 && fEta < 0.8");
    }

    if (isApplyEffWeights)
    {
      dcent = dcent.Define("denW", [&par0_Eff, &par1_Eff](double Nch, double v2, float y, float pt)
                           {
                            int ptBin = par0_Eff->GetXaxis()->FindBin(pt);
                            return std::exp(par0_Eff->GetBinContent(ptBin) + par1_Eff->GetBinContent(ptBin) * Nch * (1 + 2 * v2 * cos(y))); }, {"Nch", "v2Pub", "f2PsiDiffCorr", "fPt"});
      dcent = dcent.Define("numW", [&par0_Eff, &par1_Eff](double Nch, float pt)
                           {
                           int ptBin = par0_Eff->GetXaxis()->FindBin(pt);
                           return std::exp(par0_Eff->GetBinContent(ptBin) + par1_Eff->GetBinContent(ptBin) * Nch); }, {"Nch", "fPt"});
      dcent = dcent.Define("fEffWeight", "numW/denW");
      dcent = dcent.Define("NchVar", [](double v2, float y)
                           { return (1 + 2 * v2 * cos(y)); }, {"v2Pub", "f2PsiDiffCorr"});
      dcent = dcent.Define("NchxV2", "Nch * v2Pub");
    }
    else
    {
      dcent = dcent.Define("fEffWeight", "1");
      dcent = dcent.Define("NchVar", "1");
      dcent = dcent.Define("NchxV2", "1");
    }
    dcentPzs2 = dcentPzs2.Define("fAlphaLambda", Form("if (fSign < 0) return %.4f; else return %.4f;", AlphaLambda[2], AlphaLambda[3]));
    dcentPzs2 = dcentPzs2.Define("Pzs2LambdaFromC", "fCosThetaStarProton * sin(2*(fPhi-fPsiT0C))");
    dcentPzs2 = dcentPzs2.Define("fCos2ThetaStarProton", "fCosThetaStarProton * fCosThetaStarProton");
    dcentPzs2 = dcentPzs2.Define("Pzs2LambdaFromCAlpha", "Pzs2LambdaFromC / fAlphaLambda");
    dcentPzs2 = dcentPzs2.Define("fCosThetaStarProtonAlpha", "fCosThetaStarProton / fAlphaLambda");

    auto NchVar = dcent.Histo2D({Form("NchVar_cent%i-%i", CentFT0CMin, CentFT0CMax), "Rel. variation of Nch vs 2*(Psi-Phi)", 100, 0, 2 * TMath::Pi(), 100, 0.7, 1.3}, "f2PsiDiffCorr", "NchVar");
    NchVarVector.push_back(NchVar);
    auto NchTimesV2 = dcent.Histo2D({Form("NchTimesV2_cent%i-%i", CentFT0CMin, CentFT0CMax), "Nch * v2 vs 2*(Psi-Phi)", 100, 0, 2 * TMath::Pi(), 1000, 0, 1000}, "f2PsiDiffCorr", "NchxV2");
    NchTimesV2Vector.push_back(NchTimesV2);
    auto effWeight = dcent.Histo2D({Form("EffWeight_cent%i-%i", CentFT0CMin, CentFT0CMax), "Efficiency weight vs 2*(Psi-Phi)", 20, 0, 2 * TMath::Pi(), 100, 0.7, 1.3}, "f2PsiDiffCorr", "fEffWeight");
    effWeightVector.push_back(effWeight);
    auto effWeight3D = dcent.Histo3D({Form("EffWeight3D_cent%i-%i", CentFT0CMin, CentFT0CMax), "Efficiency weight vs Pt vs 2*(Psi-Phi)", 100, 0, 10, 20, 0, 2 * TMath::Pi(), 100, 0.7, 1.3}, "fPt", "f2PsiDiffCorr", "fEffWeight");
    effWeight3DVector.push_back(effWeight3D);

    auto PsiDiff = dcentPzs2.Histo1D({Form("2PsiDiffCorr_cent%i-%i", CentFT0CMin, CentFT0CMax), "2PsiDiffCorr", 100, 0, 2 * TMath::Pi()}, "f2PsiDiffCorr");
    PsiDiffVector.push_back(PsiDiff);
    auto PsiDiff2 = dcentPzs2.Histo1D({Form("2PsiDiff_cent%i-%i", CentFT0CMin, CentFT0CMax), "2PsiDiff", 100, -2 * TMath::Pi(), 2 * TMath::Pi()}, "f2PsiDiff");
    PsiDiff2Vector.push_back(PsiDiff2);
    auto hrapidity = dcentPzs2.Histo1D({Form("rapidity_cent%i-%i", CentFT0CMin, CentFT0CMax), "Rapidity distribution of selected candidates", 200, -2, 2}, "fRapidity");
    hrapidityVector.push_back(hrapidity);
    auto hEta = dcent.Histo1D({Form("Eta_cent%i-%i", CentFT0CMin, CentFT0CMax), "Eta distribution of selected candidates", 200, -2, 2}, "fEta");
    hEtaVector.push_back(hEta);
    auto hEtaPzs2 = dcentPzs2.Histo1D({Form("EtaPzs2_cent%i-%i", CentFT0CMin, CentFT0CMax), "Eta distribution of selected candidates for Pzs2 analysis", 200, -2, 2}, "fEta");
    hEtaPzs2Vector.push_back(hEtaPzs2);

    if (isXi)
    {
      dcentPzs2 = dcentPzs2.Define("fAlphaXi", Form("if (fSign < 0) return %.4f; else return %.4f;", AlphaH[2], AlphaH[3]));
      dcentPzs2 = dcentPzs2.Define("Pzs2Xi", "fCosThetaStarLambdaFromXi * sin(2*(fPhi-fPsiT0C))");
      dcentPzs2 = dcentPzs2.Define("fCos2ThetaStarLambda", "fCosThetaStarLambdaFromXi * fCosThetaStarLambdaFromXi");
      dcentPzs2 = dcentPzs2.Define("Pzs2XiAlpha", "fCosThetaStarLambdaFromXi * sin(2*(fPhi-fPsiT0C)) / fAlphaXi");
      dcentPzs2 = dcentPzs2.Define("fCosThetaStarLambdaFromXiAlpha", "fCosThetaStarLambdaFromXi / fAlphaXi");

      auto hMassCut = dmasscut.Histo1D({Form("massCut_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass of #Lambda#pi", 100, 1.28, 1.36}, "fMassXi");
      hMassCutVector.push_back(hMassCut);

      auto massVsPtVsV2C = dcent.Histo3D({Form("massVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs V2C", 80, 1.28, 1.36, 100, 0, 10, Nv2, Minv2, Maxv2}, "fMassXi", "fPt", v2Chosen);
      massVsPtVsV2CVector.push_back(massVsPtVsV2C);
      auto massVsPtVsV2CWeighted = dcent.Histo3D({Form("massVsPtVsV2CWeighted_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs V2C", 80, 1.28, 1.36, 100, 0, 10, Nv2, Minv2, Maxv2}, "fMassXi", "fPt", v2Chosen, "fEffWeight");
      massVsPtVsV2CWeightedVector.push_back(massVsPtVsV2CWeighted);
      auto massVsPtVsPzs2 = dcentPzs2.Histo3D({Form("massVsPtVsPzs2_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Pzs2", 80, 1.28, 1.36, 100, 0, 10, NPzs2, MinPzs2, MaxPzs2}, "fMassXi", "fPt", "Pzs2Xi");
      massVsPtVsPzs2Vector.push_back(massVsPtVsPzs2);
      auto massVsPtVsPzs2LambdaFromC = dcentPzs2.Histo3D({Form("massVsPtVsPzs2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Pzs2", 80, 1.28, 1.36, 100, 0, 10, NPzs2, MinPzs2, MaxPzs2}, "fMassXi", "fPt", "Pzs2LambdaFromC");
      massVsPtVsPzs2LambdaFromCVector.push_back(massVsPtVsPzs2LambdaFromC);
      auto massVsPsiVsPz = dcentPzs2.Histo3D({Form("massVsPsiVsPz_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Pz", 80, 1.28, 1.36, 20, 0, 2 * TMath::Pi(), NPz, MinPz, MaxPz}, "fMassXi", "f2PsiDiffCorr", "fCosThetaStarLambdaFromXi");
      massVsPsiVsPzVector.push_back(massVsPsiVsPz);
      auto massVsPsiVsPzLambdaFromC = dcentPzs2.Histo3D({Form("massVsPsiVsPzLambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Pz", 80, 1.28, 1.36, 20, 0, 2 * TMath::Pi(), NPz, MinPz, MaxPz}, "fMassXi", "f2PsiDiffCorr", "fCosThetaStarProton");
      massVsPsiVsPzLambdaFromCVector.push_back(massVsPsiVsPzLambdaFromC);

      auto massVsPtVsPzs2WithAlpha = dcentPzs2.Histo3D({Form("massVsPtVsPzs2_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Pzs2", 80, 1.28, 1.36, 100, 0, 10, NPzs2, MinPzs2WithAlphaXi, MaxPzs2WithAlphaXi}, "fMassXi", "fPt", "Pzs2XiAlpha");
      massVsPtVsPzs2VectorWithAlpha.push_back(massVsPtVsPzs2WithAlpha);
      auto massVsPtVsPzs2LambdaFromCWithAlpha = dcentPzs2.Histo3D({Form("massVsPtVsPzs2LambdaFromC_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Pzs2", 80, 1.28, 1.36, 100, 0, 10, NPzs2, MinPzs2WithAlphaXi, MaxPzs2WithAlphaXi}, "fMassXi", "fPt", "Pzs2LambdaFromCAlpha");
      massVsPtVsPzs2LambdaFromCVectorWithAlpha.push_back(massVsPtVsPzs2LambdaFromCWithAlpha);
      auto massVsPsiVsPzWithAlpha = dcentPzs2.Histo3D({Form("massVsPsiVsPz_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Pz", 80, 1.28, 1.36, 20, 0, 2 * TMath::Pi(), NPz, MinPzWithAlphaXi, MaxPzWithAlphaXi}, "fMassXi", "f2PsiDiffCorr", "fCosThetaStarLambdaFromXiAlpha");
      massVsPsiVsPzVectorWithAlpha.push_back(massVsPsiVsPzWithAlpha);
      auto massVsPsiVsPzLambdaFromCWithAlpha = dcentPzs2.Histo3D({Form("massVsPsiVsPzLambdaFromC_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Pz", 80, 1.28, 1.36, 20, 0, 2 * TMath::Pi(), NPz, MinPzWithAlphaXi, MaxPzWithAlphaXi}, "fMassXi", "f2PsiDiffCorr", "fCosThetaStarProtonAlpha");
      massVsPsiVsPzLambdaFromCVectorWithAlpha.push_back(massVsPsiVsPzLambdaFromCWithAlpha);

      auto massVsPtVsCos2 = dcentPzs2.Histo3D({Form("massVsPtVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Cos2", 80, 1.28, 1.36, 100, 0, 10, 100, 0, 1}, "fMassXi", "fPt", "fCos2ThetaStarLambda");
      massVsPtVsCos2Vector.push_back(massVsPtVsCos2);
      auto massVsPtVsCos2LambdaFromC = dcentPzs2.Histo3D({Form("massVsPtVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Cos2", 80, 1.28, 1.36, 100, 0, 10, 100, 0, 1}, "fMassXi", "fPt", "fCos2ThetaStarProton");
      massVsPtVsCos2LambdaFromCVector.push_back(massVsPtVsCos2LambdaFromC);
      auto massVsPsiVsCos2 = dcentPzs2.Histo3D({Form("massVsPsiVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Cos2", 80, 1.28, 1.36, 20, 0, 2 * TMath::Pi(), 100, 0, 1}, "fMassXi", "f2PsiDiffCorr", "fCos2ThetaStarLambda");
      massVsPsiVsCos2Vector.push_back(massVsPsiVsCos2);
      auto massVsPsiVsCos2LambdaFromC = dcentPzs2.Histo3D({Form("massVsPsiVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Cos2", 80, 1.28, 1.36, 20, 0, 2 * TMath::Pi(), 100, 0, 1}, "fMassXi", "f2PsiDiffCorr", "fCos2ThetaStarProton");
      massVsPsiVsCos2LambdaFromCVector.push_back(massVsPsiVsCos2LambdaFromC);

      // profile: mean value of v2 vs mass and pt in centrality classes
      auto profile = dcent.Profile2D({Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax), "Mean invariant mass vs Pt vs V2C", 80, 1.28, 1.36, numPtBins, PtBins}, "fMassXi", "fPt", v2Chosen);
      profileVector.push_back(profile);
    }
    else
    {
      dcentPzs2 = dcentPzs2.Define("fAlphaOmega", Form("if (fSign < 0) return %.4f; else return %.4f;", AlphaH[4], AlphaH[5]));
      dcentPzs2 = dcentPzs2.Define("Pzs2Omega", "fCosThetaStarLambdaFromOmega * sin(2*(fPhi-fPsiT0C))");
      dcentPzs2 = dcentPzs2.Define("fCos2ThetaStarLambda", "fCosThetaStarLambdaFromOmega * fCosThetaStarLambdaFromOmega");
      dcentPzs2 = dcentPzs2.Define("Pzs2OmegaAlpha", "fCosThetaStarLambdaFromOmega * sin(2*(fPhi-fPsiT0C)) / fAlphaOmega");
      dcentPzs2 = dcentPzs2.Define("fCosThetaStarLambdaFromOmegaAlpha", "fCosThetaStarLambdaFromOmega / fAlphaOmega");

      auto hMassCut = dmasscut.Histo1D({Form("massCut_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass of #LambdaK", 100, 1.6, 1.73}, "fMassOmega");
      hMassCutVector.push_back(hMassCut);

      auto massVsPtVsV2C = dcent.Histo3D({Form("massVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs V2C", 80, 1.63, 1.726, 100, 0, 10, Nv2, Minv2, Maxv2}, "fMassOmega", "fPt", v2Chosen);
      massVsPtVsV2CVector.push_back(massVsPtVsV2C);
      auto massVsPtVsV2CWeighted = dcent.Histo3D({Form("massVsPtVsV2CWeighted_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs V2C", 80, 1.63, 1.726, 100, 0, 10, Nv2, Minv2, Maxv2}, "fMassOmega", "fPt", v2Chosen, "fEffWeight");
      massVsPtVsV2CWeightedVector.push_back(massVsPtVsV2CWeighted);
      auto massVsPtVsPzs2 = dcentPzs2.Histo3D({Form("massVsPtVsPzs2_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Pzs2", 80, 1.28, 1.36, 100, 0, 10, NPzs2, MinPzs2, MaxPzs2}, "fMassOmega", "fPt", "Pzs2Omega");
      massVsPtVsPzs2Vector.push_back(massVsPtVsPzs2);
      auto massVsPtVsPzs2LambdaFromC = dcentPzs2.Histo3D({Form("massVsPtVsPzs2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Pzs2", 80, 1.63, 1.726, 100, 0, 10, NPzs2, MinPzs2, MaxPzs2}, "fMassOmega", "fPt", "Pzs2LambdaFromC");
      massVsPtVsPzs2LambdaFromCVector.push_back(massVsPtVsPzs2LambdaFromC);
      auto massVsPsiVsPz = dcentPzs2.Histo3D({Form("massVsPsiVsPz_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Pz", 80, 1.63, 1.726, 20, 0, 2 * TMath::Pi(), NPz, MinPz, MaxPz}, "fMassOmega", "f2PsiDiffCorr", "fCosThetaStarLambdaFromOmega");
      massVsPsiVsPzVector.push_back(massVsPsiVsPz);
      auto massVsPsiVsPzLambdaFromC = dcentPzs2.Histo3D({Form("massVsPsiVsPzLambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Pz", 80, 1.63, 1.726, 20, 0, 2 * TMath::Pi(), NPz, MinPz, MaxPz}, "fMassOmega", "f2PsiDiffCorr", "fCosThetaStarProton");
      massVsPsiVsPzLambdaFromCVector.push_back(massVsPsiVsPzLambdaFromC);

      auto massVsPtVsPzs2WithAlpha = dcentPzs2.Histo3D({Form("massVsPtVsPzs2_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Pzs2", 80, 1.63, 1.726, 100, 0, 10, NPzs2, MinPzs2WithAlphaOmega, MaxPzs2WithAlphaOmega}, "fMassOmega", "fPt", "Pzs2OmegaAlpha");
      massVsPtVsPzs2VectorWithAlpha.push_back(massVsPtVsPzs2WithAlpha);
      auto massVsPtVsPzs2LambdaFromCWithAlpha = dcentPzs2.Histo3D({Form("massVsPtVsPzs2LambdaFromC_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Pzs2", 80, 1.63, 1.726, 100, 0, 10, NPzs2, MinPzs2WithAlphaOmega, MaxPzs2WithAlphaOmega}, "fMassOmega", "fPt", "Pzs2LambdaFromCAlpha");
      massVsPtVsPzs2LambdaFromCVectorWithAlpha.push_back(massVsPtVsPzs2LambdaFromCWithAlpha);
      auto massVsPsiVsPzWithAlpha = dcentPzs2.Histo3D({Form("massVsPsiVsPz_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Pz", 80, 1.63, 1.726, 20, 0, 2 * TMath::Pi(), NPz, MinPzWithAlphaOmega, MaxPzWithAlphaOmega}, "fMassOmega", "f2PsiDiffCorr", "fCosThetaStarLambdaFromOmegaAlpha");
      massVsPsiVsPzVectorWithAlpha.push_back(massVsPsiVsPzWithAlpha);
      auto massVsPsiVsPzLambdaFromCWithAlpha = dcentPzs2.Histo3D({Form("massVsPsiVsPzLambdaFromC_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Pz", 80, 1.63, 1.726, 20, 0, 2 * TMath::Pi(), NPz, MinPzWithAlphaOmega, MaxPzWithAlphaOmega}, "fMassOmega", "f2PsiDiffCorr", "fCosThetaStarProtonAlpha");
      massVsPsiVsPzLambdaFromCVectorWithAlpha.push_back(massVsPsiVsPzLambdaFromCWithAlpha);

      auto massVsPtVsCos2 = dcentPzs2.Histo3D({Form("massVsPtVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Cos2", 80, 1.63, 1.726, 100, 0, 10, 100, 0, 1}, "fMassOmega", "fPt", "fCos2ThetaStarLambda");
      massVsPtVsCos2Vector.push_back(massVsPtVsCos2);
      auto massVsPtVsCos2LambdaFromC = dcentPzs2.Histo3D({Form("massVsPtVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Cos2", 80, 1.63, 1.726, 100, 0, 10, 100, 0, 1}, "fMassOmega", "fPt", "fCos2ThetaStarProton");
      massVsPtVsCos2LambdaFromCVector.push_back(massVsPtVsCos2LambdaFromC);
      auto massVsPsiVsCos2 = dcentPzs2.Histo3D({Form("massVsPsiVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Cos2", 80, 1.63, 1.726, 20, 0, 2 * TMath::Pi(), 100, 0, 1}, "fMassOmega", "f2PsiDiffCorr", "fCos2ThetaStarLambda");
      massVsPsiVsCos2Vector.push_back(massVsPsiVsCos2);
      auto massVsPsiVsCos2LambdaFromC = dcentPzs2.Histo3D({Form("massVsPsiVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Cos2", 80, 1.63, 1.726, 20, 0, 2 * TMath::Pi(), 100, 0, 1}, "fMassOmega", "f2PsiDiffCorr", "fCos2ThetaStarProton");
      massVsPsiVsCos2LambdaFromCVector.push_back(massVsPsiVsCos2LambdaFromC);

      // profile: mean value of v2 vs mass and pt in centrality classes
      auto profile = dcent.Profile2D({Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax), "Mean invariant mass vs Pt vs V2C", 80, 1.63, 1.726, 100, 0, 10, -2, 2}, "fMassOmega", "fPt", v2Chosen);
      profileVector.push_back(profile);
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
  h->Write();
  hPtvsCent_Bef->Write();
  hPtvsCent_RapCut_Bef->Write();
  hPtvsCent_Aft->Write();
  hPtvsCent_RapCut_Aft->Write();
  cMass->Write();
  MassXi->Write();
  hmass_Bef->Write();
  BDT_response_Bef->Write();
  BDT_response_Bef2->Write();
  mass_vs_BDTResponse->Write();
  hmass->Write();
  hmassvsPt->Write();
  heta->Write();
  hphi->Write();
  hEtaPhi->Write();
  BDT_response->Write();

  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    hrapidityVector[cent]->Write();
    hEtaVector[cent]->Write();
    hEtaPzs2Vector[cent]->Write();
    PsiDiffVector[cent]->Write();
    PsiDiff2Vector[cent]->Write();

    hMassCutVector[cent]->Write();
    v2CVector[cent]->Write();
    hPhiCentVector[cent]->Write();
    hPsiCentVector[cent]->Write();

    massVsPtVsV2CVector[cent]->Write();
    massVsPtVsV2CWeightedVector[cent]->Write();

    NchVarVector[cent]->Write();
    NchTimesV2Vector[cent]->Write();
    effWeightVector[cent]->Write();
    effWeight3DVector[cent]->Write();

    massVsPtVsPzs2Vector[cent]->Write();
    massVsPtVsPzs2LambdaFromCVector[cent]->Write();
    massVsPsiVsPzVector[cent]->Write();
    massVsPsiVsPzLambdaFromCVector[cent]->Write();

    massVsPtVsPzs2VectorWithAlpha[cent]->Write();
    massVsPtVsPzs2LambdaFromCVectorWithAlpha[cent]->Write();
    massVsPsiVsPzVectorWithAlpha[cent]->Write();
    massVsPsiVsPzLambdaFromCVectorWithAlpha[cent]->Write();

    massVsPtVsCos2Vector[cent]->Write();
    massVsPtVsCos2LambdaFromCVector[cent]->Write();
    massVsPsiVsCos2Vector[cent]->Write();
    massVsPsiVsCos2LambdaFromCVector[cent]->Write();

    profileVector[cent]->Write();
  }

  file->Close();

  cout << "I created the file " << file->GetName() << endl;
}
