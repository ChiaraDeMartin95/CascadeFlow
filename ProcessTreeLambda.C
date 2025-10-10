// This macro was originally written by:
// chiara.de.martin@cern.ch
#include <fstream>
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
#include "TKey.h"
#include "TChain.h"
#include <ROOT/RDataFrame.hxx>
#include "ROOT/RResultPtr.hxx"
#include "ROOT/RDFHelpers.hxx"
#include <random>
#include "ROOT/RDF/RInterface.hxx"

using namespace ROOT;
using namespace std;

void createChain(TChain& chainD, const std::string &filename, const std::string &treename) {
  TFile fileD(filename.c_str());
  TIter next(fileD.GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    std::string keyName = key->GetName();
    if (keyName.find("DF_") != std::string::npos) {
      chainD.Add((filename + "/" + keyName + "/" + chainD.GetName()).c_str());
    }
  }
}

constexpr double massSigmaParameters[4][2]{
    {4.9736e-3, 0.006815},
    {-2.39594, -2.257},
    {1.8064e-3, 0.00138},
    {1.03468e-1, 0.1898}};

Float_t Minv2 = -1;
Float_t Maxv2 = 1;
Int_t Nv2 = 200;

Float_t MinPzs2 = -1;
Float_t MinPzs2Reso[numCentLambdaOO + 1] = {-20, -23, -30, -35, -40, -60, -65, -80, -120, -160, -240};
Float_t MinPzs2WithAlphaXi = -2.8;
Float_t MinPzs2WithAlphaOmega = -65;
Float_t MaxPzs2 = 1;
Float_t MaxPzs2Reso[numCentLambdaOO + 1] = {20, 23, 30, 35, 40, 60, 65, 80, 120, 160, 240};
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

void ProcessTreeLambda(Int_t indexMultTrial = 0,
                       Bool_t isRapiditySel = ExtrisRapiditySel,
		       Bool_t isApplyResoOnTheFly = ExtrisApplyResoOnTheFly,
		       //		       Bool_t isPartialEta = ExtrisPartialEta,
                       Int_t ChosenPart = ChosenParticle,
                       TString inputFileName = SinputFileName,
                       Int_t EtaSysChoice = ExtrEtaSysChoice,
                       Bool_t isSysMultTrial = ExtrisSysMultTrial)
{

  ROOT::EnableImplicitMT();

  string v2Chosen = "fV2CEP";

  if (isSysMultTrial)
    inputFileName = SinputFileNameSyst;
  if (ChosenPart!=6 && (isLoosest || isTightest)){
    cout << "You set the loosest or the tightest topo sel, this can only be done for lambdas" << endl;
    return;
  }

  cout << "ChosenPart: " << ParticleName[ChosenPart] << endl;
  cout << "EtaSysChoice: " << EtaSysChoice << endl;
  //  cout << "IsPartialEta: only eta > 0 selected (opposite to FT0C)" << endl;

  std::vector<std::string> name;
  TString filename = "input_" +  SinputFileName + ".txt";
  std::ifstream fileIn(Form("%s", filename.Data()));

  cout << filename.Data() << endl;
  if (fileIn.is_open())  {
    std::string line;
    while (std::getline(fileIn, line)){
      name.push_back(line);
      cout << line << endl;
    }
    fileIn.close();
  }
  else  {
    std::cerr << "Unable to open file: " << std::endl;
  }
  
  const int nfiles = name.size();
  TString TreeName = "O2lambdaanalysis";

  //  TString inputFile = "TreeForAnalysis";
  //  inputFile += "/AnalysisResults_trees_" + inputFileName + "_New.root";  

  TFile *inputFile[nfiles];
  TChain chainDataMB("O2lambdaanalysis");
  for (Int_t i = 0; i < nfiles; i++){
      cout << "name " << name[i].c_str() << endl;
    inputFile[i] = TFile::Open(name[i].c_str());
    createChain(chainDataMB, name[i].c_str(), "O2lambdaanalysis");
  }
  auto originalDF = ROOT::RDataFrame(chainDataMB);
  //  RDataFrame originalDF(TreeName, inputFile);

  auto d1 = originalDF;

  TString weightCentFileName = "CentralityWeight_" + SinputFileNameCentWeight + ".root";
  TFile *weightCentFile = new TFile(weightCentFileName, "READ");
  TH3D *weights{weightCentFileName ? (TH3D *)weightCentFile->Get("hCentWeight") : nullptr};

  TString resoCentFileName = SinputFileNameResoWeight;
  TFile *resoFile = new TFile(resoCentFileName, "READ");
  TH1D *reso{resoCentFileName ? (TH1D *)resoFile->Get("hResoPerCentBins") : nullptr};
  
  auto h = d1.Histo1D("fPt");

  // invariant mass histograms
  auto hmass_Bef = d1.Histo1D({"mass_Lambda_Bef", "Invariant mass of p#pi", 100, 1.09, 1.14}, "fMassLambda");

  // apply charge selection
  string chargecut = "fSign >= 0";
  auto d2 = d1.Filter(chargecut);

  //apply eta selection
  //  if (isPartialEta)  d2 = d2.Filter("fEta > 0");

  // pt vs centrality before selections
  auto hPtvsCent_BefSel = d2.Histo2D({"PtvsCent_BefSel", "PtvsCent_BefSel", 100, 0, 100, 400, 0, 20}, "fCentFT0C", "fPt");

  auto histoBefV0Radius = d2.Histo1D({"histoBefV0Radius", "V0 Radius Distribution", 100, 0, 10}, "fV0Radius");
  auto histoBefDcaV0Daughters = d2.Histo1D({"histoBefDcaV0Daughters", "DCA V0 Daughters Distribution", 100, -2, 2}, "fDcaV0Daughters");
  auto histoBefV0CosPA = d2.Histo1D({"histoBefV0CosPA", "V0 CosPA Distribution", 100, 0.985, 1}, "fV0CosPA");
  auto histoBefDcaNegToPV = d2.Histo1D({"histoBefDcaNegToPV", "DCA Neg to PV Distribution", 100, -2, 2}, "fDcaNegToPV");
  auto histoBefDcaPosToPV = d2.Histo1D({"histoBefDcaPosToPV", "DCA Pos to PV Distribution", 100, -2, 2}, "fDcaPosToPV");

  // topological selections
  gRandom->SetSeed(0);
  //float RadiusCut = LowerlimitV0RadiusCut + (UpperlimitV0RadiusCut - LowerlimitV0RadiusCut) * gRandom->Rndm();
  //  cout << "RadiusCut " << RadiusCut << endl;

  string V0FilterString = Form("fV0Radius > %.2f", DefaultV0RadiusCut);
  if (isSysMultTrial) {
    if (isLoosest) V0FilterString = Form("fV0Radius > %.2f", LowerlimitV0RadiusCut);
    else if (isTightest)  V0FilterString =Form("fV0Radius > %.2f", UpperlimitV0RadiusCut);
    else V0FilterString = Form("fV0Radius > %.2f", LowerlimitV0RadiusCut + (UpperlimitV0RadiusCut - LowerlimitV0RadiusCut) * gRandom->Rndm());
  }

  string DcaV0DauString = Form("abs(fDcaV0Daughters) < %.2f", DefaultDcaV0DauCut);
  if (isSysMultTrial) {
    if (isLoosest) DcaV0DauString = Form("abs(fDcaV0Daughters) < %f", UpperlimitDcaV0DauCut);
    else if (isTightest) DcaV0DauString = Form("abs(fDcaV0Daughters) < %f", LowerlimitDcaV0DauCut);
    else DcaV0DauString = Form("abs(fDcaV0Daughters) < %f", LowerlimitDcaV0DauCut + (UpperlimitDcaV0DauCut - LowerlimitDcaV0DauCut) * gRandom->Rndm());
  }

  string V0CosPAString = Form("fV0CosPA > %.5f", DefaultV0CosPA);
  if (isSysMultTrial) {
    if (isLoosest) V0CosPAString = Form("fV0CosPA > %.5f", LowerlimitV0CosPA);
    else if (isTightest)  V0CosPAString = Form("fV0CosPA > %.5f", UpperlimitV0CosPA);
    else V0CosPAString = Form("fV0CosPA > %.5f", LowerlimitV0CosPA + (UpperlimitV0CosPA - LowerlimitV0CosPA) * gRandom->Rndm());
  }

  string DcaNegToPVString = Form("abs(fDcaNegToPV) > %.2f", DefaultDcaNegToPV); 
  if (isSysMultTrial) {
    if (isLoosest) DcaNegToPVString = Form("abs(fDcaNegToPV) > %f", LowerlimitDcaNegToPV);
    else if (isTightest) DcaNegToPVString =Form("abs(fDcaNegToPV) > %f", UpperlimitDcaNegToPV);
    else DcaNegToPVString = Form("abs(fDcaNegToPV) > %f", LowerlimitDcaNegToPV + (UpperlimitDcaNegToPV - LowerlimitDcaNegToPV) * gRandom->Rndm());
  }

  string DcaPosToPVString = Form("abs(fDcaPosToPV) > %.2f", DefaultDcaPosToPV); 
  if (isSysMultTrial) {
    if (isLoosest) DcaPosToPVString = Form("abs(fDcaPosToPV) > %f", LowerlimitDcaPosToPV);
    else if (isTightest) DcaPosToPVString =Form("abs(fDcaPosToPV) > %f", UpperlimitDcaPosToPV);
    else DcaPosToPVString = Form("abs(fDcaPosToPV) > %f", LowerlimitDcaPosToPV + (UpperlimitDcaPosToPV - LowerlimitDcaPosToPV) * gRandom->Rndm());
  }
  
  auto d2a = d2.Filter(V0FilterString);
  auto d2b = d2a.Filter(DcaV0DauString);
  auto d2c = d2b.Filter(V0CosPAString);
  auto d2d = d2c.Filter(DcaNegToPVString);
  auto d2e = d2d.Filter(DcaPosToPVString);

  //default cuts
  auto df_withCuts = d2.Define("cutV0Radius", [=]() { return DefaultV0RadiusCut; })
    .Define("cutDcaV0Daughters", [=]() { return DefaultDcaV0DauCut; });

  // now vary those thresholds
  auto df_varied = df_withCuts.Vary(
				    {"cutV0Radius", "cutDcaV0Daughters"},  // columns that will vary together

				    // Lambda generating 100 random variations within limits
				    [](float r, float d) {
				      const int NVAR = 2;

				      // Define the allowed limits for each variable
				      const float V0Radius_min = LowerlimitV0RadiusCut;
				      const float V0Radius_max = UpperlimitV0RadiusCut;
				      const float DcaV0Daughters_min = LowerlimitDcaV0DauCut;
				      const float DcaV0Daughters_max = UpperlimitDcaV0DauCut;

				      // Prepare the output: one RVec per column
				      RVec<RVecF> variations(2);
				      variations[0].reserve(NVAR);  // cutV0Radius
				      variations[1].reserve(NVAR);  // cutDcaV0Daughters

				      // Initialize random generator (fixed seed for reproducibility)
				      std::mt19937 gen(42);
				      std::uniform_real_distribution<float> distR(V0Radius_min, V0Radius_max);
				      std::uniform_real_distribution<float> distD(DcaV0Daughters_min, DcaV0Daughters_max);

				      // Generate random values for each variation
				      for (int i = 0; i < NVAR; ++i) {
					variations[0].push_back(distR(gen));
					variations[1].push_back(distD(gen));
				      }

				      return variations;
				    },

				    // Inputs to the Vary function (they can be the same as the varied columns)
				    {"cutV0Radius", "cutDcaV0Daughters"},

				    2,  // number of variations

				    "RandomV0Cuts"  // name of this variation set
				    );
  
  // apply the cuts using the **varied** thresholds
  auto df_selected = df_varied.Filter(
				      [](float V0Radius, float DcaV0Daughters, float cutR, float cutD){
					return V0Radius > cutR && DcaV0Daughters > cutD;
				      },
				      {"fV0Radius", "fDcaV0Daughters", "cutV0Radius", "cutDcaV0Daughters"}
				      );

  // histogram your invariant mass and pT
  auto hMassVsPt = df_selected.Histo2D(
				       {"hMassVsPt", "Mass vs pT;M_{Λ} (GeV/c^{2});p_{T} (GeV/c)", 100, 1.09, 1.14, 100, 0, 10},
				       "fMassLambda", "fPt"
				       );

  // extract all variations
  auto variationsMassVsPt = ROOT::RDF::Experimental::VariationsFor(hMassVsPt);

  auto histoV0Radius = d2e.Histo1D({"histoV0Radius", "V0 Radius Distribution", 100, 0, 10}, "fV0Radius");
  auto histoDcaV0Daughters = d2e.Histo1D({"histoDcaV0Daughters", "DCA V0 Daughters Distribution", 100, -2, 2}, "fDcaV0Daughters");
  auto histoV0CosPA = d2e.Histo1D({"histoV0CosPA", "V0 CosPA Distribution", 100, 0.985, 1}, "fV0CosPA");
  auto histoDcaNegToPV = d2e.Histo1D({"histoDcaNegToPV", "DCA Neg to PV Distribution", 100, -2, 2}, "fDcaNegToPV");
  auto histoDcaPosToPV = d2e.Histo1D({"histoDcaPosToPV", "DCA Pos to PV Distribution", 100, -2, 2}, "fDcaPosToPV");

  // pt vs centrality after selections
  auto hPtvsCent_AftSel = d2e.Histo2D({"PtvsCent_AftSel", "PtvsCent_AftSel", 100, 0, 100, 400, 0, 20}, "fCentFT0C", "fPt");

  // invariant mass histograms
  auto hmass = d2e.Histo1D({"mass_Lambda", "Invariant mass of p#pi", 100, 1.09, 1.14}, "fMassLambda");

  // invariant mass histograms vs pt
  auto hmassvsPt = d2e.Histo2D({"mass_LambdavPt", "Invariant mass of p#pi vs pT", 100, 1.09, 1.14, 100, 0, 10}, "fMassLambda", "fPt");

  // eta distributions
  //  auto heta = d2e.Histo1D({"eta", "Eta distribution of selected candidates", 200, -2, 2}, "fEta");

  // phi distributions
  auto hphi = d2e.Histo1D({"phi", "Phi distribution of selected candidates", 200, 0, 2 * TMath::Pi()}, "fPhi");

  // eta - phi distributions
  // auto hEtaPhi = d2e.Histo2D({"PhivsEta", "Phi vs Eta distribution of selected candidates", 100, -1, 1, 200, 0, 2 * TMath::Pi()}, "fEta", "fPhi");

  // create output file
  cout << "I am creating the output file " << endl;
  cout << "From the file: " << inputFileName << endl;
  TString OutputFileName = "OutputAnalysis/Output_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice];
  if (isApplyWeights)
    OutputFileName += "_Weighted";
  if (isApplyCentWeight)
    OutputFileName += "_CentWeighted";
  if (v2type == 1)
    OutputFileName += "_SP";
  if (isRun2Binning)
    OutputFileName += "_Run2Binning";
  if (!isRapiditySel){
    //    if (isPartialEta) OutputFileName += "_PositiveEta08";
    //else  OutputFileName += "_Eta08";
    OutputFileName += "_Eta08";
  }
  if (isSysMultTrial) {
    if (isLoosest)  OutputFileName += "_isLoosest";
    else if (isTightest)  OutputFileName += "_isTightest";
    else OutputFileName += Form("_SysMultTrial_%i", indexMultTrial);
  }
  if (isOOCentrality) OutputFileName += "_isOOCentrality";
  if (isApplyResoOnTheFly) OutputFileName += "_ResoOnTheFly";
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

  std::vector<ROOT::RDF::RResultPtr<TH2D>> massvsptVector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsV2CVector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsPzs2Vector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPsiVsPzVector;

  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsPzs2VectorWithAlpha; // decay parameter included in the calculation
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPsiVsPzVectorWithAlpha;  // decay parameter included in the calculation

  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPtVsCos2Vector;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> massVsPsiVsCos2Vector;

  // Containers to store histos and variations
  std::vector<decltype(ROOT::RDF::Experimental::VariationsFor(
							      std::declval<ROOT::RDF::RResultPtr<TH3D>>()))> h_variationsmassVsPtVsPzs2;
  
  df_selected = df_selected.Define("fPsiDiff", "if ((fPhi-fPsiT0C) < 0) return (fPhi-fPsiT0C+(float)TMath::Pi()); else if ((fPhi-fPsiT0C) > 2* TMath::Pi()) return (fPhi-fPsiT0C-2*(float)TMath::Pi()); else if ((fPhi-fPsiT0C) > TMath::Pi()) return (fPhi-fPsiT0C-(float)TMath::Pi()); else return (fPhi-fPsiT0C);");
  df_selected = df_selected.Define("f2PsiDiffCorr", "2*fPsiDiff");
  df_selected = df_selected.Define("f2PsiDiff", "2*fPhi-2*fPsiT0C");

  double centWeight[100];
  double resoWeight[100];
  for (Int_t cent = 0; cent < 100; cent++){
    centWeight[cent]= 1./weights->GetBinContent(cent+1);
    resoWeight[cent]= reso->GetBinContent(cent+1);
  }
  
  df_selected = df_selected.DefineSlot("fCentWeight", [&centWeight](unsigned int, float cent) {
    if (cent < 0. || cent > 100.) return 0.;
    else return centWeight[static_cast<int>(cent)];
  }, {"fCentFT0C"});
  //  df_selected.Display({"fCentFT0C", "fCentWeight"}, 128)->Print();
  //  return;

  df_selected = df_selected.DefineSlot("fResoWeight", [&resoWeight](unsigned int, float cent) {
    if (cent < 0. || cent > 100.) return 0.;
    else return resoWeight[static_cast<int>(cent)];
  }, {"fCentFT0C"});
  //  df_selected.Display({"fCentFT0C", "fResoWeight"}, 128)->Print();

  string SPzs2LambdaFinal = "fPzs2Lambda";
  if (isApplyResoOnTheFly) SPzs2LambdaFinal = "fPzs2Lambda/fResoWeight";
  //  df_selected = df_selected.Define("fCentResoWeight", Sweight);
  df_selected = df_selected.Define("fPzs2LambdaFinal", SPzs2LambdaFinal);
  //  df_selected.Display({"fCentFT0C", "fCentResoWeight"}, 512)->Print();

  cout << "I am looping over all centrality classes " << endl;
  for (Int_t cent = 0; cent < numCentLambdaOO + 1; cent++)
  {
    if (cent == numCentLambdaOO)
    { // 0-100%
      CentFT0CMin = 0;
      CentFT0CMax = 100;
    }
    else
    {
      CentFT0CMin = CentFT0CLambdaOO[cent];
      CentFT0CMax = CentFT0CLambdaOO[cent + 1];
    }

    if (ChosenPart == 6) {
      MinPzs2 = -7;
      MaxPzs2 = 7;
      if (isApplyResoOnTheFly){
	MinPzs2 = MinPzs2Reso[cent];
	MaxPzs2 = MaxPzs2Reso[cent];
      }
    }

    //    cout << "cent: " << cent << " min: " << CentFT0CMin << " max: " << CentFT0CMax << endl;
    auto dcent = df_selected.Filter(Form("fCentFT0C>=%.1f && fCentFT0C<%.1f", CentFT0CMin + 0.001, CentFT0CMax - 0.001));
    auto df_selected_cent = df_selected.Filter(Form("fCentFT0C>=%.1f && fCentFT0C<%.1f", CentFT0CMin + 0.001, CentFT0CMax - 0.001));
    auto dmasscut = dcent.Filter("fMassLambda > 1.1 && fMassLambda < 1.13");

    // histogram your invariant mass and pT
    auto massVsPtVsPzs2 = df_selected_cent.Histo3D({Form("massVsPtVsPzs2_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Pzs2", 80, 1.09, 1.14, 100, 0, 10, NPzs2, MinPzs2, MaxPzs2}, "fMassLambda", "fPt", "fPzs2LambdaFinal", "fCentWeight");
    auto variationsmassVsPtVsPzs2 = ROOT::RDF::Experimental::VariationsFor(massVsPtVsPzs2);
    massVsPtVsPzs2Vector.push_back(massVsPtVsPzs2);
    h_variationsmassVsPtVsPzs2.push_back(variationsmassVsPtVsPzs2);
    
    auto massvspt2D = dcent.Histo2D({Form("MassVsPt_cent%i-%i", CentFT0CMin, CentFT0CMax), "massvspt", 80, 1.09, 1.14, 100, 0, 10,}, "fMassLambda", "fPt");
    massvsptVector.push_back(massvspt2D);
    
    auto v2C = dcent.Histo1D({Form("v2C_cent%i-%i", CentFT0CMin, CentFT0CMax), "v2C", Nv2, Minv2, Maxv2}, v2Chosen);
    v2CVector.push_back(v2C);
    auto hPhiCent = dmasscut.Histo2D({Form("PhiHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Phi vs pt", 100, 0, 10, 100, 0, 2 * TMath::Pi()}, "fPt", "fPhi");
    hPhiCentVector.push_back(hPhiCent);
    auto hPsiCent = dmasscut.Histo1D({Form("PsiHist_cent%i-%i", CentFT0CMin, CentFT0CMax), "Psi", 100, -2 * TMath::Pi(), 2 * TMath::Pi()}, "fPsiT0C");
    hPsiCentVector.push_back(hPsiCent);

    auto PsiDiff = dcent.Histo1D({Form("2PsiDiffCorr_cent%i-%i", CentFT0CMin, CentFT0CMax), "2PsiDiffCorr", 100, 0, 2 * TMath::Pi()}, "f2PsiDiffCorr");
    PsiDiffVector.push_back(PsiDiff);
    auto PsiDiff2 = dcent.Histo1D({Form("2PsiDiff_cent%i-%i", CentFT0CMin, CentFT0CMax), "2PsiDiff", 100, -2 * TMath::Pi(), 2 * TMath::Pi()}, "f2PsiDiff");
    PsiDiff2Vector.push_back(PsiDiff2);

    auto hMassCut = dmasscut.Histo1D({Form("massCut_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass of #Lambda#pi", 100, 1.28, 1.36}, "fMassLambda");
    hMassCutVector.push_back(hMassCut);

    auto massVsPtVsV2C = dcent.Histo3D({Form("massVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs V2C", 80, 1.09, 1.14, 100, 0, 10, Nv2, Minv2, Maxv2}, "fMassLambda", "fPt", v2Chosen);
    massVsPtVsV2CVector.push_back(massVsPtVsV2C);

    auto massVsPsiVsPz = dcent.Histo3D({Form("massVsPsiVsPz_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Pz", 80, 1.09, 1.14, 20, 0, 2 * TMath::Pi(), NPz, MinPz, MaxPz}, "fMassLambda", "f2PsiDiffCorr", "fCosThetaLambda");
    massVsPsiVsPzVector.push_back(massVsPsiVsPz);

    auto massVsPtVsCos2 = dcent.Histo3D({Form("massVsPtVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs Pt vs Cos2", 80, 1.09, 1.14, 100, 0, 10, 100, 0, 1}, "fMassLambda", "fPt", "fCos2ThetaLambda");
    massVsPtVsCos2Vector.push_back(massVsPtVsCos2);
    auto massVsPsiVsCos2 = dcent.Histo3D({Form("massVsPsiVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax), "Invariant mass vs 2*(Psi-Phi) vs Cos2", 80, 1.09, 1.14, 20, 0, 2 * TMath::Pi(), 100, 0, 1}, "fMassLambda", "f2PsiDiffCorr", "fCos2ThetaLambda");
    massVsPsiVsCos2Vector.push_back(massVsPsiVsCos2);
  }

  // draw histograms
  TCanvas *cMass = new TCanvas("cMass", "cMass", 900, 600);
  TString TitleX = "M_{#Lambda#pi}";
  StyleCanvas(cMass, 0.1, 0.1, 0.03, 0.1);
  StyleHisto(*hmass_Bef, 0, 1.2 * hmass_Bef->GetMaximum(), kRed, 20, "M_{#Lambda#pi}", "Counts", "", kTRUE, 1.2, 1.4, 1.2, 1.2, 0.7);
  StyleHisto(*hmass, 0, 1.2 * hmass_Bef->GetMaximum(), kBlue, 20, "M_{#Lambda#pi}", "Counts", "", kTRUE, 1.2, 1.4, 1.2, 1.2, 0.7);
  hmass_Bef->Draw("E");
  hmass->Draw("E SAME");
  h->Write();
  hPtvsCent_BefSel->Write();
  hPtvsCent_AftSel->Write();
  cMass->Write();
  hmass_Bef->Write();
  hmass->Write();
  hmassvsPt->Write();
  hphi->Write();
  histoV0Radius->Write();
  histoDcaV0Daughters->Write();
  histoV0CosPA->Write();
  histoDcaNegToPV->Write();
  histoDcaPosToPV->Write();
  histoBefV0Radius->Write();
  histoBefDcaV0Daughters->Write();
  histoBefV0CosPA->Write();
  histoBefDcaNegToPV->Write();
  histoBefDcaPosToPV->Write();

  for (Int_t cent = 0; cent < numCentLambdaOO + 1; cent++)
  {

    auto &variationMap = h_variationsmassVsPtVsPzs2[cent];

    // Retrieve the variation tags
    auto tags = variationMap.GetKeys();

    for (auto &tag : tags) {
      std::cout << "Centrality bin " << cent
		<< " | Variation tag Pzs2: " << tag << std::endl;

      auto histo = variationMap[tag];
      TString histName = Form("massVsPtVsPzs2_cent%i_%s", cent, tag.c_str());
      histo.Write(histName);  // or histo->Draw()
    }

    massvsptVector[cent]->Write();
    massVsPtVsV2CVector[cent]->Write();
    massVsPtVsPzs2Vector[cent]->Write();
    massVsPsiVsPzVector[cent]->Write();
    massVsPtVsCos2Vector[cent]->Write();
    massVsPsiVsCos2Vector[cent]->Write();
    PsiDiffVector[cent]->Write();
    PsiDiff2Vector[cent]->Write();
    hPsiCentVector[cent]->Write();
    hPhiCentVector[cent]->Write();
  }

  file->Close();

  cout << "I created the file " << file->GetName() << endl;
}
