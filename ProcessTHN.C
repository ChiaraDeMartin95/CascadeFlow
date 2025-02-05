// This macro was originally written by:
// chiara.de.martin@cern.ch
#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
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

void ProcessTHN(Bool_t isEff = 0,
                Int_t indexMultTrial = 0,
                Int_t ChosenPart = ChosenParticle,
                TString inputFileName = SinputFileName,
                Int_t EtaSysChoice = ExtrEtaSysChoice,
                Bool_t isSysMultTrial = ExtrisSysMultTrial)
{

  // phi and efficiency weights should be applied in task

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

  inputFileName = "_TestTHN.root";
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
  THnSparseF *hV2 = (THnSparseF *)dir1->Get("h" + ParticleName[0] + "V2");
  //**********/

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;

  // acceptance correction ***
  TH3D *hmassVsPtVsCos2Theta[numCent + 1];
  TH3D *hmassVsPtVsCos2ThetaLambdaFromC[numCent + 1];
  TH3D *hmassVsPsiVsCos2Theta[numCent + 1];
  TH3D *hmassVsPsiVsCos2ThetaLambdaFromC[numCent + 1];
  TString hNameCos2Theta_3D[numCent + 1] = {""};
  TString hNameCos2ThetaLambdaFromC_3D[numCent + 1] = {""};
  TString hNameCos2ThetaVsPsi_3D[numCent + 1] = {""};
  TString hNameCos2ThetaVsPsiLambdaFromC_3D[numCent + 1] = {""};

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
  TString profName[numCent + 1] = {""};

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

    // THESE are the needed histograms
    hName[cent] = Form("massVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax);
    if (ExtrisApplyEffWeights)
      hName[cent] = Form("massVsPtVsV2CWeighted_cent%i-%i", CentFT0CMin, CentFT0CMax);
    profName[cent] = Form("ProfilemassVsPtVsV2C_cent%i-%i", CentFT0CMin, CentFT0CMax);

    hNamePzs2_3D[cent] = Form("massVsPtVsPzs2_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNamePzs2LambdaFromC_3D[cent] = Form("massVsPtVsPzs2LambdaFromC_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNamePzVsPsi_3D[cent] = Form("massVsPsiVsPz_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNamePzVsPsiLambdaFromC_3D[cent] = Form("massVsPsiVsPzLambdaFromC_WithAlpha_cent%i-%i", CentFT0CMin, CentFT0CMax);

    hNameCos2Theta_3D[cent] = Form("massVsPtVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNameCos2ThetaLambdaFromC_3D[cent] = Form("massVsPtVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNameCos2ThetaVsPsi_3D[cent] = Form("massVsPsiVsCos2_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hNameCos2ThetaVsPsiLambdaFromC_3D[cent] = Form("massVsPsiVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax);

    //->Projection(1,2,3);
  }

  // create output file
  TString SBDT = "";
  if (BDTscoreCut != DefaultBDTscoreCut)
    SBDT = Form("_BDT%.3f", BDTscoreCut);
  TString OutputFileName = "OutputAnalysis/OutputFromTHN_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (isApplyWeights)
    OutputFileName += "_Weighted";
  if (v2type == 1)
    OutputFileName += "_SP";
  if (!useCommonBDTValue)
    OutputFileName += "_BDTCentDep";
  if (isRun2Binning)
    OutputFileName += "_Run2Binning";
  OutputFileName += ".root";
  TFile *file = new TFile(OutputFileName, "RECREATE");

  // h->Write();
  // hPtvsCent_Bef->Write();
  // hPtvsCent_RapCut_Bef->Write();
  // hPtvsCent_Aft->Write();
  // hPtvsCent_RapCut_Aft->Write();
  // cMass->Write();
  // MassXi->Write();
  // hmass_Bef->Write();
  // BDT_response_Bef->Write();
  // BDT_response_Bef2->Write();
  // mass_vs_BDTResponse->Write();
  // hmass->Write();
  // hmassvsPt->Write();
  // heta->Write();
  // hphi->Write();
  // hEtaPhi->Write();
  // BDT_response->Write();

  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    // hrapidityVector[cent]->Write();
    // hEtaVector[cent]->Write();
    // PsiDiffVector[cent]->Write();
    // PsiDiff2Vector[cent]->Write();
    //
    // hMassCutVector[cent]->Write();
    // v2CVector[cent]->Write();
    // hPhiCentVector[cent]->Write();
    // hPsiCentVector[cent]->Write();
    //
    // massVsPtVsV2CVector[cent]->Write();
    // massVsPtVsV2CWeightedVector[cent]->Write();
    //
    // massVsPtVsPzs2Vector[cent]->Write();
    // massVsPtVsPzs2LambdaFromCVector[cent]->Write();
    // massVsPsiVsPzVector[cent]->Write();
    // massVsPsiVsPzLambdaFromCVector[cent]->Write();
    //
    // massVsPtVsPzs2VectorWithAlpha[cent]->Write();
    // massVsPtVsPzs2LambdaFromCVectorWithAlpha[cent]->Write();
    // massVsPsiVsPzVectorWithAlpha[cent]->Write();
    // massVsPsiVsPzLambdaFromCVectorWithAlpha[cent]->Write();
    //
    // massVsPtVsCos2Vector[cent]->Write();
    // massVsPtVsCos2LambdaFromCVector[cent]->Write();
    // massVsPsiVsCos2Vector[cent]->Write();
    // massVsPsiVsCos2LambdaFromCVector[cent]->Write();
    //
    // profileVector[cent]->Write();
  }

  file->Close();

  cout << "I created the file " << file->GetName() << endl;
}
