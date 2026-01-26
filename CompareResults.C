#include "Riostream.h"
#include "TTimer.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TH3F.h>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TProfile.h"
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TFile.h>
#include <TLine.h>
#include <TSpline.h>
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
#include "CommonVar.h"
#include "ErrRatioCorr.C"

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(1.5);
  histo->GetXaxis()->SetTitle(TitleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(TitleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.3);
  histo->SetTitle(title);
}

void StyleHistoYield(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(TitleX);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(TitleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset); // 1.2
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->SetTitle(title);
}

void SetFont(TH1F *histo)
{
  histo->GetXaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetYaxis()->SetLabelFont(43);
}
void SetTickLength(TH1F *histo, Float_t TickLengthX, Float_t TickLengthY)
{
  histo->GetXaxis()->SetTickLength(TickLengthX);
  histo->GetYaxis()->SetTickLength(TickLengthY);
}

void SetHistoTextSize(TH1F *histo, Float_t XSize, Float_t XLabelSize, Float_t XOffset, Float_t XLabelOffset, Float_t YSize, Float_t YLabelSize, Float_t YOffset, Float_t YLabelOffset)
{
  histo->GetXaxis()->SetTitleSize(XSize);
  histo->GetXaxis()->SetLabelSize(XLabelSize);
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetXaxis()->SetLabelOffset(XLabelOffset);
  histo->GetYaxis()->SetTitleSize(YSize);
  histo->GetYaxis()->SetLabelSize(YLabelSize);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetYaxis()->SetLabelOffset(YLabelOffset);
}

void StyleCanvas(TCanvas *canvas, Float_t TopMargin, Float_t BottomMargin, Float_t LeftMargin, Float_t RightMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  gPad->SetTopMargin(TopMargin);
  gPad->SetLeftMargin(LeftMargin);
  gPad->SetBottomMargin(BottomMargin);
  gPad->SetRightMargin(RightMargin);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
}

void StylePad(TPad *pad, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  pad->SetFillColor(0);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetLeftMargin(LMargin);
  pad->SetRightMargin(RMargin);
  pad->SetTopMargin(TMargin);
  pad->SetBottomMargin(BMargin);
}

Float_t YLow = 0;
Float_t YUp = 0;
Float_t YLowRatio = 0;
Float_t YUpRatio = 0;

TString hTitleX = "";
TString hTitleY = "";

Float_t xTitle = 15;
Float_t xOffset = 4;
Float_t yTitle = 30;
Float_t yOffset = 2;

Float_t xLabel = 30;
Float_t yLabel = 30;
Float_t xLabelOffset = 0.05;
Float_t yLabelOffset = 0.01;

Float_t tickX = 0.03;
Float_t tickY = 0.042;

Float_t LimSupMultRatio = 5.1;
Float_t LimInfMultRatio = 1e-2;
Float_t YoffsetSpectraRatio = 1.1;
Float_t xTitleR = 35;
Float_t xOffsetR = 1;
Float_t yTitleR = 30;
Float_t yOffsetR = 2;

Float_t xLabelR = 25;
Float_t yLabelR = 25;
Float_t xLabelOffsetR = 0.02;
Float_t yLabelOffsetR = 0.04;

TString Sinputfile = "";
TString namehisto[100] = {""};
TString CommonFileName = "";

Int_t numOptions = 0;
TH1F *hDef;
TString fileName[100] = {""};
TH1F *h[100];
TH1F *hRatio[100];
TString sleg[100] = {""};
TString Smolt[numCent + 1];
const Int_t numOptionsConst = 91;
Int_t SignRun[numOptionsConst] = {0};
Int_t nrun[numOptionsConst] = {545312, 545296, 545295, 545294, 545291, 545249, 545223, 545210, 545185, 545066, 545064, 545063, 545047, 545009,
                               544992, 544968, 544964, 544963, 544961, 544947, 544931, 544917, 544914, 544913, 544896, 544887, 544886, 544868, 544813, 544797, 544795,
                               544794, 544767, 544754, 544742, 544739, 544696, 544694, 544693, 544692, 544674, 544672, 544653, 544652, 544640, 544614, 544585, 544583,
                               544582, 544580, 544568, 544567, 544565, 544564, 544551, 544550, 544549, 544548, 544518, 544515, 544514, 544512, 544511, 544510, 544508,
                               544492, 544491, 544490, 544477, 544476, 544475, 544474, 544454, 544451, 544392, 544391, 544390, 544389, 544185, 544184, 544124, 544123,
                               544122, 544121, 544116, 544098, 544095, 544091, 544032, 544028, 544013};

void CompareResults(Int_t TypeComp = 0,
                    Int_t mult = 0,
                    Bool_t isPtAnalysis = 1,
                    Int_t ChosenPart = ChosenParticle,
                    Int_t BkgType = ExtrBkgType,
                    Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  Bool_t isRatio = 1;    // 1 for ratio, 0 for difference
  Bool_t isFullCorr = 0; // 1 for full correlation, 0 for partial correlation
  Bool_t isFitRatio = 0;
  Bool_t isStoreSyst = 0;
  TString TypeSyst = "";

  gStyle->SetOptStat(0);

  if (mult > numCent)
  {
    cout << "Multiplciity out of range" << endl;
    return;
  }

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  if (mult == numCent)
  { // 0-80%
    CentFT0CMin = 0;
    CentFT0CMax = 80;
  }
  else
  {
    CentFT0CMin = CentFT0C[mult];
    CentFT0CMax = CentFT0C[mult + 1];
  }

  Float_t MinHistoX = MinPt[ChosenPart];
  Float_t MaxHistoX = MaxPt[ChosenPart];

  // TypeComp = 0 --> weighted vs unweighted v2
  // TypeComp = 1 --> LHC23 vs LHC23zzh
  // TypeComp = 2 --> SP vs EP
  // TypeComp = 3 --> BDT cut 0.7 vs default
  // TypeComp = 4 --> v2 from fit vs v2 from histo
  // TypeComp = 5 --> resolution LF vs CFW
  // TypeComp = 6 --> v2 LF vs CFW
  // TypeComp = 7 --> v2 with and without occupancy event selection
  // TypeComp = 8 --> v2 pass4 vs pass3 (LHCzzh)
  // TypeComp = 9 --> v2 with pass4; test 5 vs test 3
  // TypeComp = 10 --> Pz from Xi vs Pz from daughter Lambda
  // TypeComp = 11 --> Pzs2 vs centrality : from Xi vs from daughter Lambda
  // TypeComp = 12 --> Pzs2 vs centrality: XiMinus vs XiPlus vs Xi
  // TypeComp = 13 --> weighted vs unweighted v2 (eff weight)
  // TypeComp = 14 --> default v2 vs corrected v2
  // TypeComp = 15 --> weighted vs corrected v2
  // TypeComp = 16 --> Pz from Xi - offline vs online acceptance
  // TypeComp = 17 --> Pz from Lambda - offline vs online acceptance
  // TypeComp = 18 --> Pz from Lambda - online vs pt vs online vs eta and pt acceptance
  // TypeComp = 19 --> Pz from Lambda vs from Xi
  // TypeComp = 20 --> Pz from Lambda fit vs no fit
  // TypeComp = 21 --> Proton acceptance 2023 pass5 vs 2023 pass4 vs eta
  // TypeComp = 22 --> Proton acceptance 2023 pass5 vs 2023 pass4 vs pt
  // TypeComp = 24, 25, 26, 27, 28 --> 2023 pass4 vs 2023 pass5 (purity, sigma, yield, Pzs2, Pzs2Error  integrated in pT vs centrality)
  // TypeComp = 29 --> Lambda Pzs2 from my code vs Preliminary one from Junlee
  // TypeComp == 31 --> Lambda Pzs2 from Tree vs THN
  // TypeComp == 32 --> Compare purity with loosest and tightest topo sel
  // TypeComp == 33 --> Compare Pzs2 with loosest and tightest topo sel
  // TypeComp == 34 --> Compare Pzs2 of Lambda with resolution applied offline or on the fly
  // TypeComp == 35 --> Resolution comparison (T0C vs T0M vs V0A vs T0A)
  // TypeComp == 36 --> 2023 pass5 with pass5 acceptance vs w/ pass4 acceptance (Pzs2, Pzs2Error  integrated in pT vs centrality)
  // TypeComp == 37 --> Compare Pzs2 of Xi with flat and non-flat event plane
  // TypeComp == 38 --> Compare Pzs2 of Lambda with different PV Z vertex selection
  // TypeComp == 39 --> Resolution comparison (T0C reso with different reference detectors)
  // TypeComp == 40 --> Resolution comparison with corrected TPC (T0C and T0A reso with different reference detectors)
  // TypeComp == 41 --> Compare Pzs2 of Lambda from THN and from Tree
  // TypeComp == 42 --> Compare Pzs2 of Lambda from Tree with less and more pt bins
  // TypeComp == 43 --> Compare Pzs2 of Lambda from two different trains
  // TypeComp == 44 --> Compare T0C event plane resolution in PbPb (normalised SP vs EP)
  // TypeComp == 45 --> Compare Pzs errors (normalised SP vs EP)
  // TypeComp == 46 --> Xi Pzs2: preliminaries vs paper proposal results
  // TypeComp == 47 --> Xi Pzs2: from fit vs assuming zero background polarization
  // TypeComp == 48 --> Compare acceptance of Lambda from different runs
  // TypeComp == 49 --> Compare resolution in PbPb with and withou occupancy cut
  // TypeComp == 50 --> Xi Pzs2: occupancy dependence in PbPb
  // TypeComp == 51 --> Resolution in OO with pt, tracks < 2 GeV/c: full sample vs partial sample
  // TypeComp == 52 --> Compare Pzs2 of Lambda with resolution from different reference detectors
  // TypeComp == 53 --> Compare Pzs2 of Lambda with Pzs,bkg = 0 vs default
  // TypeComp == 54 --> Compare Pzs2 of Lambda with different fit ranges in Pz
  // TypeComp == 55 --> Compare Pzs2 of Lambda with different fit function for inv mass bkg
  // TypeComp == 56 --> Compare different acceptance selections
  // TypeComp == 57 --> Compare purity with different selections for acceptance
  // TypeComp == 58 --> Compare Pzs2 of Lambda with |eta| < 0.8 vs |eta| < 100
  // TypeComp == 59 --> Compare Pzs2 of Lambda from THN and from tree (reso and cent weights applied)
  // TypeComp == 60 --> Compare Pzs2 of Lambda with |z| < 10 cm and |z| < 8 cm
  // TypeComp == 61 --> Compare Pzs2 of Lambda with new and old acceptance (new = higher purity)
  // TypeComp == 62 --> Compare Pzs2 of Lambda with FD correction vs without FD correction
  // TypeComp == 63 --> Compare proton acceptance with and without |eta| < 0.8 cut for daughters
  // TypeComp == 64 --> Compare proton acceptance with |eta| < 0.8 and 0 < eta < 0.8 and -0.8 < eta < 0

  // TypeComp = 0 --> weighted vs unweighted v2
  if (TypeComp == 0)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_" + SinputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[ExtrEtaSysChoice];
    CommonFileName += IsOneOrTwoGauss[ExtrUseTwoGauss];
    CommonFileName += SIsBkgParab[ExtrBkgType];
    CommonFileName += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] = "";
    fileName[1] = "_Weighted";
    namehisto[0] = "histoV2";
    namehisto[1] = "histoV2";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "Default";
    sleg[1] = "Weighted";
  }
  // TypeComp = 1 --> LHC23 vs LHC23zzh
  else if (TypeComp == 1)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23_PbPb_pass3_Train218607_";
    fileName[1] = "LHC23zzh_pass3_Train224930_";
    fileName[0] += ParticleName[ChosenPart] + SEtaSysChoice[ExtrEtaSysChoice];
    fileName[0] += IsOneOrTwoGauss[ExtrUseTwoGauss];
    fileName[0] += SIsBkgParab[ExtrBkgType];
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] += ParticleName[ChosenPart] + SEtaSysChoice[ExtrEtaSysChoice];
    fileName[1] += IsOneOrTwoGauss[ExtrUseTwoGauss];
    fileName[1] += SIsBkgParab[ExtrBkgType];
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";

    namehisto[0] = "histoV2NoFit";
    namehisto[1] = "histoV2NoFit";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "LHC23_PbPb_pass3";
    sleg[1] = "LHC23zzh_pass3";
  }
  // TypeComp = 2 --> SP vs EP
  else if (TypeComp == 2)
  {
    numOptions = 2;
    // CommonFileName = "OutputAnalysis/FitV2_" + SinputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[ExtrEtaSysChoice];
    // CommonFileName += IsOneOrTwoGauss[ExtrUseTwoGauss];
    // CommonFileName += SIsBkgParab[ExtrBkgType];
    // CommonFileName += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    // CommonFileName += "_Weighted";
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23zzh_pass3_Train226234_CFW_Xi_BkgParab";
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] = fileName[0];
    fileName[0] += "";
    fileName[1] += "_SP";
    fileName[0] += "_Run2Binning";
    fileName[1] += "_Run2Binning";
    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "EP";
    sleg[1] = "SP";
  }
  // TypeComp = 3 --> BDT cut 0.7 vs default
  else if (TypeComp == 3)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23_PbPb_pass3_Train218607_";
    fileName[1] = "LHC23zzh_pass3_Train225737_BDT0.7_";
    fileName[0] += ParticleName[ChosenPart] + SEtaSysChoice[ExtrEtaSysChoice];
    fileName[0] += IsOneOrTwoGauss[ExtrUseTwoGauss];
    fileName[0] += SIsBkgParab[ExtrBkgType];
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted_BDTCentDep";
    fileName[1] += ParticleName[ChosenPart] + SEtaSysChoice[ExtrEtaSysChoice];
    fileName[1] += IsOneOrTwoGauss[ExtrUseTwoGauss];
    fileName[1] += SIsBkgParab[ExtrBkgType];
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";

    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "";
    sleg[1] = "BDTcut > 0.97";
  }
  // v2 from fit vs v2 from histo
  else if (TypeComp == 4)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_" + SinputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[ExtrEtaSysChoice];
    CommonFileName += IsOneOrTwoGauss[ExtrUseTwoGauss];
    CommonFileName += SIsBkgParab[ExtrBkgType];
    CommonFileName += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    CommonFileName += "_Weighted";
    fileName[0] = "";
    fileName[0] = "";
    namehisto[0] = "histoV2";
    namehisto[1] = "histoV2NoFit";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "v2 from fit";
    sleg[1] = "v2 no fit";
  }
  // resolution LF vs CFW
  else if (TypeComp == 5)
  {
    numOptions = 2;
    fileName[0] = "Resolution/Resolution_EP_LF";
    fileName[1] = "Resolution/Resolution_EP_CFW";
    namehisto[0] = "hReso";
    namehisto[1] = "hReso";
    hTitleY = "Resolution";
    hTitleX = "Centrality (%)";
    YLow = 0;
    YUp = 1.8;
    YLowRatio = 0.6;
    YUpRatio = 1.2;
    sleg[0] = "LF";
    sleg[1] = "central FW";
  }
  // v2 LF vs CFW
  else if (TypeComp == 6)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23zzh_pass3_Train226234_CFW_Xi_BkgParab";
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] = "LHC23zzh_pass3_Train224930_Xi_BkgParab";
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";
    // fileName[0] += "_SP";
    // fileName[1] += "_SP";
    fileName[0] += "_Run2Binning";
    fileName[1] += "_Run2Binning";
    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "CFW";
    sleg[1] = "LF";
  }
  // v2 with and without occupancy event selection
  else if (TypeComp == 7)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23_PbPb_pass3_Train218607_Xi_BkgParab";
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] = "LHC23_PbPb_pass3_Train231308_Xi_BkgParab";
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";

    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "w/o occupancy sel.";
    sleg[1] = "w/ occupancy sel.";
  }
  // v2 with pass4 vs pass3
  else if (TypeComp == 8)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23zzh_pass3_Train224930_Xi_BkgParab";
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] = "LHC23zzh_pass4_test3_Train232412_Xi_BkgParab";
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";
    fileName[0] += "_Run2Binning";
    fileName[1] += "_Run2Binning";
    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.;
    YUpRatio = 2.;
    sleg[0] = "pass 3";
    sleg[1] = "pass 4";
  }
  // v2 with pass4; test 5 vs test 3
  else if (TypeComp == 9)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_";
    fileName[0] = "LHC23zzh_pass4_test3_Train232412_Xi_BkgParab";
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_Weighted";
    fileName[1] = "LHC23zzh_pass4_test5_Train235645_Xi_BkgParab";
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += "_Weighted";
    fileName[0] += "_Run2Binning";
    fileName[1] += "_Run2Binning";
    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2Mixed";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.;
    YUpRatio = 2.;
    sleg[0] = "test 3";
    sleg[1] = "test 5";
  }
  // TypeComp = 10 --> Pz from Xi vs Pz from daughter Lambda
  else if (TypeComp == 10)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitPzs2_" + SinputFileName + "_" + ParticleName[ChosenPart];
    CommonFileName += IsOneOrTwoGauss[UseTwoGauss];
    CommonFileName += SIsBkgParab[BkgType];
    Smolt[mult] = Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);
    CommonFileName += Smolt[mult];
    if (isApplyWeights)
      CommonFileName += "_Weighted";
    if (!useCommonBDTValue)
      CommonFileName += "_BDTCentDep";
    if (isRun2Binning)
      CommonFileName += "_Run2Binning";
    if (!isPtAnalysis)
      CommonFileName += "_vsPsi";
    fileName[0] = "";
    fileName[1] = "_PolFromLambda";
    namehisto[0] = "histoPzs2";
    namehisto[1] = "histoPzs2LambdaFromC";
    hTitleY = "P_{z}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "P_{z} from #Xi";
    sleg[1] = "P_{z} from #Lambda";
  }
  // TypeComp = 11 --> Pzs2 vs centrality : from Xi vs from daughter Lambda
  else if (TypeComp == 11)
  {
    numOptions = 2;
    isFullCorr = 0;
    CommonFileName = "Pzs2VsCentrality/Pzs2_" + SinputFileName + "_" + ParticleName[ChosenPart];
    CommonFileName += SIsBkgParab[BkgType] + "_Pzs2";
    if (isApplyWeights)
      CommonFileName += "_Weighted";
    if (!useCommonBDTValue)
      CommonFileName += "_BDTCentDep";
    if (isRun2Binning)
      CommonFileName += "_Run2Binning";
    fileName[0] = "";
    fileName[1] = "_PolFromLambda";
    fileName[0] += "_PtInt";
    fileName[1] += "_PtInt";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    MinHistoX = 0;
    MaxHistoX = 80;
    hTitleY = "P_{z,s2}";
    hTitleX = "FT0C Centrality (%)";
    YLow = -0.01;
    YUp = 0.03;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "P_{z,s2} from #Xi";
    sleg[1] = "P_{z,s2} from #Lambda";
  }
  // TypeComp = 12 --> Pzs2 vs centrality : XiPlus vs XiMinus vs Xi
  else if (TypeComp == 12)
  {
    numOptions = 3;
    CommonFileName = "Pzs2VsCentrality/Pzs2_" + SinputFileName + "_";
    // CommonFileName = "Pzs2VsCentrality/Pzs2_LHC23_pass4_QC1_Train268802_";
    fileName[0] = ParticleName[0] + SIsBkgParab[BkgType] + "_Pzs2_Weighted_PolFromLambda_PtInt";
    fileName[1] = ParticleName[2] + SIsBkgParab[BkgType] + "_Pzs2_Weighted_PolFromLambda_PtInt";
    fileName[2] = ParticleName[3] + SIsBkgParab[BkgType] + "_Pzs2_Weighted_PolFromLambda_PtInt";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    namehisto[2] = "fHistPzs";
    MinHistoX = 0;
    MaxHistoX = 80;
    hTitleY = "P_{z,s2}";
    hTitleX = "FT0C Centrality";
    YLow = -0.01;
    YUp = 0.03;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "#Xi";
    sleg[1] = "#Xi^{-}";
    sleg[2] = "#Xi^{+}";
  }
  else if (TypeComp == 13)
  {
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_" + SinputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[ExtrEtaSysChoice];
    CommonFileName += IsOneOrTwoGauss[ExtrUseTwoGauss];
    CommonFileName += SIsBkgParab[ExtrBkgType];
    CommonFileName += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    CommonFileName += "_Weighted_BDTCentDep";
    fileName[0] = "";
    fileName[1] = "_EffW";
    namehisto[0] = "histoV2";
    namehisto[1] = "histoV2";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "Default";
    sleg[1] = "Weighted";
  }
  else if (TypeComp == 14)
  {
    // TypeComp = 14 --> default vs corrected v2
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_" + SinputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[ExtrEtaSysChoice];
    CommonFileName += IsOneOrTwoGauss[ExtrUseTwoGauss];
    CommonFileName += SIsBkgParab[ExtrBkgType];
    CommonFileName += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    CommonFileName += "_Weighted_BDTCentDep";
    fileName[0] = "";
    fileName[1] = "";
    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2MixedCorr";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "Default";
    sleg[1] = "Corrected";
  }
  else if (TypeComp == 15)
  {
    // TypeComp = 15 --> weighted vs corrected v2
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitV2_" + SinputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[ExtrEtaSysChoice];
    CommonFileName += IsOneOrTwoGauss[ExtrUseTwoGauss];
    CommonFileName += SIsBkgParab[ExtrBkgType];
    CommonFileName += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    CommonFileName += "_Weighted_BDTCentDep";
    fileName[0] = "_EffW";
    fileName[1] = "";
    namehisto[0] = "histoV2Mixed";
    namehisto[1] = "histoV2MixedCorr";
    hTitleY = "v_{2}";
    hTitleX = TitleXPt;
    YLow = -0.2;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "Weighted";
    sleg[1] = "Corrected";
  }
  else if (TypeComp == 16)
  {
    // TypeComp = 16 --> Pz from Xi - offline vs online acceptance
    numOptions = 2;
    CommonFileName = "Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass4_";
    fileName[0] = "Train365784";
    fileName[1] = "Train361757";
    fileName[0] += "_Xi_BkgParab_Pzs2_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    fileName[1] += "_Xi_BkgParab_Pzs2_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleY = "P_{z,s2}";
    hTitleX = "FT0C Centrality";
    YLow = -0.02;
    YUp = 0.02;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "Offline";
    sleg[1] = "On the fly";
  }

  else if (TypeComp == 17)
  {
    // TypeComp = 17 --> Pz from Lambda - offline vs online acceptance
    numOptions = 2;
    CommonFileName = "Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass4_";
    fileName[0] = "Train354079";
    fileName[1] = "Train370610_ProtonAcc";
    fileName[0] += "_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    fileName[1] += "_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleY = "P_{z,s2}";
    hTitleX = "FT0C Centrality";
    YLow = -0.004;
    YUp = 0.012;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "Offline";
    sleg[1] = "On the fly";
  }
  else if (TypeComp == 18)
  {
    // TypeComp = 18 --> Pz from Lambda - online vs pt vs online vs eta and pt acceptance
    numOptions = 2;
    CommonFileName = "Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass4_";
    fileName[0] = "Train369742";
    fileName[1] = "Train370610_ProtonAcc";
    fileName[0] += "_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    fileName[1] += "_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleY = "P_{z,s2}";
    hTitleX = "FT0C Centrality";
    YLow = -0.004;
    YUp = 0.012;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "On the fly (vs pt)";
    sleg[1] = "On the fly (vs pt and eta)";
  }
  else if (TypeComp == 19)
  {
    // TypeComp = 19 --> Pz from Lambda vs from Xi
    numOptions = 2;
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass4_";
    fileName[0] = "Train370610_ProtonAcc";
    fileName[1] = "Train361757";
    fileName[0] += "_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    fileName[1] += "_Xi_BkgParab_Pzs2_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleY = "P_{z,s2}";
    hTitleX = "FT0C Centrality";
    YLow = -0.004;
    YUp = 0.02;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "From Lambda";
    sleg[1] = "From Xi";
  }
  else if (TypeComp == 20)
  {
    // TypeComp = 20 --> Pz from Lambda fit vs no fit
    numOptions = 2;
    CommonFileName = "Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass4_";
    fileName[0] = "Train370610_ProtonAcc";
    fileName[1] = "Train370610_ProtonAcc";
    fileName[0] += "_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1";
    fileName[1] += "_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleY = "P_{z,s2}";
    hTitleX = "FT0C Centrality";
    YLow = -0.004;
    YUp = 0.02;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    sleg[0] = "From fit";
    sleg[1] = "From Not from fit";
  }
  else if (TypeComp == 21)
  {
    // TypeComp = 21 --> Proton acceptance 2023 pass5 vs 2023 pass4 vs eta
    numOptions = 2;
    CommonFileName = "AcceptancePlots/Acceptance_LHC23_PbPb_pass";
    fileName[0] = "4_Train368064_ProtonAcc_Xi_WithAlpha_Eta08_FromTHN";
    fileName[1] = "5_Train456578_ProtonAcc_Xi_WithAlpha_Eta08_FromTHN";
    namehisto[0] = "Cos2ThetaLambdaFromCVsEta_cent60-70";
    namehisto[1] = "Cos2ThetaLambdaFromCVsEta_cent60-70";
    hTitleY = "Cos^{2}(#theta_{p})";
    hTitleX = "#eta";
    YLow = -0.004;
    YUp = 0.02;
    YLowRatio = 0.95;
    YUpRatio = 1.05;
    sleg[0] = "pass4";
    sleg[1] = "pass5";
    MinHistoX = -0.8;
    MaxHistoX = 0.8;
    yOffset = 6;
  }
  else if (TypeComp == 22)
  {
    // TypeComp = 22 --> Proton acceptance 2023 pass5 vs 2023 pass4 vs pt
    numOptions = 2;
    CommonFileName = "AcceptancePlots/Acceptance_LHC23_PbPb_pass";
    fileName[0] = "4_Train368064_ProtonAcc_Xi_WithAlpha_Eta08_FromTHN";
    fileName[1] = "5_Train456578_ProtonAcc_Xi_WithAlpha_Eta08_FromTHN";
    namehisto[0] = "Cos2ThetaLambdaFromCVsPt_cent60-70";
    namehisto[1] = "Cos2ThetaLambdaFromCVsPt_cent60-70";
    hTitleY = "Cos^{2}(#theta_{p})";
    hTitleX = "p_T (GeV/c)";
    YLow = -0.004;
    YUp = 0.02;
    YLowRatio = 0.95;
    YUpRatio = 1.05;
    sleg[0] = "pass4";
    sleg[1] = "pass5";
    yOffset = 6;
  }
  else if (TypeComp == 23)
  {
    // TypeComp = 23 --> 2023 pass4 vs 2023 pass5 (mean, sigma, purity, yield vs pT)
    numOptions = 2;
    CommonFileName = "OutputAnalysis/FitPzs2_LHC23_PbPb_pass";
    fileName[0] = "4_Train370610_ProtonAcc_Xi_BkgParab";
    fileName[1] = "5_Train456579_ProtAccFromPass4_Xi_BkgParab";
    fileName[0] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[1] += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] += "_PolFromLambda_Eta08_FromTHN_MixedBDT_TightMassCut2.1";
    fileName[1] += "_PolFromLambda_Eta08_FromTHN_MixedBDT_TightMassCut2.1";
    namehisto[0] = "histoPurity";
    namehisto[1] = "histoPurity";
    hTitleY = "Cos^{2}(#theta_{p})";
    hTitleX = "p_T (GeV/c)";
    YLow = -0.004;
    YUp = 0.02;
    YLowRatio = 0.95;
    YUpRatio = 1.05;
    sleg[0] = "pass4";
    sleg[1] = "pass5";
    yOffset = 6;
  }
  else if (TypeComp == 24 || TypeComp == 25 || TypeComp == 26 || TypeComp == 27 || TypeComp == 28)
  {
    // TypeComp = 24, 25, 26, 27, 28 --> 2023 pass4 vs 2023 pass5 (purity, sigma, yield, Pzs2, Pzs2Error  integrated in pT vs centrality)
    numOptions = 2;
    CommonFileName = "Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass";
    fileName[0] = "4_Train370610_ProtonAcc_Xi_BkgParab_Pzs2";
    fileName[1] = "5_Train456579_ProtAccFromPass4_Xi_BkgParab_Pzs2";
    fileName[0] += "_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    fileName[1] += "_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    hTitleY = "";
    hTitleX = "FT0C centrality (%)";
    // purity
    if (TypeComp == 24)
    {
      YLow = 0.8;
      YUp = 1.2;
      YLowRatio = 0.95;
      YUpRatio = 1.05;
      hTitleY = "Purity";
      namehisto[0] = "fHistPuritySummary";
      namehisto[1] = "fHistPuritySummary";
    }
    else if (TypeComp == 25)
    {
      // sigma
      YLow = 0.001;
      YUp = 0.004;
      YLowRatio = 0.8;
      YUpRatio = 1.2;
      hTitleY = "Sigma";
      namehisto[0] = "fHistSigmaSummary";
      namehisto[1] = "fHistSigmaSummary";
    }
    else if (TypeComp == 26)
    {
      // yield
      YLow = 0;
      YUp = 0.1;
      YLowRatio = 0.8;
      YUpRatio = 10;
      hTitleY = "Yield";
      namehisto[0] = "fHistYieldSummary";
      namehisto[1] = "fHistYieldSummary";
    }
    else if (TypeComp == 27)
    {
      // Pzs2
      YLow = -0.002;
      YUp = 0.015;
      YLowRatio = 0.8;
      YUpRatio = 1.2;
      hTitleY = "P_{z, s2}";
      namehisto[0] = "fHistPzs";
      namehisto[1] = "fHistPzs";
    }
    else if (TypeComp == 28)
    {
      // Pzs2 stat error
      YLow = 0.;
      YUp = 0.004;
      YLowRatio = 0.8;
      YUpRatio = 1.2;
      hTitleY = "Stat. error P_{z, s2}";
      namehisto[0] = "fHistPzsError";
      namehisto[1] = "fHistPzsError";
    }
    sleg[0] = "pass4";
    sleg[1] = "pass5";
    yOffset = 6;
    MinHistoX = 0;
    MaxHistoX = 80;
  }
  else if (TypeComp == 29 || TypeComp == 30)
  {
    // TypeComp = 29 --> Lambda Pzs2 from my code vs Preliminary one from Junlee
    numOptions = 2;
    CommonFileName = "";
    fileName[0] = "LambdaJunlee/fout_psi2_mult_WHisto";
    fileName[1] = "Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass5_Train463979_ProtAcceptanceFromSecondayLambdas_Lambda_BkgParab_Pzs2_PtInt_FromTHN_MixedBDT_TightMassCut2.1_ReducedPtBins";
    namehisto[0] = "fHistPzsLambdaJunlee";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "Preliminary";
    sleg[1] = "My result";
    YLow = -0.001;
    YUp = 0.009;
    YLowRatio = 0.5;
    YUpRatio = 2;
    // YLowRatio = -3;
    // YUpRatio = 3;
    if (TypeComp == 30)
    {
      namehisto[0] = "fHistPzsLambdaJunleeStatError";
      namehisto[1] = "fHistPzsError";
      YLow = -0.001;
      YUp = 0.004;
      YLowRatio = 0.8;
      YUpRatio = 1.5;
    }
  }
  else if (TypeComp == 31)
  {
    // TypeComp = 31 --> Lambda Pzs2 from Tree vs THN
    numOptions = 2;
    CommonFileName = "";
    fileName[0] = "Pzs2VsCentrality/Pzs2_LHC25_OO_LambdaPol_Train491711_Lambda_BkgParab_Pzs2_PtInt_TightMassCut2.1NoFit";
    fileName[1] = "Pzs2VsCentrality/Pzs2_LHC25_OO_pass2_Train500672_Lambda_BkgParab_Pzs2_PtInt_FromTHN_TightMassCut2.1NoFit";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "Tree";
    sleg[1] = "THN";
    YLow = -0.001;
    YUp = 0.009;
    YLowRatio = 0.5;
    YUpRatio = 2;
  }
  else if (TypeComp == 32)
  {
    numOptions = 3;
    /*
    // CommonFileName = "OutputAnalysis/FitPzs2_LHC25_OO_pass2_Train503805_Lambda_BkgParab";
    CommonFileName = "../OutputAnalysis/FitPzs2_LHC25_OO_pass2_Train589711_Lambda_BkgParab";
    if (mult == 8) // 70-80%
      CommonFileName += "_Cent0-100";
    else
      CommonFileName += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] = "_CentWeighted_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[1] = "_CentWeighted_Eta08_TightMassCut2.1_ReducedPtBins_isTightest_isSysLambdaMultTrial_ResoOnTheFly";
    fileName[2] = "_CentWeighted_Eta08_TightMassCut2.1_ReducedPtBins_isTightest_isSysLambdaMultTrial_ResoOnTheFly";
    // namehisto[0] = "histoPurityPtInt";
    // namehisto[1] = "histoPurityPtInt";
    // namehisto[2] = "histoPurityPtInt";
    namehisto[0] = "histoPurityPtInt";
    namehisto[1] = "histoPurityPtInt";
    namehisto[2] = "histoPurityPtInt";
    */
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2_Train589711_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[0] = "";
    fileName[1] = "_isTightest_isSysLambdaMultTrial";
    fileName[2] = "_isLoosest_isSysLambdaMultTrial";
    namehisto[0] = "fHistPuritySummary";
    namehisto[1] = "fHistPuritySummary";
    namehisto[2] = "fHistPuritySummary";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "Default";
    sleg[1] = "Tightest";
    sleg[2] = "Loosest";
    YLow = 0.8;
    YUp = 1.1;
    YLowRatio = 0.98;
    YUpRatio = 1.02;
    // namehisto[0] = "fHistYieldSummary";
    // namehisto[1] = "fHistYieldSummary";
    // namehisto[2] = "fHistYieldSummary";
    // YLow = 0;
    // YUp = 1.5;
    // YLowRatio = 0.8;
    // YUpRatio = 1.2;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 33)
  {
    numOptions = 3;
    CommonFileName = "OutputAnalysis/FitPzs2_LHC25_OO_pass2_Train503805_Lambda_BkgParab";
    if (mult == 8) // 70-80%
      CommonFileName += "_Cent0-80";
    else
      CommonFileName += Form("_Cent%i-%i", CentFT0C[mult], CentFT0C[mult + 1]);
    fileName[0] = "_TightMassCut2.1";
    fileName[1] = "_TightMassCut2.1_isTightest_SysMultTrial_0";
    fileName[2] = "_TightMassCut2.1_isLoosest_SysMultTrial_0";
    // namehisto[0] = "histoPurityPtInt";
    // namehisto[1] = "histoPurityPtInt";
    // namehisto[2] = "histoPurityPtInt";
    namehisto[0] = "histoPzs2PtInt";
    namehisto[1] = "histoPzs2PtInt";
    namehisto[2] = "histoPzs2PtInt";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "Default";
    sleg[1] = "Tightest";
    sleg[2] = "Loosest";
    // YLow = 0.;
    // YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.001;
    YUp = 0.005;
  }
  else if (TypeComp == 34)
  {
    numOptions = 2;
    CommonFileName = "Pzs2VsCentrality/Pzs2_LHC25_OO_pass2_Train510678_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins";
    fileName[0] = "_ResoOnTheFly";
    fileName[1] = "";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "Reso on the fly";
    sleg[1] = "Default";
    // YLow = 0.;
    // YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 35)
  {
    numOptions = 4;
    // fileName[0] = "Resolution/Resolution_SP_CFW_LHC25_OO_pass2_Train510916";
    // fileName[1] = "Resolution/Resolution_SP_CFW_LHC25_OO_pass2_Train515731_T0MResolution";
    // fileName[2] = "Resolution/Resolution_SP_CFW_LHC25_OO_pass2_Train515730_V0AResolution";
    fileName[0] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train549964_T0CShiftCorr";
    fileName[1] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train549965_T0MShiftCorr";
    fileName[2] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train547804_V0AShiftCorr";
    fileName[3] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train548457_T0AShiftCorr";
    namehisto[0] = "hReso";
    namehisto[1] = "hReso";
    namehisto[2] = "hReso";
    namehisto[3] = "hReso";
    hTitleY = "Resolution";
    hTitleX = "Centrality (%)";
    YLow = 0;
    YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.2;
    sleg[0] = "T0C";
    sleg[1] = "T0M";
    sleg[2] = "V0A";
    sleg[3] = "T0A";
    yOffset = 2;
  }
  else if (TypeComp == 36)
  {
    // TypeComp = 36 --> 2023 pass5 with pass5 acceptance vs w/ pass4 acceptance (Pzs2, Pzs2Error  integrated in pT vs centrality)
    numOptions = 2;
    CommonFileName = "Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass";
    fileName[0] = "5_Train456579_ProtAccFromPass4_Xi_BkgParab_Pzs2";
    fileName[1] = "5_Train534683_Xi_BkgParab_Pzs2";
    fileName[0] += "_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    fileName[1] += "_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    hTitleY = "";
    hTitleX = "FT0C centrality (%)";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
  }
  else if (TypeComp == 37)
  {
    numOptions = 2;
    CommonFileName = "Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass5";
    fileName[0] = "_Train540301_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    fileName[1] = "_Train534683_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "flat EP";
    sleg[1] = "non-flat EP";
    // YLow = 0.;
    // YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 38)
  {
    numOptions = 2;
    CommonFileName = "Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass5";
    fileName[0] = "_Train541065_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    fileName[1] = "_Train540301_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "|z| < 8 cm";
    sleg[1] = "|z| < 10 cm";
    // YLow = 0.;
    // YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 39)
  {
    numOptions = 3;
    fileName[0] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train549964_T0CShiftCorr";
    fileName[1] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train549964_T0CShiftCorr";
    fileName[2] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train549964_T0CShiftCorr";
    namehisto[0] = "hReso";
    namehisto[1] = "hResoV0ATPCA";
    namehisto[2] = "hResoV0ATPCC";
    hTitleY = "Resolution";
    hTitleX = "Centrality (%)";
    YLow = 0;
    YUp = 0.5;
    YLowRatio = 0.9;
    YUpRatio = 1.6;
    sleg[0] = "T0C with TPCA & TPCC";
    sleg[1] = "T0C with V0A & TPCA";
    sleg[2] = "T0C with V0A & TPCC";
    yOffset = 2;
  }
  else if (TypeComp == 40)
  {
    numOptions = 4;
    // fileName[0] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train556005_T0CShiftCorr_TPCCorr";
    // fileName[1] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train556005_T0CShiftCorr_TPCCorr";
    // fileName[0] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train557787_T0CShiftCorr_TPCCorr_WithT0A";
    // fileName[1] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train557787_T0CShiftCorr_TPCCorr_WithT0A";
    // fileName[2] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train557787_T0CShiftCorr_TPCCorr_WithT0A";
    // fileName[3] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train557787_T0CShiftCorr_TPCCorr_WithT0A";
    // fileName[4] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train557787_T0CShiftCorr_TPCCorr_WithT0A";
    fileName[0] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train567017";
    fileName[1] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train567017";
    fileName[2] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train567017";
    fileName[3] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train567017";
    fileName[4] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train567017";
    namehisto[0] = "hResoV0ATPCA";
    namehisto[1] = "hResoV0ATPCC";
    namehisto[2] = "hResoT0ATPCA";
    namehisto[3] = "hResoT0ATPCC";
    namehisto[4] = "hReso";
    hTitleY = "Resolution";
    hTitleX = "Centrality (%)";
    YLow = 0;
    YUp = 0.5;
    YLowRatio = 0.7;
    YUpRatio = 1.6;
    sleg[0] = "T0C with V0A & TPCA";
    sleg[1] = "T0C with V0A & TPCC";
    sleg[2] = "T0C with T0A & TPCA";
    sleg[3] = "T0C with T0A & TPCC";
    sleg[4] = "T0C with TPCA & TPCC";
    yOffset = 2;
    MinHistoX = 0;
    MaxHistoX = 100;
  }
  else if (TypeComp == 41)
  {
    numOptions = 2;
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2";
    fileName[0] = "_Train562132_wTHN_Lambda_BkgParab_Pzs2_PtInt_Eta08_FromTHN_TightMassCut2.1_ReducedPtBins";
    fileName[1] = "_Train510678_CorrectReso_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "THN";
    sleg[1] = "Tree";
    // YLow = 0.;
    // YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 42)
  {
    numOptions = 2;
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2";
    fileName[0] = "_Train510678_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly_CorrectReso_TestLeassPtBins";
    fileName[1] = "_Train510678_CorrectReso_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "less bins";
    sleg[1] = "more bins";
    // YLow = 0.;
    // YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 43)
  {
    numOptions = 2;
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2";
    fileName[0] = "_Train510678_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly_CorrectReso_TestLeassPtBins";
    fileName[1] = "_Train562850_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "Train 510678";
    sleg[1] = "Train 562850";
    // YLow = 0.;
    // YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 44)
  {
    numOptions = 2;
    fileName[0] = "../Resolution/Resolution_EP_CFW_LHC23_PbPb_pass5_Train563856";
    fileName[1] = "../Resolution/Resolution_SP_CFW_LHC23_PbPb_pass5_Train540301";
    namehisto[0] = "hReso";
    namehisto[1] = "hReso";
    hTitleY = "Resolution";
    hTitleX = "Centrality (%)";
    YLow = 0;
    YUp = 1.5;
    YLowRatio = 0.9;
    YUpRatio = 1.6;
    sleg[0] = "EP";
    sleg[1] = "Normalised SP";
    yOffset = 2;
    MinHistoX = 0;
    MaxHistoX = 80;
  }
  else if (TypeComp == 45)
  {
    numOptions = 2;
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass5_Train540301_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    fileName[0] = "_EPReso";
    fileName[1] = "";
    namehisto[0] = "fHistPzsError";
    namehisto[1] = "fHistPzsError";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "EP reso";
    sleg[1] = "SP reso";
    // YLowRatio = 0.5;
    // YUpRatio = 1.5;
    // YLow = -0.002;
    // YUp = 0.02;
    YLowRatio = 0.7;
    YUpRatio = 1.3;
    YLow = 0.;
    YUp = 0.004;
    MinHistoX = 0;
    MaxHistoX = 80;
  }
  else if (TypeComp == 46)
  {
    // TypeComp = 46 --> Xi Pzs2: preliminaries vs paper proposal results
    numOptions = 2;
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass";
    fileName[0] = "4_Train370610_ProtonAcc_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit";
    fileName[1] = "5_Train540301_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit_EPReso";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleY = "";
    hTitleX = "";
    YLow = -0.004;
    YUp = 0.02;
    YLowRatio = -6;
    YUpRatio = 6;
    sleg[0] = "pass4";
    sleg[1] = "pass5";
    yOffset = 6;
    MinHistoX = 0;
    MaxHistoX = 80;
  }
  else if (TypeComp == 47)
  {
    // TypeComp = 47 --> Xi Pzs2: from fit vs assuming zero background polarization
    numOptions = 2;
    isRatio = 0;
    isFullCorr = 0;
    isFitRatio = 1;
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass5_Train540301_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1";
    // fileName[0] = "NoFit_EPReso";
    fileName[0] = "NoFit_EPReso_NoPurityDivision";
    fileName[1] = "_EPReso";
    // fileName[2] = "_EPReso_isBkgPol0";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    namehisto[2] = "fHistPzs";
    namehisto[3] = "fHistPzs";
    hTitleY = "";
    hTitleX = "";
    YLow = -0.004;
    YUp = 0.02;
    YLowRatio = -0.001;
    YUpRatio = 0.001;
    if (isRatio == 1)
    {
      YLowRatio = 0;
      YUpRatio = 2;
    }
    // sleg[0] = "No fit";
    sleg[0] = "No fit";
    sleg[1] = "Fit";
    // sleg[3] = "Fit, Pz, bkg = 0";
    yOffset = 6;
    MinHistoX = 0;
    MaxHistoX = 80;
  }
  else if (TypeComp == 48)
  {
    // TypeComp = 48 --> Compare acceptance in different runs
    numOptions = 92;
    // Int_t nrun[numOptionsConst] = {, 545312, , 545296, , 545295, , 545294, , 545291, , 545249, 545223, 545210, 545185, 545066};
    // Int_t nrun[numOptionsConst] = {545064, 545063, 545047, 545009, 544992, 544968, 544964, 544963, 544961, 544947};
    // Int_t nrun[numOptionsConst] = {544931, 544917, 544914, 544913, 544896, 544887, 544886, 544868, 544813, 544797};
    // Int_t nrun[numOptionsConst] = {544795, 544794, 544767, 544754, 544742, 544739, 544696, 544694, 544693, 544692};
    // Int_t nrun[numOptionsConst] = {544674, 544672, 544653, 544652, 544640, 544614, 544585, 544583, 544582, 544580};
    // Int_t nrun[numOptionsConst] = {544568, 544567, 544565, 544564, 544551, 544550, 544549, 544548, 544518, 544515};
    // Int_t nrun[numOptionsConst] = {544514, 544512, 544511, 544510, 544508, 544492, 544491, 544490, 544477, 544476};
    // Int_t nrun[numOptionsConst] = {544475, 544474, 544454, 544451, 544392, 544391, 544390, 544389, 544185, 544184};
    // Int_t nrun[numOptionsConst] = {544124, 544123, 544122, 544121, 544116, 544098, 544095, 544091, 544032, 544028, 544013};

    CommonFileName = "../AcceptancePlots/Acceptance_LHC23_PbPb_pass5_";
    fileName[0] = "Train566502_Xi_WithAlpha_Eta08_FromTHN";
    for (Int_t i = 1; i < numOptions; i++)
    {
      fileName[i] = Form("Train566502_%i_Xi_WithAlpha_Eta08_FromTHN", nrun[i - 1]);
    }
    for (Int_t i = 0; i < numOptions; i++)
    {
      namehisto[i] = "Cos2ThetaLambdaFromCVsPt_cent60-70";
    }
    hTitleY = "";
    hTitleX = "";
    YLow = 0;
    YUp = 0.8;
    // YLowRatio = 0.99;
    // YUpRatio = 1.01;
    YLowRatio = 0.7;
    YUpRatio = 1.3;
    // YLowRatio = 0.;
    // YUpRatio = 5;
    sleg[0] = "Train566502";
    for (Int_t i = 1; i < numOptions; i++)
    {
      sleg[i] = Form("Run %i", nrun[i - 1]);
    }
    yOffset = 6;
    // MinHistoX = -0.795;
    // MaxHistoX = 0.795;
    MinHistoX = 0.4;
    MaxHistoX = 10;
  }
  else if (TypeComp == 49)
  {
    numOptions = 2;
    fileName[0] = "../Resolution/Resolution_EP_CFW_LHC23_PbPb_pass5_Train563856";
    fileName[1] = "../Resolution/Resolution_EP_CFW_LHC23_PbPb_pass5_Train567157_OccupancyCut";
    namehisto[0] = "hReso";
    namehisto[1] = "hReso";
    hTitleY = "Resolution";
    hTitleX = "Centrality (%)";
    YLow = 0;
    YUp = 1;
    YLowRatio = 0.9;
    YUpRatio = 1.1;
    sleg[0] = "";
    sleg[1] = "Occupancy < 3000";
    yOffset = 2;
    MinHistoX = 0;
    MaxHistoX = 80;
  }
  else if (TypeComp == 50)
  {
    // TypeComp = 50 --> Xi Pzs2: occupancy dependence
    numOptions = 2;
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC23_PbPb_pass5_Train";
    fileName[0] = "540301_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit_EPReso";
    fileName[1] = "568467_OccupancySel30000_Xi_BkgParab_Pzs2_PolFromLambda_PtInt_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit_EPReso";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleY = "";
    hTitleX = "";
    YLow = -0.004;
    YUp = 0.02;
    YLowRatio = -2;
    YUpRatio = 5;
    sleg[0] = "Default";
    sleg[1] = "Occupancy < 30000";
    yOffset = 6;
    MinHistoX = 0;
    MaxHistoX = 80;
  }
  else if (TypeComp == 51)
  {
    numOptions = 2;
    fileName[0] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train567017";                              // full OO sample
    fileName[1] = "../Resolution/Resolution_EP_CFW_LHC25_OO_pass2_Train557787_T0CShiftCorr_TPCCorr_WithT0A"; // samll OO sample
    namehisto[0] = "hResoV0ATPCA";
    namehisto[1] = "hResoV0ATPCA";
    hTitleY = "Resolution";
    hTitleX = "Centrality (%)";
    YLow = 0;
    YUp = 0.5;
    YLowRatio = 0.7;
    YUpRatio = 1.6;
    sleg[0] = "Full sample";
    sleg[1] = "Partial sample";
    yOffset = 2;
    MinHistoX = 0;
    MaxHistoX = 100;
  }
  else if (TypeComp == 52)
  {
    numOptions = 2;
    isRatio = 0;
    isFullCorr = 1;
    isStoreSyst = 1;
    TypeSyst = "Reso";
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2_Train562850_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[0] = "";
    fileName[1] = "_SystReso";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "Default reso";
    sleg[1] = "T0A + TPCc as ref. detectors";
    // YLow = 0.;
    // YUp = 0.5;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 53)
  {
    numOptions = 2;
    isRatio = 0;
    isFullCorr = 0;
    isStoreSyst = 0;
    isFitRatio = 1;
    TypeSyst = "BkgPol0";
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2_Train562850_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[0] = "";
    fileName[1] = "_isBkgPol0";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "Default";
    sleg[1] = "P_{z, s2, bkg} = 0";
    // YLow = 0.;
    // YUp = 0.5;
    YLowRatio = -0.003;
    YUpRatio = 0.003;
    // YLowRatio = 0.5;
    // YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 54)
  {
    numOptions = 2;
    isRatio = 0;
    isFullCorr = 0;
    isStoreSyst = 0;
    isFitRatio = 1;
    TypeSyst = "PzFitRange";
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2_Train562850_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[0] = "";
    fileName[1] = "_TighterPzFitRange";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "Default";
    sleg[1] = "Tighter P_{z, s2} fit range";
    // YLow = 0.;
    // YUp = 0.5;
    YLowRatio = -0.004;
    YUpRatio = 0.004;
    // YLowRatio = 0.5;
    // YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 55)
  {
    numOptions = 2;
    isRatio = 0;
    isFullCorr = 0;
    isStoreSyst = 0;
    TypeSyst = "BkgFit";
    isFitRatio = 1;
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2_Train562850_Lambda";
    fileName[0] = "_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[1] = "_BkgExpo_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "Default";
    sleg[1] = "Expo background fit";
    // YLow = 0.;
    // YUp = 0.5;
    YLowRatio = -0.004;
    YUpRatio = 0.004;
    // YLowRatio = 0.5;
    // YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 56)
  {
    // TypeComp == 56 --> Compare different acceptance selections
    numOptions = 2;
    isRatio = 1;
    isFullCorr = 0;
    isStoreSyst = 0;
    TypeSyst = "Acceptance";
    CommonFileName = "../AcceptancePlots/Acceptance_LHC25_OO_pass2_SecondaryProtonAcc_Train508938_Lambda_WithAlpha_Eta08_FromTHN_isOOCentrality";
    // fileName[0] = "";
    // fileName[1] = "_TightAcceptance";
    fileName[0] = "_TightAcceptance";
    fileName[1] = "_TightAcceptance2";
    namehisto[0] = "Cos2ThetaLambdaFromCVsPt_cent60-70";
    namehisto[1] = "Cos2ThetaLambdaFromCVsPt_cent60-70";
    // namehisto[0] = "Cos2ThetaLambdaFromCVsEta_cent60-70";
    // namehisto[1] = "Cos2ThetaLambdaFromCVsEta_cent60-70";
    hTitleY = "Cos^{2}(#theta_{p})";
    // hTitleX = "#eta";
    hTitleX = "p_{T} (GeV/c)";
    YLow = 0;
    YUp = 0.5;
    YLowRatio = 0.9;
    YUpRatio = 1.1;
    sleg[0] = "Default";
    sleg[1] = "tighter purity";
    // MinHistoX = -0.8;
    // MaxHistoX = 0.8;
    MinHistoX = 0.5;
    MaxHistoX = 10;
    yOffset = 6;
  }
  else if (TypeComp == 57)
  {
    // TypeComp == 57 --> Compare purity with different selections for acceptance
    numOptions = 2;
    isRatio = 1;
    isFullCorr = 0;
    isStoreSyst = 0;
    CommonFileName = "../OutputAnalysis/FitPzs2_LHC25_OO_pass2_Train562850_Lambda_BkgParab_Cent60-70_CentWeighted_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[0] = "";
    fileName[1] = "_isTightMassForAcceptancePurity";
    namehisto[0] = "histoPurity";
    namehisto[1] = "histoPurity";
    hTitleY = "Purity";
    hTitleX = "p_{T}";
    YLow = 0.85;
    YUp = 1.05;
    YLowRatio = 0.9;
    YUpRatio = 1.1;
    sleg[0] = "Default";
    sleg[1] = "tighter (1.114 < m < 1.117)";
    // MinHistoX = -0.8;
    // MaxHistoX = 0.8;
    MinHistoX = 0.5;
    MaxHistoX = 8;
    yOffset = 6;
  }
  else if (TypeComp == 58)
  {
    numOptions = 2;
    isRatio = 1;
    isFullCorr = 0;
    isStoreSyst = 0;
    TypeSyst = "BkgFit";
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2";
    fileName[0] = "_Train576495_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[1] = "_Train562850_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "|#eta| < 0.8";
    sleg[1] = "|#eta| < 100";
    // YLow = 0.;
    // YUp = 0.5;
    // YLowRatio = -0.0003;
    // YUpRatio = 0.0003;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 59)
  {
    // TypeComp == 59 --> Compare Pzs2 of Lambda from THN and from tree (reso and cent weights applied)
    numOptions = 2;
    isRatio = 1;
    isFullCorr = 1;
    isStoreSyst = 0;
    TypeSyst = "";
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2";
    // fileName[0] = "_Train575744_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_FromTHN_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[0] = "_Train589559_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_FromTHN_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[1] = "_Train576495_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "from THN";
    sleg[1] = "from tree";
    // YLow = 0.;
    // YUp = 0.5;
    // YLowRatio = -0.0003;
    // YUpRatio = 0.0003;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 60)
  {
    // TypeComp == 60 --> Compare Pzs2 of Lambda with |z| < 10 cm and |z| < 8 cm
    numOptions = 2;
    isRatio = 0;
    isFullCorr = 0;
    isStoreSyst = 1;
    isFitRatio = 0;
    TypeSyst = "ZVertex";
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2";
    fileName[0] = "_Train576495_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[1] = "_Train576496_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "|z| < 10 cm";
    sleg[1] = "|z| < 8 cm";
    // YLow = 0.;
    // YUp = 0.5;
    // YLowRatio = -0.0003;
    // YUpRatio = 0.0003;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 61)
  {
    // TypeComp == 61 --> Compare Pzs2 of Lambda with new and old acceptance (new = higher purity)
    numOptions = 2;
    isRatio = 0;
    isFullCorr = 0;
    isStoreSyst = 1;
    isFitRatio = 1;
    TypeSyst = "AcceptancePurity";
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2";
    fileName[0] = "_Train576495_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[1] = "_Train589711_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "Def acceptance";
    sleg[1] = "New acceptance";
    // YLow = 0.;
    // YUp = 0.5;
    // YLowRatio = -0.0003;
    // YUpRatio = 0.0003;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 62)
  {
    // TypeComp == 62 --> Compare Pzs2 of Lambda with FD correction vs without FD correction
    numOptions = 2;
    isRatio = 1;
    isFullCorr = 0;
    isStoreSyst = 0;
    isFitRatio = 0;
    TypeSyst = "AcceptancePurity";
    CommonFileName = "../Pzs2VsCentrality/Pzs2_LHC25_OO_pass2";
    fileName[0] = "_Train589711_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly";
    fileName[1] = "_Train589711_Lambda_BkgParab_Pzs2_CentWeighted_PtInt_Eta08_TightMassCut2.1_ReducedPtBins_ResoOnTheFly_FDCorrected";
    namehisto[0] = "fHistPzs";
    namehisto[1] = "fHistPzs";
    hTitleX = "FT0C centrality (%)";
    sleg[0] = "No FD correction";
    sleg[1] = "With FD correction";
    // YLow = 0.;
    // YUp = 0.5;
    // YLowRatio = -0.0003;
    // YUpRatio = 0.0003;
    YLowRatio = 0.5;
    YUpRatio = 1.5;
    YLow = -0.002;
    YUp = 0.02;
    MinHistoX = 0;
    MaxHistoX = 90;
  }
  else if (TypeComp == 63)
  {
    // TypeComp == 63 --> Compare proton acceptance with and without |eta| < 0.8 cut for daughters
    numOptions = 2;
    isRatio = 1;
    isFullCorr = 0;
    isStoreSyst = 0;
    TypeSyst = "Acceptance";
    CommonFileName = "../AcceptancePlots/Acceptance_LHC25_OO_pass2";
    fileName[0] = "_SecondaryProtonAcc_Train508938_Lambda_WithAlpha_Eta08_FromTHN_isOOCentrality";
    fileName[1] = "_Train597528_NewAcc_Lambda_WithAlpha_Eta08_FromTHN_isOOCentrality_TightAcceptance";
    // namehisto[0] = "Cos2ThetaLambdaFromCVsPt_cent60-70";
    // namehisto[1] = "Cos2ThetaLambdaFromCVsPt_cent60-70";
    namehisto[0] = "Cos2ThetaLambdaFromCVsEta_cent60-70";
    namehisto[1] = "Cos2ThetaLambdaFromCVsEta_cent60-70";
    hTitleY = "Cos^{2}(#theta_{p})";
    hTitleX = "#eta";
    // hTitleX = "p_{T} (GeV/c)";
    YLow = 0;
    YUp = 0.5;
    YLowRatio = 0.8;
    YUpRatio = 1.2;
    sleg[0] = "|#eta| < 1.5 for daughters";
    sleg[1] = "|#eta| < 0.8 for daughters";
    MinHistoX = -0.8;
    MaxHistoX = 0.8;
    // MinHistoX = 0.5;
    // MaxHistoX = 10;
    yOffset = 6;
  }
  else if (TypeComp == 64)
  {
    // TypeComp == 64 --> Compare proton acceptance with |eta| < 0.8 and 0 < eta < 0.8 and -0.8 < eta < 0
    numOptions = 3;
    isRatio = 1;
    isFullCorr = 0;
    isStoreSyst = 0;
    TypeSyst = "Acceptance";
    CommonFileName = "../AcceptancePlots/Acceptance_LHC25_OO_pass2";
    fileName[0] = "_Train597528_NewAcc_Lambda_WithAlpha_Eta08_FromTHN_isOOCentrality_TightAcceptance";
    fileName[1] = "_Train597527_NewAcc_EtaPos_Lambda_WithAlpha_Eta08_FromTHN_isOOCentrality_TightAcceptance";
    fileName[2] = "_Train597526_NewAcc_EtaNeg_Lambda_WithAlpha_Eta08_FromTHN_isOOCentrality_TightAcceptance";
    namehisto[0] = "Cos2ThetaLambdaFromCVsEta_cent60-70";
    namehisto[1] = "Cos2ThetaLambdaFromCVsEta_cent60-70";
    namehisto[2] = "Cos2ThetaLambdaFromCVsEta_cent60-70";
    hTitleY = "Cos^{2}(#theta_{p})";
    hTitleX = "#eta";
    //hTitleX = "p_{T} (GeV/c)";
    YLow = 0;
    YUp = 0.5;
    YLowRatio = 0.4;
    YUpRatio = 1.2;
    sleg[0] = "|#eta| < 0.8";
    sleg[1] = "0 < #eta < 0.8";
    sleg[2] = "-0.8 < #eta < 0";
    MinHistoX = -0.8;
    MaxHistoX = 0.8;
    //MinHistoX = 0.5;
    //MaxHistoX = 10;
    yOffset = 6;
  }
  else
  {
    cout << "TypeComp not defined" << endl;
    return;
  }
  cout << "CommonFileName: " << CommonFileName << endl;
  TString hTitleYRatio = "Ratio to " + sleg[0];
  if (!isRatio)
    hTitleYRatio = "Difference to " + sleg[0];
  if (TypeComp == 40)
    hTitleYRatio = "Ratio to default";

  double frac = 0;
  TH1F *hEvents = new TH1F("hEvents", "hEvents", 92, 0.5, 92.5);
  // TH1F *hgauss = new TH1F("hgauss", "hgauss", 100, 0.8, 1.2);
  TH1F *hgauss = new TH1F("hgauss", "hgauss", 100, -5, 5);
  TH1F *hgaussSigmaDelta = new TH1F("hgaussSigmaDelta", "hgaussSigmaDelta", 100, 0, 0.05);
  TH1F *hgaussPt[100];
  TLegend *legGauss = new TLegend(0.6, 0.7, 0.9, 0.9);
  legGauss->SetBorderSize(0);
  legGauss->SetTextSize(0.04);
  legGauss->SetFillStyle(0);
  for (Int_t j = 1; j <= 100; j++)
    hgaussPt[j - 1] = new TH1F(Form("hgaussPt_%d", j - 1), Form("hgaussPt_%d", j - 1), 50, -5, 5);

  for (Int_t i = 0; i < numOptions; i++)
  {
    cout << "CommonFileName: " << CommonFileName << endl;
    Sinputfile = CommonFileName + fileName[i] + ".root";
    cout << "Input file: " << Sinputfile << endl;
    TFile *inputFile = new TFile(Sinputfile);
    if (i == 0)
    {
      hDef = (TH1F *)inputFile->Get(namehisto[i]);
      if (!hDef)
      {
        cout << "Histogram not found" << endl;
        return;
      }
      hDef->SetName("hDefault");
    }
    else
    {
      h[i] = (TH1F *)inputFile->Get(namehisto[i]);
      if (!h[i])
      {
        cout << "Histogram not found" << endl;
        return;
      }
      h[i]->SetName(Form("hVar_%i", i));
      hRatio[i] = (TH1F *)h[i]->Clone(Form("hRatio_%i", i));
      if (isRatio)
      {
        hRatio[i]->Divide(hDef);
        if (isFullCorr)
          ErrRatioCorr(h[i], hDef, hRatio[i], 1);
        else
          ErrRatioCorr(h[i], hDef, hRatio[i], 0);
      }
      else
      {
        hRatio[i]->Add(hDef, -1);
        for (Int_t j = 1; j <= hRatio[i]->GetNbinsX(); j++)
        {
          if (isFullCorr)
            hRatio[i]->SetBinError(j, std::abs(hDef->GetBinError(j) - h[i]->GetBinError(j)));
          else
            hRatio[i]->SetBinError(j, sqrt(std::abs(pow(hDef->GetBinError(j), 2) - pow(h[i]->GetBinError(j), 2))));
        }
      }
      for (Int_t j = 1; j <= hRatio[i]->GetNbinsX(); j++)
      {
        cout << "Bin center: " << hRatio[i]->GetBinCenter(j) << endl;
        if (isRatio)
          cout << "Ratio = " << hRatio[i]->GetBinContent(j) << " +/- " << hRatio[i]->GetBinError(j) << ", |ratio -1|/sigma " << abs(hRatio[i]->GetBinContent(j) - 1) / hRatio[i]->GetBinError(j) << endl;
        else
          cout << "Difference = " << hRatio[i]->GetBinContent(j) << " +/- " << hRatio[i]->GetBinError(j) << ", |diff|/sigma " << abs(hRatio[i]->GetBinContent(j)) / hRatio[i]->GetBinError(j) << endl;
      }

      if (TypeComp == 29)
      {
        for (Int_t j = 1; j <= hRatio[i]->GetNbinsX(); j++)
        {
          hRatio[i]->SetBinContent(j, (hRatio[i]->GetBinContent(j) - hDef->GetBinContent(j)) / hDef->GetBinError(j));
        }
      }
      else if (TypeComp == 48)
      {
        ErrRatioCorr(h[i], hDef, hRatio[i], 0); // numerator is a subsample of denominator
        double weight = 0;
        double n_tot = 0;
        TFile *fileTot = new TFile("../TreeForAnalysis/AnalysisResults_LHC23_PbPb_pass5_Train566502.root");
        if (fileTot->IsOpen())
        {
          TH1F *hEventsTot = (TH1F *)fileTot->Get("lf-cascade-flow/histos/hEventVertexZ");
          if (!hEventsTot)
            return;
          n_tot = hEventsTot->GetEntries();
          fileTot->Close();
          cout << "ntot" << n_tot << endl;
        }
        else
        {
          cout << "Could not open total events file!" << endl;
        }
        double n_run = 0;
        cout << "ntot" << n_tot << endl;
        TFile *fileRun = new TFile(Form("../TreeForAnalysis/AnalysisResults_LHC23_PbPb_pass5_Train566502_%i.root", nrun[i - 1]));
        if (fileRun->IsOpen())
        {
          TH1F *hEventsPerRun = (TH1F *)fileRun->Get("lf-cascade-flow/histos/hEventVertexZ");
          n_run = hEventsPerRun->GetEntries();
          weight = (double)n_run / (double)n_tot;
          fileRun->Close();
        }
        else
        {
          cout << "Could not open run events file!" << endl;
        }
        cout << "Weight is: " << weight << " for run " << nrun[i - 1] << endl;
        for (Int_t j = 1; j <= hRatio[i]->GetNbinsX(); j++)
        {
          // OPTION 1: relative difference divided by uncertainty on the default
          float delta = h[i]->GetBinContent(j) - hDef->GetBinContent(j);
          float sigmaDelta = sqrt((pow(h[i]->GetBinError(j), 2) - pow(hDef->GetBinError(j), 2))); // numerator is a subsample of denominator
          // hRatio[i]->SetBinContent(j, delta / sigmaDelta);
          // hRatio[i]->SetBinError(j, 0);
          // END OPTION 1
          // OPTION 2: show ratio only if significantly different from zero
          if (std::abs(hRatio[i]->GetBinContent(j) - 1) < 3 * hRatio[i]->GetBinError(j))
          {
            hRatio[i]->SetBinContent(j, 1e-2);
            hRatio[i]->SetBinError(j, 0);
          }
          else
          {
            SignRun[i - 1] = 1;
            cout << "i: " << i << " Run " << nrun[i - 1] << ", bin " << j << ": delta = " << delta << ", sigmaDelta = " << sigmaDelta << ", ratio = " << hRatio[i]->GetBinContent(j) << " +/- " << hRatio[i]->GetBinError(j) << endl;
          }
          // END OPTION 2
          // hgauss->Fill(hRatio[i]->GetBinContent(j)); //FILLING WITH RATIO
          hgauss->Fill(delta / sigmaDelta, weight); // FILLING WITH PULL VALUE
          hEvents->SetBinContent(i, weight);
          hgaussSigmaDelta->Fill(sigmaDelta);
          hgaussPt[j - 1]->Fill(delta / sigmaDelta);
        }
        if (SignRun[i - 1] == 1)
        {
          cout << "Adding weights " << endl;
          frac += weight;
          cout << "frac: " << frac << endl;
        }
      }

      if (TypeComp == 30)
      {
        for (Int_t j = 1; j <= hRatio[i]->GetNbinsX(); j++)
        {
          hRatio[i]->SetBinContent(j, h[i]->GetBinContent(j) / hDef->GetBinContent(j));
          hRatio[i]->SetBinError(j, 0);
        }
      }
    }
  }

  TCanvas *cevents = new TCanvas("cevents", "cevents", 700, 500);
  hEvents->SetTitle("");
  hEvents->GetXaxis()->SetTitle("Run index");
  hEvents->GetYaxis()->SetTitle("Fraction of events");
  hEvents->Draw();

  TCanvas *cgaussSigmaDelta = new TCanvas("cgaussSigmaDelta", "cgaussSigmaDelta", 700, 500);
  hgaussSigmaDelta->SetTitle("");
  hgaussSigmaDelta->GetXaxis()->SetTitle("#sigma_{#Delta}");
  hgaussSigmaDelta->GetYaxis()->SetTitle("Entries");
  hgaussSigmaDelta->Draw();

  TCanvas *cgaus = new TCanvas("cgaus", "cgaus", 700, 500);
  legGauss->AddEntry("", Form("#mu = %.3f #pm %.3f", hgauss->GetMean(), hgauss->GetMeanError()), "");
  legGauss->AddEntry("", Form("RMS = %.3f", hgauss->GetRMS()), "");
  // TF1 *gauss = new TF1("gauss", "gaus", 0.8, 1.2);
  TF1 *gauss = new TF1("gauss", "gaus", -3, 3);
  hgauss->SetTitle("");
  // hgauss->GetXaxis()->SetTitle("hVar / hDef");
  hgauss->GetXaxis()->SetTitle("(hVar - hDef) / #sigma");
  hgauss->Fit(gauss, "R+");
  legGauss->AddEntry("", Form("#mu fit = %.3f #pm %.3f", gauss->GetParameter(1), gauss->GetParError(1)), "");
  legGauss->AddEntry("", Form("#sigma fit = %.3f #pm %.3f", gauss->GetParameter(2), gauss->GetParError(2)), "");
  legGauss->Draw();
  for (Int_t j = 1; j <= hDef->GetNbinsX(); j++)
  {
    // hgaussPt[j - 1]->Draw("same");
  }
  if (TypeComp == 48)
    cgaus->SaveAs(Form("../CompareResults/GaussianDistribution_TypeComp%i.png", TypeComp));

  TCanvas *cgausVsPt = new TCanvas("cgausVsPt", "cgausVsPt", 700, 500);
  TH1F *hGaussVsPt = (TH1F *)hDef->Clone("hGaussVsPt");
  for (Int_t j = 1; j <= hDef->GetNbinsX(); j++)
  {
    hgaussPt[j - 1]->Fit(gauss, "R+");
    hGaussVsPt->SetBinContent(j, gauss->GetParameter(1));
    // hGaussVsPt->SetBinError(j, gauss->GetParameter(2));
    hGaussVsPt->SetBinError(j, hgaussPt[j - 1]->GetRMS());
  }
  hGaussVsPt->SetTitle("");
  hGaussVsPt->GetYaxis()->SetTitle("Mean of (hVar - hDef) / #sigma");
  hGaussVsPt->GetXaxis()->SetTitle(hTitleX);
  hGaussVsPt->GetYaxis()->SetRangeUser(-3, 3);
  hGaussVsPt->Draw();
  TLine *lineZero = new TLine(hGaussVsPt->GetXaxis()->GetXmin(), 0, hGaussVsPt->GetXaxis()->GetXmax(), 0);
  lineZero->SetLineColor(kBlue);
  lineZero->SetLineStyle(2);
  lineZero->Draw("same");
  TLine *linePlus = new TLine(hGaussVsPt->GetXaxis()->GetXmin(), 1, hGaussVsPt->GetXaxis()->GetXmax(), 1);
  linePlus->SetLineColor(kGray);
  linePlus->SetLineStyle(2);
  linePlus->Draw("same");
  TLine *lineMinus = new TLine(hGaussVsPt->GetXaxis()->GetXmin(), -1, hGaussVsPt->GetXaxis()->GetXmax(), -1);
  lineMinus->SetLineColor(kGray);
  lineMinus->SetLineStyle(2);
  lineMinus->Draw("same");
  if (TypeComp == 48)
    cgausVsPt->SaveAs(Form("../CompareResults/GaussianDistributionVsEta_TypeComp%i.png", TypeComp));

  TCanvas *canvas = new TCanvas("canvas", "canvas", 700, 900);
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B
  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvas->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow, YUp, 1, 1, hTitleX, hTitleY, "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetYaxis()->SetTitleOffset(yOffset);
  hDummy->GetXaxis()->SetRangeUser(MinHistoX, MaxHistoX);
  if (TypeComp == 5 || TypeComp == 16 || TypeComp == 17 || TypeComp == 18 || TypeComp == 19 || TypeComp == 20 || TypeComp == 29)
    hDummy->GetXaxis()->SetRangeUser(0, 80);
  if (TypeComp == 34)
    hDummy->GetXaxis()->SetRangeUser(0, 90);
  if (TypeComp == 35 || TypeComp == 39 || TypeComp == 40 || TypeComp == 51)
    hDummy->GetXaxis()->SetRangeUser(0, 100);
  if (TypeComp == 43)
    hDummy->GetXaxis()->SetRangeUser(0, 90);
  if (TypeComp == 44)
    hDummy->GetXaxis()->SetRangeUser(0, 80);
  pad1->Draw();
  pad1->cd();
  // hDummy->Draw("same");
  hDummy->DrawClone("same");

  hDef->SetMarkerColor(ColorMult[5]);
  hDef->SetLineColor(ColorMult[5]);
  hDef->SetMarkerStyle(MarkerMult[0]);
  hDef->SetMarkerSize(0.6 * SizeMult[0]);
  hDef->SetTitle("");
  hDef->GetYaxis()->SetRangeUser(YLow, YUp);
  hDef->GetXaxis()->SetRangeUser(MinHistoX, MaxHistoX);
  hDef->Draw("same");
  for (Int_t i = 1; i < numOptions; i++)
  {
    int indexColor = i + 1;
    if (i > 10 && i < 20)
      indexColor = i + 1 - 10;
    else if (i >= 20 && i < 30)
      indexColor = i + 1 - 20;
    else if (i >= 30 && i < 40)
      indexColor = i + 1 - 30;
    else if (i >= 40 && i < 50)
      indexColor = i + 1 - 40;
    else if (i >= 50 && i < 60)
      indexColor = i + 1 - 50;
    else if (i >= 60 && i < 70)
      indexColor = i + 1 - 60;
    else if (i >= 70 && i < 80)
      indexColor = i + 1 - 70;
    else if (i >= 80 && i < 90)
      indexColor = i + 1 - 80;
    else if (i >= 90 && i < 100)
      indexColor = i + 1 - 90;
    h[i]->SetMarkerColor(ColorMult[indexColor]);
    h[i]->SetLineColor(ColorMult[indexColor]);
    h[i]->SetMarkerStyle(MarkerMult[indexColor]);
    h[i]->SetMarkerSize(0.6 * SizeMult[indexColor]);
    h[i]->SetTitle("");
    h[i]->Draw("same");
    if (TypeComp == 40 && i == 3)
    {
      h[i]->SetMarkerColor(kBlue + 1);
      h[i]->SetLineColor(kBlue + 1);
    }
    if (TypeComp == 40 && i == 2)
    {
      h[i]->SetMarkerColor(kCyan + 1);
      h[i]->SetLineColor(kCyan + 1);
    }
  }

  TLegend *leg = new TLegend(0.3, 0.7, 0.9, 0.9);
  if (TypeComp == 23 || TypeComp == 24 || TypeComp == 25)
    leg = new TLegend(0.3, 0.3, 0.9, 0.5);
  else if (TypeComp == 48)
  {
    leg = new TLegend(0.3, 0.5, 0.9, 0.92);
    leg->SetNColumns(2);
  }
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  if (TypeComp == 5)
    leg->AddEntry("", ParticleNameLegend[ChosenPart] + " Pb-Pb 5.36 TeV", "");
  else if (TypeComp == 11)
    leg->AddEntry("", ParticleNameLegend[ChosenPart] + " Pb-Pb 5.36 TeV", "");
  else if (TypeComp == 12 || TypeComp == 48)
    leg->AddEntry("", "Pb-Pb 5.36 TeV", "");
  else if (TypeComp == 19 || TypeComp == 24 || TypeComp == 25 || TypeComp == 26 || TypeComp == 27 || TypeComp == 28 || TypeComp == 44 || TypeComp == 45 || TypeComp == 46 || TypeComp == 47 || TypeComp == 50)
    leg->AddEntry("", "Pb-Pb 5.36 TeV", "");
  else if (TypeComp == 32 || TypeComp == 33 || TypeComp == 35 || TypeComp == 39 || TypeComp == 40 || TypeComp >= 51)
    leg->AddEntry("", "OO 5.36 TeV", "");
  else
    leg->AddEntry("", ParticleNameLegend[ChosenPart] + Form(" Pb-Pb 5.36 TeV, FT0C %i-%i", CentFT0C[mult], CentFT0C[mult + 1]) + "%", "");
  leg->AddEntry(hDef, sleg[0], "lp");
  for (Int_t i = 1; i < numOptions; i++)
  {
    leg->AddEntry(h[i], sleg[i], "lp");
  }
  leg->Draw("same");

  TF1 *FitPol0[100];
  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, -100, 100);
  for (Int_t i = 1; i <= hDummyRatio->GetNbinsX(); i++)
    hDummyRatio->SetBinContent(i, 1e-12);
  SetFont(hDummyRatio);
  StyleHistoYield(hDummyRatio, YLowRatio, YUpRatio, 1, 1, hTitleX, hTitleYRatio, "", 1, 1.15, YoffsetSpectraRatio);
  //  StyleHistoYield(hDummyRatio, 0, 0.4, 1, 1, hTitleX, hTitleYRatio, "", 1, 1.15, YoffsetSpectraRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetTickLength(hDummyRatio, tickX, tickY);
  hDummyRatio->GetXaxis()->SetRangeUser(MinHistoX, MaxHistoX);
  if (TypeComp == 5 || TypeComp == 29 || TypeComp == 30 || TypeComp == 44)
    hDummyRatio->GetXaxis()->SetRangeUser(0, 80);
  if (TypeComp == 34)
    hDummyRatio->GetXaxis()->SetRangeUser(0, 90);
  if (TypeComp == 35 || TypeComp == 39 || TypeComp == 40 || TypeComp == 51)
    hDummyRatio->GetXaxis()->SetRangeUser(0, 100);
  canvas->cd();
  padL1->Draw();
  padL1->cd();
  hDummyRatio->Draw("");
  TLegend *legFit = new TLegend(0.3, 0.7, 0.9, 0.9);
  if (isFitRatio)
  {
    legFit->SetBorderSize(0);
    legFit->SetFillStyle(0);
    legFit->SetTextSize(0.06);
    legFit->Draw("same");
  }
  for (Int_t i = 1; i < numOptions; i++)
  {
    int indexColor = i + 1;
    if (i > 10 && i < 20)
      indexColor = i + 1 - 10;
    else if (i >= 20 && i < 30)
      indexColor = i + 1 - 20;
    else if (i >= 30 && i < 40)
      indexColor = i + 1 - 30;
    else if (i >= 40 && i < 50)
      indexColor = i + 1 - 40;
    else if (i >= 50 && i < 60)
      indexColor = i + 1 - 50;
    else if (i >= 60 && i < 70)
      indexColor = i + 1 - 60;
    else if (i >= 70 && i < 80)
      indexColor = i + 1 - 70;
    else if (i >= 80 && i < 90)
      indexColor = i + 1 - 80;
    else if (i >= 90 && i < 100)
      indexColor = i + 1 - 90;

    FitPol0[i] = new TF1(Form("FitPol0_%d", i), "pol0", MinHistoX, MaxHistoX);
    FitPol0[i]->SetLineColor(ColorMult[indexColor]);
    hRatio[i]->SetMarkerColor(ColorMult[indexColor]);
    hRatio[i]->SetLineColor(ColorMult[indexColor]);
    hRatio[i]->SetMarkerStyle(MarkerMult[indexColor]);
    hRatio[i]->SetMarkerSize(0.6 * SizeMult[indexColor]);
    hRatio[i]->GetYaxis()->SetRangeUser(YLowRatio, YUpRatio);
    hRatio[i]->Draw("same");
    if (TypeComp == 40 && i == 3)
    {
      hRatio[i]->SetMarkerColor(kBlue + 1);
      hRatio[i]->SetLineColor(kBlue + 1);
    }
    if (TypeComp == 40 && i == 2)
    {
      hRatio[i]->SetMarkerColor(kCyan + 1);
      hRatio[i]->SetLineColor(kCyan + 1);
    }
    if (isFitRatio)
    {
      hRatio[i]->Fit(FitPol0[i], "R+");
      legFit->AddEntry(FitPol0[i], Form("p0 = %.5f #pm %.5f", FitPol0[i]->GetParameter(0), FitPol0[i]->GetParError(0)), "l");
    }
  }
  if (isFitRatio)
  {
    legFit->Draw("same");
  }
  TF1 *line = new TF1("line", "1", MinHistoX, MaxHistoX);
  if (TypeComp == 5 || TypeComp == 29 || TypeComp == 30)
    line = new TF1("line", "1", 0, 80);
  else if (TypeComp == 35 || TypeComp == 39 || TypeComp == 40 || TypeComp == 51)
    line = new TF1("line", "1", 0, 100);
  line->SetLineColor(kBlack);
  line->SetLineStyle(9);
  line->Draw("same");
  padL1->Modified();
  padL1->Update();

  canvas->SaveAs(Form("../CompareResults/Canvas_CompareResults%i.png", TypeComp));
  canvas->SaveAs(Form("../CompareResults/Canvas_CompareResults%i.pdf", TypeComp));

  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 900, 700); // without ratio
  StyleCanvas(canvas2, 0.05, 0.15, 0.15, 0.05);
  hDummy->GetXaxis()->SetLabelOffset(0.02);
  hDummy->GetXaxis()->SetTitleSize(30);
  hDummy->GetXaxis()->SetTitleOffset(1.5);
  hDummy->Draw("same");
  hDef->Draw("same");
  for (Int_t i = 1; i < numOptions; i++)
  {
    h[i]->Draw("same");
  }
  leg->Draw("same");
  canvas2->SaveAs(Form("../CompareResults/Canvas_CompareResults%i_NoRatio.png", TypeComp));
  canvas2->SaveAs(Form("../CompareResults/Canvas_CompareResults%i_NoRatio.pdf", TypeComp));

  TFile *outputForSyst;
  TString SoutputForSyst = "../CompareResults/SystUncertainty_" + TypeSyst + ".root";
  if (isStoreSyst && isRatio == 0)
  {
    outputForSyst = new TFile(SoutputForSyst, "RECREATE");
    for (Int_t i = 1; i < numOptions; i++)
    {
      TH1F *hRatioClone = (TH1F *)hRatio[i]->Clone(Form("hRatioClone_%i", i));
      hRatioClone->Scale(1. / TMath::Sqrt(12));
      hRatioClone->Write();
    }
    outputForSyst->Close();
    cout << "\nSyst uncertainty histograms stored in file: " << SoutputForSyst << endl;
  }
}
