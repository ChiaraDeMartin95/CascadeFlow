#include "Riostream.h"
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
#include "TGraphErrors.h"
#include "CommonVarPub.h"
// #include "CommonVarXi.h"
#include "CommonVarLambda.h"
// #include "CommonVarOmega.h"
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

Float_t YLow[numPart] = {-0.001};
Float_t YUp[numPart] = {0.011};

// Int_t ColorJunlee = kAzure - 3;
Int_t ColorJunlee = kGreen + 2;
Int_t ColorOO = kMagenta + 1;

void PzsVsPt(Int_t ChosenPart = ChosenParticle,
             Bool_t isRebin = 1,
             Bool_t isBkgPol = 1,
             Bool_t isPolFromLambda = 0,
             Bool_t isFromFit = 1,
             Bool_t isGlobalFit = 0,
             Bool_t isFitDSCB = ExtrisFitDSCB,
             Bool_t isFDCorrected = 0,
             Bool_t isTighterPzFitRange = 0,
             Bool_t isRapiditySel = ExtrisRapiditySel,
             Int_t BkgType = ExtrBkgType,
             Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  Int_t part = 0;
  if (ChosenPart == 1 || ChosenPart == 4 || ChosenPart == 5)
  {
    part = 1; // Omega
  }

  Float_t UpperRangeParticle = PtBins[numPtBins - 1];
  Float_t LowerRangeParticle = PtBins[0];

  // fileinLambda
  TString PathInLambda = "../Run2Results/HEPData-ins1891389-PzVsPt.root";
  TFile *fileInLambda = TFile::Open(PathInLambda);
  if (!fileInLambda)
  {
    cout << "No file found" << endl;
    return;
  }
  TDirectoryFile *dirLambda = (TDirectoryFile *)fileInLambda->Get("P_z vs. pT (5.02TeV, 30-50%)");
  if (!dirLambda)
  {
    cout << "No directory found" << endl;
    return;
  }
  TH1F *fHistPzsLambda = (TH1F *)dirLambda->Get("Hist1D_y1");
  if (!fHistPzsLambda)
  {
    cout << "No hist found" << endl;
    return;
  }
  TH1F *fHistPzsLambdaSist = (TH1F *)fHistPzsLambda->Clone("fHistPzsLambdaSist");
  TH1F *fHistPzsLambda_StatErr = (TH1F *)dirLambda->Get("Hist1D_y1_e1");
  if (!fHistPzsLambda_StatErr)
  {
    cout << "No hist Stat found" << endl;
    return;
  }
  TH1F *fHistPzsLambda_SystErr = (TH1F *)dirLambda->Get("Hist1D_y1_e2");
  if (!fHistPzsLambda_SystErr)
  {
    cout << "No hist Stat found" << endl;
    return;
  }
  for (Int_t i = 1; i <= fHistPzsLambda->GetNbinsX(); i++)
  {
    fHistPzsLambda->SetBinError(i, fHistPzsLambda_StatErr->GetBinContent(i));
    fHistPzsLambdaSist->SetBinError(i, fHistPzsLambda_SystErr->GetBinContent(i));
    fHistPzsLambda_StatErr->SetBinError(i, 0);
  }
  fHistPzsLambda->SetMarkerStyle(20);
  fHistPzsLambda->SetMarkerSize(1.5);
  fHistPzsLambda->SetMarkerColor(kBlue);
  fHistPzsLambda->SetLineColor(kBlue);
  fHistPzsLambdaSist->SetMarkerColor(kBlue);
  fHistPzsLambdaSist->SetLineColor(kBlue);
  fHistPzsLambda_StatErr->SetMarkerStyle(20);
  fHistPzsLambda_StatErr->SetMarkerSize(1.5);
  fHistPzsLambda_StatErr->SetMarkerColor(kBlue);
  fHistPzsLambda_StatErr->SetLineColor(kBlue);

  TString SfilePbPbJunlee = "../LambdaJunlee/HEPData-ins3169170-v1-Table_4.root";
  TFile *fileJunlee = new TFile(SfilePbPbJunlee);
  TDirectoryFile *dirJunlee = (TDirectoryFile *)fileJunlee->Get("Table 4");
  TH1F *hPzsJunlee = (TH1F *)dirJunlee->Get("Hist1D_y1");
  hPzsJunlee->SetName("hPzsJunlee");
  TH1F *hPzsJunleeStatErr = (TH1F *)dirJunlee->Get("Hist1D_y1_e1");
  hPzsJunleeStatErr->SetName("hPzsJunleeStatErr");
  TH1F *hPzsJunleeSystErr = (TH1F *)dirJunlee->Get("Hist1D_y1_e2");
  hPzsJunleeSystErr->SetName("hPzsJunleeSystErr");
  TH1F *hPzsJunleeSyst = (TH1F *)hPzsJunlee->Clone("hPzsJunleeSyst");
  for (Int_t i = 1; i <= hPzsJunlee->GetNbinsX(); i++)
  {
    hPzsJunlee->SetBinError(i, hPzsJunleeStatErr->GetBinContent(i));
    hPzsJunleeSyst->SetBinError(i, hPzsJunleeSystErr->GetBinContent(i));
  }
  hPzsJunlee->SetMarkerStyle(20);
  hPzsJunlee->SetMarkerSize(1.5);
  hPzsJunlee->SetMarkerColor(ColorJunlee);
  hPzsJunlee->SetLineColor(ColorJunlee);
  hPzsJunleeSyst->SetMarkerStyle(20);
  hPzsJunleeSyst->SetMarkerSize(1.5);
  hPzsJunleeSyst->SetMarkerColor(ColorJunlee);
  hPzsJunleeSyst->SetLineColor(ColorJunlee);

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = "../Pzs2VsPt/" + NameAnalysis[!isV2] + "_";
  stringout += SinputFileName;
  stringout += "_" + ParticleName[ChosenPart];
  if (isFitDSCB)
  {
    stringout += "_DSCB";
    if (isFixParamDSCBFromMC)
      stringout += "_FixParamFromMC";
  }
  else
    stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[BkgType];
  if (isGaussConv)
    stringout += "_GaussConv";
  if (isCombinedFit)
    stringout += "_CombinedFit";
  stringout += "_Pzs2";
  if (isApplyWeights)
    stringout += "_Weighted";
  if (isApplyCentWeight)
    stringout += "_CentWeighted";
  if (!useCommonBDTValue)
    stringout += "_BDTPtDep";
  if (isRun2Binning)
    stringout += "_Run2Binning";
  if (isPolFromLambda)
    stringout += "_PolFromLambda";
  if (!isRapiditySel)
    stringout += "_Eta08";
  stringout += STHN[ExtrisFromTHN];
  if (useMixedBDTValueInFitMacro)
    stringout += "_MixedBDT";
  if (isTightMassCut)
    stringout += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
  stringout += V2FromFit[isFromFit];
  if (isReducedPtBins)
    stringout += "_ReducedPtBins";
  if (isApplyAcceptanceInMacro)
    stringout += "_AcceptanceInMacro";
  if (ExtrisApplyResoOnTheFly)
    stringout += "_ResoOnTheFly";
  if (ChosenPart == 0)
    stringout += "_EPReso";
  if (isBkgPol == 0)
    stringout += "_isBkgPol0";
  // stringout += "_SystReso";
  if (isTighterPzFitRange)
    stringout += "_TighterPzFitRange";
  if (isFDCorrected)
    stringout += "_FDCorrected";
  if (ChosenPart >= 6 && ExtrisSysLambdaMultTrial)
  {
    if (isLoosest)
      stringout += "_isLoosest";
    else if (isTightest)
      stringout += "_isTightest";
    stringout += "_isSysLambdaMultTrial";
  }
  if (ExtrisApplyEffWeights)
    stringout += "_EffW";
  if (isRebin)
    stringout += "_Rebin";
  // stringout += "_TestMoreBins";
  if (isGlobalFit)
    stringout += "_GlobalFit";
  // stringout += "_NegativeC";
  if (ExtrisCentOmegaRed && part == 1)
    stringout += "_OmegaRedCent";
  // stringout += "_PtFix2Bins";
  if (ApplyK0sMassRej)
    stringout += "_K0sMassRej";
  if (!isT0CCentrality)
    stringout += "_T0MCentrality";
  stringoutpdf = stringout;
  stringout += ".root";

  // canvases
  gStyle->SetOptStat(0);
  TCanvas *canvasPzs = new TCanvas("canvasPzs", "canvasPzs", 900, 700);
  StyleCanvas(canvasPzs, 0.05, 0.15, 0.15, 0.05);

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *legendLambda = new TLegend(0.2, 0.63, 0.6, 0.83);
  legendLambda->SetFillStyle(0);
  legendLambda->SetTextSize(0.045);
  // legendLambda->AddEntry(fHistPzsLambda, "#Lambda + #bar{#Lambda}, Phys. Rev. Lett. 128.17 (2022)", "pl");

  TLegend *LegendTitle;
  LegendTitle = new TLegend(0.54, 0.72, 0.93, 0.9);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(33);
  LegendTitle->SetTextSize(0.04);
  LegendTitle->AddEntry("", "#bf{ALICE}", "");
  /*
  LegendTitle->AddEntry("", "PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
  if (isPolFromLambda)
  {
    if (isRapiditySel)
      LegendTitle->AddEntry("", ParticleNameLegend[ChosenPart] + " from daughter #Lambda, |#it{y}| < 0.5", "");
    else
      LegendTitle->AddEntry("", ParticleNameLegend[ChosenPart] + " from daughter #Lambda, |#it{#eta}| < 0.8", "");
  }
  else
  {
    if (isRapiditySel)
      LegendTitle->AddEntry("", ParticleNameLegend[ChosenPart] + " |#it{y}| < 0.5", "");
    else
    {
      LegendTitle->AddEntry("", ParticleNameLegend[ChosenPart] + " |#it{#eta}| < 0.8", "");
    }
  }
    */

  TH1F *fHistPzs;
  TH1F *fHistPzsPtInt;
  TH1F *fHistPzsBkg;
  TH1F *fHistPurity;
  TH1F *fHistSignificance;
  TH1F *fHistYield;
  TH1F *fHistMean;
  TH1F *fHistSigma;
  TH1F *fHistB;
  TH1F *fHistTot;
  Int_t CentFT0CMin = 0;
  Int_t CentFT0CMax = CentFT0CMaxLambdaOO;

  // filein
  TString PathIn = "../OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_";
  TString PathIn1Bin = PathIn;
  PathIn += SinputFileName;
  if (SinputFileName == "LHC25_OO_pass2_Train742311_TwoBins")
    PathIn1Bin += "LHC25_OO_pass2_Train742311";
  TString PathInAdd = "_" + ParticleName[ChosenPart];
  if (isFitDSCB)
  {
    PathInAdd += "_DSCB";
    if (isFixParamDSCBFromMC)
      PathInAdd += "_FixParamFromMC";
  }
  else
    PathInAdd += IsOneOrTwoGauss[UseTwoGauss];
  PathInAdd += SIsBkgParab[BkgType];
  if (isGaussConv)
    PathInAdd += "_GaussConv";
  if (isCombinedFit)
    PathInAdd += "_CombinedFit";
  TString Smolt = Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);
  TString SmoltBis = Form("%i#minus%i", CentFT0CMin, CentFT0CMax);
  PathInAdd += Smolt;
  if (isApplyWeights)
    PathInAdd += "_Weighted";
  if (isApplyCentWeight)
    PathInAdd += "_CentWeighted";
  if (!useCommonBDTValue)
    PathInAdd += "_BDTPtDep";
  if (isRun2Binning)
    PathInAdd += "_Run2Binning";
  if (isPolFromLambda)
    PathInAdd += "_PolFromLambda";
  if (ExtrisApplyEffWeights)
    PathInAdd += "_EffW";
  if (!isRapiditySel)
    PathInAdd += "_Eta08";
  PathInAdd += STHN[ExtrisFromTHN];
  if (useMixedBDTValueInFitMacro)
    PathInAdd += "_MixedBDT";
  if (isTightMassCut)
    PathInAdd += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
  if (isReducedPtBins)
    PathInAdd += "_ReducedPtBins";
  if (ChosenPart >= 6 && ExtrisSysLambdaMultTrial)
  {
    if (isLoosest)
      PathInAdd += "_isLoosest";
    else if (isTightest)
      PathInAdd += "_isTightest";
    PathInAdd += "_isSysLambdaMultTrial";
  }
  if (isApplyAcceptanceInMacro)
    // PathInAdd += "_AcceptanceInMacroAllLambdas";
    PathInAdd += "_AcceptanceInMacro";
  if (ExtrisApplyResoOnTheFly)
    PathInAdd += "_ResoOnTheFly";
  if (ChosenPart == 0)
    PathInAdd += "_EPReso";
  if (isBkgPol == 0)
    PathInAdd += "_isBkgPol0";
  // PathInAdd += "_SystReso";
  if (isTighterPzFitRange)
    PathInAdd += "_TighterPzFitRange";
  if (isRebin)
    PathInAdd += "_Rebin";
  // PathInAdd += "_NegativeC";
  if (ExtrisCentOmegaRed && part == 1)
    PathInAdd += "_OmegaRedPt";
  if (ChosenPart >= 6 && SinputFileName == "LHC25_OO_pass2_Train598890_MyEff")
    PathInAdd += "_050PtCut";
  // PathInAdd += "_PtFix2Bins";
  if (ApplyK0sMassRej)
    PathInAdd += "_K0sMassRej";
  if (!isT0CCentrality)
    PathInAdd += "_PtFix_T0MCentrality";
  PathInAdd += ".root";
  PathIn += PathInAdd;
  PathIn1Bin += PathInAdd;
  cout << "Path in : " << PathIn << endl;

  TFile *fileIn = TFile::Open(PathIn);
  TFile *fileIn1Bin = TFile::Open(PathIn1Bin);

  fHistPzs = (TH1F *)fileIn->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + V2FromFit[isFromFit]);
  if (isGlobalFit)
    fHistPzs = (TH1F *)fileIn->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + "Global" + V2FromFit[isFromFit]);
  TH1F *fHistPzs1Bin = nullptr;
  if (SinputFileName == "LHC25_OO_pass2_Train742311_TwoBins")
    fHistPzs1Bin = (TH1F *)fileIn1Bin->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + V2FromFit[isFromFit]);
  else
    fHistPzs1Bin = (TH1F *)fHistPzs->Clone("fHistPzs1Bin");
  fHistPzsPtInt = (TH1F *)fileIn->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + "PtInt");
  fHistPzsBkg = (TH1F *)fileIn->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + "Bkg");
  fHistPurity = (TH1F *)fileIn->Get("histoPurity");
  fHistSignificance = (TH1F *)fileIn->Get("histoSignificance");
  fHistYield = (TH1F *)fileIn->Get("histoYield");
  fHistMean = (TH1F *)fileIn->Get("histoMean");
  fHistSigma = (TH1F *)fileIn->Get("histoSigma");
  fHistB = (TH1F *)fileIn->Get("histoB");
  fHistTot = (TH1F *)fileIn->Get("histoTot");
  if (!fHistPzs)
  {
    cout << " no hist v2 / Pzs" << endl;
    return;
  }
  if (!fHistPzsBkg)
  {
    cout << " no hist v2 / Pzs Bkg" << endl;
    return;
  }
  if (!fHistPurity)
  {
    cout << " no hist purity" << endl;
    return;
  }
  if (!fHistSignificance)
  {
    cout << " no hist significance " << endl;
    return;
  }
  if (!fHistYield)
  {
    cout << " no hist yield " << endl;
    return;
  }
  if (!fHistMean)
  {
    cout << " no hist mean" << endl;
    return;
  }
  if (!fHistSigma)
  {
    cout << " no hist sigma" << endl;
    return;
  }
  if (!fHistB)
  {
    cout << " no hist B" << endl;
    return;
  }
  if (!fHistTot)
  {
    cout << " no hist Tot" << endl;
    return;
  }

  // Get histogram with syst uncertainty
  TString PathInSyst = "../Systematics/SystVsPt_" + NameAnalysis[!isV2] + "_";
  PathInSyst += SinputFileNameSyst;
  if (ChosenPart == 7 || ChosenPart == 8)
    PathInSyst += "_" + ParticleName[6];
  else
    PathInSyst += "_" + ParticleName[ChosenPart];
  if (ExtrisFitDSCB)
  {
    PathInSyst += "_DSCB";
    if (isFixParamDSCBFromMC)
      PathInSyst += "_FixParamFromMC";
  }
  else
    PathInSyst += IsOneOrTwoGauss[UseTwoGauss];
  PathInSyst += SIsBkgParab[ExtrBkgTypeSyst];
  PathInSyst += "_Pzs2";
  if (isApplyWeights)
    PathInSyst += "_Weighted";
  if (isApplyCentWeight || ChosenPart >= 6)
    PathInSyst += "_CentWeighted";
  if (!useCommonBDTValue)
    PathInSyst += "_BDTCentDep";
  if (isRun2Binning)
    PathInSyst += "_Run2Binning";
  if (isPolFromLambda)
    PathInSyst += "_PolFromLambda";
  if (!isRapiditySel)
    PathInSyst += "_Eta08";
  if (ChosenPart < 6)
    PathInSyst += STHN[ExtrisFromTHN];
  if (useMixedBDTValueInFitMacro)
    PathInSyst += "_MixedBDT";
  if (isTightMassCut)
    // PathInSyst += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
    PathInSyst += Form("_TightMassCut%.1f", 2.1);
  // PathInSyst += V2FromFit[isFromFit];
  if (isReducedPtBins)
    PathInSyst += "_ReducedPtBins";
  if (ExtrisApplyResoOnTheFly || ChosenPart >= 6)
    PathInSyst += "_ResoOnTheFly";
  if (ChosenPart == 0)
    PathInSyst += "_New";
  PathInSyst += ".root";
  cout << "Opening systematic uncertainty file: " << PathInSyst << endl;
  TFile *fileInSyst = TFile::Open(PathInSyst);
  if (!fileInSyst)
  {
    cout << "No file found" << endl;
    return;
  }
  TH1F *fHistPzsSistError = (TH1F *)fileInSyst->Get("fHistTotalError");
  if (!fHistPzsSistError)
  {
    cout << "No hist syst error found" << endl;
    return;
  }
  TH1F *fHistPzsSist = (TH1F *)fHistPzs->Clone("fHistPzsSist");
  TH1F *fHistPzsSist1Bin = (TH1F *)fHistPzs1Bin->Clone("fHistPzsSist1Bin");
  if (!isFromFit)
  {
    for (Int_t pt = 1; pt < fHistPzs->GetNbinsX(); pt++)
    {
      fHistPzs->SetBinContent(pt, fHistPzs->GetBinContent(pt) / fHistPurity->GetBinContent(pt));
      fHistPzs->SetBinError(pt, fHistPzs->GetBinError(pt) / fHistPurity->GetBinContent(pt));
      fHistPzsSist->SetBinContent(pt, fHistPzs->GetBinContent(pt) / fHistPurity->GetBinContent(pt));
      fHistPzsSist->SetBinError(pt, fHistPzs->GetBinError(pt) / fHistPurity->GetBinContent(pt));
    }
  }
  for (Int_t b = 1; b <= fHistPzsSist1Bin->GetNbinsX(); b++)
  {
    fHistPzsSist1Bin->SetBinContent(b, fHistPzs1Bin->GetBinContent(b));
    fHistPzsSist1Bin->SetBinError(b, fHistPzsSistError->GetBinContent(b));
  }
  for (Int_t b = 1; b <= fHistPzs->GetNbinsX(); b++)
  {
    // fHistPzsSistError->SetBinContent(b, fHistPzsSistError->GetBinContent(b) / fHistPuritySummary->GetBinContent(b));
    fHistPzsSist->SetBinContent(b, fHistPzs->GetBinContent(b));
    if (SinputFileName == "LHC25_OO_pass2_Train742311_TwoBins")
    {
      int binValue = fHistPzsSistError->FindBin(fHistPzsSist->GetBinCenter(b)); // bin corresponding to the center of the current bin in fHistPzsSist
      if (b == 1)
        fHistPzsSist->SetBinError(b, fHistPzsSistError->GetBinContent(binValue) / fHistPzs1Bin->GetBinContent(binValue) * fHistPzsSist->GetBinContent(b));
      else
        fHistPzsSist->SetBinError(b, fHistPzsSistError->GetBinContent(binValue));
    }
    else
    {
      fHistPzsSist->SetBinError(b, fHistPzsSistError->GetBinContent(b));
    }
  }

  TH1F *fHistPzsError = (TH1F *)fHistPzs->Clone("fHistPzsError");
  TH1F *fHistPzsSignif = (TH1F *)fHistPzs->Clone("fHistPzs");
  TH1F *fHistPzsSignifStat = (TH1F *)fHistPzs->Clone("fHistPzs");
  TH1F *fHistPzsSignifLambda = (TH1F *)fHistPzs->Clone("fHistPzs");
  for (Int_t b = 1; b <= fHistPzs->GetNbinsX(); b++)
  {
    fHistPzsError->SetBinContent(b, fHistPzs->GetBinError(b));
    fHistPzsError->SetBinError(b, 0);
    fHistPzsSistError->SetBinContent(b, fHistPzsSist->GetBinError(b));
    fHistPzsSistError->SetBinError(b, 0);
    fHistPzsSignifStat->SetBinContent(b, abs(fHistPzs->GetBinContent(b)) / fHistPzs->GetBinError(b));
    fHistPzsSignif->SetBinContent(b, abs(fHistPzs->GetBinContent(b)) / TMath::Sqrt(fHistPzs->GetBinError(b) * fHistPzs->GetBinError(b) + fHistPzsSistError->GetBinContent(b) * fHistPzsSistError->GetBinContent(b)));
    fHistPzsSignifLambda->SetBinContent(b, abs(fHistPzsLambda->GetBinContent(b)) / TMath::Sqrt(fHistPzsLambda->GetBinError(b) * fHistPzsLambda->GetBinError(b) + fHistPzsLambda_SystErr->GetBinContent(b) * fHistPzsLambda_SystErr->GetBinContent(b)));
    if (fHistPzsSignifLambda->GetBinCenter(b) > 80)
    {
      fHistPzsSignifLambda->SetBinContent(b, -1000);
    }
    // cout << "Signif. Lambda PbPb " << fHistPzsSignifLambda->GetBinContent(b) << endl;
  }

  Float_t xTitle = 30;
  Float_t xOffset = 1.3;
  Float_t yTitle = 38;   // 30
  Float_t yOffset = 1.1; // 2.2 if setmaxdigits not set

  Float_t xLabel = 30;
  Float_t yLabel = 30;
  Float_t xLabelOffset = 0.015;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.03;
  Float_t tickY = 0.025;

  TLegend *LegendPreliminary2;
  LegendPreliminary2 = new TLegend(0.06, 0.77, 0.45, 0.914);
  LegendPreliminary2->SetFillStyle(0);
  LegendPreliminary2->SetTextAlign(11);
  LegendPreliminary2->SetTextSize(0.048);
  LegendPreliminary2->AddEntry("", "#bf{ALICE}", "");
  // LegendPreliminary2->AddEntry("", "#bf{ALICE Work In Progress}", "");
  // LegendPreliminary2->AddEntry("", "Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
  if (ChosenPart >= 6)
    LegendPreliminary2->AddEntry("", "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
  else
    LegendPreliminary2->AddEntry("", "Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");

  TLegend *LegendPreliminary3;
  LegendPreliminary3 = new TLegend(0.05, 0.85, 0.45, 0.93);
  LegendPreliminary3->SetFillStyle(0);
  LegendPreliminary3->SetTextAlign(11);
  LegendPreliminary3->SetTextSize(0.048);
  if (ChosenPart >= 6)
    LegendPreliminary3->AddEntry("", "#bf{ALICE}", "");
  else
    LegendPreliminary3->AddEntry("", "ALICE, Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");

  TLegend *legendXi = new TLegend(0.06, 0.64, 0.45, 0.78);
  legendXi->SetFillStyle(0);
  legendXi->SetTextAlign(12);
  legendXi->SetTextSize(0.048);
  if (ChosenPart == 6)
    legendXi->AddEntry("", Form("#Lambda + #bar{#Lambda}, |#it{#eta} | < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "");
  else if (ChosenPart == 7)
    legendXi->AddEntry("", Form("#Lambda, |#it{#eta}| < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "");
  else if (ChosenPart == 8)
    legendXi->AddEntry("", Form("#bar{#Lambda}, |#it{#eta}| < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "");
  else if (part == 0)
    legendXi->AddEntry("", Form("#Xi^{#minus} + #bar{#Xi}^{+}, |#it{#eta} | < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "");
  else if (part == 1)
    legendXi->AddEntry("", Form("#Omega^{#minus} + #bar{#Omega}^{+}, |#it{#eta}| < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "");

  // Theory /////////////////////////////////////////
  TGraphErrors *gPzsPalermo = new TGraphErrors(20);
  for (Int_t i = 0; i < gPzsPalermo->GetN(); i++)
  {
    // cout << "Palermo: i " << i << " " << dNdEtaPalermo[gPzsPalermo->GetN() - i - 1] << " " << gPzsPalermo->GetY()[gPzsPalermo->GetN() - i - 1] << endl;
    gPzsPalermo->SetPoint(i, pTPalermo[i], Pzs2VsPtPalermo[i]);
    gPzsPalermo->SetPointError(i, 0, 0);
  }
  gPzsPalermo->SetLineColor(kBlue);
  gPzsPalermo->SetLineWidth(2);

  const Int_t npoints = 55;
  Double_t PtValues[npoints];
  Double_t PzsVsPt502[npoints];

  ifstream file("../Fig_5_Prediction_1-1-2_model_PbPb_5020_30-50_LocalSpinPol_z_s2_vs_pT.dat");

  if (!file)
  {
    cout << "Error opening file." << endl;
    return;
  }

  int count = 0;

  while (file >> PtValues[count] >> PzsVsPt502[count])
  {
    count++;
  }

  file.close();
  TGraph *gPzsVsPt502 = new TGraph(npoints);
  for (Int_t i = 0; i < gPzsVsPt502->GetN(); i++)
  {
    gPzsVsPt502->SetPoint(i, PtValues[i], PzsVsPt502[i] / 100);
  }
  gPzsVsPt502->SetLineColor(kBlack);
  gPzsVsPt502->SetLineWidth(2);

  Double_t PtValues536[npoints];
  Double_t PzsVsPt536[npoints];

  ifstream file536("../Fig_5_Prediction_1-1-2_model_PbPb_5360_30-50_LocalSpinPol_z_s2_vs_pT.dat");

  if (!file536)
  {
    cout << "Error opening file." << endl;
    return;
  }

  int count536 = 0;

  while (file536 >> PtValues536[count536] >> PzsVsPt536[count536])
  {
    count536++;
  }

  file536.close();
  TGraph *gPzsVsPt536 = new TGraph(npoints);
  for (Int_t i = 0; i < gPzsVsPt536->GetN(); i++)
  {
    gPzsVsPt536->SetPoint(i, PtValues536[i], PzsVsPt536[i] / 100);
  }
  gPzsVsPt536->SetLineColor(kGreen + 2);
  gPzsVsPt536->SetLineWidth(2);

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, -1000);
  canvasPzs->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleXPt, TitleYPzs, "", 1, 1.15, 1.6);
  StyleHistoYield(fHistPzs, YLow[part], YUp[part], ColorPart[part], MarkerPart[part], TitleXPt, TitleYPzs, "", MarkerPartSize[part], 1.15, 1.6);
  StyleHistoYield(fHistPzsSist, YLow[part], YUp[part], ColorPart[part], MarkerPart[part], TitleXPt, TitleYPzs, "", 1.5, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  if (ChosenPart >= 6)
  {
    fHistPzs->SetLineColor(ColorOO);
    fHistPzs->SetMarkerColor(ColorOO);
    fHistPzsSist->SetLineColor(ColorOO);
    fHistPzsSist->SetMarkerColor(ColorOO);
  }
  hDummy->GetXaxis()->SetRangeUser(0, UpperRangeParticle);
  hDummy->GetYaxis()->SetTitleOffset(1.4);
  hDummy->Draw("");
  fHistPzs->Draw("same ex0");
  fHistPzsSist->SetFillStyle(0);
  fHistPzsSist->Draw("same e2");
  fHistPzsLambda->Draw("same ex0");
  fHistPzsLambdaSist->SetFillStyle(0);
  fHistPzsLambdaSist->Draw("same e2");
  hPzsJunlee->Draw("same ex0");
  hPzsJunleeSyst->SetFillStyle(0);
  hPzsJunleeSyst->Draw("same e2");
  gPzsPalermo->Draw("same l");
  gPzsVsPt536->Draw("same l");
  gPzsVsPt502->Draw("same l");
  // LegendTitle->Draw("");
  // legendLambda->AddEntry(hPzsJunlee, Form("Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV, 30-50%, |#it{y}| < 0.5, #it{p}_{T} > %1.1f GeV/#it{c}", 0.5), "pl");
  // legendLambda->AddEntry(fHistPzs, Form("OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 0-50%, |#it{#eta}| < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "pl");
  legendLambda->AddEntry(fHistPzs, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 0-50\%, |#it{#eta} | < 0.8", "pl");
  // legendLambda->AddEntry(fHistPzsLambda, "Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.02, 30-50\%, |#it{y}  < 0.5", "pl");
  legendLambda->AddEntry(hPzsJunlee, "Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV, 30-50\%, |#it{y} | < 0.5", "pl");
  legendLambda->Draw("");
  canvasPzs->SaveAs("../PzsVsPt_PbPbComp.pdf");
  canvasPzs->SaveAs("../PzsVsPt_PbPbComp.png");
  canvasPzs->SaveAs("../PzsVsPt_PbPbComp.eps");

  // Pzs vs pT in OO and PbPb + ratio
  TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 800, 900);
  TPad *pTop = new TPad("pTop", "", 0, 0.35, 1, 1);
  TPad *pBot = new TPad("pBot", "", 0, 0, 1, 0.35);
  pTop->SetBottomMargin(0.02);
  pTop->SetLeftMargin(0.15);
  pTop->SetRightMargin(0.03);
  pTop->SetTopMargin(0.05);
  pBot->SetTopMargin(0.03);
  pBot->SetBottomMargin(0.3);
  pBot->SetLeftMargin(0.15);
  pBot->SetRightMargin(0.03);
  pTop->Draw();
  pBot->Draw();

  pTop->cd();
  TH1F *hTop = new TH1F("hTop", ";;#it{P}_{z,s2}", 100, 0, 6);
  hTop->GetYaxis()->SetRangeUser(-0.001, 0.0105);
  hTop->GetYaxis()->SetTitleSize(0.065);
  hTop->GetYaxis()->SetTitleOffset(1.05);
  hTop->GetYaxis()->SetLabelSize(0.05);
  hTop->GetYaxis()->SetMaxDigits(3);
  hTop->GetXaxis()->SetLabelSize(0);
  hTop->SetLineWidth(0);
  hTop->Draw("axis");
  TLine *l0 = new TLine(0, 0, 6.2, 0);
  l0->SetLineStyle(2);
  l0->Draw();
  fHistPzs->Draw("same ex0");
  fHistPzsSist->SetFillStyle(0);
  fHistPzsSist->Draw("same e2");
  //fHistPzsLambda->Draw("same ex0");
  //fHistPzsLambdaSist->SetFillStyle(0);
  //fHistPzsLambdaSist->Draw("same e2");
  hPzsJunlee->Draw("same ex0");
  hPzsJunleeSyst->Draw("same e2");
  // gPzsPalermo->Draw("same l");
  gPzsVsPt536->Draw("same l");
  // gPzsVsPt502->Draw("same l");
  legendLambda->AddEntry(gPzsVsPt536, "Model, Pb#minusPb #sqrt{s_{NN}} = 5.36 TeV", "l");
  legendLambda->Draw("");

  TLegend *LegendTitleRatio;
  LegendTitleRatio = new TLegend(0.17, 0.79, 0.3, 0.95);
  LegendTitleRatio->SetFillStyle(0);
  LegendTitleRatio->SetTextAlign(12);
  LegendTitleRatio->SetTextSize(0.05);
  LegendTitleRatio->AddEntry("", "#bf{ALICE}", "");
  LegendTitleRatio->Draw("");

  pBot->cd();
  TH1F *hBot = new TH1F("hBot", ";#it{p}_{T} (GeV/#it{c});Pb-Pb / OO", 100, 0, 6);
  hBot->GetYaxis()->SetRangeUser(-0.5, 2.5);
  //  hBot->GetYaxis()->SetRangeUser(-2, 11);
  hBot->GetYaxis()->SetNdivisions(505);
  hBot->GetYaxis()->SetTitleSize(0.1);
  hBot->GetYaxis()->SetTitleOffset(0.62);
  hBot->GetYaxis()->SetLabelSize(0.09);
  hBot->GetXaxis()->SetTitleSize(0.11);
  hBot->GetXaxis()->SetLabelSize(0.09);
  hBot->SetLineWidth(0);
  hBot->Draw("axis");
  TH1F *hRatioPbPbToOO = (TH1F *)fHistPzs1Bin->Clone("hRatioPbPbToOO");
  TH1F *hRatioPbPbToOOSist = (TH1F *)fHistPzsSist1Bin->Clone("hRatioPbPbToOOSist");
  StyleHistoYield(hRatioPbPbToOO, YLow[part], YUp[part], ColorPart[part], MarkerPart[part], TitleXPt, TitleYPzs, "", MarkerPartSize[part], 1.15, 1.6);
  StyleHistoYield(hRatioPbPbToOOSist, YLow[part], YUp[part], ColorPart[part], MarkerPart[part], TitleXPt, TitleYPzs, "", 1.5, 1.15, 1.6);
  for (Int_t i = 1; i <= hRatioPbPbToOO->GetNbinsX(); i++)
  {
    cout << "Bin center : " << hRatioPbPbToOO->GetBinCenter(i) << " Bin content : " << hRatioPbPbToOO->GetBinContent(i) << endl;
    cout << "junlee bin center: " << hPzsJunlee->GetBinCenter(i) << " Bin content : " << hPzsJunlee->GetBinContent(i) << endl;
    hRatioPbPbToOO->SetBinContent(i, hPzsJunlee->GetBinContent(i) / fHistPzs1Bin->GetBinContent(i));
    hRatioPbPbToOO->SetBinError(i, hRatioPbPbToOO->GetBinContent(i) * TMath::Sqrt(TMath::Power(hPzsJunlee->GetBinError(i) / hPzsJunlee->GetBinContent(i), 2) + TMath::Power(fHistPzs1Bin->GetBinError(i) / fHistPzs1Bin->GetBinContent(i), 2)));
    hRatioPbPbToOOSist->SetBinContent(i, hRatioPbPbToOO->GetBinContent(i));
    hRatioPbPbToOOSist->SetBinError(i, hRatioPbPbToOO->GetBinContent(i) * TMath::Sqrt(TMath::Power(hPzsJunleeSyst->GetBinError(i) / hPzsJunlee->GetBinContent(i), 2) + TMath::Power(fHistPzsSist1Bin->GetBinError(i) / fHistPzsSist1Bin->GetBinContent(i), 2)));
  }
  hRatioPbPbToOO->SetLineColor(ColorJunlee);
  hRatioPbPbToOO->SetMarkerColor(ColorJunlee);
  hRatioPbPbToOO->Draw("same ex0");
  hRatioPbPbToOOSist->SetLineColor(ColorJunlee);
  hRatioPbPbToOOSist->SetMarkerColor(ColorJunlee);
  hRatioPbPbToOOSist->SetFillStyle(0);
  hRatioPbPbToOOSist->Draw("same e2");

  // pT integrated ratio
  TString SfilePbPbJunleeVsCent = "../LambdaJunlee/HEPData-ins3169170-v1-Table_3.root";
  TFile *fileJunleeVsCent = new TFile(SfilePbPbJunleeVsCent);
  TDirectoryFile *dirJunleeVsCent = (TDirectoryFile *)fileJunleeVsCent->Get("Table 3");
  TH1F *hPzsJunleeVsCent = (TH1F *)dirJunleeVsCent->Get("Hist1D_y1");
  Double_t PzsPbPb3050 = (hPzsJunleeVsCent->GetBinContent(hPzsJunleeVsCent->FindBin(35)) * LambdaYield[0] +
                          hPzsJunleeVsCent->GetBinContent(hPzsJunleeVsCent->FindBin(45)) * LambdaYield[1]) /
                         (LambdaYield[0] + LambdaYield[1]);
  Double_t ratioPbPbToOO = PzsPbPb3050 / fHistPzsPtInt->GetBinContent(1);
  cout << "ratioPbPbToOO = " << ratioPbPbToOO << endl;
  TF1 *fRatioPbPbToOO = new TF1("fRatioPbPbToOO", "pol0", 0, 100);
  fRatioPbPbToOO->FixParameter(0, ratioPbPbToOO);
  fRatioPbPbToOO->SetLineColor(ColorJunlee);
  fRatioPbPbToOO->SetMarkerColor(ColorJunlee);
  fRatioPbPbToOO->SetLineStyle(2);
  fRatioPbPbToOO->SetRange(0.5, 6);
  fRatioPbPbToOO->Draw("same");
  TLegend *legendRatio = new TLegend(0.19, 0.72, 0.58, 0.98);
  legendRatio->SetFillStyle(0);
  legendRatio->SetTextAlign(12);
  legendRatio->SetTextSize(0.068);
  legendRatio->AddEntry("fRatioPbPbToOO", "Pb#minusPb 30-50\%/ OO 0-50\%, #it{p}_{T}-integrated", "l");
  legendRatio->Draw("");
  cRatio->Modified();
  cRatio->Update();
  cRatio->SaveAs("../PzsVsPt_RatioPbPbToOO.eps");
  cRatio->SaveAs("../PzsVsPt_RatioPbPbToOO.png");
  cRatio->SaveAs("../PzsVsPt_RatioPbPbToOO.pdf");

  // Relative stat. uncertainty
  TCanvas *canvasPzsError = new TCanvas("canvasPzsError", "canvasPzsError", 900, 700);
  StyleCanvas(canvasPzsError, 0.03, 0.15, 0.15, 0.05);
  TH1F *hDummyError = new TH1F("hDummyError", "hDummyError", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummyError->GetNbinsX(); i++)
    hDummyError->SetBinContent(i, 1e-12);
  canvasPzsError->cd();
  SetFont(hDummyError);
  StyleHistoYield(hDummyError, 0, 0.004, 1, 1, TitleXPt, "Absolute uncertainty", "", 1, 1.15, 1.6);
  if (part == 1)
    hDummyError->GetYaxis()->SetRangeUser(0, 0.015);
  StyleHistoYield(fHistPzsError, 0, 0.01, ColorPart[part], MarkerPart[part], TitleXPt, "", "", MarkerPartSize[part], 1.15, 1.6);
  StyleHistoYield(fHistPzsSistError, 0, 0.01, kGray + 2, 22, TitleXPt, "", "", MarkerPartSize[part], 1.15, 1.6);
  SetHistoTextSize(hDummyError, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  hDummyError->GetYaxis()->SetTitleOffset(1.7);
  SetTickLength(hDummyError, tickX, tickY);
  hDummyError->GetXaxis()->SetRangeUser(LowerRangeParticle, UpperRangeParticle);
  hDummyError->Draw("");
  fHistPzsError->Draw("same");
  fHistPzsSistError->Draw("same e2");
  fHistPzsLambda_StatErr->SetLineColor(kRed);
  fHistPzsLambda_StatErr->SetMarkerColor(kRed);
  fHistPzsLambda_StatErr->Draw("same e0x0");
  // LegendTitle->Draw("");
  // legendLambda->Draw("");
  // LegendPreliminary3->Draw("");

  TLegend *legendError = new TLegend(0.19, 0.57, 0.58, 0.83);
  legendError->SetFillStyle(0);
  legendError->SetTextAlign(12);
  legendError->SetTextSize(0.048);
  if (ChosenPart >= 6)
  {
    legendError->AddEntry(fHistPzsError, "stat. #Lambda + #bar{#Lambda} Run 3", "pl");
    legendError->AddEntry(fHistPzsSistError, "syst. #Lambda + #bar{#Lambda} Run 3", "pl");
  }
  else if (part == 0)
  {
    legendError->AddEntry(fHistPzsError, "stat. #Xi^{#minus} + #bar{#Xi}^{+} Run 3", "pl");
    legendError->AddEntry(fHistPzsSistError, "syst. #Xi^{#minus} + #bar{#Xi}^{+} Run 3", "pl");
  }
  else if (part == 1)
  {
    legendError->AddEntry(fHistPzsError, "stat. #Omega^{#minus} + #bar{#Omega}^{+} Run 3", "pl");
    legendError->AddEntry(fHistPzsSistError, "syst. #Omega^{#minus} + #bar{#Omega}^{+} Run 3", "pl");
  }
  legendError->AddEntry(fHistPzsLambda_StatErr, "stat. #Lambda + #bar{#Lambda} Run 2", "pl");
  legendError->Draw("");
  canvasPzsError->SaveAs(stringoutpdf + "_Error.pdf");
  canvasPzsError->SaveAs(stringoutpdf + "_Error.png");

  // Significance
  TCanvas *canvasPzsSignif = new TCanvas("canvasPzsSignif", "canvasPzsSignif", 900, 700);
  StyleCanvas(canvasPzsSignif, 0.05, 0.15, 0.15, 0.05);
  TH1F *hDummySignif = new TH1F("hDummySignif", "hDummySignif", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummySignif->GetNbinsX(); i++)
    hDummySignif->SetBinContent(i, 1e-12);
  canvasPzsSignif->cd();
  SetFont(hDummySignif);
  StyleHistoYield(hDummySignif, 0, 1.5 * fHistPzsSignif->GetBinContent(fHistPzsSignif->GetMaximumBin()), 1, 1, TitleXPt, "S / #sigma_{S}", "", 1, 1.15, 1.6);
  StyleHistoYield(fHistPzsSignif, 0, 1.5 * fHistPzsSignif->GetBinContent(fHistPzsSignif->GetMaximumBin()), ColorPart[part], MarkerPart[part], TitleXPt, "S / #sigma_{S}", "", MarkerPartSize[part], 1.15, 1.6);
  StyleHistoYield(fHistPzsSignifStat, 0, 1.5 * fHistPzsSignif->GetBinContent(fHistPzsSignif->GetMaximumBin()), kGray + 2, 22, TitleXPt, "S / #sigma_{S}", "", MarkerPartSize[part], 1.15, 1.6);
  StyleHistoYield(fHistPzsSignifLambda, 0, 1.5 * fHistPzsSignif->GetBinContent(fHistPzsSignif->GetMaximumBin()), kBlue, 20, TitleXPt, "S / #sigma_{S}", "", MarkerPartSize[part], 1.15, 1.6);
  SetHistoTextSize(hDummySignif, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummySignif, tickX, tickY);
  hDummySignif->GetXaxis()->SetRangeUser(LowerRangeParticle, UpperRangeParticle);
  hDummySignif->GetYaxis()->SetRangeUser(0, 10);
  TH1F *fHistPzsSignifUpTo50 = (TH1F *)fHistPzsSignif->Clone("fHistPzsSignifUpTo50");
  TH1F *fHistPzsSignifStatUpTo50 = (TH1F *)fHistPzsSignifStat->Clone("fHistPzsSignifStatUpTo50");
  for (Int_t b = 1; b <= fHistPzsSignifUpTo50->GetNbinsX(); b++)
  {
    if (fHistPzsSignifUpTo50->GetXaxis()->GetBinCenter(b) > UpperRangeParticle)
      fHistPzsSignifUpTo50->SetBinContent(b, -1000);
    if (fHistPzsSignifStatUpTo50->GetXaxis()->GetBinCenter(b) > UpperRangeParticle)
      fHistPzsSignifStatUpTo50->SetBinContent(b, -1000);
  }
  hDummySignif->Draw("");
  fHistPzsSignifUpTo50->Draw("same");
  fHistPzsSignifStatUpTo50->Draw("same");

  fHistPzsSignifLambda->Draw("same e0x0");
  //  LegendTitle->Draw("");

  TString titleLambda = "#Lambda + #bar{#Lambda}";
  if (ChosenParticle == 7)
    titleLambda = "#Lambda";
  if (ChosenParticle == 8)
    titleLambda = "#bar{#Lambda}";
  TLegend *legendSignif = new TLegend(0.19, 0.66, 0.58, 0.92);
  legendSignif->SetFillStyle(0);
  legendSignif->SetTextAlign(12);
  legendSignif->SetTextSize(0.048);
  if (ChosenParticle >= 6)
  {
    legendSignif->AddEntry(fHistPzsSignif, "stat. + syst. " + titleLambda + " Run 3", "pl");
    legendSignif->AddEntry(fHistPzsSignifStat, "stat. " + titleLambda + " Run 3", "pl");
  }
  else if (part == 0)
  {
    legendSignif->AddEntry(fHistPzsSignif, "stat. + syst. #Xi^{#minus} + #bar{#Xi}^{+} Run 3", "pl");
    legendSignif->AddEntry(fHistPzsSignifStat, "stat. #Xi^{#minus} + #bar{#Xi}^{+} Run 3", "pl");
    legendSignif->AddEntry(fHistPzsSignifLambda, "stat. #Lambda + #bar{#Lambda} Run 2", "pl");
  }
  else if (part == 1)
  {
    legendSignif->AddEntry(fHistPzsSignif, "stat. + syst. #Omega^{#minus} + #bar{#Omega}^{+} Run 3", "pl");
    legendSignif->AddEntry(fHistPzsSignifStat, "stat. #Omega^{#minus} + #bar{#Omega}^{+} Run 3", "pl");
    legendSignif->AddEntry(fHistPzsSignifLambda, "stat. #Lambda + #bar{#Lambda} Run 2", "pl");
  }
  legendSignif->Draw("");
  canvasPzsSignif->SaveAs("../PzsVsPt_Signif.pdf");
  canvasPzsSignif->SaveAs("../PzsVsPt_Signif.png");

  TCanvas *canvasPurity = new TCanvas("canvasPurity", "canvasPurity", 900, 700);
  StyleCanvas(canvasPurity, 0.05, 0.15, 0.15, 0.05);
  canvasPurity->cd();
  TH1F *hDummyPurity = (TH1F *)hDummy->Clone("hDummyPurity");
  hDummyPurity->GetYaxis()->SetRangeUser(0.9, 1);
  if (part == 1)
    hDummyPurity->GetYaxis()->SetRangeUser(0.5, 1);
  hDummyPurity->GetYaxis()->SetTitle("S / (S+B)");
  hDummyPurity->Draw("");
  StyleHistoYield(fHistPurity, 0.9, 1, ColorPart[part], MarkerPart[part], TitleXPt, "S / (S+B)", "", MarkerPartSize[part], 1.15, 1.6);
  hDummyPurity->GetYaxis()->SetTitleOffset(1.3);
  fHistPurity->Draw("same");
  canvasPurity->SaveAs(stringoutpdf + "_Purity.pdf");
  canvasPurity->SaveAs(stringoutpdf + "_Purity.png");

  TCanvas *canvasSignificance = new TCanvas("canvasSignificance", "canvasSignificance", 900, 700);
  StyleCanvas(canvasSignificance, 0.05, 0.15, 0.15, 0.05);
  canvasSignificance->cd();
  TH1F *hDummySignificance = (TH1F *)hDummy->Clone("hDummySignificance");
  hDummySignificance->GetYaxis()->SetTitle("S / #sqrt{S+B}");
  hDummySignificance->GetYaxis()->SetTitleOffset(1.8);
  hDummySignificance->GetYaxis()->SetRangeUser(0, 1.2 * fHistSignificance->GetBinContent(1));
  hDummySignificance->Draw("");
  StyleHistoYield(fHistSignificance, 0, 1.2 * fHistSignificance->GetBinContent(1), ColorPart[part], MarkerPart[part], TitleXPt, "S / #sqrt{S+B}", "", MarkerPartSize[part], 1.15, 1.6);
  fHistSignificance->Draw("same");
  // canvasSignificance->SaveAs("../PzsVsPt_Significance.pdf");
  // canvasSignificance->SaveAs("../PzsVsPt_Significance.png");

  TCanvas *canvasPzBkg = new TCanvas("canvasPzBkg", "canvasPzBkg", 900, 700);
  StyleCanvas(canvasPzBkg, 0.05, 0.15, 0.15, 0.05);
  canvasPzBkg->cd();
  TH1F *hDummyPzBkg = (TH1F *)hDummy->Clone("hDummyPzBkg");
  hDummyPzBkg->GetYaxis()->SetRangeUser(-0.02, 0.02);
  hDummyPzBkg->GetXaxis()->SetRangeUser(LowerRangeParticle, UpperRangeParticle);
  hDummyPzBkg->GetYaxis()->SetTitle("P_{z,s2} (Bkg)");
  hDummyPzBkg->GetYaxis()->SetTitleOffset(1.5);
  hDummyPzBkg->Draw("");
  StyleHistoYield(fHistPzsBkg, 0, 0.1, ColorPart[part], MarkerPart[part], TitleXPt, "B / N_{events}", "", MarkerPartSize[part], 1.15, 1.6);
  fHistPzsBkg->Draw("same");
  TF1 *lineAtZero = new TF1("lineAtZero", "0", LowerRangeParticle, UpperRangeParticle);
  lineAtZero->SetLineColor(kBlack);
  lineAtZero->SetLineStyle(2);
  lineAtZero->Draw("same");
  canvasPzBkg->SaveAs("../PzsVsPt_PzBkg.pdf");
  canvasPzBkg->SaveAs("../PzsVsPt_PzBkg.png");

  TCanvas *canvasYield = new TCanvas("canvasYield", "canvasYield", 900, 700);
  StyleCanvas(canvasYield, 0.05, 0.15, 0.15, 0.05);
  canvasYield->cd();
  TH1F *hDummyYield = (TH1F *)hDummy->Clone("hDummyYield");
  hDummyYield->GetYaxis()->SetRangeUser(0, 1);
  if (part == 1)
    hDummyYield->GetYaxis()->SetRangeUser(0, 0.06);
  if (ChosenPart >= 6)
    hDummyYield->GetYaxis()->SetTitle("N_{#Lambda + #bar{#Lambda}} / N_{events}");
  else if (part == 0)
    hDummyYield->GetYaxis()->SetTitle("N_{#Xi^{#minus} + #bar{#Xi}^{+}} / N_{events}");
  else if (part == 1)
    hDummyYield->GetYaxis()->SetTitle("N_{#Omega^{#minus} + #bar{#Omega}^{+}} / N_{events}");
  hDummyYield->GetYaxis()->SetTitleOffset(1.3);
  hDummyYield->Draw("");
  StyleHistoYield(fHistYield, 0, 0.015, ColorPart[part], MarkerPart[part], TitleXPt, "Yield", "", MarkerPartSize[part], 1.15, 1.6);
  fHistYield->Draw("same");
  canvasYield->SaveAs(stringoutpdf + "_Yield.pdf");
  canvasYield->SaveAs(stringoutpdf + "_Yield.png");

  TCanvas *canvasB = new TCanvas("canvasB", "canvasB", 900, 700);
  StyleCanvas(canvasB, 0.05, 0.15, 0.15, 0.05);
  canvasB->cd();
  TH1F *hDummyB = (TH1F *)hDummy->Clone("hDummyB");
  hDummyB->GetYaxis()->SetRangeUser(0, 0.1);
  if (part == 1)
    hDummyB->GetYaxis()->SetRangeUser(0, 0.01);
  hDummyB->GetYaxis()->SetTitle("B / N_{events}");
  hDummyB->GetYaxis()->SetTitleOffset(1.5);
  hDummyB->Draw("");
  StyleHistoYield(fHistB, 0, 0.1, ColorPart[part], MarkerPart[part], TitleXPt, "B / N_{events}", "", MarkerPartSize[part], 1.15, 1.6);
  fHistB->Draw("same");
  canvasB->SaveAs(stringoutpdf + "_B.pdf");
  canvasB->SaveAs(stringoutpdf + "_B.png");

  TCanvas *canvasTot = new TCanvas("canvasTot", "canvasTot", 900, 700);
  StyleCanvas(canvasTot, 0.05, 0.15, 0.15, 0.05);
  canvasTot->cd();
  TH1F *hDummyTot = (TH1F *)hDummy->Clone("hDummyTot");
  hDummyTot->GetYaxis()->SetRangeUser(0, 1);
  hDummyTot->GetYaxis()->SetTitle("(S+B) / N_{events}");
  hDummyTot->GetYaxis()->SetTitleOffset(1.5);
  hDummyTot->Draw("");
  StyleHistoYield(fHistTot, 0, 1, ColorPart[part], MarkerPart[part], TitleXPt, "(S+B) / N_{events}", "", MarkerPartSize[part], 1.15, 1.6);
  fHistTot->Draw("same");
  canvasTot->SaveAs(stringoutpdf + "_Tot.pdf");
  canvasTot->SaveAs(stringoutpdf + "_Tot.png");

  TCanvas *canvasMeanSigma = new TCanvas("canvasMeanSigma", "canvasMeanSigma", 900, 700);
  StyleCanvas(canvasMeanSigma, 0.05, 0.15, 0.15, 0.05);
  canvasMeanSigma->cd();
  TH1F *hDummySigma = (TH1F *)hDummy->Clone("hDummySigma");
  if (ChosenPart >= 6)
    hDummySigma->GetYaxis()->SetRangeUser(1.1, 1.13);
  else if (part == 0)
  {
    hDummySigma->GetYaxis()->SetRangeUser(1.31, 1.33);
    if (ExtrisFitDSCB)
      hDummySigma->GetYaxis()->SetRangeUser(ExtrLowLimitDSCB[ChosenPart] - 0.004, ExtrUpLimitDSCB[ChosenPart] + 0.004);
  }
  else if (part == 1)
    hDummySigma->GetYaxis()->SetRangeUser(1.66, 1.685);
  hDummySigma->GetYaxis()->SetTitle("#mu");
  hDummySigma->GetYaxis()->SetTitleOffset(1.6);
  hDummySigma->Draw("");
  StyleHistoYield(fHistMean, 1.31, 1.33, ColorPart[part], MarkerPart[part], TitleXPt, "#mu", "", MarkerPartSize[part], 1.15, 2);
  fHistMean->Draw("same");
  TH1F *fHistMeanMinus2Sigma = (TH1F *)fHistMean->Clone("fHistMeanMinus2Sigma");
  TH1F *fHistMeanPlus2Sigma = (TH1F *)fHistMean->Clone("fHistMeanPlus2Sigma");
  for (Int_t m = 0; m < fHistMean->GetNbinsX(); m++)
  {
    fHistMeanMinus2Sigma->SetBinContent(m + 1, fHistMean->GetBinContent(m + 1) - 2 * fHistSigma->GetBinContent(m + 1));
    fHistMeanMinus2Sigma->SetBinError(m + 1, 0);
    fHistMeanPlus2Sigma->SetBinContent(m + 1, fHistMean->GetBinContent(m + 1) + 2 * fHistSigma->GetBinContent(m + 1));
    fHistMeanPlus2Sigma->SetBinError(m + 1, 0);
  }
  fHistMeanPlus2Sigma->SetLineColor(kBlack);
  fHistMeanPlus2Sigma->Draw("same");
  fHistMeanMinus2Sigma->SetLineColor(kBlack);
  fHistMeanMinus2Sigma->Draw("same");
  canvasMeanSigma->Modified();
  canvasMeanSigma->Update();
  canvasMeanSigma->SaveAs(stringoutpdf + "_Mean.pdf");
  canvasMeanSigma->SaveAs(stringoutpdf + "_Mean.png");

  TCanvas *canvasfitPol0 = new TCanvas("canvasfitPol0", "canvasfitPol0", 900, 700);
  StyleCanvas(canvasfitPol0, 0.05, 0.15, 0.15, 0.05);
  TH1F *fHistPzsTotError = (TH1F *)fHistPzs->Clone("fHistPzsTotError");
  for (Int_t b = 1; b <= fHistPzs->GetNbinsX(); b++)
  {
    fHistPzsTotError->SetBinError(b, TMath::Sqrt(fHistPzs->GetBinError(b) * fHistPzs->GetBinError(b) +
                                                 fHistPzsSistError->GetBinContent(b) * fHistPzsSistError->GetBinContent(b)));
  }
  canvasfitPol0->cd();
  hDummy->GetYaxis()->SetMaxDigits(1);
  hDummy->Draw("");
  fHistPzsTotError->Draw("same");
  TF1 *fpol1;
  TF1 *fpol0;
  if (ChosenPart >= 6)
  {
    fpol1 = new TF1("fpol1", "pol1", 0, 6);
    fpol0 = new TF1("fpol0", "pol0", 0, 6);
  }
  else
  {
    fpol1 = new TF1("fpol1", "pol1", 0, 80);
    fpol0 = new TF1("fpol0", "pol0", 0, 80);
  }
  fpol1->SetLineColor(kBlue + 2);
  fpol0->SetLineColor(kAzure + 1);

  TH1F *fHistPzsFit = (TH1F *)fHistPzs->Clone("fHistPzsFit");
  // TH1F * fHistPzsFit = (TH1F*)fHistPzsTotError->Clone("fHistPzsFit");
  fHistPzsFit->GetXaxis()->SetRangeUser(PtBins[0], PtBins[numPtBins - 1]);
  // fHistPzsTotError->Fit("fpol1", "R+");
  // fHistPzsTotError->Fit("fpol0", "R+");
  fHistPzsFit->Fit("fpol1", "R+");
  fHistPzsFit->Fit("fpol0", "R+");
  // LegendTitle->Draw("");
  TLegend *legendMainFit = new TLegend(0.2, 0.73, 0.5, 0.88);
  legendMainFit->SetFillStyle(0);
  legendMainFit->SetTextSize(0.05);
  if (ChosenPart >= 6)
    // legendMainFit->AddEntry(fHistPzsTotError, "stat. + syst. " + titleLambda + " Run 3", "p");
    legendMainFit->AddEntry(fHistPzsTotError, "stat. " + titleLambda + " Run 3", "p");
  else if (part == 0)
    legendMainFit->AddEntry(fHistPzsTotError, "stat. + syst. #Xi^{#minus} + #bar{#Xi}^{+} Run 3", "p");
  else if (part == 1)
    legendMainFit->AddEntry(fHistPzsTotError, "stat. + syst. #Omega^{#minus} + #bar{#Omega}^{+} Run 3", "p");
  legendMainFit->Draw("");
  TLegend *legendfit = new TLegend(0.2, 0.58, 0.5, 0.73);
  legendfit->SetFillStyle(0);
  legendfit->SetTextSize(0.04);
  if (ChosenPart >= 6)
  {
    legendfit->AddEntry(fpol0, Form("pol0, Chi2/NDF = %.2f/%i, p0 = %.4f +- %.4f", fpol0->GetChisquare(), fpol0->GetNDF(), fpol0->GetParameter(0), fpol0->GetParError(0)), "l");
    legendfit->AddEntry(fpol1, Form("pol1, Chi2/NDF = %.2f/%i, m = %.5f +- %.5f", fpol1->GetChisquare(), fpol1->GetNDF(), fpol1->GetParameter(1), fpol1->GetParError(1)), "l");
  }
  else
  {
    legendfit->AddEntry(fpol0, Form("pol0, Chi2/NDF = %.2f/%i", fpol0->GetChisquare(), fpol0->GetNDF()), "l");
    legendfit->AddEntry(fpol1, Form("pol1, Chi2/NDF = %.2f/%i", fpol1->GetChisquare(), fpol1->GetNDF()), "l");
  }
  legendfit->Draw("");
  canvasfitPol0->SaveAs("../PzsVsPt_fitPol0.pdf");
  canvasfitPol0->SaveAs("../PzsVsPt_fitPol0.png");

  TF1 *lineatZero = new TF1("lineatZero", "0", 0, 10);
  lineatZero->SetLineColor(kBlack);
  lineatZero->SetLineStyle(2);
  TCanvas *canvasPzsXi = new TCanvas("canvasPzsXi", "canvasPzsXi", 900, 700);
  StyleCanvas(canvasPzsXi, 0.06, 0.12, 0.1, 0.03); // first 0.03
  canvasPzsXi->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleXPt, TitleYPzs, "", 1, 1.15, 1.8);
  StyleHistoYield(fHistPzs, YLow[part], YUp[part], kOrange + 10, 20, TitleXPt, TitleYPzs, "", 1.9, 1.15, 1.8);
  StyleHistoYield(fHistPzsSist, YLow[part], YUp[part], kOrange + 10, 20, TitleXPt, TitleYPzs, "", 1.9, 1.15, 1.8);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(LowerRangeParticle, UpperRangeParticle);
  hDummy->GetYaxis()->SetRangeUser(-0.0004, 0.0065);
  if (part == 1) // Omega
    hDummy->GetYaxis()->SetRangeUser(-0.005, 0.02);
  else if (part == 0) // Xi
    hDummy->GetYaxis()->SetRangeUser(-0.001, 0.011);
  hDummy->Draw("");
  lineatZero->Draw("same");
  if (ChosenPart >= 6)
  {
    fHistPzs->SetLineColor(ColorOO);
    fHistPzs->SetMarkerColor(ColorOO);
    fHistPzsSist->SetLineColor(ColorOO);
    fHistPzsSist->SetMarkerColor(ColorOO);
  }
  fHistPzs->DrawClone("same ex0");
  fHistPzsSist->SetFillStyle(0);
  fHistPzsSist->DrawClone("same e2");
  LegendPreliminary2->Draw("");
  legendXi->Draw("");
  TLegend *legendData = new TLegend(0.06, 0.536, 0.42, 0.736);
  legendData->SetFillStyle(0);
  legendData->SetTextAlign(12);
  legendData->SetTextSize(0.035);
  // legendData->AddEntry(fHistPzs, "Data", "pef");
  legendData->AddEntry("", "Uncertainties: stat. (bar), total sys. (open box)", "");
  legendData->AddEntry(fHistPzs, "Stat. error", "pe");
  legendData->AddEntry(fHistPzsSist, "Syst. error", "f");
  // legendData->Draw("");
  canvasPzsXi->SaveAs("../" + ParticleName[ChosenPart] + "PolVsPt.pdf");
  canvasPzsXi->SaveAs("../" + ParticleName[ChosenPart] + "PolVsPt.png");
  canvasPzsXi->SaveAs("../" + ParticleName[ChosenPart] + "PolVsPt.eps");

  cout << "\n\nSignificance in the 0-50% class" << endl;
  Float_t Pzs0To50 = 0;
  Float_t ErrPzs0To50 = 0;
  Float_t ErrPzs0To50Stat = 0;
  Float_t ErrPzs0To50Sist = 0;
  for (Int_t b = 1; b <= 5; b++)
  {
    Pzs0To50 += fHistPzs->GetBinContent(b) / (pow(fHistPzs->GetBinError(b), 2) + pow(fHistPzsSist->GetBinError(b), 2));
    ErrPzs0To50 += 1. / (pow(fHistPzs->GetBinError(b), 2) + pow(fHistPzsSist->GetBinError(b), 2));
    ErrPzs0To50Stat += 1. / (pow(fHistPzs->GetBinError(b), 2));
    ErrPzs0To50Sist += 1. / (pow(fHistPzsSist->GetBinError(b), 2));
  }
  Pzs0To50 = Pzs0To50 / ErrPzs0To50;
  ErrPzs0To50 = 1. / sqrt(ErrPzs0To50);
  ErrPzs0To50Stat = 1. / sqrt(ErrPzs0To50Stat);
  ErrPzs0To50Sist = 1. / sqrt(ErrPzs0To50Sist);
  cout << "Pzs0To50: " << Pzs0To50 << " +/- " << ErrPzs0To50 << endl;
  cout << "Significance: " << Pzs0To50 / ErrPzs0To50 << endl;
  TH1F *fHistPzs0To50 = new TH1F("fHistPzs0To50", "fHistPzs0To50", 1, 0, 50);
  TH1F *fHistPzs0To50Sist = new TH1F("fHistPzs0To50Sist", "fHistPzs0To50Sist", 1, 0, 50);
  fHistPzs0To50->SetBinContent(1, Pzs0To50);
  fHistPzs0To50Sist->SetBinContent(1, Pzs0To50);
  fHistPzs0To50->SetBinError(1, ErrPzs0To50Stat);
  fHistPzs0To50Sist->SetBinError(1, ErrPzs0To50Sist);
  fHistPzs0To50->SetLineColor(kPink + 6);
  fHistPzs0To50->SetMarkerColor(kPink + 6);
  fHistPzs0To50->SetMarkerStyle(20);
  fHistPzs0To50->SetMarkerSize(1.5);
  fHistPzs0To50Sist->SetLineColor(kPink + 6);
  fHistPzs0To50Sist->SetMarkerColor(kPink + 6);
  fHistPzs0To50Sist->SetMarkerStyle(20);
  fHistPzs0To50Sist->SetMarkerSize(1.5);

  TLegend *legendPalermo = new TLegend(0.14, 0.51, 0.5, 0.65);
  legendPalermo->SetFillStyle(0);
  legendPalermo->SetTextAlign(12);
  legendPalermo->SetTextSize(0.048);

  TLegend *legendParticles = new TLegend(0.14, 0.7, 0.51, 0.83);
  legendParticles->SetFillStyle(0);
  legendParticles->SetTextAlign(12);
  legendParticles->SetTextSize(0.042);
  legendParticles->SetMargin(0.18);
  if (ChosenPart >= 6) // Lambda
    legendParticles->AddEntry(fHistPzs, Form("%s, |#it{#eta} | < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}, OO #sqrt{#it{s}_{NN}} = 5.36 TeV", titleLambda.Data(), MinPt[ChosenPart]), "pef");
  else if (part == 0) // Xi
    legendParticles->AddEntry(fHistPzs, Form("#Xi^{#minus} + #bar{#Xi}^{+}, |#it{#eta} | < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "pl");
  else if (part == 1) // Omega
    legendParticles->AddEntry(fHistPzs, Form("#Omega^{#minus} + #bar{#Omega}^{+}, |#it{#eta} | < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "pl");
  // if (ChosenPart >= 6)
  //   legendParticles->AddEntry(fHistPzsLambdaJunlee, Form("#Lambda + #bar{#Lambda}, |#it{y} | < 0.5, #it{p}_{T} > %1.1f GeV/#it{c}, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV", 0.5), "pef");
  // else
  //   legendParticles->AddEntry(fHistPzsLambdaJunlee, Form("#Lambda + #bar{#Lambda}, |#it{y} | < 0.5, #it{p}_{T} > %1.1f GeV/#it{c}", 0.5), "pl");
  TCanvas *canvasPzsXiLambda = new TCanvas("canvasPzsXiLambda", "canvasPzsXiLambda", 900, 700);
  StyleCanvas(canvasPzsXiLambda, 0.06, 0.12, 0.1, 0.03);
  canvasPzsXiLambda->cd();
  hDummy->Draw("");
  lineatZero->Draw("same");
  TH1F *fHistPzsUpTo50 = (TH1F *)fHistPzs->Clone("fHistPzsUpTo50");
  TH1F *fHistPzsSistUpTo50 = (TH1F *)fHistPzsSist->Clone("fHistPzsSistUpTo50");
  for (Int_t b = 1; b <= fHistPzsUpTo50->GetNbinsX(); b++)
  {
    if (fHistPzsUpTo50->GetBinCenter(b) > UpperRangeParticle)
    {
      fHistPzsUpTo50->SetBinContent(b, -1000);
      fHistPzsUpTo50->SetBinError(b, 0);
      fHistPzsSistUpTo50->SetBinContent(b, -1000);
      fHistPzsSistUpTo50->SetBinError(b, 0);
    }
  }
  fHistPzsUpTo50->Draw("same ex0");
  fHistPzsSistUpTo50->SetFillStyle(0);
  fHistPzsSistUpTo50->Draw("same e2");
  legendPalermo->AddEntry(gPzsPalermo, "#Lambda + #bar{#Lambda}, Pb-Pb 5.02 TeV, #zeta/s par III", "l");
  legendPalermo->AddEntry("", "Eur. Phys. J.C 84 (2024) 9, 920", "");
  // gPzsPalermo->Draw("same l");
  LegendPreliminary3->Draw("");
  legendParticles->Draw("");
  // legendPalermo->Draw("");
  canvasPzsXiLambda->SaveAs("../XiLambdaPolVsPt.pdf");
  canvasPzsXiLambda->SaveAs("../XiLambdaPolVsPt.png");
  canvasPzsXiLambda->SaveAs("../XiLambdaPolVsPt.eps");

  TFile *fileout = new TFile(stringout, "RECREATE");
  fHistMean->Write();
  fHistSigma->Write();
  fHistPzs->Write();
  fHistPzsSist->Write();
  fHistPzsBkg->Write();
  fHistPurity->Write();
  fHistYield->Write();
  fHistSignificance->Write();
  fHistPzsError->Write();

  fileout->Close();

  cout << "Nsigma between PbPb result in first pT bin and pt-inetrgated result " << endl;
  cout << TMath::Abs(ratioPbPbToOO - hRatioPbPbToOO->GetBinContent(1)) / TMath::Sqrt(pow(hRatioPbPbToOO->GetBinError(1), 2) + pow(hRatioPbPbToOOSist->GetBinError(1), 2)) << endl;
  cout << "\nStarting from the files: " << PathIn << endl;
  cout << "and the file: " << PathInSyst << " for syst. uncertainties,\n";
  cout << "\nI have created the file:\n " << stringout << endl;
}
