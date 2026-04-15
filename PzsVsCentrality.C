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
// #include "CommonVarLambda.h"
#include "CommonVarOmega.h"
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

Int_t colorJunlee = kAzure - 3;
Int_t ColorOO = kMagenta + 1;

void PzsVsCentrality(Int_t ChosenPart = ChosenParticle,
                     Bool_t isPolFromLambda = 0,
                     Bool_t isFromFit = 1,
                     Bool_t isFDCorrected = 0,
                     Bool_t isBkgPol = 1,
                     Bool_t isTighterPzFitRange = 0,
                     Bool_t isRapiditySel = ExtrisRapiditySel,
                     Int_t BkgType = ExtrBkgType,
                     Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  if (ChosenPart < 6 && isFDCorrected)
  {
    cout << "FD correction is only available for Lambda. Please change the settings." << endl;
    return;
  }
  if (isReducedPtBins && numPtBins != numPtBinsReduced)
  {
    cout << "Reduced pt bins are selected, but numPtBins is not set to numPtBinsReduced. Please check the settings." << endl;
    return;
  }
  Int_t ChosenPt = -999;
  cout << "Type 100 if you want to analyse Pz (integrated in pT), or the number of the pt interval you want" << endl;
  cin >> ChosenPt;

  Int_t part = 0;
  if (ChosenPart == 1 || ChosenPart == 4 || ChosenPart == 5)
  {
    part = 1; // Omega
  }
  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  Int_t UpperRangeParticle = 80;
  if (ChosenPart >= 6)
    UpperRangeParticle = 50;

  if (ChosenPart >= 6)
  { // for Lambda in OO
    YLow[part] = {-0.0005};
    // YLow[part] = {-0.0004};
    //  YUp[part] = {0.035};
    YUp[part] = {0.0075};
  }
  if (part == 1) // for Omega
  {
    YLow[part] = {-0.02};
    YUp[part] = {0.025};
  }

  // filein
  TString PathIn;
  TFile *fileIn[commonNumCent + 1];

  // fileinLambda
  TString PathInLambda = "../Run2Results/HEPData-ins1891389-v1-P_z_vsCent.root";
  TFile *fileInLambda = TFile::Open(PathInLambda);
  if (!fileInLambda)
  {
    cout << "No file found" << endl;
    return;
  }
  TDirectoryFile *dirLambda = (TDirectoryFile *)fileInLambda->Get("P_z vs. centrality (5.02 TeV)");
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

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = "../Pzs2VsCentrality/" + NameAnalysis[!isV2] + "_";
  stringout += SinputFileName;
  stringout += "_" + ParticleName[ChosenPart];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[BkgType];
  stringout += "_Pzs2";
  if (isApplyWeights)
    stringout += "_Weighted";
  if (isApplyCentWeight)
    stringout += "_CentWeighted";
  if (!useCommonBDTValue)
    stringout += "_BDTCentDep";
  if (isRun2Binning)
    stringout += "_Run2Binning";
  if (isPolFromLambda)
    stringout += "_PolFromLambda";
  if (ChosenPt == 100)
    stringout += "_PtInt";
  else
    stringout += Form("_Pt%.1f-%.1f", PtBins[ChosenPt], PtBins[ChosenPt + 1]);
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
  if (ExtrisApplyResoOnTheFly)
    stringout += "_ResoOnTheFly";
  if (ChosenPart == 0)
    stringout += "_EPReso";
  if (!isFromFit)
    stringout += "_NoPurityDivision";
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
  // stringout += "_TestMoreBins";
  stringoutpdf = stringout;
  stringout += ".root";

  // FD corrected histogram
  TH1F *hFDFractionAllLambdavsCent = nullptr;
  if (isFDCorrected)
  {
    TString PathInFD = "../LambdaFDFraction" + SinputFileNameFDFraction + ".root";
    TFile *fileInFD = TFile::Open(PathInFD);
    if (!fileInFD)
    {
      cout << "No FD correction file found" << endl;
      return;
    }
    hFDFractionAllLambdavsCent = (TH1F *)fileInFD->Get("hFDFractionAllLambdavsCent");
    if (!hFDFractionAllLambdavsCent)
    {
      cout << "No FD correction hist found" << endl;
      return;
    }
  }

  // canvases
  gStyle->SetOptStat(0);
  TCanvas *canvasPzs = new TCanvas("canvasPzs", "canvasPzs", 900, 700);
  StyleCanvas(canvasPzs, 0.05, 0.15, 0.15, 0.05);

  int nBinsCent;
  double *binsCent;
  if (isOOCentrality)
  {
    nBinsCent = numCentLambdaOO;
    binsCent = fCentFT0CLambdaOO;
  }
  else if (part == 1)
  { // Omega in Pb-Pb
    nBinsCent = numCentOmega;
    binsCent = fCentFT0COmega;
  }
  else
  {
    nBinsCent = numCent;
    binsCent = fCentFT0C;
  }
  TH1F *fHistPzs = new TH1F("fHistPzs", "fHistPzs", nBinsCent, binsCent);
  TH1F *fHistPzsSist = new TH1F("fHistPzsSist", "fHistPzsSist", nBinsCent, binsCent);
  TH1F *fHistPzsError = new TH1F("fHistPzsError", "fHistPzsError", nBinsCent, binsCent);
  TH1F *fHistPzsBkgUnderPeak = new TH1F("fHistPzsBkgUnderPeak", "fHistPzsBkgUnderPeak", nBinsCent, binsCent);
  TH1F *fHistPuritySummary = new TH1F("fHistPuritySummary", "fHistPuritySummary", nBinsCent, binsCent);
  TH1F *fHistSignificanceSummary = new TH1F("fHistSignificanceSummary", "fHistSignificanceSummary", nBinsCent, binsCent);
  TH1F *fHistYieldSummary = new TH1F("fHistYieldSummary", "fHistYieldSummary", nBinsCent, binsCent);
  TH1F *fHistMeanSummary = new TH1F("fHistMeanSummary", "fHistMeanSummary", nBinsCent, binsCent);
  TH1F *fHistSigmaSummary = new TH1F("fHistSigmaSummary", "fHistSigmaSummary", nBinsCent, binsCent);
  TH1F *fHistMeanMinus2Sigma = new TH1F("fHistMeanMinus2Sigma", "fHistMeanMinus2Sigma", nBinsCent, binsCent);
  TH1F *fHistMeanPlus2Sigma = new TH1F("fHistMeanPlus2Sigma", "fHistMeanPlus2Sigma", nBinsCent, binsCent);
  TH1F *fHistBSummary = new TH1F("fHistBSummary", "fHistBSummary", nBinsCent, binsCent);
  TH1F *fHistTotSummary = new TH1F("fHistTotSummary", "fHistTotSummary", nBinsCent, binsCent);

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *legendLambda = new TLegend(0.5, 0.53, 0.9, 0.73);
  legendLambda->SetFillStyle(0);
  legendLambda->SetTextSize(0.03);
  legendLambda->AddEntry(fHistPzsLambda, "#Lambda + #bar{#Lambda}, Phys. Rev. Lett. 128.17 (2022)", "pl");

  TLegend *LegendTitle;
  LegendTitle = new TLegend(0.54, 0.72, 0.93, 0.9);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(33);
  LegendTitle->SetTextSize(0.04);
  LegendTitle->AddEntry("", "#bf{ALICE Preliminary}", "");
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
  if (ChosenPt == 100)
    LegendTitle->AddEntry("", Form("#it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "");
  else
    LegendTitle->AddEntry("", Form("%1.1f < #it{p}_{T} < %1.1f GeV/#it{c}", PtBins[ChosenPt], PtBins[ChosenPt + 1]), "");

  TLine *lineat0 = new TLine(0, 0, 100, 0);
  lineat0->SetLineColor(1);
  lineat0->SetLineStyle(2);

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  TH1F *fHistSpectrum[commonNumCent + 1];
  TH1F *fHistPzsBkg[commonNumCent + 1];
  TH1F *fHistPurity[commonNumCent + 1];
  TH1F *fHistSignificance[commonNumCent + 1];
  TH1F *fHistYield[commonNumCent + 1];
  TH1F *fHistMean[commonNumCent + 1];
  TH1F *fHistSigma[commonNumCent + 1];
  TH1F *fHistB[commonNumCent + 1];
  TH1F *fHistTot[commonNumCent + 1];
  TString Smolt[commonNumCent + 1];
  TString SmoltBis[commonNumCent + 1];
  Float_t FDFraction[commonNumCent + 1];
  // get spectra in multiplicity classes
  for (Int_t m = 0; m < commonNumCent; m++)
  {
    if (m == numCent)
    {
      CentFT0CMin = 0;
      CentFT0CMax = CentFT0CMaxPbPb;
    }
    else
    {
      CentFT0CMin = CentFT0C[m];
      CentFT0CMax = CentFT0C[m + 1];
    }
    if (part == 1)
    { // Omega in Pb-Pb
      if (m == numCentOmega)
      {
        CentFT0CMin = 0;
        CentFT0CMax = CentFT0CMaxOmega;
      }
      else
      {
        CentFT0CMin = CentFT0COmega[m];
        CentFT0CMax = CentFT0COmega[m + 1];
      }
    }
    if (isOOCentrality)
    {
      if (m == (numCentLambdaOO))
      {
        CentFT0CMin = 0;
        CentFT0CMax = CentFT0CMaxLambdaOO;
      }
      else
      {
        CentFT0CMin = CentFT0CLambdaOO[m];
        CentFT0CMax = CentFT0CLambdaOO[m + 1];
      }
    }
    PathIn = "../OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_";
    PathIn += SinputFileName;
    PathIn += "_" + ParticleName[ChosenPart];
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[BkgType];
    Smolt[m] += Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);
    SmoltBis[m] += Form("%i#minus%i", CentFT0CMin, CentFT0CMax);
    PathIn += Smolt[m];
    if (isApplyWeights)
      PathIn += "_Weighted";
    if (isApplyCentWeight)
      PathIn += "_CentWeighted";
    if (!useCommonBDTValue)
      PathIn += "_BDTCentDep";
    if (isRun2Binning)
      PathIn += "_Run2Binning";
    if (isPolFromLambda)
      PathIn += "_PolFromLambda";
    if (ExtrisApplyEffWeights)
      PathIn += "_EffW";
    if (!isRapiditySel)
      PathIn += "_Eta08";
    PathIn += STHN[ExtrisFromTHN];
    if (useMixedBDTValueInFitMacro)
      PathIn += "_MixedBDT";
    if (isTightMassCut)
      PathIn += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
    if (isReducedPtBins)
      PathIn += "_ReducedPtBins";
    // if (ChosenPart == 6 && !ExtrisFromTHN && ExtrisSysLambdaMultTrial)
    //{
    //   PathIn += "_SysMultTrial_0_isSysLambdaMultTrial";
    // }
    if (ChosenPart >= 6 && ExtrisSysLambdaMultTrial)
    {
      if (isLoosest)
        PathIn += "_isLoosest";
      else if (isTightest)
        PathIn += "_isTightest";
      PathIn += "_isSysLambdaMultTrial";
    }
    if (ExtrisApplyResoOnTheFly)
      PathIn += "_ResoOnTheFly";
    if (ChosenPart == 0)
      PathIn += "_EPReso";
    if (isBkgPol == 0)
      PathIn += "_isBkgPol0";
    // PathIn += "_SystReso";
    if (isTighterPzFitRange)
      PathIn += "_TighterPzFitRange";
    PathIn += ".root";
    cout << "Path in : " << PathIn << endl;
    fileIn[m] = TFile::Open(PathIn);
    if (ChosenPt == 100)
    {
      fHistSpectrum[m] = (TH1F *)fileIn[m]->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + "PtInt" + V2FromFit[isFromFit]);
      fHistPzsBkg[m] = (TH1F *)fileIn[m]->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + "BkgPtInt" + V2FromFit[isFromFit]);
      fHistPurity[m] = (TH1F *)fileIn[m]->Get("histoPurityPtInt");
      fHistSignificance[m] = (TH1F *)fileIn[m]->Get("histoSignificancePtInt");
      fHistYield[m] = (TH1F *)fileIn[m]->Get("histoYieldPtInt");
      fHistMean[m] = (TH1F *)fileIn[m]->Get("histoMeanPtInt");
      // fHistSigma[m] = (TH1F *)fileIn[m]->Get("histoSigmaPtIntWeighted");
      fHistSigma[m] = (TH1F *)fileIn[m]->Get("histoSigmaPtInt");
      fHistB[m] = (TH1F *)fileIn[m]->Get("histoBPtInt");
      fHistTot[m] = (TH1F *)fileIn[m]->Get("histoTotPtInt");
    }
    else
    {
      fHistSpectrum[m] = (TH1F *)fileIn[m]->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + V2FromFit[isFromFit]);
      fHistPzsBkg[m] = (TH1F *)fileIn[m]->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + "Bkg" + V2FromFit[isFromFit]);
      fHistPurity[m] = (TH1F *)fileIn[m]->Get("histoPurity");
      fHistSignificance[m] = (TH1F *)fileIn[m]->Get("histoSignificance");
      fHistYield[m] = (TH1F *)fileIn[m]->Get("histoYield");
      fHistMean[m] = (TH1F *)fileIn[m]->Get("histoMean");
      // fHistSigma[m] = (TH1F *)fileIn[m]->Get("histoSigmaWeighted");
      fHistSigma[m] = (TH1F *)fileIn[m]->Get("histoSigma");
      fHistB[m] = (TH1F *)fileIn[m]->Get("histoB");
      fHistTot[m] = (TH1F *)fileIn[m]->Get("histoTot");
    }
    if (!fHistSpectrum[m])
    {
      cout << " no hist v2 / Pzs" << endl;
      return;
    }
    if (!fHistPzsBkg[m])
    {
      cout << " no hist v2 / Pzs Bkg" << endl;
      return;
    }
    if (!fHistPurity[m])
    {
      cout << " no hist purity" << endl;
      return;
    }
    if (!fHistSignificance[m])
    {
      cout << " no hist significance " << endl;
      return;
    }
    if (!fHistYield[m])
    {
      cout << " no hist yield " << endl;
      return;
    }
    if (!fHistMean[m])
    {
      cout << " no hist mean" << endl;
      return;
    }
    if (!fHistSigma[m])
    {
      cout << " no hist sigma" << endl;
      return;
    }
    if (!fHistB[m])
    {
      cout << " no hist B" << endl;
      return;
    }
    if (!fHistTot[m])
    {
      cout << " no hist Tot" << endl;
      return;
    }
    fHistSpectrum[m]->SetName("histoPzs2_" + Smolt[m]);
    fHistPzsBkg[m]->SetName("histoPzs2Bkg_" + Smolt[m]);
    fHistPurity[m]->SetName("histoPurity_" + Smolt[m]);
    fHistSignificance[m]->SetName("histoSignificance_" + Smolt[m]);
    fHistYield[m]->SetName("histoYield_" + Smolt[m]);
    fHistMean[m]->SetName("histoMean_" + Smolt[m]);
    fHistSigma[m]->SetName("histoSigma_" + Smolt[m]);
    fHistB[m]->SetName("histoB_" + Smolt[m]);
    fHistTot[m]->SetName("histoTot_" + Smolt[m]);

    if (!isFromFit)
    {
      // fHistPzs->SetBinContent(m + 1, fHistSpectrum[m]->GetBinContent(1) / fHistPurity[m]->GetBinContent(1));
      fHistPzs->SetBinContent(m + 1, fHistSpectrum[m]->GetBinContent(1));
      // fHistPzs->SetBinError(m + 1, fHistSpectrum[m]->GetBinError(1) / fHistPurity[m]->GetBinContent(1));
      fHistPzs->SetBinError(m + 1, fHistSpectrum[m]->GetBinError(1));
      // fHistPzsError->SetBinContent(m + 1, fHistSpectrum[m]->GetBinError(1) / fHistPurity[m]->GetBinContent(1));
      fHistPzsError->SetBinContent(m + 1, fHistSpectrum[m]->GetBinError(1));
      fHistPzsError->SetBinError(m + 1, 0);
      fHistPzsBkgUnderPeak->SetBinContent(m + 1, fHistPzsBkg[m]->GetBinContent(1));
      fHistPzsBkgUnderPeak->SetBinError(m + 1, fHistPzsBkg[m]->GetBinError(1));
    }
    else
    {
      fHistPzs->SetBinContent(m + 1, fHistSpectrum[m]->GetBinContent(1));
      fHistPzs->SetBinError(m + 1, fHistSpectrum[m]->GetBinError(1));
      fHistPzsError->SetBinContent(m + 1, fHistSpectrum[m]->GetBinError(1));
      fHistPzsError->SetBinError(m + 1, 0);
      fHistPzsBkgUnderPeak->SetBinContent(m + 1, fHistPzsBkg[m]->GetBinContent(1));
      fHistPzsBkgUnderPeak->SetBinError(m + 1, fHistPzsBkg[m]->GetBinError(1));
    }
    if (isFDCorrected)
    {
      FDFraction[m] = hFDFractionAllLambdavsCent->GetBinContent(m + 1);
      cout << "FD fraction for cent " << CentFT0CMin << "-" << CentFT0CMax << " : " << FDFraction[m] << endl;
      cout << "Correction = " << 1. / (CXiToLambda * FDFraction[m] + (1 - FDFraction[m])) << endl;
      fHistPzs->SetBinContent(m + 1, fHistSpectrum[m]->GetBinContent(1) / (CXiToLambda * FDFraction[m] + (1 - FDFraction[m])));
      fHistPzs->SetBinError(m + 1, fHistSpectrum[m]->GetBinError(1) / (CXiToLambda * FDFraction[m] + (1 - FDFraction[m])));
      fHistPzsError->SetBinContent(m + 1, fHistSpectrum[m]->GetBinError(1) / (CXiToLambda * FDFraction[m] + (1 - FDFraction[m])));
      fHistPzsError->SetBinError(m + 1, 0);
    }
    cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << " ";
    cout << "Pzs2: " << fHistSpectrum[m]->GetBinContent(1) << " +- " << fHistSpectrum[m]->GetBinError(1) << endl;
    // cout << fHistPzs->GetBinContent(m + 1) << " +- " << fHistPzs->GetBinError(m + 1) << endl;

    fHistPuritySummary->SetBinContent(m + 1, fHistPurity[m]->GetBinContent(1));
    fHistPuritySummary->SetBinError(m + 1, fHistPurity[m]->GetBinError(1));

    fHistSignificanceSummary->SetBinContent(m + 1, fHistSignificance[m]->GetBinContent(1));
    fHistSignificanceSummary->SetBinError(m + 1, fHistSignificance[m]->GetBinError(1));
    // cout << fHistSignificance[m]->GetBinContent(1) << endl;

    if (ChosenPt == 100)
    {
      fHistYieldSummary->SetBinContent(m + 1, fHistYield[m]->GetBinContent(1) * (MaxPt[ChosenPart] - MinPt[ChosenPart]));
      fHistYieldSummary->SetBinError(m + 1, fHistYield[m]->GetBinError(1) * (MaxPt[ChosenPart] - MinPt[ChosenPart]));
      fHistBSummary->SetBinContent(m + 1, fHistB[m]->GetBinContent(1) * (MaxPt[ChosenPart] - MinPt[ChosenPart]));
      fHistBSummary->SetBinError(m + 1, fHistB[m]->GetBinError(1) * (MaxPt[ChosenPart] - MinPt[ChosenPart]));
      fHistTotSummary->SetBinContent(m + 1, fHistTot[m]->GetBinContent(1) * (MaxPt[ChosenPart] - MinPt[ChosenPart]));
      fHistTotSummary->SetBinError(m + 1, fHistTot[m]->GetBinError(1) * (MaxPt[ChosenPart] - MinPt[ChosenPart]));
    }
    else
    {
      fHistYieldSummary->SetBinContent(m + 1, fHistYield[m]->GetBinContent(1) * (PtBins[ChosenPt + 1] - PtBins[ChosenPt]));
      fHistYieldSummary->SetBinError(m + 1, fHistYield[m]->GetBinError(1) * (PtBins[ChosenPt + 1] - PtBins[ChosenPt]));
      fHistBSummary->SetBinContent(m + 1, fHistB[m]->GetBinContent(1) * (PtBins[ChosenPt + 1] - PtBins[ChosenPt]));
      fHistBSummary->SetBinError(m + 1, fHistB[m]->GetBinError(1) * (PtBins[ChosenPt + 1] - PtBins[ChosenPt]));
      fHistTotSummary->SetBinContent(m + 1, fHistTot[m]->GetBinContent(1) * (PtBins[ChosenPt + 1] - PtBins[ChosenPt]));
      fHistTotSummary->SetBinError(m + 1, fHistTot[m]->GetBinError(1) * (PtBins[ChosenPt + 1] - PtBins[ChosenPt]));
    }

    fHistMeanSummary->SetBinContent(m + 1, fHistMean[m]->GetBinContent(1));
    fHistMeanSummary->SetBinError(m + 1, fHistMean[m]->GetBinError(1));

    fHistSigmaSummary->SetBinContent(m + 1, fHistSigma[m]->GetBinContent(1));
    fHistSigmaSummary->SetBinError(m + 1, fHistSigma[m]->GetBinError(1));

    fHistMeanMinus2Sigma->SetBinContent(m + 1, fHistMean[m]->GetBinContent(1) - 2 * fHistSigma[m]->GetBinContent(1));
    fHistMeanMinus2Sigma->SetBinError(m + 1, 0);

    fHistMeanPlus2Sigma->SetBinContent(m + 1, fHistMean[m]->GetBinContent(1) + 2 * fHistSigma[m]->GetBinContent(1));
    fHistMeanPlus2Sigma->SetBinError(m + 1, 0);

    // cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << " Pzs2: " << fHistSpectrum[m]->GetBinContent(1) << endl;
  } // end loop on mult

  // Get histogram with syst uncertainty
  TString PathInSyst = "../Systematics/SystVsCentrality_" + NameAnalysis[!isV2] + "_";
  PathInSyst += SinputFileNameSyst;
  if (ChosenPart == 7 || ChosenPart == 8)
    PathInSyst += "_" + ParticleName[6];
  else
    PathInSyst += "_" + ParticleName[ChosenPart];
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
    PathInSyst += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
  // PathInSyst += V2FromFit[isFromFit];
  if (isReducedPtBins)
    PathInSyst += "_ReducedPtBins";
  if (ExtrisApplyResoOnTheFly || ChosenPart >= 6)
    PathInSyst += "_ResoOnTheFly";
  PathInSyst += ".root";
  // if (ChosenPart == 6)
  if (part == 1) // Omega does not have syst uncertainty yet
    PathInSyst = "../Systematics/SystVsCentrality_Pzs2_LHC23_PbPb_pass4_Train370610_ProtonAcc_Xi_BkgParab_Pzs2_PolFromLambda_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit.root";
  TFile *fileInSyst = TFile::Open(PathInSyst);
  if (!fileInSyst)
  {
    cout << "No file found" << endl;
    return;
  }
  TH1F *fHistPzsSistError = (TH1F *)fileInSyst->Get("fHistTotalErrorVsCent");
  if (!fHistPzsSistError)
  {
    cout << "No hist found" << endl;
    return;
  }
  for (Int_t b = 1; b <= fHistPzs->GetNbinsX(); b++)
  {
    fHistPzsSistError->SetBinContent(b, fHistPzsSistError->GetBinContent(b) / fHistPuritySummary->GetBinContent(b));
    fHistPzsSist->SetBinContent(b, fHistPzs->GetBinContent(b));
    fHistPzsSist->SetBinError(b, fHistPzsSistError->GetBinContent(b));
  }
  TH1F *fHistPzsSignif = (TH1F *)fHistPzs->Clone("fHistPzs");
  TH1F *fHistPzsSignifStat = (TH1F *)fHistPzs->Clone("fHistPzs");
  TH1F *fHistPzsSignifLambda = (TH1F *)fHistPzs->Clone("fHistPzs");
  for (Int_t b = 1; b <= fHistPzs->GetNbinsX(); b++)
  {
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
  LegendPreliminary2->AddEntry("", "#bf{ALICE Preliminary}", "");
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
    LegendPreliminary3->AddEntry("", "#bf{ALICE Preliminary}", "");
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

  TString SfileNeNeJunlee = "../LambdaJunlee/resOut_psi2_NeNe.root";
  TFile *fileNeNeJunlee = new TFile(SfileNeNeJunlee);
  TGraphErrors *gPzsNeNeJunlee = (TGraphErrors *)fileNeNeJunlee->Get("gpolMult_2_0");

  TString SfileLambdaJunlee = "../LambdaJunlee/fout_psi2_mult.root";
  TFile *fileLambdaJunlee = new TFile(SfileLambdaJunlee);
  TGraphErrors *gPzsLambdaJunlee = (TGraphErrors *)fileLambdaJunlee->Get("gMultStat");
  TGraphErrors *gPzsLambdaJunleeSist = (TGraphErrors *)fileLambdaJunlee->Get("gMultSyst");
  TH1F *hDummy10 = new TH1F("hDummy10", "hDummy10", 8, 0, 80);
  Double_t NeNeBinEdges[3] = {0, 40, 100};
  TH1F *hDummyNeNe = new TH1F("hDummyNeNe", "hDummyNeNe", 2, NeNeBinEdges);
  TH1F *fHistPzsLambdaJunlee = (TH1F *)hDummy10->Clone("fHistPzsLambdaJunlee");
  TH1F *fHistPzsLambdaJunleeStatError = (TH1F *)hDummy10->Clone("fHistPzsLambdaJunleeStatError");
  TH1F *fHistPzsLambdaJunleeSist = (TH1F *)hDummy10->Clone("fHistPzsLambdaJunleeSist");
  fHistPzsLambdaJunlee->Reset();
  fHistPzsLambdaJunleeSist->Reset();
  TH1F *fHistPzsLambdaNeNeJunlee = (TH1F *)hDummyNeNe->Clone("fHistPzsLambdaNeNeJunlee");
  TH1F *fHistPzsLambdaNeNeJunleeStatError = (TH1F *)hDummyNeNe->Clone("fHistPzsLambdaNeNeJunleeStatError");
  fHistPzsLambdaNeNeJunlee->Reset();
  fHistPzsLambdaNeNeJunleeStatError->Reset();
  for (Int_t b = 1; b <= gPzsLambdaJunlee->GetN(); b++)
  {
    fHistPzsLambdaJunlee->SetBinContent(b, gPzsLambdaJunlee->GetY()[b - 1]);
    fHistPzsLambdaJunlee->SetBinError(b, gPzsLambdaJunlee->GetEY()[b - 1]);
    fHistPzsLambdaJunleeSist->SetBinContent(b, gPzsLambdaJunleeSist->GetY()[b - 1]);
    fHistPzsLambdaJunleeSist->SetBinError(b, gPzsLambdaJunleeSist->GetEY()[b - 1]);
    fHistPzsLambdaJunleeStatError->SetBinContent(b, gPzsLambdaJunlee->GetEY()[b - 1]);
    fHistPzsLambdaJunleeStatError->SetBinError(b, 0);
  }
  for (Int_t b = 1; b <= gPzsNeNeJunlee->GetN(); b++)
  {
    fHistPzsLambdaNeNeJunlee->SetBinContent(b, gPzsNeNeJunlee->GetY()[b - 1]);
    fHistPzsLambdaNeNeJunlee->SetBinError(b, gPzsNeNeJunlee->GetEY()[b - 1]);
    fHistPzsLambdaNeNeJunleeStatError->SetBinContent(b, gPzsNeNeJunlee->GetEY()[b - 1]);
    fHistPzsLambdaNeNeJunleeStatError->SetBinError(b, 0);
  }
  TFile *fileoutLambdaJunlee = new TFile("../LambdaJunlee/fout_psi2_mult_WHisto.root", "RECREATE");
  fHistPzsLambdaJunlee->Write();
  fHistPzsLambdaJunleeSist->Write();
  fHistPzsLambdaJunleeStatError->Write();
  fHistPzsLambdaNeNeJunlee->Write();
  fHistPzsLambdaNeNeJunleeStatError->Write();
  fileoutLambdaJunlee->Close();

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, -1000);
  canvasPzs->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleXCent, TitleYPzs, "", 1, 1.15, 1.6);
  StyleHistoYield(fHistPzs, YLow[part], YUp[part], ColorPart[part], MarkerPart[part], TitleXCent, TitleYPzs, "", MarkerPartSize[part], 1.15, 1.6);
  StyleHistoYield(fHistPzsSist, YLow[part], YUp[part], ColorPart[part], MarkerPart[part], TitleXCent, TitleYPzs, "", 1.5, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(0, UpperRangeParticle);
  hDummy->Draw("");
  fHistPzs->Draw("same");
  fHistPzsLambda->Draw("same e0x0");
  fHistPzsLambdaSist->SetFillStyle(0);
  fHistPzsLambdaSist->Draw("same e2");
  fHistPzsSist->SetFillStyle(0);
  fHistPzsSist->Draw("same e2");
  LegendTitle->Draw("");
  legendLambda->Draw("");
  canvasPzs->SaveAs(stringoutpdf + ".pdf");
  canvasPzs->SaveAs(stringoutpdf + ".png");
  canvasPzs->SaveAs(stringoutpdf + ".eps");

  // Relative stat. uncertainty
  TCanvas *canvasPzsError = new TCanvas("canvasPzsError", "canvasPzsError", 900, 700);
  StyleCanvas(canvasPzsError, 0.03, 0.15, 0.15, 0.05);
  TH1F *hDummyError = new TH1F("hDummyError", "hDummyError", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummyError->GetNbinsX(); i++)
    hDummyError->SetBinContent(i, 1e-12);
  canvasPzsError->cd();
  SetFont(hDummyError);
  StyleHistoYield(hDummyError, 0, 0.004, 1, 1, TitleXCent, "Absolute uncertainty", "", 1, 1.15, 1.6);
  if (part == 1)
    hDummyError->GetYaxis()->SetRangeUser(0, 0.015);
  StyleHistoYield(fHistPzsError, 0, 0.01, ColorPart[part], MarkerPart[part], TitleXCent, "", "", MarkerPartSize[part], 1.15, 1.6);
  StyleHistoYield(fHistPzsSistError, 0, 0.01, kGray + 2, 22, TitleXCent, "", "", MarkerPartSize[part], 1.15, 1.6);
  SetHistoTextSize(hDummyError, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  hDummyError->GetYaxis()->SetTitleOffset(1.7);
  SetTickLength(hDummyError, tickX, tickY);
  hDummyError->GetXaxis()->SetRangeUser(0, UpperRangeParticle);
  hDummyError->Draw("");
  fHistPzsError->Draw("same");
  fHistPzsSistError->Draw("same e2");
  fHistPzsLambda_StatErr->SetLineColor(kRed);
  fHistPzsLambda_StatErr->SetMarkerColor(kRed);
  fHistPzsLambda_StatErr->Draw("same e0x0");
  // LegendTitle->Draw("");
  // legendLambda->Draw("");
  // LegendPreliminary3->Draw("");

  TH1F *fHistPzsLambdaJunleeSistError = (TH1F *)fHistPzsLambdaJunleeSist->Clone("fHistPzsLambdaJunleeSistError");
  for (Int_t i = 1; i <= fHistPzsLambdaJunleeSistError->GetNbinsX(); i++)
  {
    fHistPzsLambdaJunleeSistError->SetBinContent(i, fHistPzsLambdaJunleeSist->GetBinError(i));
    fHistPzsLambdaJunleeSistError->SetBinError(i, 0);
  }

  TH1F *fHistPzsLambdaJunleeStatErrorA = (TH1F *)fHistPzsLambdaJunlee->Clone("fHistPzsLambdaJunleeStatErrorA");
  for (Int_t i = 1; i <= fHistPzsLambdaJunleeStatErrorA->GetNbinsX(); i++)
  {
    fHistPzsLambdaJunleeStatErrorA->SetBinContent(i, fHistPzsLambdaJunlee->GetBinError(i));
    fHistPzsLambdaJunleeStatErrorA->SetBinError(i, 0);
  }
  fHistPzsLambdaJunleeStatErrorA->SetMarkerStyle(20);
  fHistPzsLambdaJunleeStatErrorA->SetMarkerSize(1.5);
  fHistPzsLambdaJunleeStatErrorA->SetMarkerColor(kBlue);
  fHistPzsLambdaJunleeStatErrorA->SetLineColor(kBlue);
  fHistPzsLambdaJunleeStatErrorA->Draw("same ");
  fHistPzsLambdaJunleeSistError->SetMarkerColor(kGreen + 1);
  fHistPzsLambdaJunleeSistError->SetLineColor(kGreen + 1);
  fHistPzsLambdaJunleeSistError->Draw("same ");

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
  legendError->AddEntry(fHistPzsLambdaJunleeStatErrorA, "stat. #Lambda + #bar{#Lambda} Run 3", "pl");
  legendError->AddEntry(fHistPzsLambdaJunleeSistError, "syst. #Lambda + #bar{#Lambda} Run 3", "pl");
  legendError->AddEntry(fHistPzsLambda_StatErr, "stat. #Lambda + #bar{#Lambda} Run 2", "pl");
  legendError->Draw("");
  canvasPzsError->SaveAs(stringoutpdf + "_Error.pdf");
  canvasPzsError->SaveAs(stringoutpdf + "_Error.png");

  // Significativity
  TCanvas *canvasPzsSignif = new TCanvas("canvasPzsSignif", "canvasPzsSignif", 900, 700);
  StyleCanvas(canvasPzsSignif, 0.05, 0.15, 0.15, 0.05);
  TH1F *hDummySignif = new TH1F("hDummySignif", "hDummySignif", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummySignif->GetNbinsX(); i++)
    hDummySignif->SetBinContent(i, 1e-12);
  canvasPzsSignif->cd();
  SetFont(hDummySignif);
  StyleHistoYield(hDummySignif, 0, 1.2 * fHistPzsSignif->GetBinContent(fHistPzsSignif->GetMaximumBin()), 1, 1, TitleXCent, "S / #sigma_{S}", "", 1, 1.15, 1.6);
  StyleHistoYield(fHistPzsSignif, 0, 1.2 * fHistPzsSignif->GetBinContent(fHistPzsSignif->GetMaximumBin()), ColorPart[part], MarkerPart[part], TitleXCent, "S / #sigma_{S}", "", MarkerPartSize[part], 1.15, 1.6);
  StyleHistoYield(fHistPzsSignifStat, 0, 1.2 * fHistPzsSignif->GetBinContent(fHistPzsSignif->GetMaximumBin()), kGray + 2, 22, TitleXCent, "S / #sigma_{S}", "", MarkerPartSize[part], 1.15, 1.6);
  StyleHistoYield(fHistPzsSignifLambda, 0, 1.2 * fHistPzsSignif->GetBinContent(fHistPzsSignif->GetMaximumBin()), kBlue, 20, TitleXCent, "S / #sigma_{S}", "", MarkerPartSize[part], 1.15, 1.6);
  SetHistoTextSize(hDummySignif, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummySignif, tickX, tickY);
  hDummySignif->GetXaxis()->SetRangeUser(0, UpperRangeParticle);
  hDummySignif->GetYaxis()->SetRangeUser(0, 7);
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

  // fHistPzsSignifLambda->Draw("same e0x0");
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
  canvasPzsSignif->SaveAs(stringoutpdf + "_Signif.pdf");
  canvasPzsSignif->SaveAs(stringoutpdf + "_Signif.png");

  TCanvas *canvasPurity = new TCanvas("canvasPurity", "canvasPurity", 900, 700);
  StyleCanvas(canvasPurity, 0.05, 0.15, 0.15, 0.05);
  canvasPurity->cd();
  TH1F *hDummyPurity = (TH1F *)hDummy->Clone("hDummyPurity");
  hDummyPurity->GetYaxis()->SetRangeUser(0.9, 1);
  if (part == 1)
    hDummyPurity->GetYaxis()->SetRangeUser(0.5, 1);
  hDummyPurity->GetYaxis()->SetTitle("S / (S+B)");
  hDummyPurity->Draw("");
  StyleHistoYield(fHistPuritySummary, 0.9, 1, ColorPart[part], MarkerPart[part], TitleXCent, "S / (S+B)", "", MarkerPartSize[part], 1.15, 1.6);
  hDummyPurity->GetYaxis()->SetTitleOffset(1.3);
  fHistPuritySummary->Draw("same");
  canvasPurity->SaveAs(stringoutpdf + "_Purity.pdf");
  canvasPurity->SaveAs(stringoutpdf + "_Purity.png");

  TCanvas *canvasSignificance = new TCanvas("canvasSignificance", "canvasSignificance", 900, 700);
  StyleCanvas(canvasSignificance, 0.05, 0.15, 0.15, 0.05);
  canvasSignificance->cd();
  TH1F *hDummySignificance = (TH1F *)hDummy->Clone("hDummySignificance");
  hDummySignificance->GetYaxis()->SetTitle("S / #sqrt{S+B}");
  hDummySignificance->GetYaxis()->SetTitleOffset(1.8);
  hDummySignificance->GetYaxis()->SetRangeUser(0, 1.2 * fHistSignificanceSummary->GetBinContent(1));
  hDummySignificance->Draw("");
  StyleHistoYield(fHistSignificanceSummary, 0, 1.2 * fHistSignificanceSummary->GetBinContent(1), ColorPart[part], MarkerPart[part], TitleXCent, "S / #sqrt{S+B}", "", MarkerPartSize[part], 1.15, 1.6);
  fHistSignificanceSummary->Draw("same");
  canvasSignificance->SaveAs(stringoutpdf + "_Significance.pdf");
  canvasSignificance->SaveAs(stringoutpdf + "_Significance.png");

  TCanvas *canvasPzBkg = new TCanvas("canvasPzBkg", "canvasPzBkg", 900, 700);
  StyleCanvas(canvasPzBkg, 0.05, 0.15, 0.15, 0.05);
  canvasPzBkg->cd();
  TH1F *hDummyPzBkg = (TH1F *)hDummy->Clone("hDummyPzBkg");
  hDummyPzBkg->GetYaxis()->SetRangeUser(-0.01, 0.01);
  hDummyPzBkg->GetXaxis()->SetRangeUser(0, UpperRangeParticle);
  hDummyPzBkg->GetYaxis()->SetTitle("P_{z,s2} (Bkg)");
  hDummyPzBkg->GetYaxis()->SetTitleOffset(1.5);
  hDummyPzBkg->Draw("");
  StyleHistoYield(fHistPzsBkgUnderPeak, 0, 0.1, ColorPart[part], MarkerPart[part], TitleXCent, "B / N_{events}", "", MarkerPartSize[part], 1.15, 1.6);
  fHistPzsBkgUnderPeak->Draw("same");
  TF1 *lineAtZero = new TF1("lineAtZero", "0", 0, UpperRangeParticle);
  lineAtZero->SetLineColor(kBlack);
  lineAtZero->SetLineStyle(2);
  lineAtZero->Draw("same");
  canvasPzBkg->SaveAs(stringoutpdf + "_PzBkg.pdf");
  canvasPzBkg->SaveAs(stringoutpdf + "_PzBkg.png");

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
  StyleHistoYield(fHistYieldSummary, 0, 0.015, ColorPart[part], MarkerPart[part], TitleXCent, "Yield", "", MarkerPartSize[part], 1.15, 1.6);
  fHistYieldSummary->Draw("same");
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
  StyleHistoYield(fHistBSummary, 0, 0.1, ColorPart[part], MarkerPart[part], TitleXCent, "B / N_{events}", "", MarkerPartSize[part], 1.15, 1.6);
  fHistBSummary->Draw("same");
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
  StyleHistoYield(fHistTotSummary, 0, 1, ColorPart[part], MarkerPart[part], TitleXCent, "(S+B) / N_{events}", "", MarkerPartSize[part], 1.15, 1.6);
  fHistTotSummary->Draw("same");
  canvasTot->SaveAs(stringoutpdf + "_Tot.pdf");
  canvasTot->SaveAs(stringoutpdf + "_Tot.png");

  TCanvas *canvasMeanSigma = new TCanvas("canvasMeanSigma", "canvasMeanSigma", 900, 700);
  StyleCanvas(canvasMeanSigma, 0.05, 0.15, 0.15, 0.05);
  canvasMeanSigma->cd();
  TH1F *hDummySigma = (TH1F *)hDummy->Clone("hDummySigma");
  if (ChosenPart >= 6)
    hDummySigma->GetYaxis()->SetRangeUser(1.1, 1.13);
  else if (part == 0)
    hDummySigma->GetYaxis()->SetRangeUser(1.31, 1.33);
  else if (part == 1)
    hDummySigma->GetYaxis()->SetRangeUser(1.66, 1.685);
  hDummySigma->GetYaxis()->SetTitle("#mu");
  hDummySigma->GetYaxis()->SetTitleOffset(1.6);
  hDummySigma->Draw("");
  StyleHistoYield(fHistMeanSummary, 1.31, 1.33, ColorPart[part], MarkerPart[part], TitleXCent, "#mu", "", MarkerPartSize[part], 1.15, 2);
  fHistMeanSummary->Draw("same");
  fHistMeanPlus2Sigma->SetLineColor(kBlack);
  fHistMeanPlus2Sigma->Draw("same");
  fHistMeanMinus2Sigma->SetLineColor(kBlack);
  fHistMeanMinus2Sigma->Draw("same");
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
    fpol1 = new TF1("fpol1", "pol1", 0, 90);
    fpol0 = new TF1("fpol0", "pol0", 0, 90);
  }
  else
  {
    fpol1 = new TF1("fpol1", "pol1", 0, 80);
    fpol0 = new TF1("fpol0", "pol0", 0, 80);
  }
  fpol1->SetLineColor(kBlue + 2);
  fpol0->SetLineColor(kAzure + 1);
  fHistPzsTotError->Fit("fpol1", "R+");
  fHistPzsTotError->Fit("fpol0", "R+");
  // LegendTitle->Draw("");
  TLegend *legendMainFit = new TLegend(0.2, 0.73, 0.5, 0.88);
  legendMainFit->SetFillStyle(0);
  legendMainFit->SetTextSize(0.05);
  if (ChosenPart >= 6)
    legendMainFit->AddEntry(fHistPzsTotError, "stat. + syst. " + titleLambda + " Run 3", "p");
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
  canvasfitPol0->SaveAs(stringoutpdf + "_fitPol0.pdf");
  canvasfitPol0->SaveAs(stringoutpdf + "_fitPol0.png");

  TF1 *lineatZero = new TF1("lineatZero", "0", 0, 100);
  lineatZero->SetLineColor(kBlack);
  lineatZero->SetLineStyle(2);
  TCanvas *canvasPzsXi = new TCanvas("canvasPzsXi", "canvasPzsXi", 900, 700);
  StyleCanvas(canvasPzsXi, 0.06, 0.12, 0.1, 0.03); // first 0.03
  canvasPzsXi->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleXCent, TitleYPzs, "", 1, 1.15, 1.8);
  StyleHistoYield(fHistPzs, YLow[part], YUp[part], kOrange + 10, 20, TitleXCent, TitleYPzs, "", 1.9, 1.15, 1.8);
  StyleHistoYield(fHistPzsSist, YLow[part], YUp[part], kOrange + 10, 20, TitleXCent, TitleYPzs, "", 1.9, 1.15, 1.8);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(0, UpperRangeParticle);
  hDummy->GetYaxis()->SetRangeUser(-0.0004, 0.0065);
  if (part == 1) // Omega
    hDummy->GetYaxis()->SetRangeUser(-0.005, 0.02);
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
  // legendData->AddEntry(fHistPzs, "Stat. error", "pe");
  // legendData->AddEntry(fHistPzsSist, "Syst. error", "f");
  legendData->Draw("");
  canvasPzsXi->SaveAs("../" + ParticleName[ChosenPart] + "PolVsCent.pdf");
  canvasPzsXi->SaveAs("../" + ParticleName[ChosenPart] + "PolVsCent.png");
  canvasPzsXi->SaveAs("../" + ParticleName[ChosenPart] + "PolVsCent.eps");

  TGraphErrors *gPzsPalermo = new TGraphErrors(9);
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

  for (Int_t i = 0; i < gPzsPalermo->GetN(); i++)
  {
    // cout << "Palermo: i " << i << " " << dNdEtaPalermo[gPzsPalermo->GetN() - i - 1] << " " << gPzsPalermo->GetY()[gPzsPalermo->GetN() - i - 1] << endl;
    gPzsPalermo->SetPoint(i, CentPalermo[i], Pzs2Palermo[i]);
    gPzsPalermo->SetPointError(i, 0, 0);
  }
  TLegend *legendPalermo = new TLegend(0.14, 0.51, 0.5, 0.65);
  legendPalermo->SetFillStyle(0);
  legendPalermo->SetTextAlign(12);
  legendPalermo->SetTextSize(0.048);
  StyleHistoYield(fHistPzsLambdaJunlee, YLow[part], YUp[part], colorJunlee, 47, TitleXCent, TitleYPzs, "", 2.1, 1.15, 1.8);
  StyleHistoYield(fHistPzsLambdaJunleeSist, YLow[part], YUp[part], colorJunlee, 47, TitleXCent, TitleYPzs, "", 2.1, 1.15, 1.8);

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
  if (ChosenPart >= 6)
    legendParticles->AddEntry(fHistPzsLambdaJunlee, Form("#Lambda + #bar{#Lambda}, |#it{y} | < 0.5, #it{p}_{T} > %1.1f GeV/#it{c}, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV", 0.5), "pef");
  else
    legendParticles->AddEntry(fHistPzsLambdaJunlee, Form("#Lambda + #bar{#Lambda}, |#it{y} | < 0.5, #it{p}_{T} > %1.1f GeV/#it{c}", 0.5), "pl");
  TCanvas *canvasPzsXiLambda = new TCanvas("canvasPzsXiLambda", "canvasPzsXiLambda", 900, 700);
  StyleCanvas(canvasPzsXiLambda, 0.06, 0.12, 0.1, 0.03);
  canvasPzsXiLambda->cd();
  hDummy->Draw("");
  lineatZero->Draw("same");
  fHistPzsLambdaJunlee->Draw("same ex0");
  fHistPzsLambdaJunleeSist->SetFillStyle(0);
  fHistPzsLambdaJunleeSist->Draw("same e2");
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
  fHistPzsLambdaNeNeJunlee->SetLineColor(kGreen + 2);
  fHistPzsLambdaNeNeJunlee->SetMarkerColor(kGreen + 2);
  fHistPzsLambdaNeNeJunlee->SetMarkerStyle(29);
  fHistPzsLambdaNeNeJunlee->SetMarkerSize(2.1);
  gPzsPalermo->SetLineColor(kBlue + 1);
  gPzsPalermo->SetMarkerColor(kBlue + 1);
  gPzsPalermo->SetLineWidth(3);
  legendPalermo->AddEntry(gPzsPalermo, "#Lambda + #bar{#Lambda}, Pb-Pb 5.02 TeV, #zeta/s par III", "l");
  legendPalermo->AddEntry("", "Eur. Phys. J.C 84 (2024) 9, 920", "");
  if (ChosenPart < 6)
    gPzsPalermo->Draw("same l");
  // fHistPzsLambdaNeNeJunlee->Draw("same ex0");
  // gPzsLambdaJunlee->Draw("same p");
  // gPzsLambdaJunleeSist->Draw("same e2");
  // fHistPzs0To50->Draw("same ex0");
  // fHistPzs0To50Sist->SetFillStyle(0);
  // fHistPzs0To50Sist->Draw("same e2");
  LegendPreliminary3->Draw("");
  legendParticles->Draw("");
  if (ChosenPart < 6)
    legendPalermo->Draw("");
  canvasPzsXiLambda->SaveAs("../XiLambdaPolVsCent.pdf");
  canvasPzsXiLambda->SaveAs("../XiLambdaPolVsCent.png");
  canvasPzsXiLambda->SaveAs("../XiLambdaPolVsCent.eps");

  Float_t xLabelMult = 35;
  Float_t xTitleMult = 35;
  Float_t xLabelOffsetMult = 0.003;

  if (ChosenPart >= 6)
  {
    TCanvas *canvasPzsVsMultiplicity = new TCanvas("canvasPzsVsMultiplicity", "canvasPzsVsMultiplicity", 900, 700);
    StyleCanvas(canvasPzsVsMultiplicity, 0.06, 0.15, 0.1, 0.03);
    canvasPzsVsMultiplicity->cd();
    gPad->SetLogx();
    TH1F *hDummyVsMultiplicity = new TH1F("hDummyVsMultiplicity", "hDummyVsMultiplicity", 2000, 0, 10000);
    for (Int_t i = 1; i <= hDummyVsMultiplicity->GetNbinsX(); i++)
      hDummyVsMultiplicity->SetBinContent(i, -1000);
    hDummyVsMultiplicity->GetYaxis()->SetMaxDigits(1);
    canvasPzsVsMultiplicity->cd();
    SetFont(hDummyVsMultiplicity);
    TString titledNdeta = "#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";
    StyleHistoYield(hDummyVsMultiplicity, YLow[part], YUp[part], 1, 1, titledNdeta, TitleYPzs, "", 1, 1.15, 1.8);
    SetHistoTextSize(hDummyVsMultiplicity, xTitleMult, xLabelMult, xOffset, xLabelOffsetMult, yTitle, yLabel, yOffset, yLabelOffset);
    SetTickLength(hDummyVsMultiplicity, tickX, tickY);
    hDummyVsMultiplicity->GetXaxis()->SetRangeUser(25, 2500);
    // hDummyVsMultiplicity->GetYaxis()->SetTitleOffset(1.7);
    hDummyVsMultiplicity->GetYaxis()->SetRangeUser(YLow[part], 0.0085);
    hDummyVsMultiplicity->Draw("");
    TF1 *lineatZeroVsMult = new TF1("lineatZeroVsMult", "0", 25, 2500);
    lineatZeroVsMult->SetLineColor(kBlack);
    lineatZeroVsMult->SetLineStyle(2);
    lineatZeroVsMult->Draw("same");
    TString SfilePzspPb = "../HEPData-ins2879229-v1-Table_1.root";
    TFile *filePzspPb = new TFile(SfilePzspPb);
    TDirectoryFile *dir = (TDirectoryFile *)filePzspPb->Get("Table 1");
    TGraphAsymmErrors *gPzspPb = (TGraphAsymmErrors *)dir->Get("Graph1D_y3");

    gStyle->SetEndErrorSize(0);
    TGraphErrors *gPzsVsMult;
    TGraphErrors *gPzsVsMultJunlee = new TGraphErrors(8);
    TGraphErrors *gPzsVsMultNeNeJunlee = new TGraphErrors(2);
    TGraphErrors *gPzsVsMultSist;
    TGraphErrors *gPzsVsMultSistJunlee = new TGraphErrors(8);
    TGraphErrors *gPzsVsMultSistNeNeJunlee = new TGraphErrors(2); // not available yet

    for (Int_t i = 1; i <= gPzsLambdaJunlee->GetN(); i++)
    {
      // cout << "i " << i << " " << dNdEtaAbhi[gPzsLambdaJunlee->GetN() - i] << " " << gPzsLambdaJunlee->GetY()[gPzsLambdaJunlee->GetN() - i] << endl;
      gPzsVsMultJunlee->SetPoint(i, dNdEtaAbhi[gPzsLambdaJunlee->GetN() - i], gPzsLambdaJunlee->GetY()[gPzsLambdaJunlee->GetN() - i]);
      gPzsVsMultJunlee->SetPointError(i, 0, gPzsLambdaJunlee->GetEY()[gPzsLambdaJunlee->GetN() - i]);
      gPzsVsMultSistJunlee->SetPoint(i, dNdEtaAbhi[gPzsLambdaJunlee->GetN() - i], gPzsLambdaJunleeSist->GetY()[gPzsLambdaJunlee->GetN() - i]);
      gPzsVsMultSistJunlee->SetPointError(i, dNdEtaAbhiErr[gPzsLambdaJunlee->GetN() - i], gPzsLambdaJunleeSist->GetEY()[gPzsLambdaJunlee->GetN() - i]);
    }
    for (Int_t i = 1; i <= gPzsNeNeJunlee->GetN(); i++)
    {
      // cout << "NeNe: i " << i << " " << dNdEtaNeNe[gPzsNeNeJunlee->GetN() - i] << " " << gPzsNeNeJunlee->GetY()[gPzsNeNeJunlee->GetN() - i] << endl;
      gPzsVsMultNeNeJunlee->SetPoint(i, dNdEtaNeNe[gPzsVsMultNeNeJunlee->GetN() - i], gPzsNeNeJunlee->GetY()[gPzsNeNeJunlee->GetN() - i]);
      gPzsVsMultNeNeJunlee->SetPointError(i, dNdEtaNeNeErr[gPzsVsMultNeNeJunlee->GetN() - i], gPzsNeNeJunlee->GetEY()[gPzsNeNeJunlee->GetN() - i]);
    }
    gPzsVsMult = new TGraphErrors(nBinsCent);
    gPzsVsMultSist = new TGraphErrors(nBinsCent);

    for (Int_t b = 1; b <= commonNumCent; b++)
    {
      if (b < 6)
        continue;
      // cout << "b " << b << " " << dNdEtaOO[commonNumCent - b] << endl;
      if (isOOCentrality)
      {
        gPzsVsMult->SetPoint(b - 1, dNdEtaOO[commonNumCent - b], fHistPzs->GetBinContent(commonNumCent - b + 1));
        gPzsVsMult->SetPointError(b - 1, dNdEtaOOErr[commonNumCent - b], fHistPzs->GetBinError(commonNumCent - b + 1));
        gPzsVsMultSist->SetPoint(b - 1, dNdEtaOO[commonNumCent - b], fHistPzsSist->GetBinContent(commonNumCent - b + 1));
        gPzsVsMultSist->SetPointError(b - 1, dNdEtaOOErrSyst[commonNumCent - b], fHistPzsSist->GetBinError(commonNumCent - b + 1));
      }
      else
      {
        gPzsVsMult->SetPoint(b - 1, dNdEtaAbhi[commonNumCent - b], fHistPzs->GetBinContent(commonNumCent - b + 1));
        gPzsVsMult->SetPointError(b - 1, dNdEtaAbhiErr[commonNumCent - b], fHistPzs->GetBinError(commonNumCent - b + 1));
        gPzsVsMultSist->SetPoint(b - 1, dNdEtaAbhi[commonNumCent - b], fHistPzsSist->GetBinContent(commonNumCent - b + 1));
        gPzsVsMultSist->SetPointError(b - 1, dNdEtaAbhiErr[commonNumCent - b], fHistPzsSist->GetBinError(commonNumCent - b + 1));
      }
    }
    for (Int_t i = 0; i < gPzspPb->GetN(); i++)
    {
      gPzspPb->SetPoint(i, gPzspPb->GetX()[i] / 4.8, gPzspPb->GetY()[i] / 100);
      gPzspPb->SetPointError(i, 0, 0, gPzspPb->GetErrorYlow(i) / 100, gPzspPb->GetErrorYhigh(i) / 100);
    }

    cout << "\n\nNsigma between OO and PbPb results at multiplicity of about 100: " << endl;
    cout << "Mult (OO, PbPb): " << gPzsVsMult->GetX()[9] << " " << gPzsVsMultJunlee->GetX()[2] << endl;
    cout << "Pzs,2 (OO, PbPb): " << gPzsVsMult->GetY()[9] << " " << gPzsVsMultJunlee->GetY()[2] << endl;
    Double_t ErrorOO = sqrt(pow(gPzsVsMult->GetErrorYlow(9), 2) + pow(gPzsVsMultSist->GetErrorYlow(9), 2));
    Double_t ErrorPbPb = sqrt(pow(gPzsVsMultJunlee->GetErrorYlow(2), 2) + pow(gPzsVsMultSistJunlee->GetErrorYlow(2), 2));
    cout << "ErrorOO: " << ErrorOO << " ErrorPbPb: " << ErrorPbPb << endl;
    cout << "Nsigma: " << fabs(gPzsVsMult->GetY()[9] - gPzsVsMultJunlee->GetY()[2]) / sqrt(ErrorOO * ErrorOO + ErrorPbPb * ErrorPbPb) << endl;

    cout << "\n\nNsigma between OO and PbPb results at multiplicity of about 50: " << endl;
    cout << "Mult (OO, PbPb): " << gPzsVsMult->GetX()[6] << " " << gPzsVsMultJunlee->GetX()[1] << endl;
    cout << "Pzs,2 (OO, PbPb): " << gPzsVsMult->GetY()[6] << " " << gPzsVsMultJunlee->GetY()[1] << endl;
    ErrorOO = sqrt(pow(gPzsVsMult->GetErrorYlow(6), 2) + pow(gPzsVsMultSist->GetErrorYlow(6), 2));
    ErrorPbPb = sqrt(pow(gPzsVsMultJunlee->GetErrorYlow(1), 2) + pow(gPzsVsMultSistJunlee->GetErrorYlow(1), 2));
    cout << "ErrorOO: " << ErrorOO << " ErrorPbPb: " << ErrorPbPb << endl;
    cout << "Nsigma: " << fabs(gPzsVsMult->GetY()[6] - gPzsVsMultJunlee->GetY()[1]) / sqrt(ErrorOO * ErrorOO + ErrorPbPb * ErrorPbPb) << endl;

    gPzspPb->SetLineColor(kBlack);
    gPzspPb->SetMarkerColor(kBlack);
    gPzspPb->SetMarkerStyle(20);
    // gPzspPb->Draw("same p");
    gPzsVsMult->SetMarkerStyle(20);
    gPzsVsMult->SetMarkerSize(1.45);
    gPzsVsMult->SetMarkerColor(ColorOO);
    gPzsVsMult->SetLineColor(ColorOO);
    gPzsVsMult->Draw("same p");
    gPzsVsMultSist->SetMarkerStyle(20);
    gPzsVsMultSist->SetMarkerSize(1.45);
    gPzsVsMultSist->SetFillStyle(0);
    gPzsVsMultSist->SetMarkerColor(ColorOO);
    gPzsVsMultSist->SetFillColor(ColorOO);
    gPzsVsMultSist->SetLineColor(ColorOO);
    gPzsVsMultSist->Draw("same p2");
    gPzsVsMultJunlee->SetMarkerStyle(47);
    gPzsVsMultJunlee->SetMarkerSize(1.7);
    gPzsVsMultJunlee->SetMarkerColor(colorJunlee);
    gPzsVsMultJunlee->SetLineColor(colorJunlee);
    gPzsVsMultJunlee->Draw("same p");
    gPzsVsMultSistJunlee->SetMarkerStyle(47);
    gPzsVsMultSistJunlee->SetMarkerSize(1.7);
    gPzsVsMultSistJunlee->SetFillStyle(0);
    gPzsVsMultSistJunlee->SetMarkerColor(colorJunlee);
    gPzsVsMultSistJunlee->SetLineColor(colorJunlee);
    gPzsVsMultSistJunlee->Draw("same p2");
    gPzsVsMultNeNeJunlee->SetMarkerStyle(20);
    gPzsVsMultNeNeJunlee->SetMarkerColor(kGreen + 2);
    gPzsVsMultNeNeJunlee->SetLineColor(kGreen + 2);
    // gPzsVsMultNeNeJunlee->Draw("same p");
    TLegend *legendSystem = new TLegend(0.2, 0.7, 0.60, 0.9);
    legendSystem->SetFillStyle(0);
    legendSystem->SetTextAlign(12);
    legendSystem->SetTextSize(0.035);
    if (ChosenPart >= 6) // Lambda
      legendSystem->AddEntry(gPzsVsMult, Form("%s, |#it{#eta} | < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}, OO #sqrt{#it{s}_{NN}} = 5.36 TeV", titleLambda.Data(), MinPt[ChosenPart]), "pl");
    else if (part == 0) // Xi
      legendSystem->AddEntry(gPzsVsMult, Form("#Xi^{#minus} + #bar{#Xi}^{+}, |#it{#eta} | < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}, Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV", MinPt[ChosenPart]), "pl");
    else if (part == 1) // Omega
      legendSystem->AddEntry(gPzsVsMult, Form("#Omega^{#minus} + #bar{#Omega}^{+}, |#it{#eta} | < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}, Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV", MinPt[ChosenPart]), "pl");

    legendSystem->AddEntry(gPzsVsMultJunlee, Form("#Lambda + #bar{#Lambda}, |#it{y} | < 0.5, #it{p}_{T} > %1.1f GeV/#it{c}, Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV", 0.5), "pl");
    LegendPreliminary3->Draw("");
    legendParticles->Draw("");
    // legendSystem->Draw("");
    //  gPzsVsMultJunleeSist->SetFillStyle(0);
    //  gPzsVsMultJunleeSist->Draw("same e2");
    TLegend *legendData2 = new TLegend(0.15, 0.557, 0.50, 0.757);
    legendData2->SetMargin(0.18);
    legendData2->SetFillStyle(0);
    legendData2->SetTextAlign(12);
    legendData2->SetTextSize(0.035);
    legendData2->AddEntry("", "Uncertainties: stat. (bar), total sys. (open box)", "");
    legendData2->Draw("");
    canvasPzsVsMultiplicity->SaveAs("../PzsVsMultiplicity.pdf");
    canvasPzsVsMultiplicity->SaveAs("../PzsVsMultiplicity.png");
    canvasPzsVsMultiplicity->SaveAs("../PzsVsMultiplicity.eps");

    TCanvas *canvasV2vsPzs = new TCanvas("canvasV2vsPzs", "canvasV2vsPzs", 900, 700);
    StyleCanvas(canvasV2vsPzs, 0.06, 0.15, 0.15, 0.03);
    canvasV2vsPzs->cd();
    TH1F *hDummyV2vsPzs = new TH1F("hDummyV2vsPzs", "hDummyV2vsPzs", 100, 0, 0.12);
    for (Int_t i = 1; i <= hDummyV2vsPzs->GetNbinsX(); i++)
      hDummyV2vsPzs->SetBinContent(i, -1000);
    canvasV2vsPzs->cd();
    SetFont(hDummyV2vsPzs);
    StyleHistoYield(hDummyV2vsPzs, -0.001, 0.006, 1, 1, "v_{2}", "P_{z,s2}", "", 1, 1.15, 1.6);
    SetHistoTextSize(hDummyV2vsPzs, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
    SetTickLength(hDummyV2vsPzs, tickX, tickY);
    hDummyV2vsPzs->GetXaxis()->SetRangeUser(0.03, 0.117);
    hDummyV2vsPzs->GetYaxis()->SetTitleOffset(1.7);
    hDummyV2vsPzs->Draw("");
    TGraphErrors *gV2vsPzsOO = new TGraphErrors(numCentLambdaOO);
    TGraphErrors *gV2vsPzsOOSyst = new TGraphErrors(numCentLambdaOO);
    TGraphErrors *gV2vsEccentricityOO = new TGraphErrors(numCentLambdaOO);
    TGraphErrors *gV2vsEccentricityOOSyst = new TGraphErrors(numCentLambdaOO);
    TGraphErrors *gV2vsMultiplicityOO = new TGraphErrors(numCentLambdaOO);
    TGraphErrors *gV2vsMultiplicityOOSyst = new TGraphErrors(numCentLambdaOO);
    TGraphErrors *gPzsVsEccentricityOO = new TGraphErrors(numCentLambdaOO);
    TGraphErrors *gPzsVsEccentricityOOSyst = new TGraphErrors(numCentLambdaOO);
    TGraphErrors *gPzs2VsEccMinusV2OO = new TGraphErrors(numCentLambdaOO);
    TGraphErrors *gPzs2VsEccMinusV2OOSyst = new TGraphErrors(numCentLambdaOO);

    Float_t v2PubOO = 0;
    Float_t v2PubOOErrStat = 0;
    Float_t v2PubOOErrSyst = 0;
    for (Int_t b = 1; b <= 5; b++)
    {
      if (b == 1)
      {
        v2PubOO = (V2OOPub[0] + V2OOPub[1] + V2OOPub[2] + V2OOPub[3] + V2OOPub[4]) / 5.0; // average 0-5%
        v2PubOO = (v2PubOO + V2OOPub[5]) / 2.0;                                           // average 0-10%
      }
      if (b == 2)
        v2PubOO = (V2OOPub[6] + V2OOPub[7]) / 2.0; // average 10-20%
      if (b == 3)
        v2PubOO = (V2OOPub[8] + V2OOPub[9]) / 2.0; // average 20-30%
      if (b == 4)
        v2PubOO = (V2OOPub[10] + V2OOPub[11]) / 2.0; // average 30-40%
      if (b == 5)
        v2PubOO = (V2OOPub[12] + V2OOPub[13]) / 2.0; // average 40-50%
      if (b == 6)
        v2PubOO = (V2OOPub[14] + V2OOPub[15]) / 2.0; // average 50-60%
      if (b == 1)
      {
        v2PubOOErrStat = TMath::Sqrt(V2OOPubErrStat[0] * V2OOPubErrStat[0] + V2OOPubErrStat[1] * V2OOPubErrStat[1] + V2OOPubErrStat[2] * V2OOPubErrStat[2] + V2OOPubErrStat[3] * V2OOPubErrStat[3] + V2OOPubErrStat[4] * V2OOPubErrStat[4]) / 5.0; // average 0-5%
        v2PubOOErrStat = TMath::Sqrt(v2PubOOErrStat * v2PubOOErrStat + V2OOPubErrStat[5] * V2OOPubErrStat[5]) / 2.0;                                                                                                                               // average 0-10%
      }
      if (b == 2)
        v2PubOOErrStat = TMath::Sqrt(V2OOPubErrStat[6] * V2OOPubErrStat[6] + V2OOPubErrStat[7] * V2OOPubErrStat[7]) / 2.0; // average 10-20%
      if (b == 3)
        v2PubOOErrStat = TMath::Sqrt(V2OOPubErrStat[8] * V2OOPubErrStat[8] + V2OOPubErrStat[9] * V2OOPubErrStat[9]) / 2.0; // average 20-30%
      if (b == 4)
        v2PubOOErrStat = TMath::Sqrt(V2OOPubErrStat[10] * V2OOPubErrStat[10] + V2OOPubErrStat[11] * V2OOPubErrStat[11]) / 2.0; // average 30-40%
      if (b == 5)
        v2PubOOErrStat = TMath::Sqrt(V2OOPubErrStat[12] * V2OOPubErrStat[12] + V2OOPubErrStat[13] * V2OOPubErrStat[13]) / 2.0; // average 40-50%
      if (b == 6)
        v2PubOOErrStat = TMath::Sqrt(V2OOPubErrStat[14] * V2OOPubErrStat[14] + V2OOPubErrStat[15] * V2OOPubErrStat[15]) / 2.0; // average 50-60%
      if (b == 1)
      {
        v2PubOOErrSyst = TMath::Sqrt(V2OOPubErrSyst[0] * V2OOPubErrSyst[0] + V2OOPubErrSyst[1] * V2OOPubErrSyst[1] + V2OOPubErrSyst[2] * V2OOPubErrSyst[2] + V2OOPubErrSyst[3] * V2OOPubErrSyst[3] + V2OOPubErrSyst[4] * V2OOPubErrSyst[4]) / 5.0; // average 0-5%
        v2PubOOErrSyst = TMath::Sqrt(v2PubOOErrSyst * v2PubOOErrSyst + V2OOPubErrSyst[5] * V2OOPubErrSyst[5]) / 2.0;                                                                                                                               // average 0-10%
      }
      if (b == 2)
        v2PubOOErrSyst = TMath::Sqrt(V2OOPubErrSyst[6] * V2OOPubErrSyst[6] + V2OOPubErrSyst[7] * V2OOPubErrSyst[7]) / 2.0; // average 10-20%
      if (b == 3)
        v2PubOOErrSyst = TMath::Sqrt(V2OOPubErrSyst[8] * V2OOPubErrSyst[8] + V2OOPubErrSyst[9] * V2OOPubErrSyst[9]) / 2.0; // average 20-30%
      if (b == 4)
        v2PubOOErrSyst = TMath::Sqrt(V2OOPubErrSyst[10] * V2OOPubErrSyst[10] + V2OOPubErrSyst[11] * V2OOPubErrSyst[11]) / 2.0; // average 30-40%
      if (b == 5)
        v2PubOOErrSyst = TMath::Sqrt(V2OOPubErrSyst[12] * V2OOPubErrSyst[12] + V2OOPubErrSyst[13] * V2OOPubErrSyst[13]) / 2.0; // average 40-50%
      if (b == 6)
        v2PubOOErrSyst = TMath::Sqrt(V2OOPubErrSyst[14] * V2OOPubErrSyst[14] + V2OOPubErrSyst[15] * V2OOPubErrSyst[15]) / 2.0; // average 50-60%
      gV2vsPzsOO->SetPoint(b - 1, v2PubOO, fHistPzs->GetBinContent(b));
      gV2vsPzsOO->SetPointError(b - 1, v2PubOOErrStat, fHistPzs->GetBinError(b));
      gV2vsPzsOOSyst->SetPoint(b - 1, v2PubOO, fHistPzs->GetBinContent(b));
      gV2vsPzsOOSyst->SetPointError(b - 1, v2PubOOErrSyst, fHistPzsSist->GetBinError(b));

      gV2vsEccentricityOO->SetPoint(b - 1, EccOO[b - 1], v2PubOO);
      gV2vsEccentricityOO->SetPointError(b - 1, 0, v2PubOOErrStat);
      gV2vsEccentricityOOSyst->SetPoint(b - 1, EccOO[b - 1], v2PubOO);
      gV2vsEccentricityOOSyst->SetPointError(b - 1, 0, v2PubOOErrSyst);

      gPzsVsEccentricityOO->SetPoint(b - 1, EccOO[b - 1], fHistPzs->GetBinContent(b));
      gPzsVsEccentricityOO->SetPointError(b - 1, 0, fHistPzs->GetBinError(b));
      gPzsVsEccentricityOOSyst->SetPoint(b - 1, EccOO[b - 1], fHistPzsSist->GetBinContent(b));
      gPzsVsEccentricityOOSyst->SetPointError(b - 1, 0, fHistPzsSist->GetBinError(b));

      gPzs2VsEccMinusV2OO->SetPoint(b - 1, EccOO[b - 1] - v2PubOO, fHistPzs->GetBinContent(b));
      gPzs2VsEccMinusV2OO->SetPointError(b - 1, 0, fHistPzs->GetBinError(b));
      gPzs2VsEccMinusV2OOSyst->SetPoint(b - 1, EccOO[b - 1] - v2PubOO, fHistPzsSist->GetBinContent(b));
      gPzs2VsEccMinusV2OOSyst->SetPointError(b - 1, 0, fHistPzsSist->GetBinError(b));

      gV2vsMultiplicityOO->SetPoint(b - 1, dNdEtaOO[b - 1], v2PubOO);
      gV2vsMultiplicityOO->SetPointError(b - 1, dNdEtaOOErr[b - 1], v2PubOOErrStat);
      gV2vsMultiplicityOOSyst->SetPoint(b - 1, dNdEtaOO[b - 1], v2PubOO);
      gV2vsMultiplicityOOSyst->SetPointError(b - 1, dNdEtaOOErrSyst[b - 1], v2PubOOErrSyst);

      // cout << "OO Cent bin " << b << " Pzs " << fHistPzs->GetBinContent(b) << " v2 " << v2PubOO << " +/- " << v2PubOOErrStat << endl;
    }
    TGraphErrors *gV2vsPzsPbPb = new TGraphErrors(numV2PbPbPubCent);
    TGraphErrors *gV2vsPzsPbPbSyst = new TGraphErrors(numV2PbPbPubCent);
    TGraphErrors *gV2vsEccentricityPbPb = new TGraphErrors(numV2PbPbPubCent);
    TGraphErrors *gV2vsEccentricityPbPbSyst = new TGraphErrors(numV2PbPbPubCent);
    TGraphErrors *gPzsVsEccentricityPbPb = new TGraphErrors(numV2PbPbPubCent);
    TGraphErrors *gPzsVsEccentricityPbPbSyst = new TGraphErrors(numV2PbPbPubCent);
    TGraphErrors *gPzs2VsEccMinusV2PbPb = new TGraphErrors(numV2PbPbPubCent);
    TGraphErrors *gPzs2VsEccMinusV2PbPbSyst = new TGraphErrors(numV2PbPbPubCent);
    TGraphErrors *gV2vsMultiplicityPbPb = new TGraphErrors(numV2PbPbPubCent);
    TGraphErrors *gV2vsMultiplicityPbPbSyst = new TGraphErrors(numV2PbPbPubCent);
    for (Int_t b = 1; b <= numV2PbPbPubCent; b++)
    {
      gV2vsPzsPbPb->SetPoint(b - 1, V2PbPbPub[b - 1], fHistPzsLambdaJunlee->GetBinContent(b));
      gV2vsPzsPbPb->SetPointError(b - 1, V2PbPbPubErrStat[b - 1], fHistPzsLambdaJunlee->GetBinError(b));
      gV2vsPzsPbPbSyst->SetPoint(b - 1, V2PbPbPub[b - 1], fHistPzsLambdaJunlee->GetBinContent(b));
      gV2vsPzsPbPbSyst->SetPointError(b - 1, V2PbPbPubErrSys[b - 1], fHistPzsLambdaJunleeSist->GetBinError(b));
      // cout << "PbPb Cent bin " << b << " Pzs " << fHistPzsLambdaJunlee->GetBinContent(b) << " v2 " << V2PbPbPub[b - 1] << " +/- " << V2PbPbPubErrStat[b - 1] << endl;

      gV2vsEccentricityPbPb->SetPoint(b - 1, EccPbPb[b - 1], V2PbPbPub[b - 1]);
      gV2vsEccentricityPbPb->SetPointError(b - 1, 0, V2PbPbPubErrStat[b - 1]);
      gV2vsEccentricityPbPbSyst->SetPoint(b - 1, EccPbPb[b - 1], V2PbPbPub[b - 1]);
      gV2vsEccentricityPbPbSyst->SetPointError(b - 1, 0, V2PbPbPubErrSys[b - 1]);

      gPzsVsEccentricityPbPb->SetPoint(b - 1, EccPbPb[b - 1], fHistPzsLambdaJunlee->GetBinContent(b));
      gPzsVsEccentricityPbPb->SetPointError(b - 1, 0, fHistPzsLambdaJunlee->GetBinError(b));
      gPzsVsEccentricityPbPbSyst->SetPoint(b - 1, EccPbPb[b - 1], fHistPzsLambdaJunleeSist->GetBinContent(b));
      gPzsVsEccentricityPbPbSyst->SetPointError(b - 1, 0, fHistPzsLambdaJunleeSist->GetBinError(b));

      gPzs2VsEccMinusV2PbPb->SetPoint(b - 1, EccPbPb[b - 1] - V2PbPbPub[b - 1], fHistPzsLambdaJunlee->GetBinContent(b));
      gPzs2VsEccMinusV2PbPb->SetPointError(b - 1, 0, fHistPzsLambdaJunlee->GetBinError(b));
      gPzs2VsEccMinusV2PbPbSyst->SetPoint(b - 1, EccPbPb[b - 1] - V2PbPbPub[b - 1], fHistPzsLambdaJunleeSist->GetBinContent(b));
      gPzs2VsEccMinusV2PbPbSyst->SetPointError(b - 1, 0, fHistPzsLambdaJunleeSist->GetBinError(b));

      gV2vsMultiplicityPbPb->SetPoint(b - 1, dNdEtaAbhi[b - 1], V2PbPbPub[b - 1]);
      gV2vsMultiplicityPbPbSyst->SetPoint(b - 1, dNdEtaAbhi[b - 1], V2PbPbPub[b - 1]);
      gV2vsMultiplicityPbPb->SetPointError(b - 1, dNdEtaAbhiErr[b - 1], V2PbPbPubErrStat[b - 1]);
      gV2vsMultiplicityPbPbSyst->SetPointError(b - 1, dNdEtaAbhiErr[b - 1], V2PbPbPubErrSys[b - 1]);
    }
    gV2vsPzsOO->SetMarkerStyle(20);
    gV2vsPzsOO->SetMarkerColor(ColorOO);
    gV2vsPzsOO->SetLineColor(ColorOO);
    gV2vsPzsOO->Draw("same p");
    gV2vsPzsOOSyst->SetMarkerStyle(20);
    gV2vsPzsOOSyst->SetFillStyle(0);
    gV2vsPzsOOSyst->SetMarkerColor(ColorOO);
    gV2vsPzsOOSyst->SetFillColor(ColorOO);
    gV2vsPzsOOSyst->SetLineColor(ColorOO);
    gV2vsPzsOOSyst->Draw("same p2");
    gV2vsPzsPbPb->SetMarkerStyle(21);
    gV2vsPzsPbPb->SetMarkerColor(colorJunlee);
    gV2vsPzsPbPb->SetLineColor(colorJunlee);
    gV2vsPzsPbPb->Draw("same p");
    gV2vsPzsPbPbSyst->SetMarkerStyle(21);
    gV2vsPzsPbPbSyst->SetFillStyle(0);
    gV2vsPzsPbPbSyst->SetMarkerColor(colorJunlee);
    gV2vsPzsPbPbSyst->SetLineColor(colorJunlee);
    gV2vsPzsPbPbSyst->Draw("same p2");
    canvasV2vsPzs->SaveAs("../V2vsPzs.pdf");
    canvasV2vsPzs->SaveAs("../V2vsPzs.png");
    canvasV2vsPzs->SaveAs("../V2vsPzs.eps");

    TCanvas *canvasV2vsEccentricity = new TCanvas("canvasV2vsEccentricity", "canvasV2vsEccentricity", 900, 700);
    StyleCanvas(canvasV2vsEccentricity, 0.06, 0.15, 0.15, 0.03);
    canvasV2vsEccentricity->cd();
    TH1F *hDummyV2vsEccentricity = new TH1F("hDummyV2vsEccentricity", "hDummyV2vsEccentricity", 1000, 0, 1);
    for (Int_t i = 1; i <= hDummyV2vsEccentricity->GetNbinsX(); i++)
      hDummyV2vsEccentricity->SetBinContent(i, -1000);
    canvasV2vsEccentricity->cd();
    SetFont(hDummyV2vsEccentricity);
    StyleHistoYield(hDummyV2vsEccentricity, 0.03, 0.12, 1, 1, "<e>", "v_{2}", "", 1, 1.15, 1.6);
    SetHistoTextSize(hDummyV2vsEccentricity, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
    SetTickLength(hDummyV2vsEccentricity, tickX, tickY);
    hDummyV2vsEccentricity->GetXaxis()->SetRangeUser(0, 1);
    hDummyV2vsEccentricity->GetYaxis()->SetTitleOffset(1.7);
    hDummyV2vsEccentricity->Draw("");
    gV2vsEccentricityOO->SetMarkerStyle(20);
    gV2vsEccentricityOO->SetMarkerColor(ColorOO);
    gV2vsEccentricityOO->SetLineColor(ColorOO);
    gV2vsEccentricityOO->Draw("same p");
    gV2vsEccentricityOOSyst->SetMarkerStyle(20);
    gV2vsEccentricityOOSyst->SetFillStyle(0);
    gV2vsEccentricityOOSyst->SetMarkerColor(ColorOO);
    gV2vsEccentricityOOSyst->SetFillColor(ColorOO);
    gV2vsEccentricityOOSyst->SetLineColor(ColorOO);
    gV2vsEccentricityOOSyst->Draw("same p2");
    gV2vsEccentricityPbPb->SetMarkerStyle(21);
    gV2vsEccentricityPbPb->SetMarkerColor(colorJunlee);
    gV2vsEccentricityPbPb->SetLineColor(colorJunlee);
    gV2vsEccentricityPbPb->Draw("same p");
    gV2vsEccentricityPbPbSyst->SetMarkerStyle(21);
    gV2vsEccentricityPbPbSyst->SetFillStyle(0);
    gV2vsEccentricityPbPbSyst->SetMarkerColor(colorJunlee);
    gV2vsEccentricityPbPbSyst->SetLineColor(colorJunlee);
    gV2vsEccentricityPbPbSyst->Draw("same p2");
    canvasV2vsEccentricity->SaveAs("../V2vsEccentricity.pdf");
    canvasV2vsEccentricity->SaveAs("../V2vsEccentricity.png");
    canvasV2vsEccentricity->SaveAs("../V2vsEccentricity.eps");

    TCanvas *canvasPzsVsEccentricity = new TCanvas("canvasPzsVsEccentricity", "canvasPzsVsEccentricity", 900, 700);
    StyleCanvas(canvasPzsVsEccentricity, 0.06, 0.15, 0.15, 0.03);
    canvasPzsVsEccentricity->cd();
    TH1F *hDummyPzsVsEccentricity = new TH1F("hDummyPzsVsEccentricity", "hDummyPzsVsEccentricity", 1000, 0, 1);
    for (Int_t i = 1; i <= hDummyPzsVsEccentricity->GetNbinsX(); i++)
      hDummyPzsVsEccentricity->SetBinContent(i, -1000);
    canvasPzsVsEccentricity->cd();
    SetFont(hDummyPzsVsEccentricity);
    StyleHistoYield(hDummyPzsVsEccentricity, -0.001, 0.006, 1, 1, "<e>", "P_{z,s2}", "", 1, 1.15, 1.6);
    SetHistoTextSize(hDummyPzsVsEccentricity, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
    SetTickLength(hDummyPzsVsEccentricity, tickX, tickY);
    hDummyPzsVsEccentricity->GetXaxis()->SetRangeUser(0, 1);
    hDummyPzsVsEccentricity->GetYaxis()->SetTitleOffset(1.7);
    hDummyPzsVsEccentricity->Draw("");
    gPzsVsEccentricityOO->SetMarkerStyle(20);
    gPzsVsEccentricityOO->SetMarkerColor(ColorOO);
    gPzsVsEccentricityOO->SetLineColor(ColorOO);
    gPzsVsEccentricityOO->Draw("same p");
    gPzsVsEccentricityOOSyst->SetMarkerStyle(20);
    gPzsVsEccentricityOOSyst->SetFillStyle(0);
    gPzsVsEccentricityOOSyst->SetMarkerColor(ColorOO);
    gPzsVsEccentricityOOSyst->SetFillColor(ColorOO);
    gPzsVsEccentricityOOSyst->SetLineColor(ColorOO);
    gPzsVsEccentricityOOSyst->Draw("same p2");
    gPzsVsEccentricityPbPb->SetMarkerStyle(21);
    gPzsVsEccentricityPbPb->SetMarkerColor(colorJunlee);
    gPzsVsEccentricityPbPb->SetLineColor(colorJunlee);
    gPzsVsEccentricityPbPb->Draw("same p");
    gPzsVsEccentricityPbPbSyst->SetMarkerStyle(21);
    gPzsVsEccentricityPbPbSyst->SetFillStyle(0);
    gPzsVsEccentricityPbPbSyst->SetMarkerColor(colorJunlee);
    gPzsVsEccentricityPbPbSyst->SetLineColor(colorJunlee);
    gPzsVsEccentricityPbPbSyst->Draw("same p2");
    canvasPzsVsEccentricity->SaveAs("../PzsVsEccentricity.pdf");
    canvasPzsVsEccentricity->SaveAs("../PzsVsEccentricity.png");
    canvasPzsVsEccentricity->SaveAs("../PzsVsEccentricity.eps");

    TCanvas *canvasV2vsMultiplicity = new TCanvas("canvasV2vsMultiplicity", "canvasV2vsMultiplicity", 900, 700);
    StyleCanvas(canvasV2vsMultiplicity, 0.06, 0.15, 0.12, 0.03);
    canvasV2vsMultiplicity->cd();
    gPad->SetLogx();
    TH1F *hDummyV2vsMultiplicity = new TH1F("hDummyV2vsMultiplicity", "hDummyV2vsMultiplicity", 2000, 0, 10000);
    for (Int_t i = 1; i <= hDummyV2vsMultiplicity->GetNbinsX(); i++)
      hDummyV2vsMultiplicity->SetBinContent(i, -1000);
    hDummyV2vsMultiplicity->GetYaxis()->SetMaxDigits(1);
    canvasV2vsMultiplicity->cd();
    SetFont(hDummyV2vsMultiplicity);
    StyleHistoYield(hDummyV2vsMultiplicity, 0.03, 0.12, 1, 1, titledNdeta, "v_{2}", "", 1, 1.15, 1.8);
    SetHistoTextSize(hDummyV2vsMultiplicity, xTitleMult, xLabelMult, xOffset, xLabelOffsetMult, yTitle, yLabel, yOffset, yLabelOffset);
    SetTickLength(hDummyV2vsMultiplicity, tickX, tickY);
    hDummyV2vsMultiplicity->GetXaxis()->SetRangeUser(25, 2500);
    hDummyV2vsMultiplicity->GetYaxis()->SetTitleOffset(1.4);
    hDummyV2vsMultiplicity->Draw("");
    lineatZeroVsMult->Draw("same");
    gV2vsMultiplicityPbPb->SetMarkerStyle(21);
    gV2vsMultiplicityPbPb->SetMarkerColor(colorJunlee);
    gV2vsMultiplicityPbPb->SetLineColor(colorJunlee);
    gV2vsMultiplicityPbPb->Draw("same p");
    gV2vsMultiplicityOO->SetMarkerStyle(20);
    gV2vsMultiplicityOO->SetMarkerColor(ColorOO);
    gV2vsMultiplicityOO->SetLineColor(ColorOO);
    gV2vsMultiplicityOO->Draw("same p");
    gV2vsMultiplicityPbPbSyst->SetMarkerStyle(21);
    gV2vsMultiplicityPbPbSyst->SetFillStyle(0);
    gV2vsMultiplicityPbPbSyst->SetMarkerColor(colorJunlee);
    gV2vsMultiplicityPbPbSyst->SetLineColor(colorJunlee);
    gV2vsMultiplicityPbPbSyst->Draw("same p2");
    gV2vsMultiplicityOOSyst->SetMarkerStyle(20);
    gV2vsMultiplicityOOSyst->SetFillStyle(0);
    gV2vsMultiplicityOOSyst->SetMarkerColor(ColorOO);
    gV2vsMultiplicityOOSyst->SetLineColor(ColorOO);
    gV2vsMultiplicityOOSyst->Draw("same p2");
    canvasV2vsMultiplicity->SaveAs("../V2vsMultiplicity.pdf");
    canvasV2vsMultiplicity->SaveAs("../V2vsMultiplicity.png");
    canvasV2vsMultiplicity->SaveAs("../V2vsMultiplicity.eps");

    TCanvas *canvasPzs2VsEccMinusV2 = new TCanvas("canvasPzs2VsEccMinusV2", "canvasPzs2VsEccMinusV2", 900, 700);
    StyleCanvas(canvasPzs2VsEccMinusV2, 0.06, 0.15, 0.15, 0.03);
    canvasPzs2VsEccMinusV2->cd();
    TH1F *hDummyPzs2VsEccMinusV2 = new TH1F("hDummyPzs2VsEccMinusV2", "hDummyPzs2VsEccMinusV2", 1000, 0, 1);
    for (Int_t i = 1; i <= hDummyPzs2VsEccMinusV2->GetNbinsX(); i++)
      hDummyPzs2VsEccMinusV2->SetBinContent(i, -1000);
    canvasPzs2VsEccMinusV2->cd();
    SetFont(hDummyPzs2VsEccMinusV2);
    StyleHistoYield(hDummyPzs2VsEccMinusV2, -0.001, 0.006, 1, 1, "<e> - v_{2}", "P_{z,s2}", "", 1, 1.15, 1.6);
    SetHistoTextSize(hDummyPzs2VsEccMinusV2, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
    SetTickLength(hDummyPzs2VsEccMinusV2, tickX, tickY);
    hDummyPzs2VsEccMinusV2->GetXaxis()->SetRangeUser(0, 1);
    hDummyPzs2VsEccMinusV2->GetYaxis()->SetTitleOffset(1.7);
    hDummyPzs2VsEccMinusV2->Draw("");
    gPzs2VsEccMinusV2OO->SetMarkerStyle(20);
    gPzs2VsEccMinusV2OO->SetMarkerColor(ColorOO);
    gPzs2VsEccMinusV2OO->SetLineColor(ColorOO);
    gPzs2VsEccMinusV2OO->Draw("same p");
    gPzs2VsEccMinusV2OOSyst->SetMarkerStyle(20);
    gPzs2VsEccMinusV2OOSyst->SetFillStyle(0);
    gPzs2VsEccMinusV2OOSyst->SetMarkerColor(ColorOO);
    gPzs2VsEccMinusV2OOSyst->SetFillColor(ColorOO);
    gPzs2VsEccMinusV2OOSyst->SetLineColor(ColorOO);
    gPzs2VsEccMinusV2OOSyst->Draw("same p2");
    gPzs2VsEccMinusV2PbPb->SetMarkerStyle(21);
    gPzs2VsEccMinusV2PbPb->SetMarkerColor(colorJunlee);
    gPzs2VsEccMinusV2PbPb->SetLineColor(colorJunlee);
    gPzs2VsEccMinusV2PbPb->Draw("same p");
    gPzs2VsEccMinusV2PbPbSyst->SetMarkerStyle(21);
    gPzs2VsEccMinusV2PbPbSyst->SetFillStyle(0);
    gPzs2VsEccMinusV2PbPbSyst->SetMarkerColor(colorJunlee);
    gPzs2VsEccMinusV2PbPbSyst->SetLineColor(colorJunlee);
    gPzs2VsEccMinusV2PbPbSyst->Draw("same p2");
    canvasPzs2VsEccMinusV2->SaveAs("../Pzs2VsEccMinusV2.pdf");
    canvasPzs2VsEccMinusV2->SaveAs("../Pzs2VsEccMinusV2.png");
    canvasPzs2VsEccMinusV2->SaveAs("../Pzs2VsEccMinusV2.eps");
  }

  TFile *fileout = new TFile(stringout, "RECREATE");
  fHistMeanSummary->Write();
  fHistSigmaSummary->Write();
  fHistPzs->Write();
  fHistPzsSist->Write();
  fHistPzsBkgUnderPeak->Write();
  fHistPuritySummary->Write();
  fHistYieldSummary->Write();
  fHistSignificanceSummary->Write();
  fHistPzsError->Write();

  fileout->Close();

  cout << "\nStarting from the files (for the different mult): " << PathIn << endl;
  cout << "and the file: " << PathInSyst << " for syst. uncertainties,\n";
  cout << "\nI have created the file:\n " << stringout << endl;
}
