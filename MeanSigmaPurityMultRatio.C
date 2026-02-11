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
// #include "CommonVar.h"
#include "CommonVarLambda.h"
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

// take spectra in input
// produces ratio of spectra wrt 0-100% multiplciity class

Float_t YLowMean[numPart] = {1.31, 1.66, 1.1};
Float_t YUpMean[numPart] = {1.327, 1.68, 1.13};
Float_t YLowSigma[numPart] = {0.0, 0.0, 0.0};
Float_t YUpSigma[numPart] = {0.006, 0.006, 0.004};
Float_t YLowPurity[numPart] = {0.8, 0.8, 0.6};
Float_t YLowV2[numPart] = {-0.3, -0.4, -0.3};
Float_t YUpV2[numPart] = {0.5, 0.5, 0.5};
Float_t YLowPzs2[numPart] = {-0.01, -0.01, -0.01};
Float_t YUpPzs2[numPart] = {0.01, 0.01, 0.01};
Float_t YLowCos2Theta[numPart] = {0, 0, 0};
Float_t YUpCos2Theta[numPart] = {0.4, 0.4, 0.4};
Float_t YLowCos2ThetaLambda[numPart] = {0, 0, 0};
Float_t YUpCos2ThetaLambda[numPart] = {0.5, 0.5, 0.5};

Float_t YLow[numPart] = {0};
Float_t YUp[numPart] = {0};

Float_t YLowRatio[numChoice] = {0.99, 0.2, 0.8, 0.1, 0, -10, -10, 0.9, 0.9, 0, 0.8, 0.8};
Float_t YUpRatio[numChoice] = {1.01, 1.8, 1.2, 4, 1, 10, 10, 1.1, 1.1, 1, 1.2, 1.2};

void MeanSigmaPurityMultRatio(Bool_t isPtAnalysis = 1,
                              Int_t ChosenPart = ChosenParticle,
                              Int_t Choice = 0,
                              Bool_t isTightMassForAcceptancePurity = 0,
                              Int_t ChosenMult = 7 /*commonNumCent*/ /*- 3*/,
                              Bool_t isRapiditySel = ExtrisRapiditySel,
                              TString OutputDir = "../MeanSigmaPurityMultClasses/",
                              Int_t BkgType = ExtrBkgType,
                              Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  Int_t part = 0;
  if (ChosenPart == 1 || ChosenPart == 4 || ChosenPart == 5)
  {
    part = 1;
  }
  else if (ChosenPart == 6)
  {
    part = 2;
  }
  if (isProducedAcceptancePlots && useMixedBDTValueInFitMacro)
  {
    cout << "Either you produce acceptance plots or you use mixed BDT values" << endl;
    return;
  }

  gStyle->SetOptStat(0);
  if ((ChosenMult > commonNumCent && !isV2) || (ChosenMult > (commonNumCent - 1) && isV2))
  {
    cout << "Chosen Mult outside of available range" << endl;
    return;
  }
  if (ChosenPart == 6)
  {
    if (Choice == 5)
      TypeHisto[Choice] = "Pzs2";
  }
  if (!isPtAnalysis)
  {
    if (Choice == 5)
    {
      TypeHisto[Choice] = "PzMixed";
      TitleY[Choice] = "P_{z}";
      YLowPzs2[part] = -0.2;
      YUpPzs2[part] = 0.2;
      YLowRatio[Choice] = 0;
      YUpRatio[Choice] = 2;
    }
    if (Choice == 6)
    {
      TypeHisto[Choice] = "PzLambdaFromCMixed";
      TitleY[Choice] = "P_{z}";
      YLowPzs2[part] = -0.01;
      YUpPzs2[part] = 0.01;
      YLowRatio[Choice] = 0;
      YUpRatio[Choice] = 2;
    }
  }
  //cout << Choice << " " << TypeHisto[Choice] << endl;
  if (Choice > (numChoice - 1))
  {
    cout << "Option not implemented" << endl;
    return;
  }
  if (Choice == 0)
  {
    YLow[part] = YLowMean[part];
    YUp[part] = YUpMean[part];
  }
  else if (Choice == 1)
  {
    YLow[part] = YLowSigma[part];
    YUp[part] = YUpSigma[part];
  }
  else if (Choice == 2)
  {
    YLow[part] = YLowPurity[part];
    YUp[part] = 1;
  }
  else if (Choice == 3)
  {
    YLow[part] = 1e-9;
    YUp[part] = 100;
  }
  else if (Choice == 4 || Choice == 9)
  {
    YLow[part] = YLowV2[part];
    YUp[part] = YUpV2[part];
  }
  else if (Choice == 5 || Choice == 6)
  {
    YLow[part] = YLowPzs2[part];
    YUp[part] = YUpPzs2[part];
  }
  else if (Choice == 7)
  {
    YLow[part] = YLowCos2Theta[part];
    YUp[part] = YUpCos2Theta[part];
  }
  else if (Choice == 8 || Choice == 10 || Choice == 11)
  {
    YLow[part] = YLowCos2ThetaLambda[part];
    YUp[part] = YUpCos2ThetaLambda[part];
  }

  // multiplicity related variables
  TString Smolt[commonNumCent + 1];
  TString SmoltBis[commonNumCent + 1];
  TString sScaleFactorFinal[commonNumCent + 1];
  Float_t ScaleFactorFinal[commonNumCent + 1];

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn;
  TFile *fileIn[commonNumCent + 1];

  // fileout name
  TString stringout;
  TString stringoutpdf;
  TString stringoutShort;
  stringout = OutputDir + "PlotRatios" + NameAnalysis[!isV2] + "_";
  stringout += SinputFileName;
  stringout += "_" + ParticleName[ChosenPart];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[BkgType];
  stringout += "_" + TypeHisto[Choice];
  stringoutShort = OutputDir + "PlotRatios" + "_" + TypeHisto[Choice] + "_";
  stringoutShort += SinputFileName;
  stringoutShort += "_" + ParticleName[ChosenPart];
  if (!isPtAnalysis)
    stringoutShort += "_vsPsi";
  if (isApplyWeights)
    stringout += "_Weighted";
  if (isApplyCentWeight)
    stringout += "_CentWeighted";
  if (!useCommonBDTValue)
    stringout += "_BDTCentDep";
  if (isRun2Binning)
    stringout += "_Run2Binning";
  if (!isPtAnalysis)
    stringout += "_vsPsi";
  if (ExtrisApplyEffWeights)
    stringout += "_EffW";
  if (!isRapiditySel)
    stringout += "_Eta08";
  stringout += STHN[ExtrisFromTHN];
  if (useMixedBDTValueInFitMacro)
    stringout += "_MixedBDT";
  if (isProducedAcceptancePlots)
    stringout += "_AcceptancePlots";
  if (isTightMassCut)
    stringout += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
  if (isReducedPtBins)
    stringout += "_ReducedPtBins";
  if (ExtrisApplyResoOnTheFly)
    stringout += "_ResoOnTheFly";
  if (ChosenPart == 0)
    stringout += "_EPReso";
  if (isTightMassForAcceptancePurity)
    stringout += "_isTightMassForAcceptancePurity";
  stringoutpdf = stringout;
  stringout += ".root";
  TFile *fileout = new TFile(stringout, "RECREATE");

  // canvases
  TCanvas *canvasPtSpectra = new TCanvas("canvasPtSpectra", "canvasPtSpectra", 700, 900);
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

  TH1F *fHistSpectrum[commonNumCent + 1];
  TH1F *fHistSpectrumScaled[commonNumCent + 1];
  TH1F *fHistSpectrumMultRatio[commonNumCent + 1];

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *legendAllMult;
  legendAllMult = new TLegend(0.22, 0.03, 0.73, 0.28);
  if (Choice == 2 && (ChosenPart == 1 || ChosenPart == 4 || ChosenPart == 5))
  {
    legendAllMult = new TLegend(0.44, 0.03, 0.9, 0.28);
  }
  legendAllMult->SetHeader("FT0C Centrality Percentile");
  legendAllMult->SetNColumns(3);
  legendAllMult->SetFillStyle(0);
  TLegendEntry *lheaderAllMult = (TLegendEntry *)legendAllMult->GetListOfPrimitives()->First();
  lheaderAllMult->SetTextSize(0.04);

  TLegend *legendAllMultBis;
  legendAllMultBis = new TLegend(0.22, 0.16, 0.73, 0.36);
  legendAllMultBis->SetHeader("FT0C Centrality Percentile");
  legendAllMultBis->SetNColumns(3);
  legendAllMultBis->SetFillStyle(0);
  TLegendEntry *lheaderAllMultBis = (TLegendEntry *)legendAllMultBis->GetListOfPrimitives()->First();
  lheaderAllMultBis->SetTextSize(0.04);

  TLegend *LegendTitle;
  if (Choice == 2)
    LegendTitle = new TLegend(0.54, 0.35, 0.95, 0.54);
  else
    LegendTitle = new TLegend(0.54, 0.7, 0.95, 0.92);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(33);
  LegendTitle->SetTextSize(0.05);
  LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  if (isOOCentrality)
    LegendTitle->AddEntry("", "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
  else
    LegendTitle->AddEntry("", "PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
  if (isV2)
    LegendTitle->AddEntry("", ParticleNameLegend[ChosenPart] + ", |#it{#eta}| < 0.8", "");
  else
  {
    if (Choice == 6)
    {
      if (isRapiditySel)
        LegendTitle->AddEntry("", ParticleNameLegend[ChosenPart] + " via daughter #Lambda P_{H}, |#it{y} | < 0.5", "");
      else
        LegendTitle->AddEntry("", ParticleNameLegend[ChosenPart] + " via daughter #Lambda P_{H}, |#it{#eta} | < 0.8", "");
    }
    else if (Choice == 10 || Choice == 11)
    {
      LegendTitle->AddEntry("", "#Lambda + #bar{#Lambda}, |#it{#eta} | < 0.8", "");
    }
    else
    {
      if (isRapiditySel)
        LegendTitle->AddEntry("", ParticleNameLegend[ChosenPart] + ", |#it{y} | < 0.5", "");
      else
        LegendTitle->AddEntry("", ParticleNameLegend[ChosenPart] + ", |#it{#eta} | < 0.8", "");
    }
  }

  TLine *lineat1Mult = new TLine(MinPt[ChosenPart], 1, MaxPt[ChosenPart], 1);
  if (!isPtAnalysis)
    lineat1Mult = new TLine(0, 1, 2 * TMath::Pi(), 1);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  // get spectra in multiplicity classes
  for (Int_t m = commonNumCent; m >= 0; m--)
  {
    if ((m == commonNumCent || m == (commonNumCent - 1)) && isV2)
      continue;
    if (isRun2Binning && (m == 0 || (m > (commonNumCent - 2) && m != commonNumCent)))
      continue;
    if (m == commonNumCent)
    { // 0-80%
      CentFT0CMin = 0;
      CentFT0CMax = 80;
    }
    else
    {
      CentFT0CMin = CentFT0C[m];
      CentFT0CMax = CentFT0C[m + 1];
    }
    if (isOOCentrality)
    {
      if (m == numCentLambdaOO)
      {
        CentFT0CMin = 0;
        CentFT0CMax = 100;
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
    if (!isPtAnalysis)
      PathIn += "_vsPsi";
    if ((Choice == 6 || Choice == 8) && ChosenPart != 6)
      PathIn += "_PolFromLambda";
    if (Choice <= 3 && ChosenPart != 6)
      PathIn += "_PolFromLambda";
    if (ExtrisApplyEffWeights)
      PathIn += "_EffW";
    if (!isRapiditySel)
      PathIn += "_Eta08";
    PathIn += STHN[ExtrisFromTHN];
    if (useMixedBDTValueInFitMacro)
      PathIn += "_MixedBDT";
    if (isProducedAcceptancePlots)
      PathIn += "_AcceptancePlots";
    if (isTightMassCut)
      PathIn += Form("_TightMassCut%.1f", Extrsigmacentral[1]);

    if (Choice == 10 || Choice == 11)
    {
      PathIn = "../AcceptancePlots/Acceptance_" + SinputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[0];
      if (isApplyWeights)
        PathIn += "_Weighted";
      // if (isApplyCentWeight)
      //   PathIn += "_CentWeighted";
      if (v2type == 1)
        PathIn += "_SP";
      if (!useCommonBDTValue)
        PathIn += "_BDTCentDep";
      if (isRun2Binning)
        PathIn += "_Run2Binning";
      if (ExtrisApplyEffWeights)
      {
        PathIn += "_EffW";
      }
      PathIn += "_WithAlpha";
      if (!isRapiditySel)
        PathIn += "_Eta08";
      PathIn += STHN[ExtrisFromTHN];
    }
    if (isReducedPtBins)
      PathIn += "_ReducedPtBins";
    if ((isOOCentrality) && (Choice == 10 || Choice == 11))
      PathIn += "_isOOCentrality";
    if (ExtrisApplyResoOnTheFly)
      PathIn += "_ResoOnTheFly";
    if (ChosenPart == 0)
      PathIn += "_EPReso";
    if (isTightMassForAcceptancePurity)
      PathIn += "_isTightMassForAcceptancePurity";
    if (Choice == 2 && isProducedAcceptancePlots)
      PathIn += "_NoMassCutForAcceptance";
    if ((Choice == 10 || Choice == 11) && ChosenParticle == 6)
      PathIn += "_TightAcceptance"; 
    PathIn += ".root";
    cout << "Path in: " << PathIn << endl;
    fileIn[m] = TFile::Open(PathIn);
    //cout << TypeHisto[Choice] + Form("_cent%i-%i", CentFT0CMin, CentFT0CMax) << endl;
    if (Choice == 10 || Choice == 11)
      fHistSpectrum[m] = (TH1F *)fileIn[m]->Get(TypeHisto[Choice] + Form("_cent%i-%i", CentFT0CMin, CentFT0CMax));
    else if (Choice == 1)
      fHistSpectrum[m] = (TH1F *)fileIn[m]->Get("histo" + TypeHisto[Choice]); // + "Weighted");
    else
      fHistSpectrum[m] = (TH1F *)fileIn[m]->Get("histo" + TypeHisto[Choice]);
    if (!fHistSpectrum[m])
    {
      cout << " no hist " << endl;
      return;
    }
    fHistSpectrum[m]->SetName("histoSpectrum_" + Smolt[m]);
  } // end loop on mult

  // draw spectra in multiplicity classes
  Float_t LimSupSpectra = 9.99;
  Float_t LimInfSpectra = 0.2 * 1e-5;

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

  TString titlex = "#it{p}_{T} (GeV/#it{c})";
  if (!isPtAnalysis)
    titlex = "2(#varphi-#Psi_{EP})";

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, -10, 10);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, -100000);
  canvasPtSpectra->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, titlex, TitleY[Choice], "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  if (Choice == 10)
    hDummy->GetXaxis()->SetRangeUser(PtBinsLambda[0], PtBinsLambda[numPtBinsLambda]);
  else if (Choice == 11)
  {
    hDummy->GetXaxis()->SetRangeUser(EtaBins[0], EtaBins[numEtaBins]);
    hDummy->GetXaxis()->SetTitle("#it{#eta}");
  }
  else
    hDummy->GetXaxis()->SetRangeUser(MinPt[ChosenPart], MaxPt[ChosenPart]);
  // hDummy->GetXaxis()->SetRangeUser(1, MaxPt[part]);
  if (!isPtAnalysis)
    hDummy->GetXaxis()->SetRangeUser(0, 2 * TMath::Pi());
  pad1->Draw();
  pad1->cd();
  if (Choice == 3)
    gPad->SetLogy();
  hDummy->Draw("same");

  Int_t NFixedBin = 0;
  for (Int_t m = commonNumCent; m >= 0; m--)
  {
    if ((m == commonNumCent || m == (commonNumCent - 1)) && isV2)
      continue;
    if (isRun2Binning && (m == 0 || (m > (commonNumCent - 2) && m != commonNumCent)))
      continue;

    ScaleFactorFinal[m] = ScaleFactor[m];
    for (Int_t b = 1; b <= fHistSpectrum[m]->GetNbinsX(); b++)
    {
      cout << "bin " << b << " content " << fHistSpectrum[m]->GetBinContent(b) << " error " << fHistSpectrum[m]->GetBinError(b) << endl;
      if (std::isnan(fHistSpectrum[m]->GetBinError(b)) || fHistSpectrum[m]->GetBinError(b)/fHistSpectrum[m]->GetBinContent(b) > 1)
      {
        fHistSpectrum[m]->SetBinError(b, fHistSpectrum[m]->GetBinError(b-1));
        NFixedBin++;
      }
      if (fHistSpectrum[m]->GetBinError(b) != fHistSpectrum[m]->GetBinError(b))
      {
        if (b < fHistSpectrum[m]->GetNbinsX())
        {
          fHistSpectrum[m]->SetBinError(b, fHistSpectrum[m]->GetBinError(b + 1) / fHistSpectrum[m]->GetBinContent(b + 1) * fHistSpectrum[m]->GetBinContent(b));
        }
      } // fuffa
    }
    fHistSpectrumScaled[m] = (TH1F *)fHistSpectrum[m]->Clone("fHistSpectrumScaled_" + Smolt[m]);
    if (Choice == 3)
      fHistSpectrumScaled[m]->Scale(ScaleFactorFinal[m]);
    SetFont(fHistSpectrumScaled[m]);
    fHistSpectrumScaled[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumScaled[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumScaled[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumScaled[m]->SetMarkerSize(0.6 * SizeMult[m]);
    fHistSpectrumScaled[m]->Draw("same e0x0");
    if (Choice == 3)
      sScaleFactorFinal[m] = Form(" (x2^{%i})", int(log2(ScaleFactorFinal[m])));
    else
      sScaleFactorFinal[m] = "";
    legendAllMult->AddEntry(fHistSpectrumScaled[m], SmoltBis[m] + "%" + sScaleFactorFinal[m] + " ", "pef");
    legendAllMultBis->AddEntry(fHistSpectrumScaled[m], SmoltBis[m] + "%", "pef");
  } // end loop on mult
  LegendTitle->Draw("");
  legendAllMult->Draw("");

  // PDG mass
  TF1 *fMassPDG = new TF1("fMassPDG", Form("%f", ParticleMassPDG[ChosenPart]), MinPt[ChosenPart], MaxPt[ChosenPart]);
  fMassPDG->SetLineColor(kBlack);
  fMassPDG->SetLineStyle(8);
  TLegend *legendMassPDG = new TLegend(0.25, 0.76, 0.5, 0.85);
  if (Choice == 0)
  {
    fMassPDG->Draw("same");
    legendMassPDG->AddEntry(fMassPDG, "PDG mass", "l");
    legendMassPDG->SetBorderSize(0);
    legendMassPDG->SetFillStyle(0);
    legendMassPDG->SetTextSize(0.035);
    legendMassPDG->Draw();
  }

  // Compute and draw spectra ratios
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

  TString TitleYSpectraRatio = "Ratio to " + SmoltBis[ChosenMult] + "%";
  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, -10, 10);
  for (Int_t i = 1; i <= hDummyRatio->GetNbinsX(); i++)
    hDummyRatio->SetBinContent(i, 1e-12);
  SetFont(hDummyRatio);
  StyleHistoYield(hDummyRatio, YLowRatio[Choice], YUpRatio[Choice], 1, 1, titlex, TitleYSpectraRatio, "", 1, 1.15, YoffsetSpectraRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetTickLength(hDummyRatio, tickX, tickY);
  if (Choice == 10)
    hDummyRatio->GetXaxis()->SetRangeUser(PtBinsLambda[0], PtBinsLambda[numPtBinsLambda]);
  else if (Choice == 11)
  {
    hDummyRatio->GetXaxis()->SetRangeUser(EtaBins[0], EtaBins[numEtaBins]);
    hDummyRatio->GetXaxis()->SetTitle("#it{#eta}");
  }
  else
  {
    hDummyRatio->GetXaxis()->SetRangeUser(MinPt[ChosenPart], MaxPt[ChosenPart]);
  }
  // hDummyRatio->GetXaxis()->SetRangeUser(1, MaxPt[part]);
  if (!isPtAnalysis)
    hDummyRatio->GetXaxis()->SetRangeUser(0, 2 * TMath::Pi());
  canvasPtSpectra->cd();
  padL1->Draw();
  padL1->cd();
  hDummyRatio->Draw("same");

  for (Int_t m = commonNumCent; m >= 0; m--)
  {
    if ((m == commonNumCent || m == (commonNumCent - 1)) && isV2)
      continue;
    if (isRun2Binning && (m == 0 || (m > (commonNumCent - 2) && m != commonNumCent)))
      continue;

    fHistSpectrumMultRatio[m] = (TH1F *)fHistSpectrum[m]->Clone("fHistSpectrumMultRatio_" + Smolt[m]);
    fHistSpectrumMultRatio[m]->Divide(fHistSpectrum[ChosenMult]);
    ErrRatioCorr(fHistSpectrum[m], fHistSpectrum[ChosenMult], fHistSpectrumMultRatio[m], 0);
    for (Int_t b = 1; b <= fHistSpectrum[m]->GetNbinsX(); b++)
    {
      // cout << "bin " << b << " " << fHistSpectrum[m]->GetBinContent(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrum[ChosenMult]->GetBinContent(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrumMultRatio[m]->GetBinContent(b) << endl;
      // fHistSpectrumMultRatio[m]->SetBinError(b, fHistSpectrumMultRatio[commonNumCent]->GetBinError(b));
    }
    fHistSpectrumMultRatio[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumMultRatio[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumMultRatio[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumMultRatio[m]->SetMarkerSize(SizeMultRatio[m]);

    if (m != ChosenMult)
    {
      fHistSpectrumMultRatio[m]->Draw("same e0x0");
    }
    lineat1Mult->Draw("same");

  } // end loop on mult

  TCanvas *canvas = new TCanvas("canvas", "canvas", 700, 900);
  StyleCanvas(canvas, 0.03, 0.1, 0.15, 0.02);
  canvas->cd();
  if (Choice == 3)
    gPad->SetLogy();
  hDummy->GetXaxis()->SetLabelSize(0.05);
  SetHistoTextSize(hDummy, xTitleR, xLabelR, 1, 0.015, yTitle, yLabel, 1.8, yLabelOffset);
  hDummy->Draw("same");
  for (Int_t m = commonNumCent; m >= 0; m--)
  {
    if ((m == commonNumCent || m == (commonNumCent - 1)) && isV2)
      continue;
    if (isRun2Binning && (m == 0 || (m > (commonNumCent - 2) && m != commonNumCent)))
      continue;

    fHistSpectrumScaled[m]->Draw("same e0x0");
  }
  LegendTitle->Draw("");
  legendAllMultBis->Draw("");

  fileout->Close();
  canvasPtSpectra->SaveAs(stringoutShort + ".pdf");
  canvasPtSpectra->SaveAs(stringoutShort + ".png");
  canvas->SaveAs(stringoutShort + "_Top.pdf");
  canvas->SaveAs(stringoutShort + "_Top.png");
  cout << "\nI have created the file:\n " << stringout << endl;

  if (NFixedBin > 0)
  {
    cout << "\nI had to fix " << NFixedBin << " bins with very large errors" << endl;
  }
}
