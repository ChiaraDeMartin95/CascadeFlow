#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "CommonVarLambda.h"

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

Float_t YLow[numPart] = {0};
Float_t YUp[numPart] = {0.0015};

Float_t AccRelError[numCent + 1] = {0.05, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.02};
Float_t LambdaDecayParameterRelError = 0.01;
Float_t TransferCoefficienctRelError = 0.0043;
Float_t RunByRunAccRelError = 0.01;
// Float_t ResoRelError[numCentLambdaOO + 1] = {0.045, 0.05, 0.06, 0.065, 0.07, 0.09, 0.127, 0.1915, 0.2865, 0.39};
Float_t ResoRelError[numCentLambdaOO + 1] = {0};
Float_t PrimaryLambdaFraction = 0.03;
Float_t SecondaryLambdaFraction = 0.1;
Float_t ZVertexErrorLambdaOO = 0.00009;

void SystematicErrorVsCent(Int_t ChosenPart = ChosenParticle,
                           Bool_t isPolFromLambda = 0,
                           Float_t nsigmaBarlowMassCut = 1.0,
                           Bool_t isFromFit = 0,
                           Bool_t isRapiditySel = ExtrisRapiditySel,
                           Int_t BkgType = ExtrBkgType,
                           Bool_t UseTwoGauss = ExtrUseTwoGauss,
                           Bool_t isPtAnalysis = 1,  // 1 for V2 vs pt and Pzs2 vs pt, 0 for Pz vs 2(phi-Psi)
                           Bool_t isPtIntegrated = 1 // 1 for results integrated in pt / phi
)
{

  Int_t ChosenPt = -999;
  cout << "Type 100 if you want to analyse Pz (integrated in pT), or the number of the pt interval you want" << endl;
  cin >> ChosenPt;

  Int_t part = 0;
  if (ChosenPart == 1 || ChosenPart == 4 || ChosenPart == 5)
  {
    part = 1;
  }
  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // fileIn Default
  TString PathIn = "";
  TFile *fileIn[commonNumCent + 1];

  TString PathIn1 = "";
  TString PathIn2 = "";
  // fileInBDT
  TString PathInBDT;
  TFile *fileInBDT[commonNumCent + 1];
  // fileIn MassCut
  TString PathInMassCut;
  TFile *fileInMassCut[commonNumCent + 1];
  // fileIn MassCutAndBDT
  TString PathInMassCutAndBDT;
  TFile *fileInMassCutAndBDT[commonNumCent + 1];

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = "../Systematics/SystVsCentrality_" + NameAnalysis[!isV2] + "_";
  stringout += SinputFileNameSyst;
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
  if (!isRapiditySel || ExtrisFromTHN)
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
  stringoutpdf = stringout;
  stringout += ".root";

  // canvases
  gStyle->SetOptStat(0);
  TCanvas *canvasError = new TCanvas("canvasError", "canvasError", 800, 600);
  StyleCanvas(canvasError, 0.05, 0.2, 0.2, 0.03);

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  TH1F *fHistPzs2[commonNumCent + 1];
  TH1F *fHistBDTError[commonNumCent + 1];
  TH1F *fHistMassCutError[commonNumCent + 1];
  TH1F *fHistMassCutAndBDTError[commonNumCent + 1];
  TH1F *fHistMeanSystMultiTrialMean[commonNumCent + 1];
  TH1F *fHistBDTErrorVsCent;
  TH1F *fHistMassCutErrorVsCent;
  TH1F *fHistMassCutAndBDTErrorVsCent;
  TH1F *fHistPrimaryLambdaErrorVsCent;
  TH1F *fHistAccErrorVsCent;
  TH1F *fHistResoErrorVsCent;
  TH1F *fHistDecayParErrorVsCent;
  TH1F *fHistTransferCoeffErrorVsCent;
  TH1F *fHistRunByRunAccErrorVsCent;
  TH1F *fHistPolBkg0ErrorVsCent;
  TH1F *fHistPzFitRangeErrorVsCent;
  TH1F *fHistBkgExpoErrorVsCent;
  TH1F *fHistZVertexErrorVsCent;
  TH1F *fHistTotalErrorVsCent;
  TH1F *fHistMeanVsCent;
  if (isOOCentrality)
  {
    fHistBDTErrorVsCent = new TH1F("fHistBDTErrorVsCent", "fHistBDTErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistMassCutErrorVsCent = new TH1F("fHistMassCutErrorVsCent", "fHistMassCutErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistMassCutAndBDTErrorVsCent = new TH1F("fHistMassCutAndBDTErrorVsCent", "fHistMassCutAndBDTErrorVsCent  ", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistPrimaryLambdaErrorVsCent = new TH1F("fHistPrimaryLambdaErrorVsCent", "fHistPrimaryLambdaErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistAccErrorVsCent = new TH1F("fHistAccErrorVsCent", "fHistAccErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistResoErrorVsCent = new TH1F("fHistResoErrorVsCent", "fHistResoErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistDecayParErrorVsCent = new TH1F("fHistDecayParErrorVsCent", "fHistDecayParErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistTransferCoeffErrorVsCent = new TH1F("fHistTransferCoeffErrorVsCent", "fHistTransferCoeffErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistRunByRunAccErrorVsCent = new TH1F("fHistRunByRunAccErrorVsCent", "fHistRunByRunAccErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistPolBkg0ErrorVsCent = new TH1F("fHistPolBkg0ErrorVsCent", "fHistPolBkg0ErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistPzFitRangeErrorVsCent = new TH1F("fHistPzFitRangeErrorVsCent", "fHistPzFitRangeErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistBkgExpoErrorVsCent = new TH1F("fHistBkgExpoErrorVsCent", "fHistBkgExpoErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistZVertexErrorVsCent = new TH1F("fHistZVertexErrorVsCent", "fHistZVertexErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistTotalErrorVsCent = new TH1F("fHistTotalErrorVsCent", "fHistTotalErrorVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
    fHistMeanVsCent = new TH1F("fHistMeanVsCent", "fHistMeanVsCent", numCentLambdaOO, fCentFT0CLambdaOO);
  }
  else
  {
    fHistBDTErrorVsCent = new TH1F("fHistBDTErrorVsCent", "fHistBDTErrorVsCent", numCent, fCentFT0C);
    fHistMassCutErrorVsCent = new TH1F("fHistMassCutErrorVsCent", "fHistMassCutErrorVsCent", numCent, fCentFT0C);
    fHistMassCutAndBDTErrorVsCent = new TH1F("fHistMassCutAndBDTErrorVsCent", "fHistMassCutAndBDTErrorVsCent", numCent, fCentFT0C);
    fHistPrimaryLambdaErrorVsCent = new TH1F("fHistPrimaryLambdaErrorVsCent", "fHistPrimaryLambdaErrorVsCent", numCent, fCentFT0C);
    fHistAccErrorVsCent = new TH1F("fHistAccErrorVsCent", "fHistAccErrorVsCent", numCent, fCentFT0C);
    fHistResoErrorVsCent = new TH1F("fHistResoErrorVsCent", "fHistResoErrorVsCent", numCent, fCentFT0C);
    fHistDecayParErrorVsCent = new TH1F("fHistDecayParErrorVsCent", "fHistDecayParErrorVsCent", numCent, fCentFT0C);
    fHistTransferCoeffErrorVsCent = new TH1F("fHistTransferCoeffErrorVsCent", "fHistTransferCoeffErrorVsCent", numCent, fCentFT0C);
    fHistRunByRunAccErrorVsCent = new TH1F("fHistRunByRunAccErrorVsCent", "fHistRunByRunAccErrorVsCent", numCent, fCentFT0C);
    fHistPolBkg0ErrorVsCent = new TH1F("fHistPolBkg0ErrorVsCent", "fHistPolBkg0ErrorVsCent", numCent, fCentFT0C);
    fHistPzFitRangeErrorVsCent = new TH1F("fHistPzFitRangeErrorVsCent", "fHistPzFitRangeErrorVsCent", numCent, fCentFT0C);
    fHistBkgExpoErrorVsCent = new TH1F("fHistBkgExpoErrorVsCent", "fHistBkgExpoErrorVsCent", numCent, fCentFT0C);
    fHistZVertexErrorVsCent = new TH1F("fHistZVertexErrorVsCent", "fHistZVertexErrorVsCent", numCent, fCentFT0C);
    fHistTotalErrorVsCent = new TH1F("fHistTotalErrorVsCent", "fHistTotalErrorVsCent", numCent, fCentFT0C);
    fHistMeanVsCent = new TH1F("fHistMeanVsCent", "fHistMeanVsCent", numCent, fCentFT0C);
  }

  TFile *fileInResoError = TFile::Open("../CompareResults/SystUncertainty_Reso1.root");
  TH1F *fHistResoError = (TH1F *)fileInResoError->Get("hRatioClone_1");
  if (!fHistResoError)
  {
    cout << "Error: histogram Reso not found" << endl;
    return;
  }
  fHistResoError->SetName("hResoSystError");

  TFile *fileInPolBkg0 = TFile::Open("../CompareResults/SystUncertainty_BkgPol0.root");
  TH1F *fHistPolBkg0Error = (TH1F *)fileInPolBkg0->Get("hRatioClone_1");
  if (!fHistPolBkg0Error)
  {
    cout << "Error: histogram PolBkg0 not found" << endl;
    return;
  }
  fHistPolBkg0Error->SetName("hPolBkg0SystError");

  TFile* fileInBkgExpo = TFile::Open("../CompareResults/SystUncertainty_BkgFit.root");
  TH1F* fHistBkgExpoError = (TH1F*)fileInBkgExpo->Get("hRatioClone_1");
  if (!fHistBkgExpoError)
  {
    cout << "Error: histogram BkgExpo not found" << endl;
    return;
  }
  fHistBkgExpoError->SetName("hBkgExpoSystError");

  TFile *fileInPzFitRange = TFile::Open("../CompareResults/SystUncertainty_PzFitRange.root");
  TH1F *fHistPzFitRangeError = (TH1F *)fileInPzFitRange->Get("hRatioClone_1");
  if (!fHistPzFitRangeError)
  {
    cout << "Error: histogram PzFitRange not found" << endl;
    return;
  }
  fHistPzFitRangeError->SetName("hPzFitRangeSystError");

  TString Smolt[commonNumCent + 1];
  // get spectra in multiplicity classes
  for (Int_t m = 0; m < commonNumCent; m++)
  {
    if (m == numCent)
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
    Smolt[m] += Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);

    PathIn = "../OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_";
    PathIn += SinputFileName;
    PathIn += "_" + ParticleName[ChosenPart];
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[BkgType];
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
    if (!isRapiditySel || ExtrisFromTHN)
      PathIn += "_Eta08";
    // if (isReducedPtBins)
    //   PathIn += "_ReducedPtBins";
    PathIn += STHN[ExtrisFromTHN];
    if (useMixedBDTValueInFitMacro)
      PathIn += "_MixedBDT";
    if (isTightMassCut)
      PathIn += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
    if (isReducedPtBins)
      PathIn += "_ReducedPtBins";
    if (ChosenPart == 6)
    {
      // PathIn += Form("_SysMultTrial_%i", 0);
      // PathIn += "_isSysLambdaMultTrial_";
      PathIn += "_ResoOnTheFly";
    }

    PathIn += ".root";
    cout << "Path in : " << PathIn << endl;
    fileIn[m] = TFile::Open(PathIn);
    fHistPzs2[m] = (TH1F *)fileIn[m]->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + "PtInt" + V2FromFit[isFromFit]);
    if (!fHistPzs2[m])
    {
      cout << " no hist v2 / Pzs" << endl;
      return;
    }

    PathIn1 = "../Systematics/SystMultiTrial_" + SinputFileNameSyst + Form("_%i-%i_", CentFT0CMin, CentFT0CMax) + ParticleName[ChosenPart] + "_";
    PathIn2 = "";
    if (isPtIntegrated)
      PathIn2 += "_PtInt";
    if (isApplyWeights)
      PathIn2 += "_Weighted";
    if (isApplyCentWeight)
      PathIn2 += "_CentWeighted";
    if (!isPtAnalysis)
      PathIn2 += "_vsPsi";
    if (!isV2 && isPolFromLambda)
      PathIn2 += "_PolFromLambda";
    if (ExtrisApplyEffWeights)
      PathIn2 += "_EffW";
    if (!isRapiditySel || ExtrisFromTHN)
      PathIn2 += "_Eta08";
    PathIn2 += STHN[ExtrisFromTHN];
    if (nsigmaBarlow != 0)
      PathIn2 += Form("_nsigmaBarlow%.1f", nsigmaBarlow);
    if (ChosenPart == 6)
    {
      if (isTightMassCut)
        PathIn2 += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
      if (isReducedPtBins)
        PathIn2 += "_ReducedPtBins";
      PathIn2 += "_ResoOnTheFly";
    }
    // PathIn2 += "_NewTest";
    PathInBDT = PathIn1 + "BDT" + PathIn2 + ".root";
    if (ChosenParticle != 6)
      cout << "Path in BDT: " << PathInBDT << endl;
    fileInBDT[m] = TFile::Open(PathInBDT);
    PathInMassCut = PathIn1 + "MassCut" + PathIn2 + Form("_nsigmaBarlow%.1f.root", nsigmaBarlowMassCut);
    if (ChosenParticle != 6)
      cout << "Path in MassCut: " << PathInMassCut << endl;
    fileInMassCut[m] = TFile::Open(PathInMassCut);
    PathInMassCutAndBDT = PathIn1;
    if (ChosenPart == 6)
      PathInMassCutAndBDT += "LambdaTopo";
    else
      PathInMassCutAndBDT += "MassAndBDTCut";
    PathInMassCutAndBDT += PathIn2 + ".root";
    cout << "Path in MassCutAndBDT: " << PathInMassCutAndBDT << endl;
    fileInMassCutAndBDT[m] = TFile::Open(PathInMassCutAndBDT);

    // fHistBDTError[m] = (TH1F *)fileInBDT[m]->Get("hAbsoluteMaxDev");
    fHistBDTError[m] = (TH1F *)fileInMassCutAndBDT[m]->Get("hAbsoluteMaxDev");
    if (!fHistBDTError[m])
    {
      cout << "Error: histogram not found" << endl;
      return;
    }
    fHistBDTError[m]->SetName("hAbsoluteSystErrorBDT_" + Smolt[m]);
    fHistBDTError[m]->Reset();

    // fHistMassCutError[m] = (TH1F *)fileInMassCut[m]->Get("hAbsoluteMaxDev");
    fHistMassCutError[m] = (TH1F *)fileInMassCutAndBDT[m]->Get("hAbsoluteMaxDev");
    if (!fHistMassCutError[m])
    {
      cout << "Error: histogram not found" << endl;
      return;
    }
    fHistMassCutError[m]->SetName("hAbsoluteSystErrorMassCut_" + Smolt[m]);
    fHistMassCutError[m]->Reset();

    // fHistMassCutAndBDTError[m] = (TH1F *)fileInMassCutAndBDT[m]->Get("hRMS");

    fHistMeanSystMultiTrialMean[m] = (TH1F *)fileInMassCutAndBDT[m]->Get("hSystMultiTrialMean");
    if (!fHistMeanSystMultiTrialMean[m])
    {
      cout << "Error: histogram not found" << endl;
      return;
    }
    fHistMeanSystMultiTrialMean[m]->SetName("hMeanSystMultiTrialMean_" + Smolt[m]);

    fHistMassCutAndBDTError[m] = (TH1F *)fileInMassCutAndBDT[m]->Get("hSystMultiTrial2");
    if (!fHistMassCutAndBDTError[m])
    {
      cout << "Error: histogram not found" << endl;
      return;
    }
    fHistMassCutAndBDTError[m]->SetName("hAbsoluteSystErrorMassCutAndBDT_" + Smolt[m]);

    if (ChosenPart == 6)
    {
      fHistPrimaryLambdaErrorVsCent->SetBinContent(m + 1, std::abs(SecondaryLambdaFraction * fHistPzs2[m]->GetBinContent(1)));
      fHistPrimaryLambdaErrorVsCent->SetBinError(m + 1, 0);
    }
    else
    {
      fHistPrimaryLambdaErrorVsCent->SetBinContent(m + 1, std::abs(PrimaryLambdaFraction * fHistPzs2[m]->GetBinContent(1)));
      fHistPrimaryLambdaErrorVsCent->SetBinError(m + 1, 0);
    }
    fHistBDTErrorVsCent->SetBinContent(m + 1, fHistBDTError[m]->GetBinContent(1));
    fHistBDTErrorVsCent->SetBinError(m + 1, 0);
    fHistMassCutErrorVsCent->SetBinContent(m + 1, fHistMassCutError[m]->GetBinContent(1));
    fHistMassCutErrorVsCent->SetBinError(m + 1, 0);
    //cout << "fHistMassCutAndBDTError[m]->GetBinContent(1): " << fHistMassCutAndBDTError[m]->GetBinContent(1) << endl;
    fHistMassCutAndBDTErrorVsCent->SetBinContent(m + 1, fHistMassCutAndBDTError[m]->GetBinContent(1));
    fHistMassCutAndBDTErrorVsCent->SetBinError(m + 1, 0);
    fHistMeanVsCent->SetBinContent(m + 1, fHistMeanSystMultiTrialMean[m]->GetBinContent(1));
    fHistMeanVsCent->SetBinError(m + 1, 0);
    if (ChosenPart == 6)
    {
      fHistAccErrorVsCent->SetBinContent(m + 1, 0);
      fHistAccErrorVsCent->SetBinError(m + 1, 0);
      fHistResoErrorVsCent->SetBinContent(m + 1, abs(fHistResoError->GetBinContent(m + 1)));
      //cout << " Reso error cent " << m << " : " << ResoRelError[m] * fHistPzs2[m]->GetBinContent(1) << endl;
      //cout << " Reso rel error cent " << m << " : " << ResoRelError[m] << endl;
      //cout << " Pzs2 cent " << m << " : " << fHistPzs2[m]->GetBinContent(1) << endl;
      fHistResoErrorVsCent->SetBinError(m + 1, 0);
      fHistDecayParErrorVsCent->SetBinContent(m + 1, std::abs(LambdaDecayParameterRelError * fHistPzs2[m]->GetBinContent(1)));
      fHistDecayParErrorVsCent->SetBinError(m + 1, 0);
      fHistTransferCoeffErrorVsCent->SetBinContent(m + 1, 0);
      fHistTransferCoeffErrorVsCent->SetBinError(m + 1, 0);
      fHistRunByRunAccErrorVsCent->SetBinContent(m + 1, 0);
      fHistRunByRunAccErrorVsCent->SetBinError(m + 1, 0);
      fHistPolBkg0ErrorVsCent->SetBinContent(m + 1, abs(fHistPolBkg0Error->GetBinContent(m + 1)));
      //fHistPolBkg0ErrorVsCent->SetBinContent(m + 1, 0); // not significant
      fHistPolBkg0ErrorVsCent->SetBinError(m + 1, 0);
      fHistPzFitRangeErrorVsCent->SetBinContent(m + 1, abs(fHistPzFitRangeError->GetBinContent(m + 1)));
      fHistPzFitRangeErrorVsCent->SetBinError(m + 1, 0);
      fHistBkgExpoErrorVsCent->SetBinContent(m + 1, abs(fHistBkgExpoError->GetBinContent(m + 1)));
      fHistBkgExpoErrorVsCent->SetBinError(m + 1, 0);
      fHistZVertexErrorVsCent->SetBinContent(m + 1, ZVertexErrorLambdaOO);
      fHistZVertexErrorVsCent->SetBinError(m + 1, 0);
    }
    else
    {
      fHistAccErrorVsCent->SetBinContent(m + 1, std::abs(AccRelError[m] * fHistPzs2[m]->GetBinContent(1)));
      fHistAccErrorVsCent->SetBinError(m + 1, 0);
      fHistResoErrorVsCent->SetBinContent(m + 1, 0);
      fHistResoErrorVsCent->SetBinError(m + 1, 0);
      fHistDecayParErrorVsCent->SetBinContent(m + 1, std::abs(LambdaDecayParameterRelError * fHistPzs2[m]->GetBinContent(1)));
      fHistDecayParErrorVsCent->SetBinError(m + 1, 0);
      fHistTransferCoeffErrorVsCent->SetBinContent(m + 1, std::abs(TransferCoefficienctRelError * fHistPzs2[m]->GetBinContent(1)));
      fHistTransferCoeffErrorVsCent->SetBinError(m + 1, 0);
      fHistRunByRunAccErrorVsCent->SetBinContent(m + 1, std::abs(RunByRunAccRelError * fHistPzs2[m]->GetBinContent(1)));
      fHistRunByRunAccErrorVsCent->SetBinError(m + 1, 0);
      fHistPolBkg0ErrorVsCent->SetBinContent(m + 1, 0);
      fHistPolBkg0ErrorVsCent->SetBinError(m + 1, 0);
      fHistPzFitRangeErrorVsCent->SetBinContent(m + 1, 0);
      fHistPzFitRangeErrorVsCent->SetBinError(m + 1, 0);
      fHistBkgExpoErrorVsCent->SetBinContent(m + 1, 0);
      fHistBkgExpoErrorVsCent->SetBinError(m + 1, 0);
      fHistZVertexErrorVsCent->SetBinContent(m + 1, 0);
      fHistZVertexErrorVsCent->SetBinError(m + 1, 0);
    }
  } // end loop on mult

  fHistBDTErrorVsCent->Smooth();
  fHistMassCutErrorVsCent->Smooth();
  fHistMassCutAndBDTErrorVsCent->Smooth();
  //fHistBkgExpoErrorVsCent->Smooth();
  //fHistPzFitRangeErrorVsCent->Smooth();
  //fHistPolBkg0ErrorVsCent->Smooth();
  // fHistResoErrorVsCent->Smooth();

  for (Int_t m = 0; m < commonNumCent; m++)
  {
    // fHistTotalErrorVsCent->SetBinContent(m + 1, TMath::Sqrt(fHistBDTErrorVsCent->GetBinContent(m + 1) * fHistBDTErrorVsCent->GetBinContent(m + 1) + fHistMassCutErrorVsCent->GetBinContent(m + 1) * fHistMassCutErrorVsCent->GetBinContent(m + 1)));
    // fHistTotalErrorVsCent->SetBinContent(m + 1, TMath::Sqrt(fHistBDTErrorVsCent->GetBinContent(m + 1) * fHistBDTErrorVsCent->GetBinContent(m + 1) +
    //                                                      fHistMassCutErrorVsCent->GetBinContent(m + 1) * fHistMassCutErrorVsCent->GetBinContent(m + 1) +
    //                                                    fHistPrimaryLambdaErrorVsCent->GetBinContent(m + 1) * fHistPrimaryLambdaErrorVsCent->GetBinContent(m + 1)));
    fHistTotalErrorVsCent->SetBinContent(m + 1, TMath::Sqrt(fHistMassCutAndBDTErrorVsCent->GetBinContent(m + 1) * fHistMassCutAndBDTErrorVsCent->GetBinContent(m + 1) +
                                                            fHistPrimaryLambdaErrorVsCent->GetBinContent(m + 1) * fHistPrimaryLambdaErrorVsCent->GetBinContent(m + 1) +
                                                            fHistAccErrorVsCent->GetBinContent(m + 1) * fHistAccErrorVsCent->GetBinContent(m + 1) +
                                                            fHistResoErrorVsCent->GetBinContent(m + 1) * fHistResoErrorVsCent->GetBinContent(m + 1) +
                                                            fHistTransferCoeffErrorVsCent->GetBinContent(m + 1) * fHistTransferCoeffErrorVsCent->GetBinContent(m + 1) +
                                                            fHistDecayParErrorVsCent->GetBinContent(m + 1) * fHistDecayParErrorVsCent->GetBinContent(m + 1) +
                                                            fHistRunByRunAccErrorVsCent->GetBinContent(m + 1) * fHistRunByRunAccErrorVsCent->GetBinContent(m + 1) +
                                                            fHistZVertexErrorVsCent->GetBinContent(m + 1) * fHistZVertexErrorVsCent->GetBinContent(m + 1) +
                                                            fHistPolBkg0ErrorVsCent->GetBinContent(m + 1) * fHistPolBkg0ErrorVsCent->GetBinContent(m + 1) +
                                                            fHistPzFitRangeErrorVsCent->GetBinContent(m + 1) * fHistPzFitRangeErrorVsCent->GetBinContent(m + 1) +
                                                            fHistBkgExpoErrorVsCent->GetBinContent(m + 1) * fHistBkgExpoErrorVsCent->GetBinContent(m + 1)));
    fHistTotalErrorVsCent->SetBinError(m + 1, 0);
  }

  Float_t xTitle = 30;
  Float_t xOffset = 1.5;
  Float_t yTitle = 30;
  Float_t yOffset = 3;

  Float_t xLabel = 30;
  Float_t yLabel = 30;
  Float_t xLabelOffset = 0.03;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.03;
  Float_t tickY = 0.042;

  if (ChosenPart == 6)
    YUp[part] = 0.0015;

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasError->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleXCent, "Absolute syst. uncertainty", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(0, 80);
  if (ChosenPart == 6)
    hDummy->GetXaxis()->SetRangeUser(0, 50);
  hDummy->Draw("");
  fHistBDTErrorVsCent->SetLineColor(kRed);
  fHistMassCutErrorVsCent->SetLineColor(kGreen + 2);
  fHistMassCutAndBDTErrorVsCent->SetLineColor(kGreen + 2);
  fHistPrimaryLambdaErrorVsCent->SetLineColor(kOrange + 2);
  fHistAccErrorVsCent->SetLineColor(kBlue);
  fHistResoErrorVsCent->SetLineColor(kMagenta);
  fHistDecayParErrorVsCent->SetLineColor(kRed + 1);
  fHistTransferCoeffErrorVsCent->SetLineColor(kViolet + 1);
  fHistRunByRunAccErrorVsCent->SetLineColor(kCyan);
  fHistPolBkg0ErrorVsCent->SetLineColor(kGray + 1);
  fHistPolBkg0ErrorVsCent->SetLineStyle(2);
  fHistPolBkg0ErrorVsCent->SetLineWidth(3);
  fHistPzFitRangeErrorVsCent->SetLineColor(kOrange - 3);
  fHistPzFitRangeErrorVsCent->SetLineStyle(2);
  fHistPzFitRangeErrorVsCent->SetLineWidth(3);
  fHistBkgExpoErrorVsCent->SetLineColor(kPink + 1);
  fHistBkgExpoErrorVsCent->SetLineStyle(2);
  fHistBkgExpoErrorVsCent->SetLineWidth(3);
  fHistZVertexErrorVsCent->SetLineColor(kGreen);
  fHistTotalErrorVsCent->SetLineColor(kBlack);
  fHistTotalErrorVsCent->SetLineWidth(4);
  // fHistBDTErrorVsCent->Draw("same");
  // fHistMassCutErrorVsCent->Draw("same");

  fHistAccErrorVsCent->Draw("same");
  fHistMassCutAndBDTErrorVsCent->Draw("same");
  fHistPrimaryLambdaErrorVsCent->Draw("same");
  fHistDecayParErrorVsCent->Draw("same");
  if (ChosenPart != 6)
  {
    fHistTransferCoeffErrorVsCent->Draw("same");
    fHistRunByRunAccErrorVsCent->Draw("same");
  }
  if (ChosenPart == 6)
  {
    fHistResoErrorVsCent->Draw("same");
    fHistPolBkg0ErrorVsCent->Draw("same");
    fHistPzFitRangeErrorVsCent->Draw("same");
    fHistBkgExpoErrorVsCent->Draw("same");
    fHistZVertexErrorVsCent->Draw("same");
  }
  fHistTotalErrorVsCent->Draw("same");
  TLegend *legend = new TLegend(0.3, 0.4, 0.9, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.05);
  if (ChosenPart == 6)
    legend->AddEntry(fHistMassCutAndBDTErrorVsCent, "Topological selections", "l");
  else
    legend->AddEntry(fHistMassCutAndBDTErrorVsCent, "Mass Cut + BDT", "l");
  legend->AddEntry(fHistPrimaryLambdaErrorVsCent, "Primary #Lambda", "l");
  // legend->AddEntry(fHistBDTErrorVsCent, "BDT", "l");
  // legend->AddEntry(fHistMassCutErrorVsCent, "Mass Cut", "l");
  if (ChosenPart != 6)
  {
    legend->AddEntry(fHistAccErrorVsCent, "Acceptance", "l");
    legend->AddEntry(fHistTransferCoeffErrorVsCent, "Transfer coefficient", "l");
    legend->AddEntry(fHistRunByRunAccErrorVsCent, "Run-by-run acceptance dep.", "l");
  }
  if (ChosenPart == 6)
  {
    legend->AddEntry(fHistResoErrorVsCent, "Resolution", "l");
    legend->AddEntry(fHistZVertexErrorVsCent, "Z_{vtx} selection", "l");
    legend->AddEntry(fHistPolBkg0ErrorVsCent, "P_{z, s2, bkg} = 0", "l");
    legend->AddEntry(fHistPzFitRangeErrorVsCent, "P_{z} fit range", "l");
    legend->AddEntry(fHistBkgExpoErrorVsCent, "Background fit function", "l");
  }
  legend->AddEntry(fHistDecayParErrorVsCent, "Decay parameter", "l");
  legend->AddEntry(fHistTotalErrorVsCent, "Total", "l");
  legend->Draw("same");
  canvasError->SaveAs(Form("../AbsoluteUncertaintySummary_%s.png", ParticleName[ChosenPart].Data()));
  canvasError->SaveAs(Form("../AbsoluteUncertaintySummary_%s.pdf", ParticleName[ChosenPart].Data()));

  // RelErrors
  TH1F *fHistBDTRelErrorVsCent = (TH1F *)fHistBDTErrorVsCent->Clone("fHistBDTRelErrorVsCent");
  TH1F *fHistMassCutRelErrorVsCent = (TH1F *)fHistMassCutErrorVsCent->Clone("fHistMassCutRelErrorVsCent");
  TH1F *fHistMassCutAndBDTRelErrorVsCent = (TH1F *)fHistMassCutAndBDTErrorVsCent->Clone("fHistMassCutAndBDTRelErrorVsCent");
  TH1F *fHistPrimaryLambdaRelErrorVsCent = (TH1F *)fHistPrimaryLambdaErrorVsCent->Clone("fHistPrimaryLambdaRelErrorVsCent");
  TH1F *fHistAccRelErrorVsCent = (TH1F *)fHistAccErrorVsCent->Clone("fHistAccRelErrorVsCent");
  TH1F *fHistResoRelErrorVsCent = (TH1F *)fHistResoErrorVsCent->Clone("fHistResoRelErrorVsCent");
  TH1F *fHistDecayParRelErrorVsCent = (TH1F *)fHistDecayParErrorVsCent->Clone("fHistDecayParRelErrorVsCent");
  TH1F *fHistTransferCoeffRelErrorVsCent = (TH1F *)fHistTransferCoeffErrorVsCent->Clone("fHistTransferCoeffRelErrorVsCent");
  TH1F *fHistRunByRunAccRelErrorVsCent = (TH1F *)fHistRunByRunAccErrorVsCent->Clone("fHistRunByRunAccRelErrorVsCent");
  TH1F *fHistPolBkg0RelErrorVsCent = (TH1F *)fHistPolBkg0ErrorVsCent->Clone("fHistPolBkg0RelErrorVsCent");
  TH1F *fHistPzFitRangeRelErrorVsCent = (TH1F *)fHistPzFitRangeErrorVsCent->Clone("fHistPzFitRangeRelErrorVsCent");
  TH1F *fHistBkgExpoRelErrorVsCent = (TH1F *)fHistBkgExpoErrorVsCent->Clone("fHistBkgExpoRelErrorVsCent");
  TH1F *fHistZVertexRelErrorVsCent = (TH1F *)fHistZVertexErrorVsCent->Clone("fHistZVertexRelErrorVsCent");
  TH1F *fHistTotalRelErrorVsCent = (TH1F *)fHistTotalErrorVsCent->Clone("fHistTotalRelErrorVsCent");
  for (Int_t m = 0; m < commonNumCent; m++)
  {
    fHistBDTRelErrorVsCent->SetBinContent(m + 1, fHistBDTErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistMassCutRelErrorVsCent->SetBinContent(m + 1, fHistMassCutErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistMassCutAndBDTRelErrorVsCent->SetBinContent(m + 1, fHistMassCutAndBDTErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistPrimaryLambdaRelErrorVsCent->SetBinContent(m + 1, fHistPrimaryLambdaErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistAccRelErrorVsCent->SetBinContent(m + 1, fHistAccErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistResoRelErrorVsCent->SetBinContent(m + 1, fHistResoErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistDecayParRelErrorVsCent->SetBinContent(m + 1, fHistDecayParErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistTransferCoeffRelErrorVsCent->SetBinContent(m + 1, fHistTransferCoeffErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistRunByRunAccRelErrorVsCent->SetBinContent(m + 1, fHistRunByRunAccErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistPolBkg0RelErrorVsCent->SetBinContent(m + 1, fHistPolBkg0ErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistPzFitRangeRelErrorVsCent->SetBinContent(m + 1, fHistPzFitRangeErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistBkgExpoRelErrorVsCent->SetBinContent(m + 1, fHistBkgExpoErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistZVertexRelErrorVsCent->SetBinContent(m + 1, fHistZVertexErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistTotalRelErrorVsCent->SetBinContent(m + 1, fHistTotalErrorVsCent->GetBinContent(m + 1) / TMath::Abs(fHistPzs2[m]->GetBinContent(1)));
    fHistBDTRelErrorVsCent->SetBinError(m + 1, 0);
    fHistMassCutRelErrorVsCent->SetBinError(m + 1, 0);
    fHistMassCutAndBDTRelErrorVsCent->SetBinError(m + 1, 0);
    fHistPrimaryLambdaRelErrorVsCent->SetBinError(m + 1, 0);
    fHistAccRelErrorVsCent->SetBinError(m + 1, 0);
    fHistResoRelErrorVsCent->SetBinError(m + 1, 0);
    fHistDecayParRelErrorVsCent->SetBinError(m + 1, 0);
    fHistTransferCoeffRelErrorVsCent->SetBinError(m + 1, 0);
    fHistRunByRunAccRelErrorVsCent->SetBinError(m + 1, 0);
    fHistPolBkg0RelErrorVsCent->SetBinError(m + 1, 0);
    fHistPzFitRangeRelErrorVsCent->SetBinError(m + 1, 0);
    fHistBkgExpoRelErrorVsCent->SetBinError(m + 1, 0);
    fHistZVertexRelErrorVsCent->SetBinError(m + 1, 0);
    fHistTotalRelErrorVsCent->SetBinError(m + 1, 0);
  }

  TCanvas *canvasRelError = new TCanvas("canvasRelError", "canvasRelError", 800, 600);
  StyleCanvas(canvasRelError, 0.05, 0.2, 0.2, 0.03);
  canvasRelError->cd();
  TH1F *hDummyRelError = new TH1F("hDummyRelError", "hDummyRelError", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummyRelError->GetNbinsX(); i++)
    hDummyRelError->SetBinContent(i, 1e-12);
  canvasRelError->cd();
  SetFont(hDummyRelError);
  StyleHistoYield(hDummyRelError, 0, 0.2, 1, 1, TitleXCent, "Relative syst. uncertainty", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummyRelError, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, 2, yLabelOffset);
  SetTickLength(hDummyRelError, tickX, tickY);
  hDummyRelError->GetXaxis()->SetRangeUser(0, 80);
  if (ChosenPart == 6)
    hDummyRelError->GetXaxis()->SetRangeUser(0, 50);
  hDummyRelError->Draw("");
  fHistBDTRelErrorVsCent->Draw("same");
  fHistMassCutRelErrorVsCent->Draw("same");
  fHistAccRelErrorVsCent->Draw("same");
  fHistMassCutAndBDTRelErrorVsCent->Draw("same"); 
  fHistPrimaryLambdaRelErrorVsCent->Draw("same");
  fHistDecayParRelErrorVsCent->Draw("same");
  if (ChosenPart != 6)
  {
    fHistTransferCoeffRelErrorVsCent->Draw("same");
    fHistRunByRunAccRelErrorVsCent->Draw("same");
  }
  if (ChosenPart == 6)
  {    
    fHistResoRelErrorVsCent->Draw("same");
    fHistPolBkg0RelErrorVsCent->Draw("same");
    fHistPzFitRangeRelErrorVsCent->Draw("same");
    fHistBkgExpoRelErrorVsCent->Draw("same");
    fHistZVertexRelErrorVsCent->Draw("same");
  }
  fHistTotalRelErrorVsCent->Draw("same");
  //legend->Draw("same");
  canvasRelError->SaveAs(Form("../RelativeUncertaintySummary_%s.png", ParticleName[ChosenPart].Data()));
  canvasRelError->SaveAs(Form("../RelativeUncertaintySummary_%s.pdf", ParticleName[ChosenPart].Data()));

  TFile *fileout = new TFile(stringout, "RECREATE");
  fHistMassCutAndBDTErrorVsCent->Write();
  fHistDecayParErrorVsCent->Write();
  if (ChosenPart != 6)
  {
    fHistBDTErrorVsCent->Write();
    fHistMassCutErrorVsCent->Write();
    fHistTransferCoeffErrorVsCent->Write();
    fHistRunByRunAccErrorVsCent->Write();
  }
  if (ChosenPart == 6)
  {
    fHistPrimaryLambdaErrorVsCent->Write();
    fHistResoErrorVsCent->Write();
    // fHistPolBkg0ErrorVsCent->Write();
    fHistZVertexErrorVsCent->Write();
  }
  fHistTotalErrorVsCent->Write();
  fHistMeanVsCent->Write();
  fileout->Close();

  cout << "\nStarting from the files (for the different mult): " << PathInBDT << endl;
  cout << "\nI have created the file:\n " << stringout << endl;
}