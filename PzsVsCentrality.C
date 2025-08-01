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
#include "TGraphErrors.h"
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

Float_t YLow[numPart] = {-0.001};
Float_t YUp[numPart] = {0.011};

void PzsVsCentrality(Int_t ChosenPart = ChosenParticle,
                     Bool_t isPolFromLambda = 0,
                     Bool_t isFromFit = 0,
                     Bool_t isRapiditySel = ExtrisRapiditySel,
                     Int_t BkgType = ExtrBkgType,
                     Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  if (isReducedPtBins && numPtBins != numPtBinsReduced)
  {
    cout << "Reduced pt bins are selected, but numPtBins is not set to numPtBinsReduced. Please check the settings." << endl;
    return;
  }
  Int_t ChosenPt = -999;
  cout << "Type 100 if you want to analyse Pz (integrated in pT), or the number of the pt interval you want" << endl;
  cin >> ChosenPt;

  if (ChosenPart == 6 && !isFromFit)
  {
    cout << "You have chosen the #Lambda particle and are not using the fit. Please select a different option." << endl;
    return;
  }

  Int_t part = 0;
  if (ChosenPart == 1 || ChosenPart == 4 || ChosenPart == 5)
  {
    part = 1;
  }
  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn;
  TFile *fileIn[numCent + 1];

  // fileinLambda
  TString PathInLambda = "Run2Results/HEPData-ins1891389-v1-P_z_vsCent.root";
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
  stringout = "Pzs2VsCentrality/" + NameAnalysis[!isV2] + "_";
  stringout += SinputFileName;
  stringout += "_" + ParticleName[ChosenPart];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[BkgType];
  stringout += "_Pzs2";
  if (isApplyWeights)
    stringout += "_Weighted";
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
  stringoutpdf = stringout;
  stringout += ".root";

  // canvases
  gStyle->SetOptStat(0);
  TCanvas *canvasPzs = new TCanvas("canvasPzs", "canvasPzs", 900, 700);
  StyleCanvas(canvasPzs, 0.05, 0.15, 0.15, 0.05);
  TH1F *fHistPzs = new TH1F("fHistPzs", "fHistPzs", numCent, fCentFT0C);
  TH1F *fHistPzsSist = new TH1F("fHistPzsSist", "fHistPzsSist", numCent, fCentFT0C);
  TH1F *fHistPzsError = new TH1F("fHistPzsError", "fHistPzsError", numCent, fCentFT0C);
  TH1F *fHistPuritySummary = new TH1F("fHistPuritySummary", "fHistPuritySummary", numCent, fCentFT0C);
  TH1F *fHistSignificanceSummary = new TH1F("fHistSignificanceSummary", "fHistSignificanceSummary", numCent, fCentFT0C);
  TH1F *fHistYieldSummary = new TH1F("fHistYieldSummary", "fHistYieldSummary", numCent, fCentFT0C);
  TH1F *fHistMeanSummary = new TH1F("fHistMeanSummary", "fHistMeanSummary", numCent, fCentFT0C);
  TH1F *fHistSigmaSummary = new TH1F("fHistSigmaSummary", "fHistSigmaSummary", numCent, fCentFT0C);
  TH1F *fHistMeanMinus2Sigma = new TH1F("fHistMeanMinus2Sigma", "fHistMeanMinus2Sigma", numCent, fCentFT0C);
  TH1F *fHistMeanPlus2Sigma = new TH1F("fHistMeanPlus2Sigma", "fHistMeanPlus2Sigma", numCent, fCentFT0C);
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
  LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
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
      LegendTitle->AddEntry("", ParticleNameLegend[ChosenPart] + " |#it{#eta}| < 0.8", "");
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
  TH1F *fHistSpectrum[numCent + 1];
  TH1F *fHistPurity[numCent + 1];
  TH1F *fHistSignificance[numCent + 1];
  TH1F *fHistYield[numCent + 1];
  TH1F *fHistMean[numCent + 1];
  TH1F *fHistSigma[numCent + 1];
  TString Smolt[numCent + 1];
  TString SmoltBis[numCent + 1];
  // get spectra in multiplicity classes
  for (Int_t m = 0; m < numCent; m++)
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

    PathIn = "OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_";
    PathIn += SinputFileName;
    PathIn += "_" + ParticleName[ChosenPart];
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[BkgType];
    Smolt[m] += Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);
    SmoltBis[m] += Form("%i#minus%i", CentFT0CMin, CentFT0CMax);
    PathIn += Smolt[m];
    if (isApplyWeights)
      PathIn += "_Weighted";
    if (!useCommonBDTValue)
      PathIn += "_BDTCentDep";
    if (isRun2Binning)
      PathIn += "_Run2Binning";
    if (isPolFromLambda)
      PathIn += "_PolFromLambda";
    if (!isRapiditySel)
      PathIn += "_Eta08";
    PathIn += STHN[ExtrisFromTHN];
    if (useMixedBDTValueInFitMacro)
      PathIn += "_MixedBDT";
    if (isTightMassCut)
      PathIn += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
    if (isReducedPtBins)
      PathIn += "_ReducedPtBins";
    PathIn += ".root";
    cout << "Path in : " << PathIn << endl;
    fileIn[m] = TFile::Open(PathIn);
    if (ChosenPt == 100)
    {
      fHistSpectrum[m] = (TH1F *)fileIn[m]->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + "PtInt" + V2FromFit[isFromFit]);
      fHistPurity[m] = (TH1F *)fileIn[m]->Get("histoPurityPtInt");
      fHistSignificance[m] = (TH1F *)fileIn[m]->Get("histoSignificancePtInt");
      fHistYield[m] = (TH1F *)fileIn[m]->Get("histoYieldPtInt");
      fHistMean[m] = (TH1F *)fileIn[m]->Get("histoMeanPtInt");
      // fHistSigma[m] = (TH1F *)fileIn[m]->Get("histoSigmaPtIntWeighted");
      fHistSigma[m] = (TH1F *)fileIn[m]->Get("histoSigmaPtInt");
    }
    else
    {
      fHistSpectrum[m] = (TH1F *)fileIn[m]->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + V2FromFit[isFromFit]);
      fHistPurity[m] = (TH1F *)fileIn[m]->Get("histoPurity");
      fHistSignificance[m] = (TH1F *)fileIn[m]->Get("histoSignificance");
      fHistYield[m] = (TH1F *)fileIn[m]->Get("histoYield");
      fHistMean[m] = (TH1F *)fileIn[m]->Get("histoMean");
      // fHistSigma[m] = (TH1F *)fileIn[m]->Get("histoSigmaWeighted");
      fHistSigma[m] = (TH1F *)fileIn[m]->Get("histoSigma");
    }
    if (!fHistSpectrum[m])
    {
      cout << " no hist v2 / Pzs" << endl;
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
    fHistSpectrum[m]->SetName("histoPzs2_" + Smolt[m]);
    fHistPurity[m]->SetName("histoPurity_" + Smolt[m]);
    fHistSignificance[m]->SetName("histoSignificance_" + Smolt[m]);
    fHistYield[m]->SetName("histoYield_" + Smolt[m]);
    fHistMean[m]->SetName("histoMean_" + Smolt[m]);
    fHistSigma[m]->SetName("histoSigma_" + Smolt[m]);

    fHistPzs->SetBinContent(m + 1, fHistSpectrum[m]->GetBinContent(1) / fHistPurity[m]->GetBinContent(1));
    fHistPzs->SetBinError(m + 1, fHistSpectrum[m]->GetBinError(1) / fHistPurity[m]->GetBinContent(1));

    fHistPzsError->SetBinContent(m + 1, fHistSpectrum[m]->GetBinError(1) / fHistPurity[m]->GetBinContent(1));
    fHistPzsError->SetBinError(m + 1, 0);

    fHistPuritySummary->SetBinContent(m + 1, fHistPurity[m]->GetBinContent(1));
    fHistPuritySummary->SetBinError(m + 1, fHistPurity[m]->GetBinError(1));

    fHistSignificanceSummary->SetBinContent(m + 1, fHistSignificance[m]->GetBinContent(1));
    fHistSignificanceSummary->SetBinError(m + 1, fHistSignificance[m]->GetBinError(1));
    cout << fHistSignificance[m]->GetBinContent(1) << endl;

    fHistYieldSummary->SetBinContent(m + 1, fHistYield[m]->GetBinContent(1));
    fHistYieldSummary->SetBinError(m + 1, fHistYield[m]->GetBinError(1));

    fHistMeanSummary->SetBinContent(m + 1, fHistMean[m]->GetBinContent(1));
    fHistMeanSummary->SetBinError(m + 1, fHistMean[m]->GetBinError(1));

    fHistSigmaSummary->SetBinContent(m + 1, fHistSigma[m]->GetBinContent(1));
    fHistSigmaSummary->SetBinError(m + 1, fHistSigma[m]->GetBinError(1));

    fHistMeanMinus2Sigma->SetBinContent(m + 1, fHistMean[m]->GetBinContent(1) - 2 * fHistSigma[m]->GetBinContent(1));
    fHistMeanMinus2Sigma->SetBinError(m + 1, 0);

    fHistMeanPlus2Sigma->SetBinContent(m + 1, fHistMean[m]->GetBinContent(1) + 2 * fHistSigma[m]->GetBinContent(1));
    fHistMeanPlus2Sigma->SetBinError(m + 1, 0);

    cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << " Pzs2: " << fHistSpectrum[m]->GetBinContent(1) << endl;
  } // end loop on mult

  // Get histogram with syst uncertainty
  TString PathInSyst = "Systematics/SystVsCentrality_" + NameAnalysis[!isV2] + "_";
  PathInSyst += SinputFileNameSyst;
  PathInSyst += "_" + ParticleName[ChosenPart];
  PathInSyst += IsOneOrTwoGauss[UseTwoGauss];
  PathInSyst += SIsBkgParab[BkgType];
  PathInSyst += "_Pzs2";
  if (isApplyWeights)
    PathInSyst += "_Weighted";
  if (!useCommonBDTValue)
    PathInSyst += "_BDTCentDep";
  if (isRun2Binning)
    PathInSyst += "_Run2Binning";
  // if (isPolFromLambda)
  PathInSyst += "_PolFromLambda";
  if (!isRapiditySel)
    PathInSyst += "_Eta08";
  PathInSyst += STHN[ExtrisFromTHN];
  if (useMixedBDTValueInFitMacro)
    PathInSyst += "_MixedBDT";
  if (isTightMassCut)
    PathInSyst += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
  PathInSyst += V2FromFit[isFromFit];
  PathInSyst += ".root";
  if (ChosenPart == 6)
    PathInSyst = "Systematics/SystVsCentrality_Pzs2_LHC23_PbPb_pass4_Train370610_ProtonAcc_Xi_BkgParab_Pzs2_PolFromLambda_Eta08_FromTHN_MixedBDT_TightMassCut2.1NoFit.root";
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

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 8000, 0, 80);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, -1000);
  canvasPzs->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleXCent, TitleYPzs, "", 1, 1.15, 1.6);
  StyleHistoYield(fHistPzs, YLow[part], YUp[part], ColorPart[part], MarkerPart[part], TitleXCent, TitleYPzs, "", MarkerPartSize[part], 1.15, 1.6);
  StyleHistoYield(fHistPzsSist, YLow[part], YUp[part], ColorPart[part], MarkerPart[part], TitleXCent, TitleYPzs, "", 1.5, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(0, 80);
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
  StyleCanvas(canvasPzsError, 0.05, 0.15, 0.15, 0.05);
  TH1F *hDummyError = new TH1F("hDummyError", "hDummyError", 8000, 0, 80);
  for (Int_t i = 1; i <= hDummyError->GetNbinsX(); i++)
    hDummyError->SetBinContent(i, 1e-12);
  canvasPzsError->cd();
  SetFont(hDummyError);
  StyleHistoYield(hDummyError, 0, 0.01, 1, 1, TitleXCent, "Absolute stat. uncertainty", "", 1, 1.15, 1.6);
  StyleHistoYield(fHistPzsError, 0, 0.01, ColorPart[part], MarkerPart[part], TitleXCent, "Absolute stat. uncertainty", "", MarkerPartSize[part], 1.15, 1.6);
  StyleHistoYield(fHistPzsSistError, 0, 0.01, kGray + 2, 22, TitleXCent, "Absolute syst. uncertainty", "", MarkerPartSize[part], 1.15, 1.6);
  SetHistoTextSize(hDummyError, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummyError, tickX, tickY);
  hDummyError->GetXaxis()->SetRangeUser(0, 80);
  hDummyError->Draw("");
  fHistPzsError->Draw("same");
  fHistPzsSistError->Draw("same e2");
  fHistPzsLambda_StatErr->Draw("same e0x0");
  LegendTitle->Draw("");
  legendLambda->Draw("");
  canvasPzsError->SaveAs(stringoutpdf + "_Error.pdf");
  canvasPzsError->SaveAs(stringoutpdf + "_Error.png");

  // Significativity
  TCanvas *canvasPzsSignif = new TCanvas("canvasPzsSignif", "canvasPzsSignif", 900, 700);
  StyleCanvas(canvasPzsSignif, 0.05, 0.15, 0.15, 0.05);
  TH1F *hDummySignif = new TH1F("hDummySignif", "hDummySignif", 8000, 0, 80);
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
  hDummySignif->GetXaxis()->SetRangeUser(0, 80);
  hDummySignif->GetYaxis()->SetRangeUser(0, 7);
  hDummySignif->Draw("");
  fHistPzsSignif->Draw("same");
  fHistPzsSignifStat->Draw("same");
  fHistPzsSignifLambda->Draw("same e0x0");
  LegendTitle->Draw("");
  canvasPzsSignif->SaveAs(stringoutpdf + "_Signif.pdf");
  canvasPzsSignif->SaveAs(stringoutpdf + "_Signif.png");

  TCanvas *canvasPurity = new TCanvas("canvasPurity", "canvasPurity", 900, 700);
  StyleCanvas(canvasPurity, 0.05, 0.15, 0.15, 0.05);
  canvasPurity->cd();
  TH1F *hDummyPurity = (TH1F *)hDummy->Clone("hDummyPurity");
  hDummyPurity->GetYaxis()->SetRangeUser(0.9, 1);
  hDummyPurity->GetYaxis()->SetTitle("S / (S+B)");
  hDummyPurity->Draw("");
  StyleHistoYield(fHistPuritySummary, 0.9, 1, ColorPart[part], MarkerPart[part], TitleXCent, "S / (S+B)", "", MarkerPartSize[part], 1.15, 1.6);
  fHistPuritySummary->Draw("same");
  canvasPurity->SaveAs(stringoutpdf + "_Purity.pdf");
  canvasPurity->SaveAs(stringoutpdf + "_Purity.png");

  TCanvas *canvasSignificance = new TCanvas("canvasSignificance", "canvasSignificance", 900, 700);
  StyleCanvas(canvasSignificance, 0.05, 0.15, 0.15, 0.05);
  canvasSignificance->cd();
  TH1F *hDummySignificance = (TH1F *)hDummy->Clone("hDummySignificance");
  hDummySignificance->GetYaxis()->SetTitle("S / #sqrt{S+B}");
  hDummySignificance->Draw("");
  StyleHistoYield(fHistSignificanceSummary, 0, 1.2 * fHistSignificanceSummary->GetBinContent(fHistSignificanceSummary->GetMaximum()), ColorPart[part], MarkerPart[part], TitleXCent, "S / #sqrt{S+B}", "", MarkerPartSize[part], 1.15, 1.6);
  for (Int_t b = 1; b <= fHistSignificanceSummary->GetNbinsX(); b++)
  {
    cout << fHistSignificanceSummary->GetBinContent(b) << endl;
  }
  fHistSignificanceSummary->Draw("same");
  canvasSignificance->SaveAs(stringoutpdf + "_Significance.pdf");
  canvasSignificance->SaveAs(stringoutpdf + "_Significance.png");

  TCanvas *canvasYield = new TCanvas("canvasYield", "canvasYield", 900, 700);
  StyleCanvas(canvasYield, 0.05, 0.15, 0.15, 0.05);
  canvasYield->cd();
  TH1F *hDummyYield = (TH1F *)hDummy->Clone("hDummyYield");
  hDummyYield->GetYaxis()->SetRangeUser(0, 0.015);
  hDummyYield->GetYaxis()->SetTitle("Yield");
  hDummyYield->Draw("");
  StyleHistoYield(fHistYieldSummary, 0, 0.015, ColorPart[part], MarkerPart[part], TitleXCent, "Yield", "", MarkerPartSize[part], 1.15, 1.6);
  fHistYieldSummary->Draw("same");

  TCanvas *canvasMeanSigma = new TCanvas("canvasMeanSigma", "canvasMeanSigma", 900, 700);
  StyleCanvas(canvasMeanSigma, 0.05, 0.15, 0.15, 0.05);
  canvasMeanSigma->cd();
  TH1F *hDummySigma = (TH1F *)hDummy->Clone("hDummySigma");
  hDummySigma->GetYaxis()->SetRangeUser(1.31, 1.33);
  hDummySigma->GetYaxis()->SetTitle("#mu");
  hDummySigma->Draw("");
  StyleHistoYield(fHistMeanSummary, 1.31, 1.33, ColorPart[part], MarkerPart[part], TitleXCent, "#mu", "", MarkerPartSize[part], 1.15, 1.6);
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
  TF1 *fpol1 = new TF1("fpol1", "pol1", 0, 80);
  TF1 *fpol0 = new TF1("fpol0", "pol0", 0, 80);
  fpol1->SetLineColor(kBlue + 2);
  fpol0->SetLineColor(kAzure + 1);
  fHistPzsTotError->Fit("fpol1", "R+");
  fHistPzsTotError->Fit("fpol0", "R+");
  LegendTitle->Draw("");
  TLegend *legendfit = new TLegend(0.2, 0.58, 0.5, 0.73);
  legendfit->SetFillStyle(0);
  legendfit->SetTextSize(0.04);
  legendfit->AddEntry(fpol0, Form("pol0, Chi2/NDF = %.2f/%i", fpol0->GetChisquare(), fpol0->GetNDF()), "l");
  legendfit->AddEntry(fpol1, Form("pol1, Chi2/NDF = %.2f/%i", fpol1->GetChisquare(), fpol1->GetNDF()), "l");
  legendfit->Draw("");
  canvasfitPol0->SaveAs(stringoutpdf + "_fitPol0.pdf");
  canvasfitPol0->SaveAs(stringoutpdf + "_fitPol0.png");

  TLegend *LegendPreliminary;
  LegendPreliminary = new TLegend(0.12, 0.70, 0.51, 0.94);
  LegendPreliminary->SetFillStyle(0);
  LegendPreliminary->SetTextAlign(11);
  LegendPreliminary->SetTextSize(0.04);
  LegendPreliminary->AddEntry("", "#bf{ALICE Preliminary}", "");
  LegendPreliminary->AddEntry("", "Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
  LegendPreliminary->AddEntry("", "#Xi^{#minus} + #bar{#Xi}^{+}, |#it{#eta} | < 0.5", "");
  LegendPreliminary->AddEntry("", Form("%1.1f < #it{p}_{T} < %1.1f GeV/#it{c}", MinPt[ChosenPart], 8.), "");

  TLegend *LegendPreliminary2;
  LegendPreliminary2 = new TLegend(0.06, 0.777, 0.45, 0.917);
  LegendPreliminary2->SetFillStyle(0);
  LegendPreliminary2->SetTextAlign(11);
  LegendPreliminary2->SetTextSize(0.048);
  LegendPreliminary2->AddEntry("", "#bf{ALICE Preliminary}", "");
  LegendPreliminary2->AddEntry("", "Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");

  TLegend *legendXi = new TLegend(0.06, 0.643, 0.45, 0.784);
  legendXi->SetFillStyle(0);
  legendXi->SetTextAlign(12);
  legendXi->SetTextSize(0.048);
  if (ChosenPart == 6)
    legendXi->AddEntry("", Form("#Lambda + #bar{#Lambda}, |#it{y} | < 0.5, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "");
  else
    legendXi->AddEntry("", Form("#Xi^{#minus} + #bar{#Xi}^{+}, |#it{#eta} | < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "");

  TF1 *lineatZero = new TF1("lineatZero", "0", 0, 80);
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
  hDummy->GetXaxis()->SetRangeUser(0, 80);
  hDummy->Draw("");
  lineatZero->Draw("same");
  fHistPzs->DrawClone("same ex0");
  fHistPzsSist->SetFillStyle(0);
  fHistPzsSist->DrawClone("same e2");
  // LegendPreliminary->Draw("");
  LegendPreliminary2->Draw("");
  legendXi->Draw("");
  canvasPzsXi->SaveAs("XiPolVsCent.pdf");
  canvasPzsXi->SaveAs("XiPolVsCent.png");
  canvasPzsXi->SaveAs("XiPolVsCent.eps");

  TString SfileLambdaJunlee = "LambdaJunlee/fout_psi2_mult.root";
  TFile *fileLambdaJunlee = new TFile(SfileLambdaJunlee);
  TGraphErrors *gPzsLambdaJunlee = (TGraphErrors *)fileLambdaJunlee->Get("gMultStat");
  TGraphErrors *gPzsLambdaJunleeSist = (TGraphErrors *)fileLambdaJunlee->Get("gMultSyst");
  TH1F *fHistPzsLambdaJunlee = (TH1F *)fHistPzs->Clone("fHistPzsLambdaJunlee");
  TH1F *fHistPzsLambdaJunleeSist = (TH1F *)fHistPzsSist->Clone("fHistPzsLambdaJunleeSist");
  fHistPzsLambdaJunlee->Reset();
  fHistPzsLambdaJunleeSist->Reset();
  for (Int_t b = 1; b <= fHistPzs->GetNbinsX(); b++)
  {
    fHistPzsLambdaJunlee->SetBinContent(b, gPzsLambdaJunlee->GetY()[b - 1]);
    fHistPzsLambdaJunlee->SetBinError(b, gPzsLambdaJunlee->GetEY()[b - 1]);
    fHistPzsLambdaJunleeSist->SetBinContent(b, gPzsLambdaJunleeSist->GetY()[b - 1]);
    fHistPzsLambdaJunleeSist->SetBinError(b, gPzsLambdaJunleeSist->GetEY()[b - 1]);
  }

  // Int_t colorJunlee = kGreen + 2;
  Int_t colorJunlee = kAzure - 3;
  StyleHistoYield(fHistPzsLambdaJunlee, YLow[part], YUp[part], colorJunlee, 47, TitleXCent, TitleYPzs, "", 2.1, 1.15, 1.8);
  StyleHistoYield(fHistPzsLambdaJunleeSist, YLow[part], YUp[part], colorJunlee, 47, TitleXCent, TitleYPzs, "", 2.1, 1.15, 1.8);

  // TLegend *legendParticles = new TLegend(0.20, 0.63, 0.60, 0.77);
  TLegend *legendParticles = new TLegend(0.136, 0.61, 0.526, 0.75);
  legendParticles->SetFillStyle(0);
  legendParticles->SetTextAlign(12);
  legendParticles->SetTextSize(0.048);
  if (ChosenParticle == 6)
    legendParticles->AddEntry(fHistPzs, Form("#Lambda + #bar{#Lambda}, |#it{y} | < 0.5, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "pl");
  else
    legendParticles->AddEntry(fHistPzs, Form("#Xi^{#minus} + #bar{#Xi}^{+}, |#it{#eta} | < 0.8, #it{p}_{T} > %1.1f GeV/#it{c}", MinPt[ChosenPart]), "pl");
  legendParticles->AddEntry(fHistPzsLambdaJunlee, Form("#Lambda + #bar{#Lambda}, |#it{y} | < 0.5, #it{p}_{T} > %1.1f GeV/#it{c}", 0.5), "pl");
  TCanvas *canvasPzsXiLambda = new TCanvas("canvasPzsXiLambda", "canvasPzsXiLambda", 900, 700);
  StyleCanvas(canvasPzsXiLambda, 0.06, 0.12, 0.1, 0.03);
  canvasPzsXiLambda->cd();
  hDummy->Draw("");
  lineatZero->Draw("same");
  fHistPzsLambdaJunlee->Draw("same ex0");
  fHistPzsLambdaJunleeSist->SetFillStyle(0);
  fHistPzsLambdaJunleeSist->Draw("same e2");
  fHistPzs->Draw("same ex0");
  fHistPzsSist->SetFillStyle(0);
  fHistPzsSist->Draw("same e2");
  // gPzsLambdaJunlee->Draw("same p");
  // gPzsLambdaJunleeSist->Draw("same e2");
  LegendPreliminary2->Draw("");
  legendParticles->Draw("");
  canvasPzsXiLambda->SaveAs("XiLambdaPolVsCent.pdf");
  canvasPzsXiLambda->SaveAs("XiLambdaPolVsCent.png");
  canvasPzsXiLambda->SaveAs("XiLambdaPolVsCent.eps");

  TFile *fileout = new TFile(stringout, "RECREATE");
  fHistMeanSummary->Write();
  fHistSigmaSummary->Write();
  fHistPzs->Write();
  fHistPzsSist->Write();
  fHistPuritySummary->Write();
  fHistYieldSummary->Write();
  fHistSignificanceSummary->Write();
  fHistPzsError->Write();

  fileout->Close();

  cout << "\nStarting from the files (for the different mult): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;

  if (ChosenPart == 6)
    cout << "\n\nWARNING: Syst path for Lambda might be dummy!! " << endl;
}
