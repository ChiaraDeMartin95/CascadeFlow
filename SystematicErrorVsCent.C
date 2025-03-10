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

Float_t YLow[numPart] = {0};
Float_t YUp[numPart] = {0.004};

void SystematicErrorVsCent(Int_t ChosenPart = ChosenParticle,
                           Bool_t isPolFromLambda = 0,
                           Float_t nsigmaBarlowMassCut = 1.0,
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
  TFile *fileIn[numCent + 1];

  TString PathIn1 = "";
  TString PathIn2 = "";
  // fileInBDT
  TString PathInBDT;
  TFile *fileInBDT[numCent + 1];
  // fileIn MassCut
  TString PathInMassCut;
  TFile *fileInMassCut[numCent + 1];

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = "Systematics/SystVsCentrality_" + NameAnalysis[!isV2] + "_";
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
  if (!isRapiditySel || ExtrisFromTHN)
    stringout += "_Eta08";
  stringout += STHN[ExtrisFromTHN];
  if (useMixedBDTValueInFitMacro)
    stringout += "_MixedBDT";
  if (isTightMassCut)
    stringout += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
  stringoutpdf = stringout;
  stringout += ".root";

  // canvases
  gStyle->SetOptStat(0);
  TCanvas *canvasError = new TCanvas("canvasError", "canvasError", 800, 600);
  StyleCanvas(canvasError, 0.05, 0.2, 0.2, 0.03);

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  TH1F *fHistPzs2[numCent + 1];
  TH1F *fHistBDTError[numCent + 1];
  TH1F *fHistMassCutError[numCent + 1];
  TH1F *fHistBDTErrorVsCent = new TH1F("fHistBDTErrorVsCent", "fHistBDTErrorVsCent", numCent, fCentFT0C);
  TH1F *fHistMassCutErrorVsCent = new TH1F("fHistMassCutErrorVsCent", "fHistMassCutErrorVsCent", numCent, fCentFT0C);
  TH1F *fHistPrimaryLambdaErrorVsCent = new TH1F("fHistPrimaryLambdaErrorVsCent", "fHistPrimaryLambdaErrorVsCent", numCent, fCentFT0C);
  TH1F *fHistTotalErrorVsCent = new TH1F("fHistTotalErrorVsCent", "fHistTotalErrorVsCent", numCent, fCentFT0C);
  TString Smolt[numCent + 1];
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
    Smolt[m] += Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);

    PathIn = "OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_";
    PathIn += SinputFileName;
    PathIn += "_" + ParticleName[ChosenPart];
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[BkgType];
    PathIn += Smolt[m];
    if (isApplyWeights)
      PathIn += "_Weighted";
    if (!useCommonBDTValue)
      PathIn += "_BDTCentDep";
    if (isRun2Binning)
      PathIn += "_Run2Binning";
    if (isPolFromLambda)
      PathIn += "_PolFromLambda";
    if (!isRapiditySel || ExtrisFromTHN)
      PathIn += "_Eta08";
    PathIn += STHN[ExtrisFromTHN];
    if (useMixedBDTValueInFitMacro)
      PathIn += "_MixedBDT";
    if (isTightMassCut)
      PathIn += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
    PathIn += ".root";
    cout << "Path in : " << PathIn << endl;
    fileIn[m] = TFile::Open(PathIn);
    fHistPzs2[m] = (TH1F *)fileIn[m]->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + "PtInt");
    if (!fHistPzs2[m])
    {
      cout << " no hist v2 / Pzs" << endl;
      return;
    }

    PathIn1 = "Systematics/SystMultiTrial_" + SinputFileName + Form("_%i-%i_", CentFT0CMin, CentFT0CMax) + ParticleName[ChosenPart] + "_";
    PathIn2 = "";
    if (isPtIntegrated)
      PathIn2 += "_PtInt";
    if (isApplyWeights)
      PathIn2 += "_Weighted";
    if (!isPtAnalysis)
      PathIn2 += "_vsPsi";
    if (!isV2 && isPolFromLambda)
      PathIn2 += "_PolFromLambda";
    if (ExtrisApplyEffWeights)
      PathIn2 += "_EffW";
    if (!isRapiditySel || ExtrisFromTHN)
      PathIn2 += "_Eta08";
    PathIn2 += STHN[ExtrisFromTHN];

    PathInBDT = PathIn1 + "BDT" + PathIn2 + ".root";
    cout << "Path in BDT: " << PathInBDT << endl;
    fileInBDT[m] = TFile::Open(PathInBDT);
    PathInMassCut = PathIn1 + "MassCut" + PathIn2 + Form("_nsigmaBarlow%.1f.root", nsigmaBarlowMassCut);
    cout << "Path in MassCut: " << PathInMassCut << endl;
    fileInMassCut[m] = TFile::Open(PathInMassCut);

    fHistBDTError[m] = (TH1F *)fileInBDT[m]->Get("hAbsoluteMaxDev");
    if (!fHistBDTError[m])
    {
      cout << "Error: histogram not found" << endl;
      return;
    }
    fHistBDTError[m]->SetName("hAbsoluteSystErrorBDT_" + Smolt[m]);

    fHistMassCutError[m] = (TH1F *)fileInMassCut[m]->Get("hAbsoluteMaxDev");
    if (!fHistMassCutError[m])
    {
      cout << "Error: histogram not found" << endl;
      return;
    }
    fHistMassCutError[m]->SetName("hAbsoluteSystErrorMassCut_" + Smolt[m]);

    fHistPrimaryLambdaErrorVsCent->SetBinContent(m + 1, 0.05 * fHistPzs2[m]->GetBinContent(1));
    fHistPrimaryLambdaErrorVsCent->SetBinError(m + 1, 0);
    fHistBDTErrorVsCent->SetBinContent(m + 1, fHistBDTError[m]->GetBinContent(1));
    fHistBDTErrorVsCent->SetBinError(m + 1, 0);
    fHistMassCutErrorVsCent->SetBinContent(m + 1, fHistMassCutError[m]->GetBinContent(1));
    fHistMassCutErrorVsCent->SetBinError(m + 1, 0);
  } // end loop on mult
  fHistBDTErrorVsCent->Smooth();
  fHistMassCutErrorVsCent->Smooth();

  for (Int_t m = 0; m < numCent; m++)
  {
    // fHistTotalErrorVsCent->SetBinContent(m + 1, TMath::Sqrt(fHistBDTErrorVsCent->GetBinContent(m + 1) * fHistBDTErrorVsCent->GetBinContent(m + 1) + fHistMassCutErrorVsCent->GetBinContent(m + 1) * fHistMassCutErrorVsCent->GetBinContent(m + 1)));
    fHistTotalErrorVsCent->SetBinContent(m + 1, TMath::Sqrt(fHistBDTErrorVsCent->GetBinContent(m + 1) * fHistBDTErrorVsCent->GetBinContent(m + 1) +
                                                            fHistMassCutErrorVsCent->GetBinContent(m + 1) * fHistMassCutErrorVsCent->GetBinContent(m + 1) +
                                                            fHistPrimaryLambdaErrorVsCent->GetBinContent(m + 1) * fHistPrimaryLambdaErrorVsCent->GetBinContent(m + 1)));
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

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 8000, 0, 80);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasError->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleXCent, "Absolute syst. uncertainty", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(0, 80);
  hDummy->Draw("");
  fHistBDTErrorVsCent->SetLineColor(kRed);
  fHistMassCutErrorVsCent->SetLineColor(kGreen + 2);
  fHistPrimaryLambdaErrorVsCent->SetLineColor(kOrange+2);
  fHistTotalErrorVsCent->SetLineColor(kBlack);
  fHistTotalErrorVsCent->SetLineWidth(2);
  fHistBDTErrorVsCent->Draw("same");
  fHistMassCutErrorVsCent->Draw("same");
  fHistPrimaryLambdaErrorVsCent->Draw("same");
  fHistTotalErrorVsCent->Draw("same");
  TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.05);
  legend->AddEntry(fHistPrimaryLambdaErrorVsCent, "Primary #Lambda", "l");
  legend->AddEntry(fHistBDTErrorVsCent, "BDT", "l");
  legend->AddEntry(fHistMassCutErrorVsCent, "Mass Cut", "l");
  legend->AddEntry(fHistTotalErrorVsCent, "Total", "l");
  legend->Draw("same");
  canvasError->SaveAs("AbsoluteUncertaintySummary.png");

  TFile *fileout = new TFile(stringout, "RECREATE");
  fHistBDTErrorVsCent->Write();
  fHistMassCutErrorVsCent->Write();
  fHistTotalErrorVsCent->Write();

  cout << "\nStarting from the files (for the different mult): " << PathInBDT << endl;
  cout << "\nI have created the file:\n " << stringout << endl;
}
