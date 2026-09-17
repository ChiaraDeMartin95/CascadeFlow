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
#include "CommonVarPub.h"
// #include "CommonVarXi.h"
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
Float_t ResoRelError[numCentLambdaOO + 1] = {0};
Float_t PrimaryLambdaFraction = 0;
Float_t SecondaryLambdaFraction = 0.008; // only considering feed-down from Xi decays, not from higher mass resonances decaying into Lambdas
Float_t ZVertexErrorLambdaOO = 0.00009;

void SystematicErrorVsPt(Bool_t isPtAnalysis = 1, // 1 for V2 vs pt and Pzs2 vs pt, 0 for Pz vs 2(phi-Psi)
                         Int_t ChosenPart = ChosenParticle,
                         Bool_t isPolFromLambda = 0,
                         Float_t nsigmaBarlowMassCut = 1.0,
                         Bool_t isFromFit = 1,
                         Bool_t isRapiditySel = ExtrisRapiditySel,
                         Int_t BkgType = ExtrBkgType,
                         Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  TString SVariableDep = isPtAnalysis ? "Pt" : "Pz";

  Float_t UpperRangeParticle = PtBins[numPtBins - 1];
  Float_t LowerRangeParticle = PtBins[0];
  if (!isPtAnalysis)
  {
    UpperRangeParticle = 2 * TMath::Pi();
    LowerRangeParticle = 0;
  }

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
  // fileIn TopoSel
  TString PathInMassCutAndBDT;
  TFile *fileInMassCutAndBDT[commonNumCent + 1];

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = "../Systematics/SystVsPt_" + NameAnalysis[!isV2] + "_";
  stringout += SinputFileNameSyst;
  stringout += "_" + ParticleName[ChosenPart];
  if (ExtrisFitDSCB)
  {
    stringout += "_DSCB";
    if (isFixParamDSCBFromMC)
      stringout += "_FixParamFromMC";
  }
  else
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
  if (ChosenPart == 0)
    stringout += "_New";
  stringoutpdf = stringout;
  stringout += ".root";

  // canvases
  gStyle->SetOptStat(0);
  TCanvas *canvasError = new TCanvas("canvasError", "canvasError", 800, 600);
  StyleCanvas(canvasError, 0.05, 0.2, 0.2, 0.03);

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;

  TH1F *fHistPzs2 = nullptr;
  TH1F *fHistMassCutAndBDTError = nullptr;
  TH1F *fHistMeanSystMultiTrialMean = nullptr;
  TH1F *fHistPrimaryLambdaError = nullptr;
  TH1F *fHistDecayParError = nullptr;
  TH1F *fHistZVertexError = nullptr;

  TH1F *fHistTotalError = nullptr;
  TH1F *fHistMean = nullptr;

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

  TFile *fileInBkgExpo = TFile::Open("../CompareResults/SystUncertainty_2GaussPlusPol2.root");
  TH1F *fHistBkgExpoError = (TH1F *)fileInBkgExpo->Get("hRatioClone_1");
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
  for (Int_t m = 0; m <= numCentLambdaOO; m++)
  {
    if (m != numCentLambdaOO)
      continue;

    CentFT0CMin = 0;
    CentFT0CMax = CentFT0CMaxLambdaOO;

    Smolt[m] += Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);

    PathIn = "../OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_";
    PathIn += SinputFileName;
    PathIn += "_" + ParticleName[ChosenPart];
    if (ExtrisFitDSCB)
    {
      PathIn += "_DSCB";
      if (isFixParamDSCBFromMC)
        PathIn += "_FixParamFromMC";
    }
    else
      PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[BkgType];
    PathIn += Smolt[m];
    if (isApplyWeights)
      PathIn += "_Weighted";
    if (isApplyCentWeight)
      PathIn += "_CentWeighted";
    if (!isPtAnalysis)
      PathIn += "_vsPsi";
    if (!useCommonBDTValue)
      PathIn += "_BDTCentDep";
    if (isRun2Binning)
      PathIn += "_Run2Binning";
    if (isPolFromLambda)
      PathIn += "_PolFromLambda";
    if (ExtrisApplyEffWeights)
      PathIn += "_EffW";
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
    if (isApplyAcceptanceInMacro)
      PathIn += "_AcceptanceInMacro";
    PathIn += "_ResoOnTheFly";

    PathIn += ".root";
    cout << "Path in : " << PathIn << endl;
    fileIn[m] = TFile::Open(PathIn);
    if (isPtAnalysis)
      fHistPzs2 = (TH1F *)fileIn[m]->Get("histoPzs2" + sPolFromLambda[isPolFromLambda] + V2FromFit[isFromFit]);
    else
      fHistPzs2 = (TH1F *)fileIn[m]->Get("histoPz" + sPolFromLambda[isPolFromLambda] + V2FromFit[isFromFit]);
    if (!fHistPzs2)
    {
      cout << " no hist v2 / Pzs" << endl;
      return;
    }

    PathIn1 = "../Systematics/SystMultiTrial_" + SinputFileNameSyst + Form("_%i-%i_", CentFT0CMin, CentFT0CMax) + ParticleName[ChosenPart] + "_";
    PathIn2 = "";
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

    PathInMassCutAndBDT = PathIn1;
    PathInMassCutAndBDT += "LambdaTopo";
    PathInMassCutAndBDT += PathIn2;
    PathInMassCutAndBDT += ".root";
    cout << "Path in MassCutAndBDT: " << PathInMassCutAndBDT << endl;
    fileInMassCutAndBDT[m] = TFile::Open(PathInMassCutAndBDT);

    // mean value of distribution with different topo variations
    fHistMeanSystMultiTrialMean = (TH1F *)fileInMassCutAndBDT[m]->Get("hSystMultiTrialMean");
    if (!fHistMeanSystMultiTrialMean)
    {
      cout << "Error: histogram not found" << endl;
      return;
    }
    fHistMeanSystMultiTrialMean->SetName("hMeanSystMultiTrialMean_" + Smolt[m]);

    // topological variations
    fHistMassCutAndBDTError = (TH1F *)fileInMassCutAndBDT[m]->Get("hSystMultiTrial2");
    if (!fHistMassCutAndBDTError)
    {
      cout << "Error: histogram not found" << endl;
      return;
    }
    fHistMassCutAndBDTError->SetName("hAbsoluteSystErrorMassCutAndBDT_" + Smolt[m]);

    // primary Lambdas
    fHistPrimaryLambdaError = (TH1F *)fHistPzs2->Clone("fHistPrimaryLambdaError");
    fHistPrimaryLambdaError->Reset();
    for (Int_t pt = 1; pt <= fHistPzs2->GetNbinsX(); pt++)
    {
      fHistPrimaryLambdaError->SetBinContent(pt, std::abs(SecondaryLambdaFraction * fHistPzs2->GetBinContent(pt)));
      fHistPrimaryLambdaError->SetBinError(pt, 0);
    }

    // decay parameter error
    fHistDecayParError = (TH1F *)fHistPzs2->Clone("fHistDecayParError");
    fHistDecayParError->Reset();
    for (Int_t pt = 1; pt <= fHistPzs2->GetNbinsX(); pt++)
    {
      fHistDecayParError->SetBinContent(pt, std::abs(LambdaDecayParameterRelError * fHistPzs2->GetBinContent(pt)));
      fHistDecayParError->SetBinError(pt, 0);
    }

    // Z vertex error
    fHistZVertexError = (TH1F *)fHistPzs2->Clone("fHistZVertexError");
    fHistZVertexError->Reset();
    for (Int_t pt = 1; pt <= fHistPzs2->GetNbinsX(); pt++)
    {
      fHistZVertexError->SetBinContent(pt, ZVertexErrorLambdaOO);
      fHistZVertexError->SetBinError(pt, 0);
    }

    fHistTotalError = (TH1F *)fHistPzs2->Clone("fHistTotalError");
    fHistTotalError->Reset();
    for (Int_t pt = 1; pt <= fHistPzs2->GetNbinsX(); pt++)
    {
      fHistTotalError->SetBinContent(pt, TMath::Sqrt(fHistMassCutAndBDTError->GetBinContent(pt) * fHistMassCutAndBDTError->GetBinContent(pt) +
                                                     fHistPrimaryLambdaError->GetBinContent(pt) * fHistPrimaryLambdaError->GetBinContent(pt) +
                                                     fHistResoError->GetBinContent(pt) * fHistResoError->GetBinContent(pt) +
                                                     fHistDecayParError->GetBinContent(pt) * fHistDecayParError->GetBinContent(pt) +
                                                     fHistZVertexError->GetBinContent(pt) * fHistZVertexError->GetBinContent(pt) +
                                                     fHistPolBkg0Error->GetBinContent(pt) * fHistPolBkg0Error->GetBinContent(pt) +
                                                     fHistPzFitRangeError->GetBinContent(pt) * fHistPzFitRangeError->GetBinContent(pt) +
                                                     fHistBkgExpoError->GetBinContent(pt) * fHistBkgExpoError->GetBinContent(pt)));
      fHistTotalError->SetBinError(pt, 0);
    }
  } // end loop on mult

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

  TString TitleX = TitleXPt;
  if (!isPtAnalysis)
    TitleX = TitleXPsi;
  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasError->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleX, "Absolute syst. uncertainty", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(LowerRangeParticle, UpperRangeParticle);
  hDummy->Draw("");
  fHistMassCutAndBDTError->SetLineColor(kGreen + 2);
  fHistPrimaryLambdaError->SetLineColor(kOrange + 2);
  fHistResoError->SetLineColor(kMagenta);
  fHistDecayParError->SetLineColor(kRed + 1);
  fHistPolBkg0Error->SetLineColor(kGray + 1);
  fHistPolBkg0Error->SetLineStyle(2);
  fHistPolBkg0Error->SetLineWidth(3);
  fHistPzFitRangeError->SetLineColor(kOrange - 3);
  fHistPzFitRangeError->SetLineStyle(2);
  fHistPzFitRangeError->SetLineWidth(3);
  fHistBkgExpoError->SetLineColor(kPink + 1);
  fHistBkgExpoError->SetLineStyle(2);
  fHistBkgExpoError->SetLineWidth(3);
  fHistZVertexError->SetLineColor(kGreen);
  fHistTotalError->SetLineColor(kBlack);
  fHistTotalError->SetLineWidth(4);

  fHistDecayParError->Draw("same");
  fHistPrimaryLambdaError->Draw("same");
  fHistResoError->Draw("same");
  fHistPolBkg0Error->Draw("same");
  fHistPzFitRangeError->Draw("same");
  fHistBkgExpoError->Draw("same");
  fHistZVertexError->Draw("same");
  fHistMassCutAndBDTError->Draw("same");

  fHistTotalError->Draw("same");
  TLegend *legend = new TLegend(0.3, 0.4, 0.9, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.05);
  legend->AddEntry(fHistMassCutAndBDTError, "Topological selections", "l");
  legend->AddEntry(fHistPrimaryLambdaError, "Primary #Lambda", "l");
  legend->AddEntry(fHistResoError, "Resolution", "l");
  legend->AddEntry(fHistZVertexError, "Z_{vtx} selection", "l");
  legend->AddEntry(fHistPolBkg0Error, "P_{z, s2, bkg} = 0", "l");
  legend->AddEntry(fHistPzFitRangeError, "P_{z} fit range", "l");
  legend->AddEntry(fHistBkgExpoError, "Background fit function", "l");
  legend->AddEntry(fHistDecayParError, "Decay parameter", "l");
  legend->AddEntry(fHistTotalError, "Total", "l");
  legend->Draw("same");
  canvasError->SaveAs(Form("../AbsoluteUncertaintySummary_%s_%s.png", SVariableDep.Data(), ParticleName[ChosenPart].Data()));
  canvasError->SaveAs(Form("../AbsoluteUncertaintySummary_%s_%s.pdf", SVariableDep.Data(), ParticleName[ChosenPart].Data()));

  // RelErrors
  TH1F *fHistMassCutAndBDTRelError = (TH1F *)fHistMassCutAndBDTError->Clone("fHistMassCutAndBDTRelError");
  TH1F *fHistPrimaryLambdaRelError = (TH1F *)fHistPrimaryLambdaError->Clone("fHistPrimaryLambdaRelError");
  TH1F *fHistResoRelError = (TH1F *)fHistResoError->Clone("fHistResoRelError");
  TH1F *fHistDecayParRelError = (TH1F *)fHistDecayParError->Clone("fHistDecayParRelError");
  TH1F *fHistPolBkg0RelError = (TH1F *)fHistPolBkg0Error->Clone("fHistPolBkg0RelError");
  TH1F *fHistPzFitRangeRelError = (TH1F *)fHistPzFitRangeError->Clone("fHistPzFitRangeRelError");
  TH1F *fHistBkgExpoRelError = (TH1F *)fHistBkgExpoError->Clone("fHistBkgExpoRelError");
  TH1F *fHistZVertexRelError = (TH1F *)fHistZVertexError->Clone("fHistZVertexRelError");
  TH1F *fHistTotalRelError = (TH1F *)fHistTotalError->Clone("fHistTotalRelError");
  for (Int_t pt = 1; pt <= fHistPzs2->GetNbinsX(); pt++)
  {
    fHistMassCutAndBDTRelError->SetBinContent(pt, fHistMassCutAndBDTError->GetBinContent(pt) / TMath::Abs(fHistPzs2->GetBinContent(pt)));
    fHistPrimaryLambdaRelError->SetBinContent(pt, fHistPrimaryLambdaError->GetBinContent(pt) / TMath::Abs(fHistPzs2->GetBinContent(pt)));
    fHistResoRelError->SetBinContent(pt, fHistResoError->GetBinContent(pt) / TMath::Abs(fHistPzs2->GetBinContent(pt)));
    fHistDecayParRelError->SetBinContent(pt, fHistDecayParError->GetBinContent(pt) / TMath::Abs(fHistPzs2->GetBinContent(pt)));
    fHistPolBkg0RelError->SetBinContent(pt, fHistPolBkg0Error->GetBinContent(pt) / TMath::Abs(fHistPzs2->GetBinContent(pt)));
    fHistPzFitRangeRelError->SetBinContent(pt, fHistPzFitRangeError->GetBinContent(pt) / TMath::Abs(fHistPzs2->GetBinContent(pt)));
    fHistBkgExpoRelError->SetBinContent(pt, fHistBkgExpoError->GetBinContent(pt) / TMath::Abs(fHistPzs2->GetBinContent(pt)));
    fHistZVertexRelError->SetBinContent(pt, fHistZVertexError->GetBinContent(pt) / TMath::Abs(fHistPzs2->GetBinContent(pt)));
    fHistTotalRelError->SetBinContent(pt, fHistTotalError->GetBinContent(pt) / TMath::Abs(fHistPzs2->GetBinContent(pt)));
    fHistMassCutAndBDTRelError->SetBinError(pt, 0);
    fHistPrimaryLambdaRelError->SetBinError(pt, 0);
    fHistResoRelError->SetBinError(pt, 0);
    fHistDecayParRelError->SetBinError(pt, 0);
    fHistPolBkg0RelError->SetBinError(pt, 0);
    fHistPzFitRangeRelError->SetBinError(pt, 0);
    fHistBkgExpoRelError->SetBinError(pt, 0);
    fHistZVertexRelError->SetBinError(pt, 0);
    fHistTotalRelError->SetBinError(pt, 0);
  }

  TCanvas *canvasRelError = new TCanvas("canvasRelError", "canvasRelError", 800, 600);
  StyleCanvas(canvasRelError, 0.05, 0.2, 0.2, 0.03);
  canvasRelError->cd();
  TH1F *hDummyRelError = new TH1F("hDummyRelError", "hDummyRelError", 10000, 0, 100);
  for (Int_t i = 1; i <= hDummyRelError->GetNbinsX(); i++)
    hDummyRelError->SetBinContent(i, 1e-12);
  canvasRelError->cd();
  SetFont(hDummyRelError);
  StyleHistoYield(hDummyRelError, 0, 0.2, 1, 1, TitleX, "Relative syst. uncertainty", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummyRelError, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, 2, yLabelOffset);
  SetTickLength(hDummyRelError, tickX, tickY);
  hDummyRelError->GetXaxis()->SetRangeUser(LowerRangeParticle, UpperRangeParticle);
  hDummyRelError->Draw("");
  fHistMassCutAndBDTRelError->Draw("same");
  fHistPrimaryLambdaRelError->Draw("same");
  fHistDecayParRelError->Draw("same");
  fHistResoRelError->Draw("same");
  fHistPolBkg0RelError->Draw("same");
  fHistPzFitRangeRelError->Draw("same");
  fHistBkgExpoRelError->Draw("same");
  fHistZVertexRelError->Draw("same");
  fHistTotalRelError->Draw("same");
  // legend->Draw("same");
  canvasRelError->SaveAs(Form("../RelativeUncertaintySummary_%s_%s.png", SVariableDep.Data(), ParticleName[ChosenPart].Data()));
  canvasRelError->SaveAs(Form("../RelativeUncertaintySummary_%s_%s.pdf", SVariableDep.Data(), ParticleName[ChosenPart].Data()));

  TFile *fileout = new TFile(stringout, "RECREATE");
  fHistMassCutAndBDTError->Write();
  fHistDecayParError->Write();
  fHistPrimaryLambdaError->Write();
  fHistResoError->Write();
  fHistPolBkg0Error->Write();
  fHistPzFitRangeError->Write();
  fHistBkgExpoError->Write();
  fHistZVertexError->Write();
  fHistTotalError->Write();
  fileout->Close();

  cout << "\nI have created the file:\n " << stringout << endl;
}