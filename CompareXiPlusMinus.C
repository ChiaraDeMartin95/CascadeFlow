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

// take spectra in input
// produces ratio of spectra wrt 0-100% multiplciity class

Float_t YLowMean[numPart] = {1.31, 1.66};
Float_t YUpMean[numPart] = {1.327, 1.68};
Float_t YLowSigma[numPart] = {0.0, 0.0};
Float_t YUpSigma[numPart] = {0.006, 0.006};
Float_t YLowPurity[numPart] = {0.8, 0};
Float_t YLowV2[numPart] = {-0.4, -0.4};
Float_t YUpV2[numPart] = {0.5, 0.5};
Float_t YLowPzs2[numPart] = {-0.1, -0.1};
Float_t YUpPzs2[numPart] = {0.1, 0.1};

Float_t YLow[numPart] = {0};
Float_t YUp[numPart] = {0};

// choices: 0: mean, 1: sigma, 2: purity, 3: yield, 4: v2, 5: Pzs2 (Pz if vs Psi), 6: Pzs2FromLambda (Pz if vs Psi)
Float_t YLowRatio[numChoice] = {0.99, 0.2, 0.8, 0.1, 0, -1.5, -1.5};
Float_t YUpRatio[numChoice] = {1.01, 1.8, 1.2, 4, 2, -0.5, -0.5};

void CompareXiPlusMinus(Bool_t isPtAnalysis = 1,
                        Int_t ChosenPart = ChosenParticle,
                        Int_t Choice = 0,
                        Int_t ChosenMult = 5, /*numCent,*/
                        TString OutputDir = "MeanSigmaPurityMultClasses/",
                        Int_t BkgType = ExtrBkgType,
                        Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  Int_t part = 0;
  if (ChosenPart == 1 || ChosenPart == 4 || ChosenPart == 5)
  {
    part = 1;
  }
  gStyle->SetOptStat(0);
  if (ChosenPart != 0 && ChosenPart != 1)
  {
    cout << "ChosenParticle must be 0 (Xi) or 1 (Omega)" << endl;
    return;
  }
  if ((ChosenMult > numCent && !isV2) || (ChosenMult > (numCent - 1) && isV2))
  {
    cout << "Chosen Mult outside of available range" << endl;
    return;
  }
  if (!isPtAnalysis)
  {
    if (Choice == 5)
    {
      TypeHisto[Choice] = "Pz";
      TitleY[Choice] = "P_{z}";
    }
    if (Choice == 6)
    {
      TypeHisto[Choice] = "PzLambdaFromC";
      TitleY[Choice] = "P_{z}";
    }
  }
  cout << Choice << " " << TypeHisto[Choice] << endl;
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
  else if (Choice == 4)
  {
    YLow[part] = YLowV2[part];
    YUp[part] = YUpV2[part];
  }
  else if (Choice == 5)
  {
    YLow[part] = YLowPzs2[part];
    YUp[part] = YUpPzs2[part];
  }
  else if (Choice == 6)
  {
    YLow[part] = YLowPzs2[part];
    YUp[part] = YUpPzs2[part];
  }

  // multiplicity related variables
  TString Smolt[numCent + 1];
  TString SmoltBis[numCent + 1];
  TString sScaleFactorFinal[numCent + 1];
  Float_t ScaleFactorFinal[numCent + 1];

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn;
  TFile *fileIn[numCent + 1];

  // cent min and max
  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  if (ChosenMult == numCent)
  { // 0-80%
    CentFT0CMin = 0;
    CentFT0CMax = 80;
  }
  else
  {
    CentFT0CMin = CentFT0C[ChosenMult];
    CentFT0CMax = CentFT0C[ChosenMult + 1];
  }

  // canvases
  TCanvas *canvasPtSpectra = new TCanvas("canvasPtSpectra", "canvasPtSpectra", 700, 900);
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

  TH1F *fHistSpectrum[numCent + 1][3];
  TH1F *fHistSpectrumScaled[numCent + 1][3];
  TH1F *fHistSpectrumMultRatio[numCent + 1][3];

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *legendAllMult;
  legendAllMult = new TLegend(0.22, 0.03, 0.73, 0.22);
  if (Choice == 2 && (ChosenPart == 1 || ChosenPart == 4 || ChosenPart == 5))
  {
    legendAllMult = new TLegend(0.44, 0.03, 0.9, 0.28);
  }
  legendAllMult->SetHeader("");
  legendAllMult->SetNColumns(3);
  legendAllMult->SetFillStyle(0);
  TLegendEntry *lheaderAllMult = (TLegendEntry *)legendAllMult->GetListOfPrimitives()->First();
  lheaderAllMult->SetTextSize(0.04);

  TLegend *LegendTitle;
  if (Choice == 2)
    LegendTitle = new TLegend(0.54, 0.5, 0.95, 0.72);
  else
    LegendTitle = new TLegend(0.54, 0.7, 0.95, 0.92);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(33);
  LegendTitle->SetTextSize(0.04);
  LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  LegendTitle->AddEntry("", "Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
  if (isV2)
    LegendTitle->AddEntry("", ParticleNameLegend[part] + ", |#it{#eta}| < 0.8", "");
  else
    LegendTitle->AddEntry("", ParticleNameLegend[part] + ", |#it{y}| < 0.5", "");
  LegendTitle->AddEntry("", Form("FT0C %i#minus%i", CentFT0CMin, CentFT0CMax), "");

  TLine *lineat1Mult = new TLine(MinPt[part], 1, MaxPt[part], 1);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  // get spectra in multiplicity classes
  TString Name[3] = {"Minus", "Plus", ""};
  for (Int_t charge = 0; charge <= 2; charge++)
  {
    if (Choice >= 5 && charge == 2)
      continue; // Pzs is studied for separate charged only
    // charge = 0: NegCharge, charge = 1: PosCharge, charge = 2: AllCharge

    PathIn = "OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_";
    PathIn += SinputFileName;
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += "_" + ParticleName[ChosenPart] + Name[charge] + SIsBkgParab[BkgType];
    Smolt[ChosenMult] = Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);
    SmoltBis[ChosenMult] = Form("%i#minus%i", CentFT0CMin, CentFT0CMax);
    PathIn += Smolt[ChosenMult];
    if (isApplyWeights)
      PathIn += "_Weighted";
    if (!useCommonBDTValue)
      PathIn += "_BDTCentDep";
    if (isRun2Binning)
      PathIn += "_Run2Binning";
    if (!isPtAnalysis)
      PathIn += "_vsPsi";
    if (Choice == 6)
      PathIn += "_PolFromLambda";
    PathIn += ".root";
    cout << "Path in : " << PathIn << endl;
    fileIn[ChosenMult] = TFile::Open(PathIn);
    if (Choice == 1)
      fHistSpectrum[ChosenMult][charge] = (TH1F *)fileIn[ChosenMult]->Get("histo" + TypeHisto[Choice] + "Weighted");
    else
      fHistSpectrum[ChosenMult][charge] = (TH1F *)fileIn[ChosenMult]->Get("histo" + TypeHisto[Choice]);
    if (!fHistSpectrum[ChosenMult][charge])
    {
      cout << " no hist " << endl;
      return;
    }
    fHistSpectrum[ChosenMult][charge]->SetName("histoSpectrum_" + Smolt[ChosenMult] + "_" + ParticleName[charge]);
  }

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

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 8);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasPtSpectra->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, titlex, TitleY[Choice], "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  if (isPtAnalysis)
    hDummy->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
  else
    hDummy->GetXaxis()->SetRangeUser(0, 2 * TMath::Pi());
  pad1->Draw();
  pad1->cd();
  if (Choice == 3)
    gPad->SetLogy();
  hDummy->Draw("same");

  for (Int_t charge = 0; charge <= 2; charge++)
  {
    if (Choice >= 5 && charge == 2)
      continue; // Pzs is studied for separate charged only
    ScaleFactorFinal[ChosenMult] = ScaleFactor[ChosenMult];
    fHistSpectrum[ChosenMult][charge] = (TH1F *)fHistSpectrum[ChosenMult][charge]->Clone("fHistSpectrumScaled_" + Smolt[ChosenMult]);
    if (Choice == 3)
      fHistSpectrum[ChosenMult][charge]->Scale(ScaleFactorFinal[ChosenMult]);
    for (Int_t b = 1; b <= fHistSpectrum[ChosenMult][charge]->GetNbinsX(); b++)
    {
      // cout << "bin " << b << " " << fHistSpectrum[ChosenMult][charge]->GetBinContent(b) << "+-" << fHistSpectrum[ChosenMult][charge]->GetBinError(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrum[ChosenMult][charge]->GetBinContent(b) << "+-" << fHistSpectrum[ChosenMult][charge]->GetBinError(b) << endl;
    }
    SetFont(fHistSpectrum[ChosenMult][charge]);
    fHistSpectrum[ChosenMult][charge]->SetMarkerColor(ColorMult[charge]);
    fHistSpectrum[ChosenMult][charge]->SetLineColor(ColorMult[charge]);
    fHistSpectrum[ChosenMult][charge]->SetMarkerStyle(MarkerMult[charge]);
    fHistSpectrum[ChosenMult][charge]->SetMarkerSize(0.6 * SizeMult[charge]);
    fHistSpectrum[ChosenMult][charge]->Draw("same e0x0");
    if (Choice == 3)
      sScaleFactorFinal[ChosenMult] = Form(" (x2^{%i})", int(log2(ScaleFactorFinal[ChosenMult])));
    else
      sScaleFactorFinal[ChosenMult] = "";
    if (part == 0)
      legendAllMult->AddEntry(fHistSpectrum[ChosenMult][charge], ParticleNameLegend[2 + charge] + sScaleFactorFinal[ChosenMult] + " ", "pef");
    else
      legendAllMult->AddEntry(fHistSpectrum[ChosenMult][charge], ParticleNameLegend[4 + charge] + sScaleFactorFinal[ChosenMult] + " ", "pef");
  } // end loop on charge
  LegendTitle->Draw("");
  legendAllMult->Draw("");

  // PDG mass
  TF1 *fMassPDG = new TF1("fMassPDG", Form("%f", ParticleMassPDG[part]), MinPt[part], MaxPt[part]);
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

  TString TitleYSpectraRatio = "Ratio to ";
  if (part == 0)
    TitleYSpectraRatio += ParticleNameLegend[3];
  else
    TitleYSpectraRatio += ParticleNameLegend[5];
  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, 0, 8);
  for (Int_t i = 1; i <= hDummyRatio->GetNbinsX(); i++)
    hDummyRatio->SetBinContent(i, 1e-12);
  SetFont(hDummyRatio);
  StyleHistoYield(hDummyRatio, YLowRatio[Choice], YUpRatio[Choice], 1, 1, titlex, TitleYSpectraRatio, "", 1, 1.15, YoffsetSpectraRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetTickLength(hDummyRatio, tickX, tickY);
  if (isPtAnalysis)
    hDummyRatio->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
  else
    hDummyRatio->GetXaxis()->SetRangeUser(0, 2 * TMath::Pi());
  canvasPtSpectra->cd();
  padL1->Draw();
  padL1->cd();
  hDummyRatio->Draw("same");

  for (Int_t charge = 0; charge <= 2; charge++)
  {
    if (charge == 2)
      continue;

    fHistSpectrumMultRatio[ChosenMult][charge] = (TH1F *)fHistSpectrum[ChosenMult][charge]->Clone("fHistSpectrumMultRatio_" + Smolt[ChosenMult]);
    fHistSpectrumMultRatio[ChosenMult][charge]->Divide(fHistSpectrum[ChosenMult][1]);
    // fHistSpectrumMultRatio[ChosenMult][charge]->Add(fHistSpectrum[ChosenMult][1], -1);
    ErrRatioCorr(fHistSpectrum[ChosenMult][charge], fHistSpectrum[ChosenMult][1], fHistSpectrumMultRatio[ChosenMult][charge], 0);
    for (Int_t b = 1; b <= fHistSpectrum[ChosenMult][charge]->GetNbinsX(); b++)
    {
      // cout << "bin " << b << " " << fHistSpectrum[ChosenMult][charge]->GetBinContent(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrum[ChosenMult]->GetBinContent(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrumMultRatio[ChosenMult][charge]->GetBinContent(b) << endl;
    }
    fHistSpectrumMultRatio[ChosenMult][charge]->SetMarkerColor(ColorMult[charge]);
    fHistSpectrumMultRatio[ChosenMult][charge]->SetLineColor(ColorMult[charge]);
    fHistSpectrumMultRatio[ChosenMult][charge]->SetMarkerStyle(MarkerMult[charge]);
    fHistSpectrumMultRatio[ChosenMult][charge]->SetMarkerSize(SizeMultRatio[charge]);
    if (charge == 0)
    {
      fHistSpectrumMultRatio[ChosenMult][charge]->Draw("same e0x0");
    }
    lineat1Mult->Draw("same");
  }

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = OutputDir + "PlotXiMinusPlusRatios" + NameAnalysis[!isV2] + "_";
  stringout += SinputFileName;
  stringout += Smolt[ChosenMult];
  stringout += "_" + ParticleName[ChosenPart];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[BkgType];
  stringout += "_" + TypeHisto[Choice];
  if (isApplyWeights)
    stringout += "_Weighted";
  if (!useCommonBDTValue)
    stringout += "_BDTCentDep";
  if (isRun2Binning)
    stringout += "_Run2Binning";
  if (!isPtAnalysis)
    stringout += "_vsPsi";
  stringoutpdf = stringout;
  stringout += "_5Cent.root";
  TFile *fileout = new TFile(stringout, "RECREATE");
  fileout->Close();

  canvasPtSpectra->SaveAs(stringoutpdf + ".pdf");
  canvasPtSpectra->SaveAs(stringoutpdf + ".png");
  cout << "\nStarting from the files (for the different mult): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
}
