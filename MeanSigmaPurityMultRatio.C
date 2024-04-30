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
Float_t YLowPurity[numPart] = {0, 0};
Float_t YLowV2[numPart] = {-0.4, -0.4};
Float_t YUpV2[numPart] = {0.5, 0.5};

Float_t YLow[numPart] = {0};
Float_t YUp[numPart] = {0};

Float_t YLowRatio[numChoice] = {0.99, 0.2, 0, 0.1, -1};
Float_t YUpRatio[numChoice] = {1.01, 1.8, 1.2, 4, 2};

TString NameRun2Dir[numPart] = {"chists_Xim_variations", "chists_Omm_variations"};

TString TypeHistoLucia[numChoice - 1] = {"fHistMean", "fHistSigma", "fHistPutityBinCount_Fit", "hFinalSpectrumCount_Fit"};
Int_t MarkerMultLucia[11] = {24, 24, 24, 25, 27, 28, 30, 24, 25, 27, 28};
TString MultChoiceLucia[11] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
TString MultClassLucia[11] = {"0-90%", "0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%"};

void MeanSigmaPurityMultRatio(Bool_t isXi = ChosenParticleXi,
                              Int_t Choice = 0,
                              Int_t part = ExtrParticle,
                              Int_t ChosenMultLucia = 0,
                              Int_t ChosenMult = numCent - 3,
                              Bool_t isDrawRun2 = 1,
                              Bool_t isDrawLuciaRun3 = 0,
                              TString SysPath = "",
                              TString OutputDir = "MeanSigmaPurityMultClasses/",
                              TString inputFileName = SinputFileName,
                              Int_t BkgType = ExtrBkgType,
                              Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  gStyle->SetOptStat(0);
  if (ChosenMult > numCent - 1)
  {
    cout << "Chosen Mult outside of available range" << endl;
    return;
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
    if (isDrawLuciaRun3)
    {
      YLow[part] = 0;
      YUp[part] = 0.05;
    }
  }
  else if (Choice == 4)
  {
    YLow[part] = YLowV2[part];
    YUp[part] = YUpV2[part];
  }

  // multiplicity related variables
  TString Smolt[numCent];
  TString SmoltBis[numCent];
  TString sScaleFactorFinal[numCent];
  Float_t ScaleFactorFinal[numCent];

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn;
  TFile *fileIn[numCent];

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = OutputDir + "PlotRatios_";
  stringout += SinputFileName;
  stringout += "_" + ParticleName[!isXi];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[BkgType];
  stringout += "_" + TypeHisto[Choice];
  stringoutpdf = stringout;
  stringout += "_5Cent.root";
  TFile *fileout = new TFile(stringout, "RECREATE");

  // canvases
  TCanvas *canvasPtSpectra = new TCanvas("canvasPtSpectra", "canvasPtSpectra", 700, 900);
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

  TH1F *fHistSpectrum[numCent];
  TH1F *fHistSpectrumScaled[numCent];
  TH1F *fHistSpectrumMultRatio[numCent];

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *legendAllMult;
  legendAllMult = new TLegend(0.22, 0.03, 0.73, 0.28);
  if (Choice == 2 && !isXi)
  {
    legendAllMult = new TLegend(0.44, 0.03, 0.9, 0.28);
  }
  else if (isDrawLuciaRun3)
  {
    legendAllMult = new TLegend(0.45, 0.52, 0.97, 0.78);
  }
  legendAllMult->SetHeader("FT0C Centrality Percentile");
  legendAllMult->SetNColumns(3);
  legendAllMult->SetFillStyle(0);
  TLegendEntry *lheaderAllMult = (TLegendEntry *)legendAllMult->GetListOfPrimitives()->First();
  lheaderAllMult->SetTextSize(0.04);

  TLegend *LegendTitle;
  if (Choice == 2)
    LegendTitle = new TLegend(0.54, 0.55, 0.95, 0.72);
  else
    LegendTitle = new TLegend(0.54, 0.75, 0.95, 0.92);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(33);
  LegendTitle->SetTextSize(0.04);
  LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  LegendTitle->AddEntry("", "PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
  LegendTitle->AddEntry("", ParticleName[part] + ", |#it{#eta}| < 0.8", "");

  TLine *lineat1Mult = new TLine(MinPt[part], 1, MaxPt[part], 1);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  // get spectra in multiplicity classes
  for (Int_t m = numCent - 1; m >= 0; m--)
  {
    if (m == 0 || m > (numCent - 3))
      continue;
    PathIn = "OutputAnalysis/FitV2_";
    PathIn += SinputFileName;
    PathIn += "_" + ParticleName[!isXi];
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[BkgType];
    Smolt[m] += Form("_Cent%i-%i", CentFT0C[m], CentFT0C[m + 1]);
    SmoltBis[m] += Form("%i#minus%i", CentFT0C[m], CentFT0C[m + 1]);
    PathIn += Smolt[m];
    PathIn += ".root";
    cout << "Path in : " << PathIn << endl;

    fileIn[m] = TFile::Open(PathIn);
    fHistSpectrum[m] = (TH1F *)fileIn[m]->Get("histo" + TypeHisto[Choice]);
    fHistSpectrum[m]->SetName("histoSpectrum_" + Smolt[m]);
    if (!fHistSpectrum[m])
    {
      cout << " no hist " << endl;
      return;
    }
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

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 8);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasPtSpectra->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleXPt, TitleY[Choice], "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
  pad1->Draw();
  pad1->cd();
  if (Choice == 3 && !isDrawLuciaRun3)
    gPad->SetLogy();
  hDummy->Draw("same");

  for (Int_t m = numCent - 1; m >= 0; m--)
  {
    if (m == 0 || m > (numCent - 3))
      continue;
    ScaleFactorFinal[m] = ScaleFactor[m];
    fHistSpectrumScaled[m] = (TH1F *)fHistSpectrum[m]->Clone("fHistSpectrumScaled_" + Smolt[m]);
    if (Choice == 3 && !isDrawLuciaRun3)
      fHistSpectrumScaled[m]->Scale(ScaleFactorFinal[m]);
    for (Int_t b = 1; b <= fHistSpectrum[m]->GetNbinsX(); b++)
    {
      // cout << "bin " << b << " " << fHistSpectrum[m]->GetBinContent(b) << "+-" << fHistSpectrum[m]->GetBinError(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrumScaled[m]->GetBinContent(b) << "+-" << fHistSpectrumScaled[m]->GetBinError(b) << endl;
    }
    SetFont(fHistSpectrumScaled[m]);
    fHistSpectrumScaled[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumScaled[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumScaled[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumScaled[m]->SetMarkerSize(0.6 * SizeMult[m]);
    fHistSpectrumScaled[m]->Draw("same e0x0");
    if (Choice == 3 && !isDrawLuciaRun3)
      sScaleFactorFinal[m] = Form(" (x2^{%i})", int(log2(ScaleFactorFinal[m])));
    else
      sScaleFactorFinal[m] = "";
    legendAllMult->AddEntry(fHistSpectrumScaled[m], SmoltBis[m] + "%" + sScaleFactorFinal[m] + " ", "pef");
  } // end loop on mult
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

  // Add Run 2 values for mean and sigma
  TFile *fileRun2 = new TFile("Run2MeanSigma/Alessandrosignal.root");
  TDirectoryFile *dirRun2;
  TDirectoryFile *dirRun2Def;
  TDirectoryFile *dirRun2Checks;
  TDirectoryFile *dirRun2Params;
  TDirectoryFile *dirRun2Choice;
  TLegend *legendRun2 = new TLegend(0.23, 0.85, 0.57, 0.94);
  legendRun2->SetBorderSize(0);
  legendRun2->SetFillStyle(0);
  legendRun2->SetTextSize(0.04);

  if (isDrawRun2 && (Choice == 0 || Choice == 1))
  {
    dirRun2 = (TDirectoryFile *)fileRun2->Get(NameRun2Dir[part]);
    dirRun2Def = (TDirectoryFile *)dirRun2->Get("Default");
    dirRun2Checks = (TDirectoryFile *)dirRun2Def->Get("Checks");
    dirRun2Params = (TDirectoryFile *)dirRun2Checks->Get("SignalParams");
    dirRun2Choice = (TDirectoryFile *)dirRun2Params->Get(TypeHisto[Choice]);
    TH1F *histoRun2 = (TH1F *)dirRun2Choice->Get(TypeHisto[Choice] + "_9"); // 0-90% class
    histoRun2->SetMarkerColor(kBlack);
    histoRun2->SetLineColor(kBlack);
    histoRun2->SetMarkerStyle(20);
    histoRun2->SetMarkerSize(1.5);
    histoRun2->Draw("same e0x0");
    legendRun2->AddEntry(histoRun2, "Run 2, 0-90%", "pe");
    legendRun2->Draw();
  }

  // Add Run 3 values computed by Lucia
  TLegend *legendRun3 = new TLegend(0.23, 0.3, 0.57, 0.4);
  if (Choice == 3)
    legendRun3 = new TLegend(0.20, 0.85, 0.54, 0.95);
  legendRun3->SetBorderSize(0);
  legendRun3->SetFillStyle(0);
  legendRun3->SetTextSize(0.04);

  if (isDrawLuciaRun3 && Choice != 4 && part == 0)
  {
    for (Int_t m = 0; m < 11; m++)
    {
      // if (m != ChosenMultLucia)
      // continue; // decomment if you want only ChosenMultLucia histo to be displayed
      if (m < 3 || m > 8)
        continue;
      ChosenMultLucia = m;
      TString fileLuciaRun3 = "Files_Chiara/";
      fileLuciaRun3 += ParticleName[part];
      fileLuciaRun3 += "Minus_LHC23_pass2_newCuts_" + MultChoiceLucia[ChosenMultLucia] + "_AfterSel_DoubleGauss.root";
      cout << "file Lucia " << fileLuciaRun3 << endl;
      TFile *fileRun3 = new TFile(fileLuciaRun3);
      TH1F *histoRun3 = (TH1F *)fileRun3->Get(TypeHistoLucia[Choice]);
      histoRun3->SetName("histoRun3_" + MultChoiceLucia[ChosenMultLucia]);
      histoRun3->SetMarkerColor(kGray + 2);
      histoRun3->SetLineColor(kGray + 2);
      histoRun3->SetMarkerStyle(24);
      if (ChosenMultLucia > 1)
      {
        histoRun3->SetMarkerColor(ColorMult[ChosenMultLucia - 2]);
        histoRun3->SetLineColor(ColorMult[ChosenMultLucia - 2]);
        histoRun3->SetMarkerStyle(MarkerMultLucia[ChosenMultLucia]);
      }
      histoRun3->SetMarkerSize(1.5);
      histoRun3->DrawCopy("same e0x0");
      if (m == ChosenMultLucia)
        legendRun3->AddEntry(histoRun3, "Run 3, Lucia " + MultClassLucia[ChosenMultLucia], "pe");
    }
    //legendRun3->Draw();
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
  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, 0, 8);
  for (Int_t i = 1; i <= hDummyRatio->GetNbinsX(); i++)
    hDummyRatio->SetBinContent(i, 1e-12);
  SetFont(hDummyRatio);
  StyleHistoYield(hDummyRatio, YLowRatio[Choice], YUpRatio[Choice], 1, 1, TitleXPt, TitleYSpectraRatio, "", 1, 1.15, YoffsetSpectraRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetTickLength(hDummyRatio, tickX, tickY);
  hDummyRatio->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
  canvasPtSpectra->cd();
  padL1->Draw();
  padL1->cd();
  hDummyRatio->Draw("same");

  for (Int_t m = numCent - 1; m >= 0; m--)
  {
    if (m == 0 || m > (numCent - 3))
      continue;
    fHistSpectrumMultRatio[m] = (TH1F *)fHistSpectrum[m]->Clone("fHistSpectrumMultRatio_" + Smolt[m]);
    fHistSpectrumMultRatio[m]->Divide(fHistSpectrum[ChosenMult]);
    ErrRatioCorr(fHistSpectrum[m], fHistSpectrum[ChosenMult], fHistSpectrumMultRatio[m], 0);
    for (Int_t b = 1; b <= fHistSpectrum[m]->GetNbinsX(); b++)
    {
      // cout << "bin " << b << " " << fHistSpectrum[m]->GetBinContent(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrum[ChosenMult]->GetBinContent(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrumMultRatio[m]->GetBinContent(b) << endl;
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
  fileout->Close();
  canvasPtSpectra->SaveAs(stringoutpdf + ".pdf");
  canvasPtSpectra->SaveAs(stringoutpdf + ".png");
  cout << "\nStarting from the files (for the different mult): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
}
