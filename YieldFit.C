#include "Riostream.h"
#include "string"
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
#include "ErrRatioCorr.C"
#include "CommonVarPub.h"
#include "CommonVarLambda.h"

TF1 *GetLevidNdptTimesPt(Double_t norm, Double_t n, Double_t temp, Double_t mass, const char *name)
{

  // Levi function, dNdpt
  char formula[500];

  snprintf(formula, 500, "( x*[0]*([1]-1)*([1]-2)  )/( [1]*[2]*( [1]*[2]+[3]*([1]-2) )  ) * pow( 1 + (sqrt([3]*[3]+x*x) -[3])/([1]*[2])  , -[1])");
  TF1 *fLastFunc = new TF1(name, formula, 0, 10);
  fLastFunc->SetParameters(norm, n, temp, mass);
  fLastFunc->SetParLimits(2, 0.01, 10);
  fLastFunc->SetParNames("norm (dN/dy)", "n", "T", "mass");
  fLastFunc->FixParameter(3, mass);
  fLastFunc->SetLineWidth(1);
  return fLastFunc;
}

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
// fits them to get pt-integrated yields

void YieldFit(Int_t typefit = 3, // mT scaling, Boltzmann, Fermi-Direc, Levi
              Int_t ChosenPart = ChosenParticle)
{

  if (typefit != 3)
  {
    cout << "Not yet implemented" << endl;
    return;
  }
  // multiplicity related variables
  const int numMult = 10; // number of multiplicity bins
  TString Smolt[numMult + 1] = {"0_5", "5_10", "10_20", "20_30", "30_40", "40_50", "50_60", "60_70", "70_80", "80_90", "90_100"};
  TString SmoltBis[numMult + 1];

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn = "../XiSpectraSara/XiSpectraOO_Sara.root";
  TFile *fileIn = TFile::Open(PathIn);
  if (!fileIn)
  {
    cout << "Error: cannot open input file " << PathIn << endl;
    return;
  }
  cout << "Input file opened successfully: " << PathIn << endl;

  // fileout name
  TString stringout = "../XiSpectraSara/fitfunctions.root";

  // canvases
  TCanvas *canvasPtSpectra = new TCanvas("canvasPtSpectra", "canvasPtSpectra", 700, 900);
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

  TH1F *fHistSpectrumStat[numMult + 1];
  TH1F *fHistSpectrumSist[numMult + 1];
  TH1F *fHistSpectrumStatScaled[numMult + 1];
  TH1F *fHistSpectrumSistScaled[numMult + 1];
  TH1F *fHistSpectrumSistScaledForLegend[numMult + 1];
  TH1F *fHistSpectrumRatioFit[numMult + 1];
  TH1F *fHistSpectrumSistMultRatio[numMult + 1];

  gStyle->SetOptStat(0);

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *legendAllMult = new TLegend(0.22, 0.03, 0.73, 0.28);
  legendAllMult->SetNColumns(3);
  legendAllMult->SetFillStyle(0);

  TLegend *LegendTitle = new TLegend(0.54, 0.75, 0.95, 0.92);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(33);
  LegendTitle->SetTextSize(0.04);

  TLegend *LegendPub = new TLegend(0.54, 0.65, 0.95, 0.73);
  LegendPub->SetFillStyle(0);
  LegendPub->SetTextAlign(23);
  LegendPub->SetTextSize(0.025);

  TLine *lineat1Mult = new TLine(0, 1, 8, 1);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  // draw spectra in multiplicity classes
  TString sScaleFactorFinal[numMult + 1] = {};
  Float_t ScaleFactorFinal[numMult + 1] = {512, 256, 128, 64, 32, 16, 8, 4, 2, 1};

  // get spectra in multiplicity classes
  for (Int_t m = numMult; m >= 0; m--)
  {
    cout << "Processing multiplicity class: " << Smolt[m] << endl;
    fHistSpectrumStat[m] = (TH1F *)fileIn->Get("hCorrectedYields_Xi_vs_pt_mult_" + Smolt[m] + "_INEL > 0");
    if (!fHistSpectrumStat[m])
    {
      cout << " no hist spectrum stat" << endl;
      return;
    }

    fHistSpectrumSist[m] = (TH1F *)fileIn->Get("hCorrectedYields_syst_Xi_vs_pt_mult_" + Smolt[m]);
    if (!fHistSpectrumSist[m])
    {
      cout << " no hist spectrum sist" << endl;
      return;
    }

    if (m == numMult)
    {
      ColorMult[m] = kBlack;
      MarkerMult[m] = 20;
      SizeMult[m] = 1.2;
      ScaleFactorFinal[m] = 1;
    }
    cout << "ScaleFactorFinal[" << m << "] = " << ScaleFactorFinal[m] << endl;

    fHistSpectrumSistScaled[m] = (TH1F *)fHistSpectrumSist[m]->Clone("fHistSpectrumSistScaled_" + Smolt[m]);
    fHistSpectrumStatScaled[m] = (TH1F *)fHistSpectrumStat[m]->Clone("fHistSpectrumStatScaled_" + Smolt[m]);

    fHistSpectrumRatioFit[m] = (TH1F *)fHistSpectrumStat[m]->Clone(Form("HistRatioFit_m%i_typefit%i", m, typefit));

  } // end loop on mult

  TLegend *legendfit = new TLegend(0.25, 0.25, 0.4, 0.45);
  legendfit->SetFillStyle(0);
  legendfit->SetTextSize(0.04);
  TString namepwgfunc[numMult + 1];

  // fit spectra
  const int numfittipo = 5;
  Int_t ColorFit[numfittipo + 1] = {860, 881, 868, 628, 419, kAzure + 7};
  TString nameFit[numfittipo + 1] = {"MTExp", "Boltzmann", "Fermi-Dirac", "Levy"};
  TFitResultPtr fFitResultPtr0[numMult + 1];
  TF1 *fit_pwgfunc[numMult + 1];
  TF1 *fit_pwgfunc_Scaled[numMult + 1];
  Float_t Chi2NDF[numMult + 1] = {0};
  Float_t Temp[numMult + 1] = {0};
  Float_t TempError[numMult + 1] = {0};

  Double_t LowRange[numMult + 1] = {0, 0, 0, 0, 0, 0};
  Double_t UpRange[numMult + 1] = {4, 4, 4, 4, 4, 4};

  TLegend *legendfitSummary = new TLegend(0.73, 0.65, 0.92, 0.75);
  legendfitSummary->SetFillStyle(0);
  legendfitSummary->SetTextSize(0.04);
  legendfitSummary->SetTextAlign(33);
  legendfitSummary->AddEntry("", nameFit[typefit] + " fit", "");

  Int_t factor = 1;

  for (Int_t m = numMult; m >= 0; m--)
  {
    cout << "Processing multiplicity class: " << m << endl;
    namepwgfunc[m] = Form("fitpwgfunc_m%i_fit%i", m, typefit);
    fit_pwgfunc[m] = GetLevidNdptTimesPt(0.04 * factor, 0.03, 0.1, ParticleMassPDG[ChosenPart], namepwgfunc[m]);
    fit_pwgfunc[m]->SetParLimits(0, 0, fHistSpectrumStat[m]->GetBinContent(fHistSpectrumStat[m]->GetMaximumBin()) * 0.5 * 10); // norm
    fit_pwgfunc[m]->SetParLimits(1, 2, 30);                                                                                    // n
    fit_pwgfunc[m]->SetParLimits(2, 0.1, 10);                                                                                  // T
    fit_pwgfunc[m]->SetParameter(2, 0.7);

    fit_pwgfunc[m]->SetLineColor(ColorMult[m]);
    fit_pwgfunc[m]->SetLineStyle(7);
    fit_pwgfunc[m]->SetLineWidth(2);
    fit_pwgfunc[m]->SetRange(LowRange[m], UpRange[m]);
    fFitResultPtr0[m] = fHistSpectrumStat[m]->Fit(fit_pwgfunc[m], "SR0I");
    fit_pwgfunc[m]->SetRange(0, 7);
    Chi2NDF[m] = fit_pwgfunc[m]->GetChisquare() / fit_pwgfunc[m]->GetNDF();
  } // end loop on multiplicity classes

  cout << "Finished to fit spectra " << endl;

  Float_t LimSupSpectra = 9999.99;
  Float_t LimInfSpectra = 0.2 * 1e-6;
  Float_t xTitle = 15;
  Float_t xOffset = 4;
  Float_t yTitle = 30;
  Float_t yOffset = 2;

  Float_t xLabel = 30;
  Float_t yLabel = 30;
  Float_t xLabelOffset = 0.05;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.035;
  Float_t tickY = 0.035;

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 1000, 0, 8);
  canvasPtSpectra->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, LimInfSpectra, LimSupSpectra, 1, 1, TitleXPt, "", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  pad1->Draw();
  pad1->cd();
  gPad->SetLogy();
  hDummy->Draw("");

  for (Int_t m = numMult; m >= 0; m--)
  {
    fHistSpectrumStatScaled[m]->Scale(ScaleFactorFinal[m]);
    fHistSpectrumSistScaled[m]->Scale(ScaleFactorFinal[m]);
    fHistSpectrumStatScaled[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumStatScaled[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumStatScaled[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumStatScaled[m]->SetMarkerSize(SizeMult[m]);
    fHistSpectrumSistScaled[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumSistScaled[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumSistScaled[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumSistScaled[m]->SetMarkerSize(SizeMult[m]);
    fHistSpectrumStatScaled[m]->Draw("same e0x0");
    fHistSpectrumSistScaled[m]->SetFillStyle(0);
    fHistSpectrumSistScaled[m]->Draw("same e2");

    fit_pwgfunc_Scaled[m] = (TF1 *)fit_pwgfunc[m]->Clone(namepwgfunc[m] + "_Scaled");
    fit_pwgfunc_Scaled[m]->SetLineColor(ColorMult[m]);
    fit_pwgfunc_Scaled[m]->SetParameter(0, fit_pwgfunc[m]->GetParameter(0) * ScaleFactorFinal[m]);
    fit_pwgfunc_Scaled[m]->Draw("same");
  }
  legendfit->Draw("");
  LegendTitle->Draw("");
  legendAllMult->Draw("");

  // Compute and draw spectra ratios
  Float_t LimSupMultRatio = 1.9;
  Float_t LimInfMultRatio = 1e-2;
  Float_t YoffsetSpectraRatio = 1.1;
  Float_t xTitleR = 25;
  Float_t xOffsetR = 4;
  Float_t yTitleR = 30;
  Float_t yOffsetR = 2;

  Float_t xLabelR = 25;
  Float_t yLabelR = 25;
  Float_t xLabelOffsetR = 0.025;
  Float_t yLabelOffsetR = 0.04;

  Float_t tickXR = 0.045;
  Float_t tickYR = 0.045;

  TF1 *rettaUno = new TF1("rettaUno", "pol0", 0, 8);
  rettaUno->FixParameter(0, 1);

  TString TitleYSpectraRatio = "Data/fit";
  canvasPtSpectra->cd();
  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 1000, 0, 8);
  SetFont(hDummyRatio);
  StyleHistoYield(hDummyRatio, LimInfMultRatio, LimSupMultRatio, 1, 1, TitleXPt, TitleYSpectraRatio, "", 1, 1.15, YoffsetSpectraRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetTickLength(hDummyRatio, tickXR, tickYR);
  padL1->Draw();
  padL1->cd();
  hDummyRatio->Draw("");

  for (Int_t m = numMult; m >= 0; m--)
  {
    for (Int_t b = 1; b <= fHistSpectrumStat[m]->GetNbinsX(); b++)
    {
      fHistSpectrumRatioFit[m]->SetBinContent(b, fHistSpectrumStat[m]->GetBinContent(b) * fHistSpectrumStat[m]->GetBinWidth(b) / fit_pwgfunc[m]->Integral(fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b), fHistSpectrumStat[m]->GetXaxis()->GetBinUpEdge(b)));
      fHistSpectrumRatioFit[m]->SetBinError(b, fHistSpectrumStat[m]->GetBinError(b) * fHistSpectrumStat[m]->GetBinWidth(b) / fit_pwgfunc[m]->Integral(fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b), fHistSpectrumStat[m]->GetXaxis()->GetBinUpEdge(b)));
    }
    rettaUno->SetLineColor(kBlack);
    fHistSpectrumRatioFit[m]->Draw("samee");
    rettaUno->Draw("same");

    SetFont(fHistSpectrumRatioFit[m]);
    StyleHistoYield(fHistSpectrumRatioFit[m], LimInfMultRatio, LimSupMultRatio, ColorMult[m], MarkerMult[m], TitleXPt, TitleYSpectraRatio, "", SizeMultRatio[m], 1.15, YoffsetSpectraRatio);
    SetHistoTextSize(fHistSpectrumRatioFit[m], xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
    SetTickLength(fHistSpectrumRatioFit[m], tickXR, tickYR);

    fHistSpectrumRatioFit[m]->Draw("same e0x0");

  } // end loop on mult

  TFile *outputFile = new TFile(stringout, "RECREATE");
  outputFile->cd();
  for (Int_t m = 0; m <= numMult; m++)
  {
    fit_pwgfunc[m]->Write();
  }
  outputFile->Close();
  cout << "Output file " << stringout << " has been created and written." << endl;
}
