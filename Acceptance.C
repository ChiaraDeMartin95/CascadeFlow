#include <Riostream.h>
#include "TProfile2D.h"
#include <string>
#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TLine.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TCutG.h>
#include "TFitResult.h"
#include "TLegend.h"
#include "CommonVarLambda.h"

void StyleCanvas(TCanvas *canvas, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(LMargin);
  canvas->SetRightMargin(RMargin);
  canvas->SetTopMargin(TMargin);
  canvas->SetBottomMargin(BMargin);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  // gStyle->SetPalette(55, 0);
}

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX,
                TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)
    histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset);
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

void Acceptance(Int_t indexMultTrial = 0,
                Int_t isTightAcceptance = 0,
                Int_t ChosenPart = ChosenParticle,
                TString inputFileName = SinputFileName,
                Bool_t isRapiditySel = ExtrisRapiditySel,
                Int_t EtaSysChoice = ExtrEtaSysChoice,
                Bool_t isSysMultTrial = ExtrisSysMultTrial)
{

  if (isSysMultTrial)
    inputFileName = SinputFileNameSyst;
  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trialsBDT)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;
  TString SBDT = "";
  if (BDTscoreCut != DefaultBDTscoreCut)
    SBDT = Form("_BDT%.3f", BDTscoreCut);

  TString SinputFile = "../OutputAnalysis/Output" + STHN[ExtrisFromTHN] + "_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice]; // + SBDT;
  if (isApplyWeights)
    SinputFile += "_Weighted";
  if (v2type == 1)
    SinputFile += "_SP";
  if (!useCommonBDTValue)
    SinputFile += "_BDTCentDep";
  if (isRun2Binning)
    SinputFile += "_Run2Binning";
  if (!isRapiditySel)
    SinputFile += "_Eta08";
  SinputFile += SBDT;
  if (isOOCentrality)
    SinputFile += "_isOOCentrality";
  if (isTightAcceptance == 1)
    SinputFile += "_TightAcceptance";
  else if (isTightAcceptance == 2)
    SinputFile += "_TighterAcceptance2";
  SinputFile += ".root";
  cout << "Input file: " << SinputFile << endl;
  TFile *inputFile = new TFile(SinputFile);
  if (!inputFile)
  {
    cout << "File not found" << endl;
    return;
  }

  TH3D *hEtaVsPtVsCos2ThetaLambdaFromC[commonNumCent + 1];
  TString hNameCos2ThetaLambdaFromC_Eta3D[commonNumCent + 1] = {""};

  TH2F *hEtaVsCos2ThetaLambdaFromC[commonNumCent + 1][numPtBinsLambda + 1];
  TString hNameEtaCos2ThetaLambdaFromC[commonNumCent + 1][numPtBinsLambda + 1] = {""};
  TH1F *hCos2ThetaLambdaFromCVsEta[commonNumCent + 1][numPtBinsLambda + 1];
  TString hNameCos2ThetaLambdaFromCVsEta[commonNumCent + 1][numPtBinsLambda + 1] = {""};
  TH1F *hEta[commonNumCent + 1][numPtBinsLambda + 1];
  TString hNameEta[commonNumCent + 1][numPtBinsLambda + 1] = {""};

  TH2F *hPtVsCos2ThetaLambdaFromC[commonNumCent + 1][numEtaBins + 1];
  TString hNamePtVsCos2ThetaLambdaFromC[commonNumCent + 1][numEtaBins + 1] = {""};
  TH1F *hCos2ThetaLambdaFromCVsPt[commonNumCent + 1][numEtaBins + 1];
  TString hNameCos2ThetaLambdaFromCVsPt[commonNumCent + 1][numEtaBins + 1] = {""};
  TH1F *hPtLambda[commonNumCent + 1][numEtaBins + 1];
  TString hNamePtLambda[commonNumCent + 1][numEtaBins + 1] = {""};

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;

  TProfile2D *pCos2ThetaLambdaFromC[commonNumCent + 1];
  TString pNameCos2ThetaLambdaFromC[commonNumCent + 1];
  TH2F *hCos2ThetaLambdaFromC2D[commonNumCent + 1];

  for (Int_t cent = 0; cent < commonNumCent + 1; cent++)
  {
    if (cent == numCent)
    { // 0-80%
      CentFT0CMin = 0;
      CentFT0CMax = 80;
    }
    else
    {
      CentFT0CMin = CentFT0C[cent];
      CentFT0CMax = CentFT0C[cent + 1];
    }
    if (isOOCentrality)
    {
      if (cent == numCentLambdaOO)
      {
        CentFT0CMin = 0;
        CentFT0CMax = 100;
      }
      else
      {
        CentFT0CMin = CentFT0CLambdaOO[cent];
        CentFT0CMax = CentFT0CLambdaOO[cent + 1];
      }
    }
    hNameCos2ThetaLambdaFromC_Eta3D[cent] = Form("etaVsPtVsCos2LambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax);
    hEtaVsPtVsCos2ThetaLambdaFromC[cent] = (TH3D *)inputFile->Get(hNameCos2ThetaLambdaFromC_Eta3D[cent]);
    pNameCos2ThetaLambdaFromC[cent] = Form("pCos2ThetaLambdaFromC_cent%i-%i", CentFT0CMin, CentFT0CMax);
    if (!hEtaVsPtVsCos2ThetaLambdaFromC[cent])
    {
      cout << "Histogram hEtaVsPtVsCos2ThetaLambdaFromC not available" << endl;
    }

    pCos2ThetaLambdaFromC[cent] = (TProfile2D *)hEtaVsPtVsCos2ThetaLambdaFromC[cent]->Project3DProfile("xy");
    pCos2ThetaLambdaFromC[cent]->SetName(pNameCos2ThetaLambdaFromC[cent]);
    hCos2ThetaLambdaFromC2D[cent] = new TH2F(Form("histoCos2ThetaLambdaFromCNoFit2D_cent%i-%i", CentFT0CMin, CentFT0CMax), Form("histoCos2ThetaLambdaFromCNoFit2D_cent%i-%i", CentFT0CMin, CentFT0CMax), numPtBinsLambda, PtBinsLambda, numEtaBins, EtaBins);

    for (Int_t pt = 0; pt < numPtBinsLambda + 1; pt++)
    {
      hNameEtaCos2ThetaLambdaFromC[cent][pt] = Form("EtaVsCos2ThetaLambdaFromC_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      hNameEta[cent][pt] = Form("Eta_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
      if (pt == numPtBinsLambda)
      {
        hNameCos2ThetaLambdaFromCVsEta[cent][pt] = Form("Cos2ThetaLambdaFromCVsEta_cent%i-%i", CentFT0CMin, CentFT0CMax);
        hEtaVsPtVsCos2ThetaLambdaFromC[cent]->GetYaxis()->SetRangeUser(PtBinsLambda[0] + 0.0001, PtBinsLambda[numPtBinsLambda] - 0.0001);
      }
      else
      {
        hNameCos2ThetaLambdaFromCVsEta[cent][pt] = Form("Cos2ThetaLambdaFromCVsEta_cent%i-%i_pt%i", CentFT0CMin, CentFT0CMax, pt);
        hEtaVsPtVsCos2ThetaLambdaFromC[cent]->GetYaxis()->SetRangeUser(PtBinsLambda[pt] + 0.0001, PtBinsLambda[pt + 1] - 0.0001);
      }
      hEtaVsCos2ThetaLambdaFromC[cent][pt] = (TH2F *)hEtaVsPtVsCos2ThetaLambdaFromC[cent]->Project3D("xze"); // eta vs cos2thetaLambdaFromC
      hEtaVsCos2ThetaLambdaFromC[cent][pt]->SetName(hNameEtaCos2ThetaLambdaFromC[cent][pt]);

      // hEta[cent][pt] = (TH1F *)hEtaVsPtVsCos2ThetaLambdaFromC[cent]->Project3D("xe");
      hEta[cent][pt] = new TH1F(hNameEta[cent][pt], hNameEta[cent][pt], numEtaBins, EtaBins);

      hCos2ThetaLambdaFromCVsEta[cent][pt] = (TH1F *)hEta[cent][pt]->Clone(hNameCos2ThetaLambdaFromCVsEta[cent][pt]);
      hCos2ThetaLambdaFromCVsEta[cent][pt]->Reset();
      for (Int_t bin = 0; bin < hEta[cent][pt]->GetNbinsX(); bin++)
      {
        TH1D *htemp = (TH1D *)hEtaVsCos2ThetaLambdaFromC[cent][pt]->ProjectionX(Form("_htemp_%i", bin), bin + 1, bin + 1);
        hCos2ThetaLambdaFromCVsEta[cent][pt]->SetBinContent(bin + 1, htemp->GetMean());
        hCos2ThetaLambdaFromCVsEta[cent][pt]->SetBinError(bin + 1, htemp->GetMeanError());
        hCos2ThetaLambdaFromC2D[cent]->SetBinContent(pt + 1, bin + 1, htemp->GetMean());
        hCos2ThetaLambdaFromC2D[cent]->SetBinError(pt + 1, bin + 1, htemp->GetMeanError());
      }
    }
    hEtaVsPtVsCos2ThetaLambdaFromC[cent]->GetYaxis()->SetRangeUser(PtBinsLambda[0] + 0.0001, PtBinsLambda[numPtBinsLambda] - 0.0001);
    for (Int_t eta = 0; eta < numEtaBins + 1; eta++)
    {
      hNamePtVsCos2ThetaLambdaFromC[cent][eta] = Form("PtVsCos2ThetaLambdaFromC_cent%i-%i_eta%i", CentFT0CMin, CentFT0CMax, eta);
      hNamePtLambda[cent][eta] = Form("PtLambda_cent%i-%i_eta%i", CentFT0CMin, CentFT0CMax, eta);

      if (eta == numEtaBins)
      {
        hEtaVsPtVsCos2ThetaLambdaFromC[cent]->GetXaxis()->SetRangeUser(EtaBins[0] + 0.0001, EtaBins[numEtaBins] - 0.0001);
        hNameCos2ThetaLambdaFromCVsPt[cent][eta] = Form("Cos2ThetaLambdaFromCVsPt_cent%i-%i", CentFT0CMin, CentFT0CMax);
      }
      else
      {
        hEtaVsPtVsCos2ThetaLambdaFromC[cent]->GetXaxis()->SetRangeUser(EtaBins[eta] + 0.0001, EtaBins[eta + 1] - 0.0001);
        hNameCos2ThetaLambdaFromCVsPt[cent][eta] = Form("Cos2ThetaLambdaFromCVsPt_cent%i-%i_eta%i", CentFT0CMin, CentFT0CMax, eta);
      }

      hPtVsCos2ThetaLambdaFromC[cent][eta] = (TH2F *)hEtaVsPtVsCos2ThetaLambdaFromC[cent]->Project3D("yze");
      hPtVsCos2ThetaLambdaFromC[cent][eta]->SetName(hNamePtVsCos2ThetaLambdaFromC[cent][eta]);

      // hPtLambda[cent][eta] = (TH1F *)hEtaVsPtVsCos2ThetaLambdaFromC[cent]->Project3D("ye");
      hPtLambda[cent][eta] = new TH1F(hNamePtLambda[cent][eta], hNamePtLambda[cent][eta], numPtBinsLambda, PtBinsLambda);

      hCos2ThetaLambdaFromCVsPt[cent][eta] = (TH1F *)hPtLambda[cent][eta]->Clone(hNameCos2ThetaLambdaFromCVsPt[cent][eta]);
      hCos2ThetaLambdaFromCVsPt[cent][eta]->Reset();

      for (Int_t bin = 0; bin < numPtBinsLambda + 1; bin++)
      {
        TH1D *htemp = (TH1D *)hPtVsCos2ThetaLambdaFromC[cent][eta]->ProjectionX(Form("_htemp_%i", bin), bin + 1, bin + 1);
        hCos2ThetaLambdaFromCVsPt[cent][eta]->SetBinContent(bin + 1, htemp->GetMean());
        hCos2ThetaLambdaFromCVsPt[cent][eta]->SetBinError(bin + 1, htemp->GetMeanError());
      }
    }
  }

  TCanvas *canvasAcceptance = new TCanvas("canvasAcceptance", "canvasAcceptance", 1200, 800);
  StyleCanvas(canvasAcceptance, 0.15, 0.1, 0.1, 0.1);
  TLegend *legendPt = new TLegend(0.6, 0.6, 0.9, 0.9);
  legendPt->SetBorderSize(0);
  legendPt->SetFillColor(0);
  legendPt->SetTextSize(0.03);
  TLegend *legendEta = new TLegend(0.6, 0.6, 0.9, 0.9);
  legendEta->SetBorderSize(0);
  legendEta->SetFillColor(0);
  legendEta->SetTextSize(0.03);
  canvasAcceptance->Divide(2, 2);
  // cos2 vs eta for different pt bins
  canvasAcceptance->cd(1);
  for (Int_t pt = 0; pt < numPtBinsLambda + 1; pt++)
  {
    Int_t cent = 6;
    StyleHisto(hCos2ThetaLambdaFromCVsEta[cent][pt], 0, 0.5, ColorMult[pt], 20, "#eta_{#Lambda}", "cos^{2}(#theta_{p})", "", 0, -0.8, 0.8, 1.2, 1., 0.7);
    hCos2ThetaLambdaFromCVsEta[cent][pt]->GetYaxis()->SetRangeUser(0, 0.5);
    hCos2ThetaLambdaFromCVsEta[cent][pt]->SetMarkerStyle(20);
    hCos2ThetaLambdaFromCVsEta[cent][pt]->SetMarkerSize(0.5);
    hCos2ThetaLambdaFromCVsEta[cent][pt]->SetMarkerColor(ColorMult[pt]);
    hCos2ThetaLambdaFromCVsEta[cent][pt]->SetLineColor(ColorMult[pt]);
    hCos2ThetaLambdaFromCVsEta[cent][pt]->Draw("same");
    legendPt->AddEntry(hCos2ThetaLambdaFromCVsEta[cent][pt], Form("%.1f < p_{T} = %.1f GeV/c", PtBinsLambda[pt], PtBinsLambda[pt + 1]), "pl");
  }
  legendPt->Draw("same");
  canvasAcceptance->cd(2);
  for (Int_t eta = 0; eta < numEtaBins + 1; eta++)
  {
    Int_t cent = 6;
    StyleHisto(hCos2ThetaLambdaFromCVsPt[cent][eta], 0, 0.5, ColorMult[eta], 20, "p_{T} (GeV/c)", "cos^{2}(#theta_{p})", "", 0, 0.5, 2.5, 1.2, 1., 0.7);
    hCos2ThetaLambdaFromCVsPt[cent][eta]->GetXaxis()->SetRangeUser(hPtLambda[cent][eta]->GetXaxis()->GetBinLowEdge(1), hPtLambda[cent][eta]->GetXaxis()->GetBinUpEdge(hPtLambda[cent][eta]->GetNbinsX()));
    hCos2ThetaLambdaFromCVsPt[cent][eta]->GetYaxis()->SetRangeUser(0, 0.5);
    hCos2ThetaLambdaFromCVsPt[cent][eta]->SetMarkerStyle(20);
    hCos2ThetaLambdaFromCVsPt[cent][eta]->SetMarkerSize(0.5);
    hCos2ThetaLambdaFromCVsPt[cent][eta]->SetMarkerColor(ColorMult[eta]);
    hCos2ThetaLambdaFromCVsPt[cent][eta]->SetLineColor(ColorMult[eta]);
    hCos2ThetaLambdaFromCVsPt[cent][eta]->Draw("same");
    legendEta->AddEntry(hCos2ThetaLambdaFromCVsPt[cent][eta], Form("%.1f < #eta = %.1f", EtaBins[eta], EtaBins[eta + 1]), "pl");
  }
  legendEta->Draw("same");
  canvasAcceptance->cd(3);
  for (Int_t pt = 0; pt < numPtBinsLambda + 1; pt++)
  {
    Int_t cent = 6;
    StyleHisto(hEta[cent][pt], 0, 0.5, ColorMult[pt], 20, "#eta_{#Lambda}", "dN/d#eta", "", 0, -0.8, 0.8, 1.2, 1., 0.7);
    hEta[cent][pt]->GetYaxis()->SetRangeUser(0, 0.5);
    hEta[cent][pt]->SetMarkerStyle(20);
    hEta[cent][pt]->SetMarkerSize(0.5);
    hEta[cent][pt]->SetMarkerColor(ColorMult[pt]);
    hEta[cent][pt]->SetLineColor(ColorMult[pt]);
    hEta[cent][pt]->Draw("same");
  }
  canvasAcceptance->cd(4);
  for (Int_t eta = 0; eta < numEtaBins + 1; eta++)
  {
    Int_t cent = 6;
    StyleHisto(hPtLambda[cent][eta], 0, 0.5, ColorMult[eta], 20, "p_{T} (GeV/c)", "dN/dp_{T}", "", 0, 0, 20, 1.2, 1., 0.7);
    hPtLambda[cent][eta]->GetXaxis()->SetRangeUser(hPtLambda[cent][eta]->GetXaxis()->GetBinLowEdge(1), hPtLambda[cent][eta]->GetXaxis()->GetBinUpEdge(hPtLambda[cent][eta]->GetNbinsX()));
    hPtLambda[cent][eta]->GetYaxis()->SetRangeUser(0, 0.5);
    hPtLambda[cent][eta]->SetMarkerStyle(20);
    hPtLambda[cent][eta]->SetMarkerSize(0.5);
    hPtLambda[cent][eta]->SetMarkerColor(ColorMult[eta]);
    hPtLambda[cent][eta]->SetLineColor(ColorMult[eta]);
    hPtLambda[cent][eta]->Draw("same");
  }

  TString SOutputFile = "../AcceptancePlots/Acceptance_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (isApplyWeights)
    SOutputFile += "_Weighted";
  if (v2type == 1)
    SOutputFile += "_SP";
  if (!useCommonBDTValue)
    SOutputFile += "_BDTCentDep";
  if (isRun2Binning)
    SOutputFile += "_Run2Binning";
  if (ExtrisApplyEffWeights)
  {
    SOutputFile += "_EffW";
  }
  SOutputFile += "_WithAlpha";
  // if (!isRapiditySel || ExtrisFromTHN)
  if (!isRapiditySel)
    SOutputFile += "_Eta08";
  SOutputFile += STHN[ExtrisFromTHN];
  if (isOOCentrality)
    SOutputFile += "_isOOCentrality";
  if (isTightAcceptance == 1)
    SOutputFile += "_TightAcceptance";
  else if (isTightAcceptance == 2)
    SOutputFile += "_TightAcceptance2";

  SOutputFile += ".root";
  TFile *file = new TFile(SOutputFile, "RECREATE");
  TList *listAcceptance = new TList();
  for (Int_t cent = 0; cent < commonNumCent + 1; cent++)
  {
    listAcceptance->Add(hCos2ThetaLambdaFromC2D[cent]);
    // pCos2ThetaLambdaFromC[cent]->Write();
    hCos2ThetaLambdaFromC2D[cent]->Write();
    hCos2ThetaLambdaFromCVsEta[cent][numPtBinsLambda]->Write();
    hCos2ThetaLambdaFromCVsPt[cent][numEtaBins]->Write();
  }
  listAcceptance->Write("ccdb_object", TObject::kSingleKey);
  file->Close();
  cout << "I produced the file: " << SOutputFile << endl;
}
