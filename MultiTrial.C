// This macro was originally written by:
// francesca.ercolessi@cern.ch
// and later modified by:
// chiara.de.martin@cern.ch

#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TPad.h"
#include "StyleFile.h"
#include "CommonVar.h"

double ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr)
{
  // Error in a Ratio
  if (B != 0)
  {
    Double_t errorfromtop = Aerr * Aerr / (B * B);
    Double_t errorfrombottom = ((A * A) / (B * B * B * B)) * Berr * Berr;
    return TMath::Sqrt(TMath::Abs(errorfromtop - errorfrombottom));
  }
  return 1.;
}

void DivideAndComputeRogerBarlow(TH1F *hVar, TH1F *hDef)
{
  // Use Roger Barlow "sigma_{delta}" as errors for ratios
  Double_t lhVarNBins = hVar->GetNbinsX();
  Double_t lhDefNBins = hDef->GetNbinsX();

  if (lhVarNBins != lhDefNBins)
  {
    cout << "Problem! Number of bins doesn't match! " << endl;
    return;
  }

  Double_t lSigmaDelta[100];
  for (Int_t i = 1; i < hVar->GetNbinsX() + 1; i++)
  {
    // Computation of roger barlow sigma_{delta}
    lSigmaDelta[i] = TMath::Sqrt(TMath::Abs(TMath::Power(hVar->GetBinError(i), 2) - TMath::Power(hDef->GetBinError(i), 2)));
    // Computation of relationship to hDef for plotting in ratio plot
    if (abs(hDef->GetBinContent(i)) > 1e-12)
      lSigmaDelta[i] /= hDef->GetBinContent(i);
    else
      lSigmaDelta[i] = 0;
  }
  // Regular Division
  hVar->Divide(hDef);
  // Replace Errors
  for (Int_t i = 1; i < hVar->GetNbinsX() + 1; i++)
  {
    hVar->SetBinError(i, abs(lSigmaDelta[i]));
  }
}

double PassRogerBarlowCriterion(int nsigmas, Double_t dev, Double_t RBsigma)
{

  if (TMath::Abs(dev) > (nsigmas * RBsigma))
  {
    return dev;
  }
  else
  {
    return 0.;
  }
}

TH1F *makeSystPlots(int num = 1, TString Sdef = "", TString Svaried = "", TString histoName = "")
{

  TFile *fdef = TFile::Open(Sdef);
  TFile *fvaried = TFile::Open(Svaried);

  TH1F *hVariedCut = (TH1F *)fvaried->Get(histoName);
  hVariedCut->SetName(Form("hVarCutSet%i", num));
  TH1F *hDefault = (TH1F *)fdef->Get(histoName);
  hDefault->SetName(Form("hDefault"));

  TH1F *hDev = (TH1F *)hDefault->Clone("hDev");
  hDev->Reset();

  DivideAndComputeRogerBarlow(hVariedCut, hDefault);
  // assigns to hVariedCut: hR = hvar/hdef as bin content, sigmaBarlow (sB) as error

  for (int i = 1; i <= hDev->GetNbinsX(); i++) // pt bins
  {
    double dev = hVariedCut->GetBinContent(i) - 1; // hR - 1
    double err = hVariedCut->GetBinError(i);       // sB

    hDev->SetBinContent(i, abs(PassRogerBarlowCriterion(nsigmaBarlow, dev, err))); // rel. syst. error = hR-1 if |hR-1| > 1*sB
  }

  hDev->GetYaxis()->SetRangeUser(0., 0.5);

  // histo with rel. syst. error (if variation passes the Barlow cut)
  return hDev;
}

TH1F *makeNSigmaBarlowPlots(int num = 1, TString Sdef = "", TString Svaried = "", TString histoName = "")
{

  TFile *fdef = TFile::Open(Sdef);
  TFile *fvaried = TFile::Open(Svaried);

  TH1F *hVariedCut = (TH1F *)fvaried->Get(histoName);
  hVariedCut->SetName(Form("hVarCutSet%i", num));
  TH1F *hDefault = (TH1F *)fdef->Get(histoName);
  hDefault->SetName(Form("hDefault"));

  TH1F *hNSigmaBarlow = (TH1F *)hDefault->Clone("hNSigmaBarlow");
  hNSigmaBarlow->Reset();

  DivideAndComputeRogerBarlow(hVariedCut, hDefault);
  // assigns to hVariedCut: hR = hvar/hdef as bin content, sigmaBarlow (sB) as error

  for (int i = 1; i <= hNSigmaBarlow->GetNbinsX(); i++) // pt bins
  {
    double dev = hVariedCut->GetBinContent(i) - 1;   // hR - 1
    double err = hVariedCut->GetBinError(i);         // sB
    hNSigmaBarlow->SetBinContent(i, abs(dev) / err); // rel. syst. error = hR-1 if |hR-1| > 1*sB
  }

  hNSigmaBarlow->GetYaxis()->SetRangeUser(0, 5);

  // histo with NSigmaBarlow
  return hNSigmaBarlow;
}

Int_t IndexNotDisplayed = 15;
TString SisPtIntegrated[2] = {"", "PtInt"};

void MultiTrial(
    Int_t mul = 0,
    Int_t Choice = 0,          // 0 = V2Mixed, 1 = Pz(s2)Mixed, 2 = Pz(s2)LambdaFromCMixed
    Bool_t isPtAnalysis = 1,   // 1 for V2 vs pt and Pzs2 vs pt, 0 for Pz vs 2(phi-Psi)
    Bool_t isPtIntegrated = 1, // 1 for results integrated in pt / phi
    TString SisSyst = "BDT",
    Int_t ChosenPart = ChosenParticle,
    Bool_t isRapiditySel = ExtrisRapiditySel,
    Int_t BkgType = ExtrBkgType,
    TString inputFileName = SinputFileNameSyst,
    Bool_t UseTwoGauss = ExtrUseTwoGauss)
{
  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  if (mul == numCent)
  { // 0-80%
    CentFT0CMin = 0;
    CentFT0CMax = 80;
  }
  else
  {
    CentFT0CMin = CentFT0C[mul];
    CentFT0CMax = CentFT0C[mul + 1];
  }

  // histoName
  Bool_t isPolFromLambda = 0;
  TString TypeHisto = "";
  if (Choice == 0)
    TypeHisto = "V2";
  else if (Choice == 1)
  {
    TypeHisto = "Pzs";
    if (!isPtAnalysis)
      TypeHisto = "Pz";
  }
  else if (Choice == 2)
  {
    isPolFromLambda = 1;
    TypeHisto = "Pzs2LambdaFromC";
    if (!isPtAnalysis)
      TypeHisto = "PzLambdaFromC";
  }
  TypeHisto += SisPtIntegrated[isPtIntegrated];
  TypeHisto += "Mixed";

  TString histoName = "histo" + TypeHisto;

  TString Suffix = inputFileName + Form("_%i-%i_", CentFT0C[mul], CentFT0C[mul + 1]) + ParticleName[ChosenPart] + "_" + SisSyst + "_" + SisPtIntegrated[isPtIntegrated];

  Int_t trials = 0;
  if (SisSyst == "BDT")
    trials = trialsBDT;
  else if (SisSyst == "eta")
    trials = 2; // eta > 0 and eta < 0
  else if (SisSyst == "IR")
    trials = 5; // different interaction rates
  TString Sdef = "OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_" + inputFileName + "_" + ParticleName[ChosenPart];
  Sdef += IsOneOrTwoGauss[UseTwoGauss];
  Sdef += SIsBkgParab[BkgType];
  Sdef += Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);
  if (isApplyWeights)
    Sdef += "_Weighted";
  if (v2type == 1)
    Sdef += "_SP";
  if (!useCommonBDTValue)
    Sdef += "_BDTCentDep";
  if (!isPtAnalysis)
    Sdef += "_vsPsi";
  if (!isV2 && isPolFromLambda)
    Sdef += "_PolFromLambda";
  if (ExtrisApplyEffWeights)
  {
    Sdef += "_EffW";
  }
  // if (!ExtrisFromTHN) Sdef += "_WithAlpha";
  TString Sdef1 = "";
  if (!isRapiditySel || ExtrisFromTHN)
    Sdef1 = "_Eta08";
  Sdef1 += STHN[ExtrisFromTHN];
  TString SdefFinal = Sdef + Sdef1;

  TString Svaried = "";

  cout << "Default input file: " << SdefFinal << endl;
  TFile *fdef = TFile::Open(SdefFinal + ".root");
  if (!fdef)
  {
    cout << "File not found: " << SdefFinal << endl;
    return;
  }

  // v2 plots
  TString Titlecv2 = "v2";
  if (Choice == 1)
  {
    Titlecv2 = "Pzs";
    if (!isPtAnalysis)
      Titlecv2 = "Pz";
  }
  else if (Choice == 2)
  {
    Titlecv2 = "Pzs2LambdaFromC";
    if (!isPtAnalysis)
      Titlecv2 = "PzLambdaFromC";
  }
  TCanvas *cv2 = new TCanvas("cv2", Titlecv2, 1000, 800);
  StyleCanvas(cv2, 0.15, 0.05, 0.05, 0.15);
  cv2->cd();
  cout << "histoName: " << histoName << endl;
  TH1F *hDefault = (TH1F *)fdef->Get(histoName);
  if (!hDefault)
  {
    cout << "histo not found in " << SdefFinal << endl;
    return;
  }
  hDefault->SetName("hDefault");
  hDefault->SetLineColor(kBlack);
  hDefault->GetYaxis()->SetTitle(Titlecv2);
  hDefault->Draw();

  TH1F *hDefaultError = (TH1F *)hDefault->Clone("hError");
  hDefaultError->Reset();
  for (int i = 1; i <= hDefault->GetNbinsX(); i++)
  {
    hDefaultError->SetBinContent(i, hDefault->GetBinError(i));
  }

  TH1F *hDefaultYield = (TH1F *)fdef->Get("histoYield" + SisPtIntegrated[isPtIntegrated]);
  if (!hDefaultYield)
  {
    cout << "histoYield not found in " << SdefFinal << endl;
    return;
  }
  hDefaultYield->SetName("hDefaultYield");
  TH1F *hDefaultPurity = (TH1F *)fdef->Get("histoPurity" + SisPtIntegrated[isPtIntegrated]);
  if (!hDefaultPurity)
  {
    cout << "histoPurity not found in " << SdefFinal << endl;
    return;
  }
  hDefaultPurity->SetName("hDefaultPurity");
  TH1F *hDefaultSignificance = (TH1F *)fdef->Get("histoSignificance" + SisPtIntegrated[isPtIntegrated]);
  if (!hDefaultSignificance)
  {
    cout << "histoSignificance not found in " << SdefFinal << endl;
    return;
  }
  hDefaultSignificance->SetName("hDefaultSignificance");

  TLegend *legTrial;
  if (SisSyst == "IR")
    legTrial = new TLegend(0.66, 0.7, 0.96, 0.9);
  else
    legTrial = new TLegend(0.66, 0.2, 0.96, 0.5);
  legTrial->SetBorderSize(0);
  legTrial->SetFillStyle(0);
  legTrial->SetTextSize(0.03);
  legTrial->SetTextFont(42);
  if (SisSyst == "BDT")
    legTrial->AddEntry(hDefault, Form("BDT score > %.3f", DefaultBDTscoreCut), "pl");
  else if (SisSyst == "eta")
    legTrial->AddEntry(hDefault, "|#eta| < 0.8", "pl");
  else if (SisSyst == "IR")
    legTrial->AddEntry(hDefault, "PbPb 2023", "pl");

  TH1F *h[trials];
  TH1F *hAbsoluteSyst[trials];
  TH1F *hNSigmaBarlow[trials];
  TH1F *hRatio[trials];
  TH1F *hError[trials];
  TH1F *hErrorRatio[trials];
  TH1F *hRawYield[trials];
  TH1F *hRawYieldRatio[trials];
  TH1F *hRawYieldRatioToTighter[trials];
  TH1F *hPurity[trials];
  TH1F *hPurityRatio[trials];
  TH1F *hSignificance[trials];
  TH1F *hSignificanceRatio[trials];
  Float_t BDTscoreCut = 0;

  for (int i = 0; i < trials; i++)
  {
    cout << "Trial: " << i << endl;
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trials * i;
    TString SBDT = Form("_BDT%.3f", BDTscoreCut);
    if (SisSyst == "BDT")
      Svaried = Sdef + SBDT + Sdef1;
    else if (SisSyst == "eta")
      Svaried = Sdef + SEtaSysChoice[i + 1];
    else if (SisSyst == "IR")
    {
      Svaried = "OutputAnalysis/FitV2_" + inputFileNameIR + SIRChoice[i + 1] + "_" + ParticleName[ChosenPart];
      Svaried += IsOneOrTwoGauss[UseTwoGauss];
      Svaried += SIsBkgParab[BkgType];
      Svaried += Form("_Cent%i-%i", CentFT0C[mul], CentFT0C[mul + 1]);
    }

    cout << "InputFile - variation: " << Svaried << endl;

    TFile *fvaried = TFile::Open(Svaried + ".root");
    if (!fvaried)
    {
      cout << "File not found: " << Svaried << endl;
      return;
    }
    TH1F *hVariedCut = (TH1F *)fvaried->Get(histoName);
    if (!hVariedCut)
    {
      cout << "histo not found in " << Svaried << endl;
      return;
    }

    hVariedCut->SetName(histoName + Form("_%i", i));
    hRatio[i] = (TH1F *)hVariedCut->Clone(histoName + Form("Ratio_%i", i));
    hRatio[i]->Divide(hDefault);

    hRawYield[i] = (TH1F *)fvaried->Get("histoYield" + SisPtIntegrated[isPtIntegrated]);
    if (!hRawYield[i])
    {
      cout << "histoYield not found in " << Svaried << endl;
      return;
    }
    hRawYield[i]->SetName(Form("histoYield_%i", i));
    hRawYieldRatio[i] = (TH1F *)hRawYield[i]->Clone(Form("histoYieldRatio_%i", i));
    hRawYieldRatio[i]->Divide(hDefaultYield);

    hPurity[i] = (TH1F *)fvaried->Get("histoPurity" + SisPtIntegrated[isPtIntegrated]);
    if (!hPurity[i])
    {
      cout << "histoPurity not found in " << Svaried << endl;
      return;
    }
    hPurity[i]->SetName(Form("histoPurity_%i", i));
    hPurityRatio[i] = (TH1F *)hPurity[i]->Clone(Form("histoPurity_%i", i));
    hPurityRatio[i]->Divide(hDefaultPurity);

    hSignificance[i] = (TH1F *)fvaried->Get("histoSignificance" + SisPtIntegrated[isPtIntegrated]);
    if (!hSignificance[i])
    {
      cout << "histoSignificance not found in " << Svaried << endl;
      return;
    }
    hSignificance[i]->SetName(Form("histoSignificance_%i", i));
    hSignificanceRatio[i] = (TH1F *)hSignificance[i]->Clone(Form("histoSignificanceRatio_%i", i));
    hSignificanceRatio[i]->Divide(hDefaultSignificance);

    hError[i] = (TH1F *)hVariedCut->Clone(Form("hError_%i", i));
    hError[i]->Reset();
    for (int j = 1; j <= hVariedCut->GetNbinsX(); j++)
    {
      hError[i]->SetBinContent(j, hVariedCut->GetBinError(j));
    }
    hErrorRatio[i] = (TH1F *)hError[i]->Clone(Form("hErrorRatio_%i", i));
    hErrorRatio[i]->Divide(hDefaultError);

    // v2 plots
    cv2->cd();
    hVariedCut->SetLineColor(ColorMult[i]);
    hVariedCut->SetMarkerColor(ColorMult[i]);
    hVariedCut->SetMarkerStyle(MarkerMult[i]);
    if (SisSyst == "BDT")
      legTrial->AddEntry(hVariedCut, Form("BDT score > %.3f", BDTscoreCut), "pl");
    else if (SisSyst == "eta")
      legTrial->AddEntry(hVariedCut, SEtaSysChoice[i + 1], "pl");
    else if (SisSyst == "IR")
      legTrial->AddEntry(hVariedCut, SIRValue[i + 1], "pl");
    hVariedCut->Draw("same");

    h[i] = makeSystPlots(i + 1, SdefFinal + ".root", Svaried + ".root", histoName);
    hAbsoluteSyst[i] = (TH1F *)h[i]->Clone(Form("hAbsoluteSyst_%i", i));
    hAbsoluteSyst[i]->Multiply(hDefault);
    hNSigmaBarlow[i] = makeNSigmaBarlowPlots(i + 1, SdefFinal + ".root", Svaried + ".root", histoName);
  }

  for (Int_t i = 0; i < trials; i++)
  {
    hRawYieldRatioToTighter[i] = (TH1F *)hRawYield[i]->Clone(Form("histoYieldRatioToTighter_%i", i));
    if (i == (trials - 1))
      continue;
    else
      hRawYieldRatioToTighter[i]->Divide(hRawYield[i + 1]);
  }

  legTrial->Draw();
  if (isV2)
  {
    cv2->SaveAs("Systematics/V2MultTrial" + Suffix + ".pdf");
    cv2->SaveAs("Systematics/V2MultTrial" + Suffix + ".png");
  }
  else
  {
    cv2->SaveAs("Systematics/Pzs2MultTrial" + Suffix + ".pdf");
    cv2->SaveAs("Systematics/Pzs2MultTrial" + Suffix + ".png");
  }

  TCanvas *cNSigmaBarlow = new TCanvas("cNSigmaBarlow", "cNSigmaBarlow", 1000, 800);
  StyleCanvas(cNSigmaBarlow, 0.15, 0.05, 0.05, 0.15);
  cNSigmaBarlow->cd();
  TH1F *hdummy = (TH1F *)hDefault->Clone("hdummy");
  hdummy->Reset();
  hdummy->GetYaxis()->SetRangeUser(0, 5.);
  hdummy->GetYaxis()->SetTitle("N_{#sigma}^{Barlow}");
  hdummy->Draw();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hNSigmaBarlow[i]->SetLineColor(ColorMult[i]);
    hNSigmaBarlow[i]->SetLineWidth(2);
    hNSigmaBarlow[i]->SetMarkerColor(ColorMult[i]);
    hNSigmaBarlow[i]->SetMarkerStyle(MarkerMult[i]);
    hNSigmaBarlow[i]->Draw("same");
  }

  TCanvas *cDev = new TCanvas("cDev", "cDev", 1000, 800);
  StyleCanvas(cDev, 0.15, 0.05, 0.05, 0.15);
  cDev->cd();
  for (Int_t i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hAbsoluteSyst[i]->SetLineColor(ColorMult[i]);
    hAbsoluteSyst[i]->SetLineWidth(2);
    hAbsoluteSyst[i]->SetMarkerColor(ColorMult[i]);
    hAbsoluteSyst[i]->SetMarkerStyle(MarkerMult[i]);
    hAbsoluteSyst[i]->GetYaxis()->SetTitle("Absolute deviation from default");
    hAbsoluteSyst[i]->GetYaxis()->SetRangeUser(0, 0.5);
    hAbsoluteSyst[i]->Draw("same");
  }

  // ratio of varied v2 to default one
  TCanvas *cv2Ratio = new TCanvas("cv2Ratio", "cv2Ratio", 1000, 800);
  StyleCanvas(cv2Ratio, 0.15, 0.05, 0.05, 0.15);
  cv2Ratio->cd();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hRatio[i]->SetLineColor(ColorMult[i]);
    hRatio[i]->SetMarkerColor(ColorMult[i]);
    hRatio[i]->SetMarkerStyle(MarkerMult[i]);
    hRatio[i]->SetTitle("");
    hRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hRatio[i]->GetYaxis()->SetRangeUser(0, 2);
    hRatio[i]->Draw("same");
  }
  TF1 *lineat1 = new TF1("lineat1", "1", 0, 8);
  lineat1->SetLineColor(kBlack);
  lineat1->SetLineStyle(7);
  legTrial->Draw();
  lineat1->Draw("same");
  if (isV2)
  {
    cv2Ratio->SaveAs("Systematics/v2RatioMultTrial" + Suffix + ".pdf");
    cv2Ratio->SaveAs("Systematics/v2RatioMultTrial" + Suffix + ".png");
  }
  else
  {
    cv2Ratio->SaveAs("Systematics/Pzs2RatioMultTrial" + Suffix + ".pdf");
    cv2Ratio->SaveAs("Systematics/Pzs2RatioMultTrial" + Suffix + ".png");
  }

  // raw yields plots
  TCanvas *cYield = new TCanvas("cYield", "cYield", 1000, 800);
  StyleCanvas(cYield, 0.15, 0.05, 0.05, 0.15);
  cYield->cd();
  hDefaultYield->SetTitle("");
  hDefaultYield->SetLineColor(kBlack);
  hDefaultYield->GetYaxis()->SetRangeUser(0, 1.8 * hDefaultYield->GetMaximum());
  hDefaultYield->Draw();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hRawYield[i]->SetLineColor(ColorMult[i]);
    hRawYield[i]->SetMarkerColor(ColorMult[i]);
    hRawYield[i]->SetMarkerStyle(MarkerMult[i]);
    hRawYield[i]->Draw("same");
  }
  legTrial->Draw();
  cYield->SaveAs("Systematics/YieldMultTrial" + Suffix + ".pdf");
  cYield->SaveAs("Systematics/YieldMultTrial" + Suffix + ".png");

  // ratio of raw yields to default one
  TCanvas *cYieldRatio = new TCanvas("cYieldRatio", "cYieldRatio", 1000, 800);
  StyleCanvas(cYieldRatio, 0.15, 0.05, 0.05, 0.15);
  cYieldRatio->cd();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hRawYieldRatio[i]->SetLineColor(ColorMult[i]);
    hRawYieldRatio[i]->SetMarkerColor(ColorMult[i]);
    hRawYieldRatio[i]->SetMarkerStyle(MarkerMult[i]);
    hRawYieldRatio[i]->SetTitle("");
    hRawYieldRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hRawYieldRatio[i]->GetYaxis()->SetRangeUser(0, 2.5);
    hRawYieldRatio[i]->Draw("same");
  }
  legTrial->Draw();
  lineat1->Draw("same");
  cYieldRatio->SaveAs("Systematics/YieldRatioMultTrial" + Suffix + ".pdf");
  cYieldRatio->SaveAs("Systematics/YieldRatioMultTrial" + Suffix + ".png");

  // ratio of raw yields to previous one
  TCanvas *cYieldRatioToTighter = new TCanvas("cYieldRatioToTighter", "cYieldRatioToTighter", 1000, 800);
  StyleCanvas(cYieldRatioToTighter, 0.15, 0.05, 0.05, 0.15);
  cYieldRatioToTighter->cd();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hRawYieldRatioToTighter[i]->SetLineColor(ColorMult[i]);
    hRawYieldRatioToTighter[i]->SetMarkerColor(ColorMult[i]);
    hRawYieldRatioToTighter[i]->SetMarkerStyle(MarkerMult[i]);
    hRawYieldRatioToTighter[i]->SetTitle("");
    hRawYieldRatioToTighter[i]->GetYaxis()->SetTitle("Ratio to tighter cut");
    hRawYieldRatioToTighter[i]->GetYaxis()->SetRangeUser(0.9, 1.5);
    hRawYieldRatioToTighter[i]->Draw("same");
  }
  legTrial->Draw();
  lineat1->Draw("same");
  cYieldRatioToTighter->SaveAs("Systematics/YieldRatioToTighterMultTrial" + Suffix + ".pdf");
  cYieldRatioToTighter->SaveAs("Systematics/YieldRatioToTighterMultTrial" + Suffix + ".png");

  // purity plots
  TCanvas *cPurity = new TCanvas("cPurity", "cPurity", 1000, 800);
  StyleCanvas(cPurity, 0.15, 0.05, 0.05, 0.15);
  cPurity->cd();
  hDefaultPurity->SetTitle("");
  hDefaultPurity->SetLineColor(kBlack);
  hDefaultPurity->GetYaxis()->SetRangeUser(0, 1);
  hDefaultPurity->Draw();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hPurity[i]->SetLineColor(ColorMult[i]);
    hPurity[i]->SetMarkerColor(ColorMult[i]);
    hPurity[i]->SetMarkerStyle(MarkerMult[i]);
    hPurity[i]->Draw("same");
  }
  legTrial->Draw();
  cPurity->SaveAs("Systematics/PurityMultTrial" + inputFileName + Suffix + ".pdf");
  cPurity->SaveAs("Systematics/PurityMultTrial" + inputFileName + Suffix + ".png");

  // ratio of purities to default one
  TCanvas *cPurityRatio = new TCanvas("cPurityRatio", "cPurityRatio", 1000, 800);
  StyleCanvas(cPurityRatio, 0.15, 0.05, 0.05, 0.15);
  cPurityRatio->cd();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hPurityRatio[i]->SetLineColor(ColorMult[i]);
    hPurityRatio[i]->SetMarkerColor(ColorMult[i]);
    hPurityRatio[i]->SetMarkerStyle(MarkerMult[i]);
    hPurityRatio[i]->SetTitle("");
    hPurityRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hPurityRatio[i]->GetYaxis()->SetRangeUser(0, 1.1);
    hPurityRatio[i]->Draw("same");
  }
  lineat1->Draw("same");
  legTrial->Draw();
  cPurityRatio->SaveAs("Systematics/PurityRatioMultTrial" + Suffix + ".pdf");
  cPurityRatio->SaveAs("Systematics/PurityRatioMultTrial" + Suffix + ".png");

  // significance plots
  TCanvas *cSignificance = new TCanvas("cSignificance", "cSignificance", 1000, 800);
  StyleCanvas(cSignificance, 0.15, 0.05, 0.05, 0.15);
  cSignificance->cd();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hSignificance[i]->SetLineColor(ColorMult[i]);
    hSignificance[i]->SetMarkerColor(ColorMult[i]);
    hSignificance[i]->SetMarkerStyle(MarkerMult[i]);
    hSignificance[i]->SetTitle("");
    hSignificance[i]->Draw("same");
  }
  legTrial->Draw();
  cSignificance->SaveAs("Systematics/SignificanceMultTrial" + Suffix + ".pdf");
  cSignificance->SaveAs("Systematics/SignificanceMultTrial" + Suffix + ".png");

  // ratio of significances to default one
  TCanvas *cSignificanceRatio = new TCanvas("cSignificanceRatio", "cSignificanceRatio", 1000, 800);
  StyleCanvas(cSignificanceRatio, 0.15, 0.05, 0.05, 0.15);
  cSignificanceRatio->cd();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hSignificanceRatio[i]->SetLineColor(ColorMult[i]);
    hSignificanceRatio[i]->SetMarkerColor(ColorMult[i]);
    hSignificanceRatio[i]->SetMarkerStyle(MarkerMult[i]);
    hSignificanceRatio[i]->SetTitle("");
    hSignificanceRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hSignificanceRatio[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
    hSignificanceRatio[i]->Draw("same");
  }
  legTrial->Draw();
  lineat1->Draw("same");
  cSignificanceRatio->SaveAs("Systematics/SignificanceRatioMultTrial" + Suffix + ".pdf");
  cSignificanceRatio->SaveAs("Systematics/SignificanceRatioMultTrial" + Suffix + ".png");

  // systematic computation taking max variation wrt to default
  const int bins = h[0]->GetNbinsX(); // pt bins
  TH1F *hMaxDev = (TH1F *)h[0]->Clone("hMaxDev");
  hMaxDev->Reset();
  TH1F *hAbsoluteMaxDev = (TH1F *)hAbsoluteSyst[0]->Clone("hAbsoluteMaxDev");
  hAbsoluteMaxDev->Reset();
  for (int pt = 0; pt < bins; pt++) // loop over pT bins
  {
    // cout << "check " << hMaxDev->GetBinContent(pt + 1) << " " << hAbsoluteMaxDev->GetBinContent(pt + 1) << endl;
    for (int i = 0; i < trials; i++)
    {
      if (i == IndexNotDisplayed)
        continue;
      if (TMath::Abs(h[i]->GetBinContent(pt + 1)) > TMath::Abs(hMaxDev->GetBinContent(pt + 1)))
      {
        hMaxDev->SetBinContent(pt + 1, h[i]->GetBinContent(pt + 1));
      }
      hMaxDev->SetBinError(pt + 1, 0.);
      // cout << "trial number: " << i << endl;
      // cout << "check " << hAbsoluteSyst[i]->GetBinContent(pt + 1) << " " << hAbsoluteMaxDev->GetBinContent(pt + 1) << endl;
      if (TMath::Abs(hAbsoluteSyst[i]->GetBinContent(pt + 1)) > TMath::Abs(hAbsoluteMaxDev->GetBinContent(pt + 1)))
      {
        hAbsoluteMaxDev->SetBinContent(pt + 1, TMath::Abs(hAbsoluteSyst[i]->GetBinContent(pt + 1)));
      }
      // cout << "check " << hAbsoluteSyst[i]->GetBinContent(pt + 1) << " " << hAbsoluteMaxDev->GetBinContent(pt + 1) << endl;
      hAbsoluteMaxDev->SetBinError(pt + 1, 0.);
    }
  }

  //  gaussian distributions + fits
  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1200);
  c2->Divide(4, 4);

  TH1F *hPtDev[bins];
  TF1 *fgaus[bins];

  for (int pt = 0; pt < bins; pt++) // loop over pT bins
  {
    hPtDev[pt] = new TH1F(Form("hPtDev%i", pt), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys}/Y_{def} - 1;Counts", h[0]->GetBinLowEdge(pt + 1), h[0]->GetBinLowEdge(pt + 2)), 15, -0.5, +0.5);
    hPtDev[pt]->SetStats(0);
    fgaus[pt] = new TF1(Form("fgaus%i", pt), "gaus", -0.5, +0.5);
  }

  for (int i = 1; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    for (int pt = 0; pt < bins; pt++)
    {
      hPtDev[pt]->Fill(h[i]->GetBinContent(pt + 1));
    }
  }

  TLatex *ltx = new TLatex();
  ltx->SetTextSize(0.05);
  ltx->SetTextColor(kRed);

  for (int pt = 0; pt < bins; pt++)
  {
    c2->cd(pt + 1);
    c2->cd(pt + 1)->SetBottomMargin(0.15);
    hPtDev[pt]->GetYaxis()->SetRangeUser(0., hPtDev[pt]->GetMaximum() * 1.2);
    hPtDev[pt]->GetXaxis()->SetTitleSize(0.06);
    hPtDev[pt]->GetXaxis()->SetTitleOffset(1.);
    hPtDev[pt]->Draw("EP");
    hPtDev[pt]->SetMarkerStyle(kFullCircle);
    hPtDev[pt]->Draw("EP SAME");

    fgaus[pt]->SetParameter(0, hPtDev[pt]->GetMaximum());
    hPtDev[pt]->Fit(fgaus[pt], "R+");
    ltx->DrawLatexNDC(0.6, 0.8, Form("#mu = %.3f", fgaus[pt]->GetParameter(1)));
    ltx->DrawLatexNDC(0.6, 0.7, Form("#sigma = %.3f", fgaus[pt]->GetParameter(2)));
    ltx->DrawLatexNDC(0.6, 0.6, Form("#chi^{2}/ndf = %.1f", fgaus[pt]->GetChisquare() / fgaus[pt]->GetNDF()));
  }

  c2->SaveAs("Systematics/MultTrial" + Suffix + ".pdf");
  c2->SaveAs("Systematics/MultTrial" + Suffix + ".png");

  // relative syst. uncertainty
  TCanvas *cSyst = new TCanvas("cSyst", "cSyst", 1000, 800);
  TH1F *hSystMultiTrial = (TH1F *)h[0]->Clone("hSystMultiTrial");
  hSystMultiTrial->Reset();
  for (int pt = 0; pt < bins; pt++)
  {
    hSystMultiTrial->SetBinContent(pt + 1, fgaus[pt]->GetParameter(2));
    hSystMultiTrial->SetBinError(pt + 1, 0);
  }
  hSystMultiTrial->GetYaxis()->SetRangeUser(0., 0.5);
  hSystMultiTrial->GetYaxis()->SetTitle("Rel. syst. error");
  hSystMultiTrial->SetLineColor(kBlack);
  if (!isPtIntegrated)
    hSystMultiTrial->Smooth();
  hSystMultiTrial->Draw();

  // relative syst. uncertainty Max Variation
  // TCanvas *cSystMax = new TCanvas("cSystMax", "cSystMax", 1000, 800);
  hMaxDev->SetLineColor(kRed);
  if (!isPtIntegrated)
    hMaxDev->Smooth();
  hMaxDev->Draw("same");

  // Absolute syst. uncertainty from maximum deviation
  TCanvas *cSystAbsMaxDev = new TCanvas("cSystAbsMaxDev", "cSystAbsMaxDev", 1000, 800);
  hAbsoluteMaxDev->SetLineColor(kGray + 1);
  hAbsoluteMaxDev->SetMarkerColor(kGray + 1);
  hAbsoluteMaxDev->SetLineWidth(1);
  hAbsoluteMaxDev->SetMarkerStyle(33);
  hAbsoluteMaxDev->GetYaxis()->SetTitle("Absolute syst. uncertainty from BDT score variation");
  hAbsoluteMaxDev->GetYaxis()->SetRangeUser(0., 1.2 * hAbsoluteMaxDev->GetBinContent(hAbsoluteMaxDev->GetMaximumBin()));
  if (!isPtIntegrated)
    hAbsoluteMaxDev->Smooth();
  hAbsoluteMaxDev->Draw("same");

  // Stat. uncertainties
  TCanvas *cStat = new TCanvas("cStat", "cStat", 1000, 800);
  StyleCanvas(cStat, 0.15, 0.05, 0.05, 0.15);
  cStat->cd();
  hDefaultError->SetLineColor(kBlack);
  hDefaultError->GetYaxis()->SetRangeUser(0, 0.5);
  if (isPtIntegrated)
    hDefaultError->GetYaxis()->SetRangeUser(0, 0.001);
  hDefaultError->GetYaxis()->SetTitle("Statistical uncertainty");
  hDefaultError->Draw();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hError[i]->SetLineColor(ColorMult[i]);
    hError[i]->SetMarkerColor(ColorMult[i]);
    hError[i]->SetMarkerStyle(MarkerMult[i]);
    hError[i]->Draw("same");
  }
  hAbsoluteMaxDev->Draw("same");
  legTrial->Draw();
  cStat->SaveAs("Systematics/StatErrorMultTrial" + Suffix + ".pdf");
  cStat->SaveAs("Systematics/StatErrorMultTrial" + Suffix + ".png");

  // ratio of stat. uncertainties to default one
  TCanvas *cStatRatio = new TCanvas("cStatRatio", "cStatRatio", 1000, 800);
  StyleCanvas(cStatRatio, 0.15, 0.05, 0.05, 0.15);
  cStatRatio->cd();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    hErrorRatio[i]->SetLineColor(ColorMult[i]);
    hErrorRatio[i]->SetMarkerColor(ColorMult[i]);
    hErrorRatio[i]->SetMarkerStyle(MarkerMult[i]);
    hErrorRatio[i]->SetTitle("");
    hErrorRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hErrorRatio[i]->GetYaxis()->SetRangeUser(0.6, 1.4);
    hErrorRatio[i]->Draw("same");
  }
  lineat1->Draw("same");
  legTrial->Draw();
  cStatRatio->SaveAs("Systematics/StatErrorRatioMultTrial" + Suffix + ".pdf");
  cStatRatio->SaveAs("Systematics/StatErrorRatioMultTrial" + Suffix + ".png");

  // save histos in output files
  TString OutputFile = "Systematics/SystMultiTrial_" + inputFileName + Form("_%i-%i_", CentFT0CMin, CentFT0CMax) + ParticleName[ChosenPart] + "_";
  OutputFile += SisSyst;
  if (isPtIntegrated)
    OutputFile += "_PtInt";
  if (isApplyWeights)
    OutputFile += "_Weighted";
  if (!isPtAnalysis)
    OutputFile += "_vsPsi";
  if (!isV2 && isPolFromLambda)
    OutputFile += "_PolFromLambda";
  if (ExtrisApplyEffWeights)
    OutputFile += "_EffW";
  if (!isRapiditySel || ExtrisFromTHN)
    OutputFile += "_Eta08";
  OutputFile += STHN[ExtrisFromTHN];
  OutputFile += ".root";

  TFile *Write = new TFile(OutputFile, "RECREATE");
  for (int pt = 0; pt < bins; pt++)
  {
    hPtDev[pt]->Write();
    fgaus[pt]->Write();
  }
  hSystMultiTrial->Write();
  hMaxDev->Write();
  hAbsoluteMaxDev->Write();

  cout << "\n\nHo creato il file: " << OutputFile << endl;
}
