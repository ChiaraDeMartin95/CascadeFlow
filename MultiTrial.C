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

double PassRogerBarlowCriterion(float nsigmas, Double_t dev, Double_t RBsigma)
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

  TH1F *hVariedCutCopy = (TH1F *)hVariedCut->Clone("hVarCutCopy");

  DivideAndComputeRogerBarlow(hVariedCut, hDefault);
  // assigns to hVariedCut: hR = hvar/hdef as bin content, sigmaBarlow (sB) as error

  for (int i = 1; i <= hDev->GetNbinsX(); i++) // pt bins
  {
    double dev = hVariedCut->GetBinContent(i) - 1; // hR - 1
    double err = hVariedCut->GetBinError(i);       // sB

    hDev->SetBinContent(i, abs(PassRogerBarlowCriterion(nsigmaBarlow, dev, err))); // rel. syst. error = hR-1 if |hR-1| > 1*sB

    cout << "default: " << hDefault->GetBinContent(i) << " +- " << hDefault->GetBinError(i) << endl;
    cout << "varied: " << hVariedCutCopy->GetBinContent(i) << " +- " << hVariedCutCopy->GetBinError(i) << endl;
    cout << "difference " << hVariedCutCopy->GetBinContent(i) - hDefault->GetBinContent(i) << endl;
    cout << "SigmaBarlow: " << err * hDefault->GetBinContent(i) << endl;
    cout << "Deviation: " << dev * hDefault->GetBinContent(i) << endl;
    cout << "nsigma " << dev / err << endl;
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

Int_t IndexNotDisplayed = 5;
TString SisPtIntegrated[2] = {"", "PtInt"};

void MultiTrial(
    Int_t mul = 0,
    Int_t Choice = 0,          // 0 = V2Mixed, 1 = Pz(s2)Mixed, 2 = Pz(s2)LambdaFromCMixed
    Bool_t isPtAnalysis = 1,   // 1 for V2 vs pt and Pzs2 vs pt, 0 for Pz vs 2(phi-Psi)
    Bool_t isPtIntegrated = 1, // 1 for results integrated in pt / phi
    TString SisSyst = /*"LambdaTopo"*/ "MassAndBDTCut",
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
  TString TypeCos2Theta = "";
  if (Choice == 0)
    TypeHisto = "V2";
  else if (Choice == 1)
  {
    TypeHisto = "Pzs";
    if (!isPtAnalysis)
      TypeHisto = "Pz";
    TypeCos2Theta = "histoCos2Theta";
  }
  else if (Choice == 2)
  {
    isPolFromLambda = 1;
    TypeHisto = "Pzs2LambdaFromC";
    if (!isPtAnalysis)
      TypeHisto = "PzLambdaFromC";
    TypeCos2Theta = "histoCos2ThetaLambdaFromC";
  }
  else if (Choice == 3)
  {
    TypeHisto = "Pzs2";
    TypeCos2Theta = "histoCos2Theta";
  }
  TypeHisto += SisPtIntegrated[isPtIntegrated];
  // TypeHisto += "Mixed";
  TypeHisto += "NoFit";
  TypeCos2Theta += SisPtIntegrated[isPtIntegrated];
  TypeCos2Theta += "NoFit";

  TString histoName = "histo" + TypeHisto;

  TString Suffix = inputFileName + Form("_%i-%i_", CentFT0C[mul], CentFT0C[mul + 1]) + ParticleName[ChosenPart] + "_" + SisSyst + "_" + SisPtIntegrated[isPtIntegrated];

  if (isPolFromLambda)
  {
    if (CentFT0CMin == 50)
    {
      IndexNotDisplayed = 18;
    }
    if (CentFT0CMin == 60)
    {
      IndexNotDisplayed = 5;
    }
  }
  IndexNotDisplayed = 1000000;

  Int_t trials = 0;
  if (SisSyst == "BDT")
    trials = trialsBDT;
  else if (SisSyst == "LambdaTopo")
    trials = trialsLambdaTopo; // different lambda topologies
  else if (SisSyst == "eta")
    trials = 2; // eta > 0 and eta < 0
  else if (SisSyst == "IR")
    trials = 5; // different interaction rates
  else if (SisSyst == "MassCut")
    trials = trialsMassCut; // different mass cuts
  else if (SisSyst == "MassAndBDTCut")
  {
    if (mul < 3)
      trials = trialsBDT * trialsMassCut;
    else if (mul < 4)
      trials = trialsBDT * 5;
    else if (mul < 5)
      trials = trialsBDT * 4;
    else
      trials = trialsBDT * 3;
  }

  cout << "The number of trials to be investigated are: " << trials << endl;

  TString Sdef = "OutputAnalysis/Fit" + NameAnalysis[!isV2] + "_" + inputFileName + "_" + ParticleName[ChosenPart];
  Sdef += IsOneOrTwoGauss[UseTwoGauss];
  Sdef += SIsBkgParab[BkgType];
  Sdef += Form("_Cent%i-%i", CentFT0CMin, CentFT0CMax);
  if (isApplyWeights)
    Sdef += "_Weighted";
  if (isApplyCentWeight)
    Sdef += "_CentWeighted";
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
  // if (useMixedBDTValueInFitMacro)
  if (ChosenPart != 6)
    SdefFinal += "_MixedBDT";
  if (isTightMassCut)
    SdefFinal += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
  if (isReducedPtBins)
    SdefFinal += "_ReducedPtBins";

  TString Svaried = "";

  cout << "Default input file: " << SdefFinal << endl;
  TFile *fdef = TFile::Open(SdefFinal + ".root");
  if (!fdef)
  {
    cout << "File not found: " << SdefFinal << endl;
    return;
  }

  // get default BDT score
  TH1F *hBDT = (TH1F *)fdef->Get("histoAppliedBDT");
  if (isPtIntegrated)
    hBDT = (TH1F *)fdef->Get("histoAppliedBDTPtInt");
  if (!hBDT)
  {
    cout << "histoAppliedBDT not found in " << SdefFinal << endl;
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
  cout << "histoName: " << histoName << endl;
  TH1F *hDefault = (TH1F *)fdef->Get(histoName);
  if (!hDefault)
  {
    cout << "histo not found in " << SdefFinal << endl;
    return;
  }
  hDefault->SetName("hDefault");
  hDefault->SetLineColor(kBlack);
  hDefault->SetLineStyle(2);
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
  TH1F *hDefaultCos2Theta = (TH1F *)fdef->Get(TypeCos2Theta);
  if (!hDefaultCos2Theta)
  {
    cout << "histoCos2Theta not found in " << SdefFinal << endl;
    return;
  }
  hDefaultCos2Theta->SetName("hDefaultCos2Theta");

  Float_t LegDefaultBDTscoreCut = hBDT->GetBinContent(1);
  if (ChosenPart != 6)
    cout << "Default BDT score: " << LegDefaultBDTscoreCut << endl;

  TLegend *legTrial;
  if (SisSyst == "IR")
    legTrial = new TLegend(0.66, 0.7, 0.96, 0.9);
  else
    // legTrial = new TLegend(0.66, 0.2, 0.96, 0.5);
    legTrial = new TLegend(0.66, 0.2, 0.96, 0.8);
  legTrial->SetBorderSize(0);
  legTrial->SetFillStyle(0);
  legTrial->SetTextSize(0.03);
  legTrial->SetTextFont(42);
  if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
  {
    // if (useMixedBDTValueInFitMacro)
    if (kTRUE)
    {
      if (isPtIntegrated)
      {
        if (SisSyst == "MassAndBDTCut")
          legTrial->AddEntry(hDefault, Form("BDT score > %.3f, Mass cut: %.1f", LegDefaultBDTscoreCut, Extrsigmacentral[1]), "pl");
        else
          legTrial->AddEntry(hDefault, Form("BDT score > %.3f", LegDefaultBDTscoreCut), "pl");
      }
      else
        legTrial->AddEntry(hDefault, "BDT selection p_{T} dependent", "pl");
    }
    else
    {
      legTrial->AddEntry(hDefault, Form("BDT score > %.3f", DefaultBDTscoreCut), "pl");
    }
  }
  else if (SisSyst == "eta")
    legTrial->AddEntry(hDefault, "|#eta| < 0.8", "pl");
  else if (SisSyst == "IR")
    legTrial->AddEntry(hDefault, "PbPb 2023", "pl");
  else if (SisSyst == "MassCut")
    legTrial->AddEntry(hDefault, Form("Mass cut: %.1f", Extrsigmacentral[1]), "pl");
  TLegend *legTrialReduced = (TLegend *)legTrial->Clone("legTrialReduced");

  TH1F *h[trials];
  TH1F *hVariedCut[trials];
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
  TH1F *hCos2Theta[trials];
  TH1F *hCos2ThetaRatio[trials];
  Float_t BDTscoreCut[trials];
  Float_t LowLimitSysXi[trials];
  Float_t UpLimitSysXi[trials];
  TFile *fvaried[trials];

  TF1 *lineatnSigmaBarlow = new TF1("lineatnSigmaBarlow", "pol0", 0, 8);
  lineatnSigmaBarlow->SetLineColor(kBlack);
  lineatnSigmaBarlow->SetLineStyle(10);
  lineatnSigmaBarlow->SetParameter(0, nsigmaBarlow);

  for (int i = 0; i < trials; i++)
  {
    BDTscoreCut[i] = 0;
    LowLimitSysXi[i] = 0;
    UpLimitSysXi[i] = 0;
    cout << "\nTrial: " << i << endl;
    if (SisSyst == "BDT")
      BDTscoreCut[i] = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trials * i;
    else if (SisSyst == "MassAndBDTCut")
    {
      BDTscoreCut[i] = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * (i % trialsBDT);
      LowLimitSysXi[i] = ExtrLowLimitSysXi[i / trialsBDT];
      UpLimitSysXi[i] = ExtrUpLimitSysXi[i / trialsBDT];
      cout << "i % trialsBDT " << i % trialsBDT << endl;
      cout << "i / trialsBDT " << i / trialsBDT << endl;
      cout << "BDT score cut: " << BDTscoreCut[i] << endl;
      cout << "Mass cut: " << LowLimitSysXi[i] << " - " << UpLimitSysXi[i] << endl;
    }
    TString SBDT = Form("_BDT%.3f", BDTscoreCut[i]);
    if (SisSyst == "BDT")
    {
      Svaried = Sdef + SBDT + Sdef1;
      if (isTightMassCut)
        Svaried += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
    }
    else if (SisSyst == "eta")
      Svaried = Sdef + SEtaSysChoice[i + 1];
    else if (SisSyst == "IR")
    {
      Svaried = "OutputAnalysis/FitV2_" + inputFileNameIR + SIRChoice[i + 1] + "_" + ParticleName[ChosenPart];
      Svaried += IsOneOrTwoGauss[UseTwoGauss];
      Svaried += SIsBkgParab[BkgType];
      Svaried += Form("_Cent%i-%i", CentFT0C[mul], CentFT0C[mul + 1]);
    }
    else if (SisSyst == "MassCut")
    {
      Svaried = Sdef + Sdef1;
      // if (useMixedBDTValueInFitMacro)
      Svaried += "_MixedBDT";
      Svaried += Form("_TightMassCutSyst%i", i);
    }
    else if (SisSyst == "MassAndBDTCut")
    {
      Svaried = Sdef + SBDT + Sdef1;
      Svaried += Form("_TightMassCutSyst%i", i / trialsBDT);
    }
    else if (SisSyst == "LambdaTopo")
    {
      Svaried = Sdef;
      if (!isRapiditySel || ExtrisFromTHN)
        Svaried += "_Eta08";
      if (isTightMassCut)
        Svaried += Form("_TightMassCut%.1f", Extrsigmacentral[1]);
      if (isReducedPtBins)
        Svaried += "_ReducedPtBins";
      Svaried += Form("_SysMultTrial_%i", i);
      Svaried += "_isSysLambdaMultTrial_ResoOnTheFly";
      cout << "Svaried " << Svaried << endl;
    }

    cout << "InputFile - variation: " << Svaried << ".root" << endl;

    fvaried[i] = TFile::Open(Svaried + ".root");
    if (!fvaried[i])
    {
      cout << "File not found: " << Svaried << ".root" << endl;
      return;
    }
    hVariedCut[i] = (TH1F *)fvaried[i]->Get(histoName);
    if (!hVariedCut[i])
    {
      cout << "histo not found in " << Svaried << endl;
      return;
    }
    hVariedCut[i]->SetName(histoName + Form("_%i", i));
    hRatio[i] = (TH1F *)hVariedCut[i]->Clone(histoName + Form("Ratio_%i", i));
    hRatio[i]->Divide(hDefault);

    hRawYield[i] = (TH1F *)fvaried[i]->Get("histoYield" + SisPtIntegrated[isPtIntegrated]);
    if (!hRawYield[i])
    {
      cout << "histoYield not found in " << Svaried << endl;
      return;
    }
    hRawYield[i]->SetName(Form("histoYield_%i", i));
    hRawYieldRatio[i] = (TH1F *)hRawYield[i]->Clone(Form("histoYieldRatio_%i", i));
    hRawYieldRatio[i]->Divide(hDefaultYield);

    hPurity[i] = (TH1F *)fvaried[i]->Get("histoPurity" + SisPtIntegrated[isPtIntegrated]);
    if (!hPurity[i])
    {
      cout << "histoPurity not found in " << Svaried << endl;
      return;
    }
    hPurity[i]->SetName(Form("histoPurity_%i", i));
    hPurityRatio[i] = (TH1F *)hPurity[i]->Clone(Form("histoPurity_%i", i));
    hPurityRatio[i]->Divide(hDefaultPurity);

    hSignificance[i] = (TH1F *)fvaried[i]->Get("histoSignificance" + SisPtIntegrated[isPtIntegrated]);
    if (!hSignificance[i])
    {
      cout << "histoSignificance not found in " << Svaried << endl;
      return;
    }
    hSignificance[i]->SetName(Form("histoSignificance_%i", i));
    hSignificanceRatio[i] = (TH1F *)hSignificance[i]->Clone(Form("histoSignificanceRatio_%i", i));
    hSignificanceRatio[i]->Divide(hDefaultSignificance);

    hCos2Theta[i] = (TH1F *)fvaried[i]->Get(TypeCos2Theta);
    if (!hCos2Theta[i])
    {
      cout << "histoCos2Theta not found in " << Svaried << endl;
      return;
    }
    hCos2Theta[i]->SetName(Form("histoCos2Theta_%i", i));
    hCos2ThetaRatio[i] = (TH1F *)hCos2Theta[i]->Clone(Form("histoCos2ThetaRatio_%i", i));
    hCos2ThetaRatio[i]->Divide(hDefaultCos2Theta);

    hError[i] = (TH1F *)hVariedCut[i]->Clone(Form("hError_%i", i));
    hError[i]->Reset();
    for (int j = 1; j <= hVariedCut[i]->GetNbinsX(); j++)
    {
      hError[i]->SetBinContent(j, hVariedCut[i]->GetBinError(j));
    }
    hErrorRatio[i] = (TH1F *)hError[i]->Clone(Form("hErrorRatio_%i", i));
    hErrorRatio[i]->Divide(hDefaultError);

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

  // v2 plots
  TCanvas *cv2 = new TCanvas("cv2", Titlecv2, 1000, 800);
  StyleCanvas(cv2, 0.15, 0.05, 0.05, 0.15);
  cv2->cd();
  Int_t ColorIndex = -1;
  Int_t NumberOfActualTrials = 0;

  for (int i = 0; i < trials; i++)
  {
    cout << "\nTrial: " << i << endl;
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    cout << "here I am " << endl;
    NumberOfActualTrials++;
    ColorIndex += 1;
    hVariedCut[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hVariedCut[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hVariedCut[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
    if (SisSyst == "BDT")
    {
      legTrial->AddEntry(hVariedCut[i], Form("BDT score > %.3f", BDTscoreCut[i]), "pl");
      legTrialReduced->AddEntry(hVariedCut[i], Form("BDT score > %.3f", BDTscoreCut[i]), "pl");
    }
    else if (SisSyst == "eta")
      legTrial->AddEntry(hVariedCut[i], SEtaSysChoice[i + 1], "pl");
    else if (SisSyst == "IR")
      legTrial->AddEntry(hVariedCut[i], SIRValue[i + 1], "pl");
    else if (SisSyst == "MassCut")
      legTrial->AddEntry(hVariedCut[i], Form("Mass cut: %.3f - %.3f", ExtrLowLimitSysXi[i], ExtrUpLimitSysXi[i]), "pl");
    else if (SisSyst == "MassAndBDTCut")
      legTrial->AddEntry(hVariedCut[i], Form("BDT score > %.3f, Mass cut: %.3f - %.3f", BDTscoreCut[i], LowLimitSysXi[i], UpLimitSysXi[i]), "pl");
    hVariedCut[i]->Draw("same");
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
  ColorIndex = -1;
  TH1F *hdummy = (TH1F *)hDefault->Clone("hdummy");
  hdummy->Reset();
  hdummy->GetYaxis()->SetRangeUser(0, 5.);
  hdummy->GetYaxis()->SetTitle("N_{#sigma}^{Barlow}");
  hdummy->Draw();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      // cout << "Selected BDT before: " << BDTscoreCut[i] << endl;
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
      // cout << "Selected BDT: " << BDTscoreCut[i] << endl;
    }
    ColorIndex += 1;
    if (isPtIntegrated)
    {
      // if (useMixedBDTValueInFitMacro && BDTscoreCut[i] == LegDefaultBDTscoreCut)
      if (kTRUE && BDTscoreCut[i] == LegDefaultBDTscoreCut)
        continue;
      // if (!useMixedBDTValueInFitMacro && BDTscoreCut[i] == DefaultBDTscoreCut)
      if (!kTRUE && BDTscoreCut[i] == DefaultBDTscoreCut)
        continue;
    }
    hNSigmaBarlow[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hNSigmaBarlow[i]->SetLineWidth(2);
    hNSigmaBarlow[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hNSigmaBarlow[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
    hNSigmaBarlow[i]->Draw("same");
  }
  lineatnSigmaBarlow->Draw("same");

  // ratio of varied v2 to default one
  TCanvas *cv2Ratio = new TCanvas("cv2Ratio", "cv2Ratio", 1000, 800);
  StyleCanvas(cv2Ratio, 0.15, 0.05, 0.05, 0.15);
  cv2Ratio->cd();
  ColorIndex = -1;
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hRatio[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hRatio[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hRatio[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
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
    // cv2Ratio->SaveAs("Systematics/Pzs2RatioMultTrial" + Suffix + ".pdf");
    // cv2Ratio->SaveAs("Systematics/Pzs2RatioMultTrial" + Suffix + ".png");
    cv2Ratio->Close();
  }

  // raw yields plots
  TCanvas *cYield = new TCanvas("cYield", "cYield", 1000, 800);
  StyleCanvas(cYield, 0.15, 0.05, 0.05, 0.15);
  cYield->cd();
  ColorIndex = -1;
  hDefaultYield->SetTitle("");
  hDefaultYield->SetLineColor(kBlack);
  hDefaultYield->SetLineStyle(2);
  hDefaultYield->GetYaxis()->SetRangeUser(0, 1.8 * hDefaultYield->GetMaximum());
  hDefaultYield->Draw();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hRawYield[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hRawYield[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hRawYield[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
    hRawYield[i]->Draw("same");
  }
  legTrial->Draw();
  cYield->SaveAs("Systematics/YieldMultTrial" + Suffix + ".pdf");
  cYield->SaveAs("Systematics/YieldMultTrial" + Suffix + ".png");

  // ratio of raw yields to default one
  TCanvas *cYieldRatio = new TCanvas("cYieldRatio", "cYieldRatio", 1000, 800);
  StyleCanvas(cYieldRatio, 0.15, 0.05, 0.05, 0.15);
  cYieldRatio->cd();
  ColorIndex = -1;
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hRawYieldRatio[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hRawYieldRatio[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hRawYieldRatio[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
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
  ColorIndex = -1;
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hRawYieldRatioToTighter[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hRawYieldRatioToTighter[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hRawYieldRatioToTighter[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
    hRawYieldRatioToTighter[i]->SetTitle("");
    hRawYieldRatioToTighter[i]->GetYaxis()->SetTitle("Ratio to tighter cut");
    hRawYieldRatioToTighter[i]->GetYaxis()->SetRangeUser(0.9, 1.5);
    hRawYieldRatioToTighter[i]->Draw("same");
  }
  // legTrial->Draw();
  lineat1->Draw("same");
  cYieldRatioToTighter->SaveAs("Systematics/YieldRatioToTighterMultTrial" + Suffix + ".pdf");
  cYieldRatioToTighter->SaveAs("Systematics/YieldRatioToTighterMultTrial" + Suffix + ".png");

  // purity plots
  TCanvas *cPurity = new TCanvas("cPurity", "cPurity", 1000, 800);
  StyleCanvas(cPurity, 0.15, 0.05, 0.05, 0.15);
  cPurity->cd();
  ColorIndex = -1;
  hDefaultPurity->SetTitle("");
  hDefaultPurity->SetLineColor(kBlack);
  hDefaultPurity->SetLineStyle(2);
  hDefaultPurity->GetYaxis()->SetRangeUser(0, 1);
  hDefaultPurity->Draw();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hPurity[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hPurity[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hPurity[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
    hPurity[i]->Draw("same");
  }
  legTrialReduced->Draw();
  cPurity->SaveAs("Systematics/PurityMultTrial" + inputFileName + Suffix + ".pdf");
  cPurity->SaveAs("Systematics/PurityMultTrial" + inputFileName + Suffix + ".png");

  // ratio of purities to default one
  TCanvas *cPurityRatio = new TCanvas("cPurityRatio", "cPurityRatio", 1000, 800);
  StyleCanvas(cPurityRatio, 0.15, 0.05, 0.05, 0.15);
  cPurityRatio->cd();
  ColorIndex = -1;
  TH1F *hDummyPurityRatio = (TH1F *)hDefaultPurity->Clone("hDummyPurityRatio");
  hDummyPurityRatio->Reset();
  hDummyPurityRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
  hDummyPurityRatio->Draw();
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hPurityRatio[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hPurityRatio[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hPurityRatio[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
    hPurityRatio[i]->SetTitle("");
    hPurityRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hPurityRatio[i]->GetYaxis()->SetRangeUser(0, 1.1);
    hPurityRatio[i]->Draw("same");
  }
  lineat1->Draw("same");
  legTrialReduced->Draw();
  cPurityRatio->SaveAs("Systematics/PurityRatioMultTrial" + Suffix + ".pdf");
  cPurityRatio->SaveAs("Systematics/PurityRatioMultTrial" + Suffix + ".png");

  // significance plots
  TCanvas *cSignificance = new TCanvas("cSignificance", "cSignificance", 1000, 800);
  StyleCanvas(cSignificance, 0.15, 0.05, 0.05, 0.15);
  cSignificance->cd();
  ColorIndex = -1;
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hSignificance[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hSignificance[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hSignificance[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
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
  ColorIndex = -1;
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hSignificanceRatio[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hSignificanceRatio[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hSignificanceRatio[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
    hSignificanceRatio[i]->SetTitle("");
    hSignificanceRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hSignificanceRatio[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
    hSignificanceRatio[i]->Draw("same");
  }
  legTrial->Draw();
  lineat1->Draw("same");
  cSignificanceRatio->SaveAs("Systematics/SignificanceRatioMultTrial" + Suffix + ".pdf");
  cSignificanceRatio->SaveAs("Systematics/SignificanceRatioMultTrial" + Suffix + ".png");

  // acceptance plots
  TCanvas *cAcceptance = new TCanvas("cAcceptance", "cAcceptance", 1000, 800);
  StyleCanvas(cAcceptance, 0.15, 0.05, 0.05, 0.15);
  cAcceptance->cd();
  ColorIndex = -1;
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hCos2Theta[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hCos2Theta[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hCos2Theta[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
    hCos2Theta[i]->GetYaxis()->SetRangeUser(0.31, 0.33);
    hCos2Theta[i]->SetTitle("");
    hCos2Theta[i]->Draw("same");
  }
  legTrial->Draw();
  cAcceptance->SaveAs("Systematics/AcceptanceMultTrial" + Suffix + ".pdf");
  cAcceptance->SaveAs("Systematics/AcceptanceMultTrial" + Suffix + ".png");

  // ratio of acceptances to default one
  TCanvas *cAcceptanceRatio = new TCanvas("cAcceptanceRatio", "cAcceptanceRatio", 1000, 800);
  StyleCanvas(cAcceptanceRatio, 0.15, 0.05, 0.05, 0.15);
  cAcceptanceRatio->cd();
  ColorIndex = -1;
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hCos2ThetaRatio[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hCos2ThetaRatio[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hCos2ThetaRatio[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
    hCos2ThetaRatio[i]->SetTitle("");
    hCos2ThetaRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hCos2ThetaRatio[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
    // hCos2ThetaRatio[i]->Draw("same");
  }
  legTrialReduced->Draw();
  lineat1->Draw("same");
  cAcceptanceRatio->SaveAs("Systematics/AcceptanceRatioMultTrial" + Suffix + ".pdf");
  cAcceptanceRatio->SaveAs("Systematics/AcceptanceRatioMultTrial" + Suffix + ".png");

  // systematic computation taking max variation wrt to default
  const int bins = h[0]->GetNbinsX(); // pt bins
  TH1F *hMaxDev = (TH1F *)h[0]->Clone("hMaxDev");
  hMaxDev->Reset();
  TH1F *hAbsoluteMaxDev = (TH1F *)hAbsoluteSyst[0]->Clone("hAbsoluteMaxDev");
  hAbsoluteMaxDev->Reset();
  TH1F *hRMS = (TH1F *)hAbsoluteSyst[0]->Clone("hRMS");
  hRMS->Reset();
  Float_t RMS = 0;
  Int_t counter = 0;

  for (int pt = 0; pt < bins; pt++) // loop over pT bins
  {
    // cout << "check " << hMaxDev->GetBinContent(pt + 1) << " " << hAbsoluteMaxDev->GetBinContent(pt + 1) << endl;
    for (int i = 0; i < trials; i++)
    {
      if (i == IndexNotDisplayed)
        continue;
      // skip those variations for which the fit is not good
      if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
      {
        if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
          continue;
      }
      if (TMath::Abs(h[i]->GetBinContent(pt + 1)) > TMath::Abs(hMaxDev->GetBinContent(pt + 1)))
      {
        hMaxDev->SetBinContent(pt + 1, h[i]->GetBinContent(pt + 1));
      }
      hMaxDev->SetBinError(pt + 1, 0.);
      if (TMath::Abs(hAbsoluteSyst[i]->GetBinContent(pt + 1)) > TMath::Abs(hAbsoluteMaxDev->GetBinContent(pt + 1)))
      {
        hAbsoluteMaxDev->SetBinContent(pt + 1, TMath::Abs(hAbsoluteSyst[i]->GetBinContent(pt + 1)));
      }
      hAbsoluteMaxDev->SetBinError(pt + 1, 0.);
      RMS += pow(hAbsoluteSyst[i]->GetBinContent(pt + 1), 2);
      counter += 1;
    }
    RMS = sqrt(RMS / counter);
    hRMS->SetBinContent(pt + 1, RMS);
  }

  //  gaussian distributions + fits
  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1200);
  // c2->Divide(4,4);

  TH1F *hPtDev[bins];
  TF1 *fgaus[bins];

  for (int pt = 0; pt < bins; pt++) // loop over pT bins
  {
    hPtDev[pt] = new TH1F(Form("hPtDev%i", pt), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys}/Y_{def} - 1;Counts", h[0]->GetBinLowEdge(pt + 1), h[0]->GetBinLowEdge(pt + 2)), 30, -2., +2.);
    if (mul == 1)
      hPtDev[pt] = new TH1F(Form("hPtDev%i", pt), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys}/Y_{def} - 1;Counts", h[0]->GetBinLowEdge(pt + 1), h[0]->GetBinLowEdge(pt + 2)), 30, -0.5, +0.5);
    else if (mul == 2)
      hPtDev[pt] = new TH1F(Form("hPtDev%i", pt), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys}/Y_{def} - 1;Counts", h[0]->GetBinLowEdge(pt + 1), h[0]->GetBinLowEdge(pt + 2)), 30, -10.0, +10.0);
    else if (mul == 3)
      hPtDev[pt] = new TH1F(Form("hPtDev%i", pt), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys}/Y_{def} - 1;Counts", h[0]->GetBinLowEdge(pt + 1), h[0]->GetBinLowEdge(pt + 2)), 30, -2.0, +2.0);
    hPtDev[pt]->SetStats(0);
    fgaus[pt] = new TF1(Form("fgaus%i", pt), "gaus", -2.5, +2.5);
  }

  for (int i = 1; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
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
    hPtDev[pt]->Fit(fgaus[pt], "LLR+");
    ltx->DrawLatexNDC(0.6, 0.8, Form("#mu = %.3f", fgaus[pt]->GetParameter(1)));
    ltx->DrawLatexNDC(0.6, 0.7, Form("#sigma = %.3f", fgaus[pt]->GetParameter(2)));
    ltx->DrawLatexNDC(0.6, 0.6, Form("#chi^{2}/ndf = %.1f", fgaus[pt]->GetChisquare() / fgaus[pt]->GetNDF()));
  }

  c2->SaveAs("Systematics/MultTrial" + Suffix + ".pdf");
  c2->SaveAs("Systematics/MultTrial" + Suffix + ".png");
  // c2->Close();

  // relative syst. uncertainty from fit
  TCanvas *cSystGaussFit = new TCanvas("cSystGaussFit", "cSystGaussFit", 1000, 800);
  TH1F *hSystMultiTrial = (TH1F *)h[0]->Clone("hSystMultiTrial");
  TH1F *hSystMultiTrialRMS = (TH1F *)h[0]->Clone("hSystMultiTrialRMS");
  hSystMultiTrial->Reset();
  hSystMultiTrialRMS->Reset();
  for (int pt = 0; pt < bins; pt++)
  {
    hSystMultiTrial->SetBinContent(pt + 1, fgaus[pt]->GetParameter(2));
    hSystMultiTrial->SetBinError(pt + 1, 0);
    hSystMultiTrialRMS->SetBinContent(pt + 1, hPtDev[pt]->GetRMS());
    hSystMultiTrialRMS->SetBinError(pt + 1, 0);
  }
  hSystMultiTrial->GetYaxis()->SetRangeUser(0., 0.5);
  hSystMultiTrial->GetYaxis()->SetTitle("Rel. syst. error");
  hSystMultiTrial->SetLineColor(kBlack);
  if (!isPtIntegrated)
    hSystMultiTrial->Smooth();
  hSystMultiTrial->Draw();
  cSystGaussFit->Close();

  // Relative syst. uncertainty from maximum deviation
  TCanvas *cSystMaxDev = new TCanvas("cSystMaxDev", "cSystMaxDev", 1000, 800);
  hMaxDev->SetLineColor(kRed);
  if (!isPtIntegrated)
    hMaxDev->Smooth();
  hMaxDev->Draw("same");
  if (!isV2)
    cSystMaxDev->Close();

  TCanvas *cSystAbsDev = new TCanvas("cSystAbsDev", "cSystAbsDev", 1000, 800);
  StyleCanvas(cSystAbsDev, 0.15, 0.05, 0.05, 0.15);
  cSystAbsDev->cd();
  ColorIndex = -1;
  for (Int_t i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hAbsoluteSyst[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hAbsoluteSyst[i]->SetLineWidth(2);
    hAbsoluteSyst[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hAbsoluteSyst[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
    hAbsoluteSyst[i]->GetYaxis()->SetTitle("Absolute deviation from default");
    hAbsoluteSyst[i]->GetYaxis()->SetRangeUser(-1.2 * hAbsoluteMaxDev->GetBinContent(hAbsoluteMaxDev->GetMaximumBin()), 1.2 * hAbsoluteMaxDev->GetBinContent(hAbsoluteMaxDev->GetMaximumBin()));
    hAbsoluteSyst[i]->Draw("same");
  }
  legTrialReduced->Draw();

  // Gaussian distribution of maximum deviations
  TCanvas *cgaus = new TCanvas("cgaus", "cgaus", 1000, 800);
  StyleCanvas(cgaus, 0.15, 0.05, 0.05, 0.15);
  TH1F *hCollectionAbsoluteSyst[bins];
  TF1 *fgaus2[bins];
  for (int pt = 0; pt < bins; pt++) // loop over pT bins
  {
    cgaus->cd(pt + 1);
    fgaus2[pt] = new TF1(Form("fgaus2%i", pt), "gaus", -2 * hAbsoluteMaxDev->GetBinContent(hAbsoluteMaxDev->GetMaximumBin()), 2 * hAbsoluteMaxDev->GetBinContent(hAbsoluteMaxDev->GetMaximumBin()));
    hCollectionAbsoluteSyst[pt] = new TH1F(Form("hCollectionAbsoluteSyst%i", pt), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys};Counts", h[0]->GetBinLowEdge(pt + 1), h[0]->GetBinLowEdge(pt + 2)), 30, -2 * hAbsoluteMaxDev->GetBinContent(hAbsoluteMaxDev->GetMaximumBin()), 2 * hAbsoluteMaxDev->GetBinContent(hAbsoluteMaxDev->GetMaximumBin()));
    for (int i = 0; i < trials; i++)
    {
      if (i == IndexNotDisplayed)
        continue;
      // skip those variations for which the fit is not good
      if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
      {
        if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
          continue;
      }
      hCollectionAbsoluteSyst[pt]->Fill(hAbsoluteSyst[i]->GetBinContent(pt + 1));
    }
    cout << "Fitting with gaus2 " << endl;
    hCollectionAbsoluteSyst[pt]->Draw("EP SAME");
    fgaus2[pt]->SetParameter(0, hCollectionAbsoluteSyst[pt]->GetMaximum());
    fgaus2[pt]->SetLineColor(kBlue);
    hCollectionAbsoluteSyst[pt]->Fit(fgaus2[pt], "LLR+");
    fgaus2[pt]->Draw("same");
  }
  cgaus->SaveAs("Systematics/MultTrial_AbsoluteSystGauss" + Suffix + ".pdf");
  cgaus->SaveAs("Systematics/MultTrial_AbsoluteSystGauss" + Suffix + ".png");
  // cgaus->Close();
  
  TH1F *hSystMultiTrial2 = (TH1F *)h[0]->Clone("hSystMultiTrial2");
  hSystMultiTrial2->Reset();
  for (int pt = 0; pt < bins; pt++)
  {
    hSystMultiTrial2->SetBinContent(pt + 1, fgaus2[pt]->GetParameter(2));
    hSystMultiTrial2->SetBinError(pt + 1, 0);
  }

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
  hRMS->SetLineColor(kRed);
  hRMS->SetMarkerColor(kRed);
  hRMS->SetLineWidth(1);
  hRMS->SetMarkerStyle(33);
  hRMS->Draw("same");
  TH1F *hSystFromGauss = (TH1F *)hSystMultiTrial->Clone("hSystFromGauss");
  hSystFromGauss->Reset();
  hSystFromGauss->SetBinContent(1, abs(hSystMultiTrial->GetBinContent(1) * hDefault->GetBinContent(1)));
  hSystFromGauss->SetBinError(1, 0);
  hSystFromGauss->SetLineColor(kBlue);
  hSystFromGauss->SetMarkerColor(kBlue);
  hSystFromGauss->SetMarkerStyle(33);
  hSystFromGauss->Draw("same");

  // Stat. uncertainties
  TCanvas *cStat = new TCanvas("cStat", "cStat", 1000, 800);
  StyleCanvas(cStat, 0.15, 0.05, 0.05, 0.15);
  cStat->cd();
  hDefaultError->SetLineColor(kBlack);
  hDefaultError->SetLineStyle(2);
  hDefaultError->GetYaxis()->SetRangeUser(0, 0.5);
  if (isPtIntegrated)
    hDefaultError->GetYaxis()->SetRangeUser(0, 0.003);
  hDefaultError->GetYaxis()->SetTitle("Statistical uncertainty");
  hDefaultError->Draw();
  ColorIndex = -1;
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hError[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hError[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hError[i]->SetMarkerStyle(MarkerMult[ColorIndex % numCent]);
    hError[i]->Draw("same");
  }
  hAbsoluteMaxDev->Draw("same");
  hRMS->Draw("same");
  hSystFromGauss->Draw("same");
  legTrialReduced->Draw();
  cStat->SaveAs("Systematics/StatErrorMultTrial" + Suffix + ".pdf");
  cStat->SaveAs("Systematics/StatErrorMultTrial" + Suffix + ".png");

  // ratio of stat. uncertainties to default one
  TCanvas *cStatRatio = new TCanvas("cStatRatio", "cStatRatio", 1000, 800);
  StyleCanvas(cStatRatio, 0.15, 0.05, 0.05, 0.15);
  cStatRatio->cd();
  ColorIndex = -1;
  for (int i = 0; i < trials; i++)
  {
    if (i == IndexNotDisplayed)
      continue;
    if (SisSyst == "BDT" || SisSyst == "MassAndBDTCut")
    {
      if (BDTscoreCut[i] < MinBDTscorePtInt[mul] || BDTscoreCut[i] > (MaxBDTscorePtInt[mul] + 0.001))
        continue;
    }
    ColorIndex += 1;
    hErrorRatio[i]->SetLineColor(ColorMult[ColorIndex % numCent]);
    hErrorRatio[i]->SetMarkerColor(ColorMult[ColorIndex % numCent]);
    hErrorRatio[i]->SetMarkerStyle(MarkerMult[i]);
    hErrorRatio[i]->SetTitle("");
    hErrorRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hErrorRatio[i]->GetYaxis()->SetRangeUser(0.2, 1.8);
    hErrorRatio[i]->Draw("same");
  }
  lineat1->Draw("same");
  legTrialReduced->Draw();
  cStatRatio->SaveAs("Systematics/StatErrorRatioMultTrial" + Suffix + ".pdf");
  cStatRatio->SaveAs("Systematics/StatErrorRatioMultTrial" + Suffix + ".png");

  // save histos in output files
  TString OutputFile = "Systematics/SystMultiTrial_" + inputFileName + Form("_%i-%i_", CentFT0CMin, CentFT0CMax) + ParticleName[ChosenPart] + "_";
  OutputFile += SisSyst;
  if (isPtIntegrated)
    OutputFile += "_PtInt";
  if (isApplyWeights)
    OutputFile += "_Weighted";
  if (isApplyCentWeight)
    OutputFile += "_CentWeighted";
  if (!isPtAnalysis)
    OutputFile += "_vsPsi";
  if (!isV2 && isPolFromLambda)
    OutputFile += "_PolFromLambda";
  if (ExtrisApplyEffWeights)
    OutputFile += "_EffW";
  if (!isRapiditySel || ExtrisFromTHN)
    OutputFile += "_Eta08";
  OutputFile += STHN[ExtrisFromTHN];
  if (nsigmaBarlow != 0)
    OutputFile += Form("_nsigmaBarlow%.1f", nsigmaBarlow);
  OutputFile += ".root";

  TFile *Write = new TFile(OutputFile, "RECREATE");
  for (int pt = 0; pt < bins; pt++)
  {
    hPtDev[pt]->Write();
    fgaus[pt]->Write();
    fgaus2[pt]->Write();
  }
  hSystMultiTrial->Write();
  hSystMultiTrialRMS->Write();
  hSystMultiTrial2->Write();
  hMaxDev->Write();
  hAbsoluteMaxDev->Write();
  hRMS->Write();
  hSystFromGauss->Write();
  Write->Close();

  cout << "\n\nHo creato il file: " << OutputFile << endl;
  cout << "The number of actual trials is: " << NumberOfActualTrials << endl;
  cout << "The default value is: " << hDefault->GetBinContent(1) << endl;
  cout << "Syst error (RMS of results): " << hRMS->GetBinContent(1) << endl;
  cout << "Syst. error (gauss fit to results distribution): " << hSystMultiTrial2->GetBinContent(1) << endl;
  cout << "Syst. error from gaus fit (from rel. deviation): " << abs(hSystMultiTrial->GetBinContent(1) * hDefault->GetBinContent(1)) << endl;
  cout << "Syst. error from RMS (of rel. deviations): " << abs(hSystMultiTrialRMS->GetBinContent(1) * hDefault->GetBinContent(1)) << endl;
  cout << "Stat error: " << hDefaultError->GetBinContent(1) << endl;
}
