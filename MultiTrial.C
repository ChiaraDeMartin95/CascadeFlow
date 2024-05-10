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

TH1F *makeSystPlots(int num = 1, TString Sdef = "", TString Svaried = "")
{

  TFile *fdef = TFile::Open(Sdef);
  TFile *fvaried = TFile::Open(Svaried);

  TH1F *hVariedCut = (TH1F *)fvaried->Get("histoV2");
  hVariedCut->SetName(Form("hVarCutSet%i", num));
  TH1F *hDefault = (TH1F *)fdef->Get("histoV2");
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

TH1F *makeNSigmaBarlowPlots(int num = 1, TString Sdef = "", TString Svaried = "")
{

  TFile *fdef = TFile::Open(Sdef);
  TFile *fvaried = TFile::Open(Svaried);

  TH1F *hVariedCut = (TH1F *)fvaried->Get("histoV2");
  hVariedCut->SetName(Form("hVarCutSet%i", num));
  TH1F *hDefault = (TH1F *)fdef->Get("histoV2");
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

void MultiTrial(
    Int_t mul = 0,
    TString SisSyst = "BDT",
    Bool_t isXi = ChosenParticleXi,
    Int_t EtaSysChoice = ExtrEtaSysChoice,
    Int_t BkgType = ExtrBkgType,
    TString inputFileName = SinputFileNameSyst,
    Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  TString Suffix = inputFileName + Form("_%i-%i_", CentFT0C[mul], CentFT0C[mul + 1]) + ParticleName[!isXi] + "_" + SisSyst + ".pdf";

  Int_t trials = 0;
  if (SisSyst == "BDT")
    trials = trialsBDT;
  else if (SisSyst == "eta")
    trials = 2; // eta > 0 and eta < 0
  else if (SisSyst == "IR")
    trials = 5; // different interaction rates
  TString Sdef = "OutputAnalysis/FitV2_" + inputFileName + "_" + ParticleName[!isXi];
  Sdef += IsOneOrTwoGauss[UseTwoGauss];
  Sdef += SIsBkgParab[BkgType];
  Sdef += Form("_Cent%i-%i", CentFT0C[mul], CentFT0C[mul + 1]);

  TString Svaried = "";

  TFile *fdef = TFile::Open(Sdef + ".root");
  if (!fdef)
  {
    cout << "File not found: " << Sdef << endl;
    return;
  }

  // v2 plots
  TCanvas *cv2 = new TCanvas("cv2", "cv2", 1000, 800);
  StyleCanvas(cv2, 0.15, 0.05, 0.05, 0.15);
  cv2->cd();
  TH1F *hDefault = (TH1F *)fdef->Get("histoV2");
  hDefault->SetName("hDefault");
  hDefault->SetLineColor(kBlack);
  hDefault->Draw();
  TH1F *hDefaultYield = (TH1F *)fdef->Get("histoYield");
  hDefaultYield->SetName("hDefaultYield");
  TH1F *hDefaultPurity = (TH1F *)fdef->Get("histoPurity");
  hDefaultPurity->SetName("hDefaultPurity");
  TLegend *legTrial;
  if (SisSyst == "IR") legTrial = new TLegend(0.66, 0.7, 0.96, 0.9);
  else legTrial = new TLegend(0.66, 0.2, 0.96, 0.5);
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
  TH1F *hNSigmaBarlow[trials];
  TH1F *hRatio[trials];
  TH1F *hRawYield[trials];
  TH1F *hRawYieldRatio[trials];
  TH1F *hPurity[trials];
  TH1F *hPurityRatio[trials];
  Float_t BDTscoreCut = 0;
  for (int i = 0; i < trials; i++)
  {
    if (i == 8)
      continue;

    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trials * i;
    TString SBDT = Form("_BDT%.3f", BDTscoreCut);
    if (SisSyst == "BDT")
      Svaried = Sdef + SBDT;
    else if (SisSyst == "eta")
      Svaried = Sdef + SEtaSysChoice[i + 1];
    else if (SisSyst == "IR")
    {
      Svaried = "OutputAnalysis/FitV2_" + inputFileNameIR + SIRChoice[i + 1] + "_" + ParticleName[!isXi];
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
    TH1F *hVariedCut = (TH1F *)fvaried->Get("histoV2");
    hVariedCut->SetName(Form("histoV2_%i", i));
    hRatio[i] = (TH1F *)hVariedCut->Clone(Form("histoV2Ratio_%i", i));
    hRatio[i]->Divide(hDefault);

    hRawYield[i] = (TH1F *)fvaried->Get("histoYield");
    hRawYield[i]->SetName(Form("histoYield_%i", i));
    hRawYieldRatio[i] = (TH1F *)hRawYield[i]->Clone(Form("histoYieldRatio_%i", i));
    hRawYieldRatio[i]->Divide(hDefaultYield);

    hPurity[i] = (TH1F *)fvaried->Get("histoPurity");
    hPurity[i]->SetName(Form("histoPurity_%i", i));
    hPurityRatio[i] = (TH1F *)hPurity[i]->Clone(Form("histoPurity_%i", i));
    hPurityRatio[i]->Divide(hDefaultPurity);

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

    h[i] = makeSystPlots(i + 1, Sdef + ".root", Svaried + ".root");
    hNSigmaBarlow[i] = makeNSigmaBarlowPlots(i + 1, Sdef + ".root", Svaried + ".root");
  }
  legTrial->Draw();
  cv2->SaveAs("Systematics/V2MultTrial" + Suffix + ".pdf");
  cv2->SaveAs("Systematics/V2MultTrial" + Suffix + ".png");

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
    if (i == 8)
      continue;
    hNSigmaBarlow[i]->SetLineColor(ColorMult[i]);
    hNSigmaBarlow[i]->SetLineWidth(2);
    hNSigmaBarlow[i]->SetMarkerColor(ColorMult[i]);
    hNSigmaBarlow[i]->SetMarkerStyle(MarkerMult[i]);
    hNSigmaBarlow[i]->Draw("same");
  }

  // ratio of varied v2 to default one
  TCanvas *cv2Ratio = new TCanvas("cv2Ratio", "cv2Ratio", 1000, 800);
  StyleCanvas(cv2Ratio, 0.15, 0.05, 0.05, 0.15);
  cv2Ratio->cd();
  for (int i = 0; i < trials; i++)
  {
    if (i == 8)
      continue;
    hRatio[i]->SetLineColor(ColorMult[i]);
    hRatio[i]->SetMarkerColor(ColorMult[i]);
    hRatio[i]->SetMarkerStyle(MarkerMult[i]);
    hRatio[i]->SetTitle("");
    hRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hRatio[i]->GetYaxis()->SetRangeUser(0, 2);
    hRatio[i]->Draw("same");
  }
  TF1 *lineat1 = new TF1("lineat1", "1", 0, 5);
  lineat1->SetLineColor(kBlack);
  lineat1->SetLineStyle(7);
  legTrial->Draw();
  lineat1->Draw("same");
  cv2Ratio->SaveAs("Systematics/v2RatioMultTrial" + Suffix + ".pdf");
  cv2Ratio->SaveAs("Systematics/v2RatioMultTrial" + Suffix + ".png");

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
    if (i == 8)
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
    if (i == 8)
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
    if (i == 8)
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
    if (i == 8)
      continue;
    hPurityRatio[i]->SetLineColor(ColorMult[i]);
    hPurityRatio[i]->SetMarkerColor(ColorMult[i]);
    hPurityRatio[i]->SetMarkerStyle(MarkerMult[i]);
    hPurityRatio[i]->SetTitle("");
    hPurityRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hPurityRatio[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
    hPurityRatio[i]->Draw("same");
  }
  lineat1->Draw("same");
  legTrial->Draw();
  cPurityRatio->SaveAs("Systematics/PurityRatioMultTrial" + Suffix + ".pdf");
  cPurityRatio->SaveAs("Systematics/PurityRatioMultTrial" + Suffix + ".png");

  // systematic computation taking 0.5 * (max - min) variation
  const int bins = h[0]->GetNbinsX(); // pt bins
  TH1F *hMaxDev = (TH1F *)h[0]->Clone("hMaxDev");
  hMaxDev->Reset();
  for (int pt = 0; pt < bins; pt++) // loop over pT bins
  {
    for (int i = 1; i < trials; i++)
    {
      if (i == 8)
        continue;
      if (TMath::Abs(h[i]->GetBinContent(pt + 1)) > TMath::Abs(hMaxDev->GetBinContent(pt + 1)))
      {
        hMaxDev->SetBinContent(pt + 1, h[i]->GetBinContent(pt + 1));
      }
      hMaxDev->SetBinError(pt + 1, 0.);
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
    if (i == 8)
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
  hSystMultiTrial->Smooth();
  hSystMultiTrial->Draw();

  // relative syst. uncertainty Max Variation
  // TCanvas *cSystMax = new TCanvas("cSystMax", "cSystMax", 1000, 800);
  hMaxDev->SetLineColor(kRed);
  hMaxDev->Smooth();
  hMaxDev->Draw("same");

  // save histos in output files
  TString OutputFile = "Systematics/SystMultiTrial_" + inputFileName + "_" + CentFT0C[mul] + "_" + ParticleName[!isXi] + "_" + SisSyst + ".root";
  TFile *Write = new TFile(OutputFile, "RECREATE");
  for (int pt = 0; pt < bins; pt++)
  {
    hPtDev[pt]->Write();
    fgaus[pt]->Write();
  }
  hSystMultiTrial->Write();
  hMaxDev->Write();
}
