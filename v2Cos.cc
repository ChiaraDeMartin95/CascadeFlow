#include "CommonVar.h"
#include <ROOT/RDataFrame.hxx>
#include "RooRealVar.h"
#include "RooCrystalBall.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "TF1.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TStyle.h"

constexpr int kNPtBins = 18;
constexpr float kPtBounds[2]{0.8, 4.4};

void v2Cos()
{
  gStyle->SetOptStat(0);
  const int colors[10] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#0033a1")};

  ROOT::EnableImplicitMT();
  ROOT::RDataFrame d1("O2cascanalysis", "TreeForAnalysis/AnalysisResults_trees_LHC23_PbPb_pass3_Train207098_merged.root");
  auto d2 = d1.Define("deltaPhi", "(fPhi - fPsiT0C) < 0 ? fPhi - fPsiT0C + 2 * TMath::Pi() : (fPhi - fPsiT0C > 2 * TMath::Pi() ? fPhi - fPsiT0C - 2 * TMath::Pi() : fPhi - fPsiT0C)").Filter("fBDTResponseXi > 0.98");

  RooRealVar m("m", "m", 1.28, 1.36);
  RooRealVar mu("mu", "mu", 1.32, 1.31, 1.33);
  RooRealVar sigma("sigma", "sigma", 0.002, 0.001, 0.005);
  RooRealVar alpha("alpha", "alpha", 1.5, 0.5, 2.5);
  RooRealVar n("n", "n", 5, 1, 10);
  RooCrystalBall cb("cb", "cb", m, mu, sigma, alpha, n, true);
  RooRealVar tau("tau", "tau", -1.0, -10.0, 0.0);
  RooExponential exp("exp", "exp", m, tau);
  RooRealVar sig("sig", "sigF", 20, 0., 1000000.);
  RooRealVar bkg("bkg", "bkgF", 20, 0., 1000000.);
  RooAddPdf model("model", "model", RooArgList(cb, exp), RooArgList(sig, bkg));

  TFile hepdata("hepdata.root");
  std::array<int,4> tables = {79, 87, 95, 103};

  TFile f("v2Cos.root", "RECREATE");

  int cent[6] = {10, 20, 30, 40, 50, 60};
  for (int iCent{0}; iCent < 5; ++iCent)
  {
    auto centDir = f.mkdir(Form("cent_%d_%d", cent[iCent], cent[iCent + 1]));
    auto h = d2.Filter(Form("fCentFT0C > %d && fCentFT0C < %d",cent[iCent], cent[iCent + 1])).Histo3D({"deltaCos", ";#it{p}_{T} (GeV/#it{c});#Delta#varphi;m (GeV/#it{c}^2)", kNPtBins, kPtBounds[0], kPtBounds[1], 14, 0, TMath::TwoPi(), 60, 1.29, 1.35}, "fPt", "deltaPhi", "fMassXi");

    centDir->mkdir("massFit");
    centDir->cd("massFit");
    TH2D h2("h2", ";#it{p}_{T} (GeV/#it{c});#Delta#varphi;d#it{N}/d#Delta#varphi", kNPtBins, kPtBounds[0], kPtBounds[1], 14, 0, TMath::TwoPi());
    for (int iPt{1}; iPt <= h->GetXaxis()->GetNbins(); ++iPt)
    {
      h->GetXaxis()->SetRange(iPt, iPt);
      auto hSlice = (TH2 *)h->Project3D("yz");
      for (int iPhi{1}; iPhi <= hSlice->GetYaxis()->GetNbins(); ++iPhi)
      {
        auto hSlicePhi = hSlice->ProjectionX(Form("hSlicePhi_%d_%d", iPt, iPhi), iPhi, iPhi);
        RooDataHist dh(Form("dh_%d_%d", iPt, iPhi), Form("dh_%d_%d", iPt, iPhi), RooArgList(m), hSlicePhi);
        model.fitTo(dh);
        auto frame = m.frame();
        frame->SetTitle(Form("p_{T} = %2.2f GeV/#it{c}, #Delta#varphi = %2.2f", h->GetXaxis()->GetBinCenter(iPt), hSlice->GetYaxis()->GetBinCenter(iPhi)));
        dh.plotOn(frame);
        model.plotOn(frame);
        model.plotOn(frame, RooFit::Name("sig"), RooFit::Components(cb), RooFit::LineStyle(kDashed));
        frame->Write(Form("massFit_%d_%d", iPt, iPhi));
        h2.SetBinContent(iPt, iPhi, sig.getVal());
        h2.SetBinError(iPt, iPhi, sig.getError());
      }
    }
    centDir->mkdir("v2Fit");
    centDir->cd("v2Fit");
    TH1D h1(Form("v2_%i", iCent), Form("%i-%i%%;#it{p}_{T} (GeV/#it{c});v_{2}", cent[iCent], cent[iCent + 1]), kNPtBins, kPtBounds[0], kPtBounds[1]);
    TF1 f1("f1", "[0] * (1 + 2*([1] * cos(2 * x) + [2] * cos(3 * x) + [3] * cos(4 * x) + [4] * cos(5 * x)))", 0, TMath::TwoPi());
    f1.SetParNames("Norm", "v2", "v3", "v4", "v5");
    // f1.FixParameter(2, 0);
    // f1.FixParameter(3, 0);
    for (int iPt{1}; iPt <= h->GetXaxis()->GetNbins(); ++iPt)
    {
      auto hSlice = h->ProjectionY(Form("vnFit_%i", iPt), iPt, iPt);
      hSlice->SetTitle(Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c};#Delta#varphi;d#it{N}/d#Delta#varphi", h->GetXaxis()->GetBinLowEdge(iPt), h->GetXaxis()->GetBinUpEdge(iPt)));
      hSlice->Fit(&f1, "Q");
      hSlice->SetDrawOption("E");
      hSlice->Write();
      h1.SetBinContent(iPt, f1.GetParameter(1));
      h1.SetBinError(iPt, f1.GetParError(1));
    }
    centDir->cd();
    h1.Scale(1. / ftcReso[iCent + 1]);
    h1.SetMarkerStyle(20);
    h1.SetMarkerSize(0.8);
    h1.SetMarkerColor(colors[iCent]);
    h1.SetLineColor(colors[iCent]);
    h1.Write("v2");
    if (iCent < tables.size()) {
      auto graph = (TGraphAsymmErrors*)hepdata.Get(Form("Table %d/Graph1D_y1", tables[iCent]));
      h1.SetMinimum(0);
      h1.SetMaximum(0.3);
      graph->SetLineColor(kBlack);
      graph->SetMarkerColor(kBlack);
      graph->SetMarkerStyle(21);
      graph->SetMarkerSize(0.8);
      TCanvas c(Form("compare%i", iCent), Form("compare%i", iCent));
      h1.Draw("E");
      graph->Draw("P");
      c.Write();
      c.SaveAs(Form("compare%i.png", iCent));
    }
    TFile fChiara(Form("ChiaraResults/FitV2_LHC23_PbPb_pass3_Train207098_Xi_BkgParab_Cent%d-%d.root", cent[iCent], cent[iCent + 1]));
    TH1* v2Chiara = (TH1*)fChiara.Get("histoV2");
    v2Chiara->SetMarkerStyle(21);
    v2Chiara->SetMarkerSize(0.8);
    v2Chiara->SetMarkerColor(kBlack);
    TCanvas cChiara("cChiara", "cChiara");
    h1.Draw("E");
    v2Chiara->Draw("P SAME");
    cChiara.Write();
    cChiara.SaveAs(Form("compareChiara%i.png", iCent));
  }
}