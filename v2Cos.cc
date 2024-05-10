#include "CommonVar.h"
#include <ROOT/RDataFrame.hxx>
#include "RooRealVar.h"
#include "RooCrystalBall.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "TF1.h"

constexpr int kNPtBins = 18;
constexpr float kPtBounds[2]{0.8, 4.4};

void v2Cos()
{
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
        frame->SetTitle(Form("p_{T} = %f, #Delta#varphi = %f", h->GetXaxis()->GetBinCenter(iPt), hSlice->GetYaxis()->GetBinCenter(iPhi)));
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
    TF1 f1("f1", "[0] * (1 + 2*[1]*cos(2 * x))", 0, TMath::TwoPi());
    for (int iPt{1}; iPt <= h->GetXaxis()->GetNbins(); ++iPt)
    {
      auto hSlice = h->ProjectionY(Form("hSlice_%d", iPt), iPt, iPt);
      hSlice->SetTitle(Form("p_{T} = %f", h->GetXaxis()->GetBinCenter(iPt)));
      hSlice->Fit(&f1, "Q");
      hSlice->SetDrawOption("E");
      hSlice->Write(Form("cosFit_%i", iPt));
      h1.SetBinContent(iPt, f1.GetParameter(1));
      h1.SetBinError(iPt, f1.GetParError(1));
    }
    centDir->cd();
    h1.Write("v2");
  }
}