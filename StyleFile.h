// This macro was originally written by:
// chiara.de.martin@cern.ch

void StyleHisto(TH1D &histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRange,
                Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize)
{
  histo.GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)
    histo.GetXaxis()->SetRangeUser(XLow, XUp);
  histo.SetLineColor(color);
  histo.SetMarkerColor(color);
  histo.SetMarkerStyle(style);
  histo.SetMarkerSize(mSize);
  histo.GetXaxis()->SetTitle(titleX);
  histo.GetXaxis()->SetLabelSize(0.05);
  histo.GetXaxis()->SetTitleSize(0.05);
  histo.GetXaxis()->SetTitleOffset(xOffset);
  histo.GetYaxis()->SetTitle(titleY);
  histo.GetYaxis()->SetTitleSize(0.05);
  histo.GetYaxis()->SetLabelSize(0.05);
  histo.GetYaxis()->SetTitleOffset(yOffset);
  histo.SetTitle(title);
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

//void StyleTGraphErrors(TGraphAsymmErrors *tgraph, Int_t color, Int_t style, Float_t mSize, Int_t linestyle){
  //tgraph->SetLineColor(color);
  //tgraph->SetLineWidth(3);
  //tgraph->SetMarkerColor(color);
  //tgraph->SetMarkerStyle(style);
  //tgraph->SetMarkerSize(mSize);
  //tgraph->SetLineStyle(linestyle);
//}

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

