import numpy as np
from xlrd import open_workbook
from ROOT import TCanvas, TH1F, TGraph, TGraphErrors, TH1D, TLegend, TFile, TDirectory, gROOT, gStyle

gStyle.SetOptFit(1111)

totalNbOfFiles = 388;
filenamelist = open("filenamelist.txt");

#prepare histograms and root file
rootfile = TFile("Histograms1.root","recreate");
hist_spec_list = []
hist_deri_list = []
hist_comptonedge = TH1F("comptonedge_stable","",40,680,720)
gComptonVsRun = TGraphErrors()
gChi2VsRun = TGraph()
filecount = 0;
for i in range(0,totalNbOfFiles,1):
    filename = filenamelist.readline().rstrip("\r\n")
    #filename = [s[0:12] for s in filename]
    print filename
    hist_spec_name = filename + "_spec" + str(i)
    hist_spec = TH1F(hist_spec_name,"",4096,0,4096)
    hist_spec_list.append(hist_spec)
    hist_deri_name = filename + "_deri" + str(i)
    hist_deri = TH1F(hist_deri_name,"",4096,0,4096)
    hist_deri_list.append(hist_deri)
    wb = open_workbook(filename)
    totalCnts = 0
    for sheet in wb.sheets():
        number_of_rows = sheet.nrows
        number_of_columns = sheet.ncols
        print number_of_rows
        firstbinvalue = sheet.cell(16,1).value
        for row in range(16, number_of_rows):
            value  = (sheet.cell(row,1).value)
            hist_spec_list[i].Fill(row-16,value)
            if row > 16:
                hist_deri_list[i].Fill(row-16+1,value-firstbinvalue)
                firstbinvalue = value
            totalCnts += value
    hist_spec_list[i].SetEntries(totalCnts)
    hist_spec_list[i].GetXaxis().SetRangeUser(0,1000)
    hist_spec_list[i].GetXaxis().SetTitle("Channel number")
    hist_spec_list[i].GetYaxis().SetRangeUser(0,4000)
    hist_spec_list[i].GetYaxis().SetTitle("Counts")
    hist_spec_list[i].Write()
    hist_deri_list[i].SetEntries(totalCnts)
    #hist_deri_list[i].Rebin(8)
    hist_deri_list[i].GetXaxis().SetRangeUser(500,1000)
    hist_deri_list[i].GetXaxis().SetTitle("Channel number")
    #hist_deri_list[i].GetYaxis().SetRangeUser(-50,10)
    hist_deri_list[i].GetYaxis().SetRangeUser(-10,10)
    hist_deri_list[i].GetYaxis().SetTitle("Counts")
    mean = hist_deri_list[i].GetMean()
    rms = hist_deri_list[i].GetRMS()
    hist_deri_list[i].Fit("pol2","Q","",mean-rms,mean+rms*1.5)
    chi2 = hist_deri_list[i].GetFunction("pol2").GetChisquare()
    ndf = hist_deri_list[i].GetFunction("pol2").GetNDF()
#    hist_deri_list[i].Fit("pol2","Q","",580,880)
    title = hist_deri_name + "_" + str(15*i) + "_minutes"
    hist_deri_list[i].SetTitle(title);
    hist_deri_list[i].Write()
    #c = TCanvas()
    #hist_deri_list[i].Draw();
    #c.Print("Animation1.gif+50");
    #c.Close()
    par1 = hist_deri_list[i].GetFunction("pol2").GetParameter(1)
    par2 = hist_deri_list[i].GetFunction("pol2").GetParameter(2)
    err1 = hist_deri_list[i].GetFunction("pol2").GetParError(1)
    err2 = hist_deri_list[i].GetFunction("pol2").GetParError(2)
    comptonEdge = -0.5*par1/par2
    err = comptonEdge*np.sqrt(((err1/par1)**2+(err2/par2)**2))
    gComptonVsRun.SetPoint(i,15.0*i,comptonEdge)
    gComptonVsRun.SetPointError(i,0,err)
    gChi2VsRun.SetPoint(i,15.0*i,chi2/ndf)
    if 2000 < 15.0*i:
        hist_comptonedge.Fill(comptonEdge)
    print chi2
canvas = TCanvas()
canvas.SetName("Results")
canvas.Divide(1,2)
canvas.cd(1)
gComptonVsRun.SetName("ComptonEdge_vs_Time")
gComptonVsRun.SetMarkerStyle(21);
gComptonVsRun.GetXaxis().SetTitle("Time [min]")
gComptonVsRun.GetYaxis().SetTitle("Channel number")
gComptonVsRun.Draw("ACP")
gComptonVsRun.Write()
canvas.cd(2)
gChi2VsRun.SetName("Fit_Chi2_NDF")
gChi2VsRun.SetMarkerStyle(21);
gChi2VsRun.GetXaxis().SetTitle("Time [min]")
gChi2VsRun.GetYaxis().SetTitle("Chi2/NDF")
gChi2VsRun.Draw("ACP")
gChi2VsRun.Write()
hist_comptonedge.GetXaxis().SetTitle("Compton edge (Channel number)")
hist_comptonedge.GetYaxis().SetTitle("Counts")
hist_comptonedge.Write()
rootfile.Close()

