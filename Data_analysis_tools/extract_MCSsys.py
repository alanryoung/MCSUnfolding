from ROOT import TCanvas, TFile, TH1D, TLegend, TText, TLatex
import sys, math
from math import sqrt

class plotgen:
    
    def __init__(self, filename):
        self.fname = filename
        self.desc = self.fname.split("_")
        self.histvarnames = ['thetaX','thetaY','thetaScatt','theta2Scatt']
        self.histstatenames = ['measdataCobb','refconv_GEANT','refconv_Cobb', 'ref', 'GEANT']
        self.histstatedesc = ['Raw Data',
                              'GEANT4 Default MCS', 'Carlisle-Cobb']
        self.histcolors = [1, 4, 2, 6]
        self.histopts  = ['lp','lp','lp','lp']
        self.RMS = {}
        self.RMSErr = {}
        self.Chi2 = {}
        self.RMSsysdiff = {}
        self.RMSsyserr  = {}
        self.sysFiles   = []
        self.eventnorm  = {}
        self.h0integral = {}
        self.useeventnorm = True

        
    def formatHists(self, hist, i):
        # hist.Sumw2()
        hist.SetMarkerColor(self.histcolors[i])
        hist.GetXaxis().SetLabelSize(0.05)
        hist.GetXaxis().SetTitleSize(0.05)
        hist.GetYaxis().SetLabelSize(0.05)
        hist.GetYaxis().SetTitleSize(0.05)
        hist.SetLineColor(self.histcolors[i])
        hist.SetMarkerColor(self.histcolors[i])
        hist.SetLineWidth(3)
        if i > 0: hist.SetMarkerStyle(21+self.histcolors[i])
        else:  hist.SetMarkerStyle(20)
        hist.SetStats(0)


    def addToSysErr(self, i, histl, histh, base, scale, histvar, sysname):
        norm = 0
        if i == 0:
            norm = self.eventnorm[histvar]
        else:
            norm = histl.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]
            # if histvar == 'thetaScatt' or histvar == 'theta2Scatt':
            #     norm /= base.GetXaxis().GetBinWidth(10)/histl.GetXaxis().GetBinWidth(10)
        histl.Scale(1./norm)
        norm = 0
        if i == 0:
            norm = self.eventnorm[histvar]
        else:
            norm = histh.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]
            # if histvar == 'thetaScatt' or histvar == 'theta2Scatt':
            #     norm /= base.GetXaxis().GetBinWidth(10)/histh.GetXaxis().GetBinWidth(10)
        histh.Scale(1./norm)
        if histvar == 'thetaScatt':
            self.RMSsysdiff[histvar][self.histstatenames[i]][sysname] \
                = 1000*(histh.GetMean() - histl.GetMean())
            self.RMSsyserr[histvar][self.histstatenames[i]][sysname] \
                = 1000*math.fabs(histh.GetMean() - histl.GetMean()) * scale
            histl.GetXaxis().SetRangeUser(-0.0,0.06)
            histh.GetXaxis().SetRangeUser(-0.0,0.06)
        elif histvar == 'theta2Scatt':
            self.RMSsysdiff[histvar][self.histstatenames[i]][sysname] \
                = 1000*(sqrt(histh.GetMean()) - \
                        sqrt(histl.GetMean()))
            self.RMSsyserr[histvar][self.histstatenames[i]][sysname] \
                = 1000*math.fabs(sqrt(histh.GetMean()) - \
                        sqrt(histl.GetMean())) * scale
            histl.GetXaxis().SetRangeUser(-0.0,0.0036)
            histh.GetXaxis().SetRangeUser(-0.0,0.0036)
        else:
            histl.GetXaxis().SetRangeUser(-0.06,0.06)
            histh.GetXaxis().SetRangeUser(-0.06,0.06)
            histl.Fit("gaus","Q0","same",-0.045,0.045)
            histh.Fit("gaus","Q0","same",-0.045,0.045)
            #print sysname
	    #print histvar
	    #print histh
	    #print histl
	    self.RMSsysdiff[histvar][self.histstatenames[i]][sysname] \
                = 1000*(histh.GetFunction('gaus').GetParameter('Sigma') - \
                        histl.GetFunction('gaus').GetParameter('Sigma'))
            self.RMSsyserr[histvar][self.histstatenames[i]][sysname] \
                = math.fabs(self.RMSsysdiff[histvar][self.histstatenames[i]][sysname]) * scale 
            
    def addToRMS(self, i, hist, base, resplot, histvar):
        norm = 0
        if i == 0:
            if self.useeventnorm:
                self.eventnorm[histvar]  = hist.GetEntries()
            self.h0integral[histvar] = hist.Integral()
            norm = self.eventnorm[histvar]
            #print histvar, i, norm, hist.Integral(), hist.Integral()/norm
        else:
            norm = hist.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]
            # if histvar == 'thetaScatt' or histvar == 'theta2Scatt':
            #     norm /= base.GetXaxis().GetBinWidth(10)/hist.GetXaxis().GetBinWidth(10)
	hist.Scale(1./norm)
        if histvar == 'thetaScatt':
            self.RMS[histvar][self.histstatenames[i]]    = (1000*hist.GetMean())
            self.RMSErr[histvar][self.histstatenames[i]] = (1000*hist.GetMeanError())
            hist.GetXaxis().SetRangeUser(-0.0,0.06)
            resplot.GetXaxis().SetRangeUser(-0.0,0.06)
        elif histvar == 'theta2Scatt':
            hist.Fit("expo","Q0","same",-0.0,0.002)
            self.RMS[histvar][self.histstatenames[i]]    = (1000/sqrt(-hist.GetFunction('expo').GetParameter('Slope')))
            self.RMSErr[histvar][self.histstatenames[i]] = (1000*hist.GetFunction('expo').GetParError(1)
                                                            /(-hist.GetFunction('expo').GetParameter('Slope'))**(3./2.))
            hist.GetXaxis().SetRangeUser(-0.0,0.0036)
            resplot.GetXaxis().SetRangeUser(-0.0,0.0036)
        else:
            hist.GetXaxis().SetRangeUser(-0.06,0.06)
            resplot.GetXaxis().SetRangeUser(-0.06,0.06)
            hist.Fit("gaus","Q0","same",-0.045,0.045)
            self.RMS[histvar][self.histstatenames[i]]    = (1000*hist.GetFunction('gaus').GetParameter('Sigma'))
            self.RMSErr[histvar][self.histstatenames[i]] = (1000*hist.GetFunction('gaus').GetParError(2))

    def calculateChi2(self, i, hist, base, resplot, histvar, pname):
        
        chi2 = 0
        ndf = 0
        # if i%2 == 0 and i > 0:
        sysLowF  = [TFile(x[0]) for x in self.sysFiles]
        sysLowH  = [f.Get(histvar + '_' + self.histstatenames[i]) for f in sysLowF]

        sysHighF = [TFile(x[1]) for x in self.sysFiles]
        sysHighH = [f.Get(histvar + '_' + self.histstatenames[i]) for f in sysHighF]
        sysScale = [x[2] for x in self.sysFiles]
        syshists = []
        for q in range(len(self.sysFiles)):
            syshists.append(resplot.Clone())
            syshists[-1].SetName(resplot.GetName()+"_"+self.sysFiles[q][3])
            syshists[-1].GetYaxis().SetTitle("Upper - Lower Systematic")
            
        self.RMSsysdiff[histvar][self.histstatenames[i]] = {}
        self.RMSsyserr[histvar][self.histstatenames[i]]  = {}
        self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum2"] = 0
        self.RMSsyserr[histvar][self.histstatenames[i]]["Sum2"] = 0
        self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum"] = 0
        self.RMSsyserr[histvar][self.histstatenames[i]]["Sum"] = 0
        binmax = 0.045
        if histvar == 'theta2Scatt':
            binmax = 0.0036
        for k in range(len(self.sysFiles)):
	    if pname=="Truth" and i==0:
		    continue
            if histvar == 'theta2Scatt':
                sysLowH[k].Rebin(8)
                sysHighH[k].Rebin(8)
            elif histvar == 'thetaScatt':
                sysLowH[k].Rebin(1)
                sysHighH[k].Rebin(1)
            elif self.sysFiles[k][-1] != "TOF":
		#print self.sysFiles[k] 
		#print sysLowH[k] 
		sysLowH[k].Rebin(4)
                sysHighH[k].Rebin(4)
            #print sysLowH[k]

	    self.formatHists(sysLowH[k], i)
            self.formatHists(sysHighH[k], i)
            self.addToSysErr(i, sysHighH[k], sysLowH[k], base, sysScale[k], histvar, self.sysFiles[k][3])
            self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum2"] \
                += self.RMSsysdiff[histvar][self.histstatenames[i]][self.sysFiles[k][3]] *\
                self.RMSsysdiff[histvar][self.histstatenames[i]][self.sysFiles[k][3]] 
            self.RMSsyserr[histvar][self.histstatenames[i]]["Sum2"] \
                += self.RMSsyserr[histvar][self.histstatenames[i]][self.sysFiles[k][3]]  *\
                self.RMSsyserr[histvar][self.histstatenames[i]][self.sysFiles[k][3]]
            self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum"] \
                = sqrt( self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum2"] )
            self.RMSsyserr[histvar][self.histstatenames[i]]["Sum"] \
                = sqrt( self.RMSsyserr[histvar][self.histstatenames[i]]["Sum2"] )
        for j in range(0, hist.GetNbinsX()+1):
	    if pname=="Truth" and i==0:
		    continue
            res = hist.GetBinContent(j) - base.GetBinContent(j)                
            #sys = hsysup[i].GetBinContent(j) - hsysdn[i].GetBinContent(j)
            #sys *= 0.06/0.8
            err2 = 0
	    err2hist = 0
	    if math.fabs(res) > 0:
                err2 = hist.GetBinError(j)**2  # + base.GetBinError(j)**2
		err2hist = base.GetBinError(j)**2
	    err2hist = base.GetBinError(j)**2
            if math.fabs(hist.GetXaxis().GetBinCenter(j)) < 0.6:
                for k in range(len(self.sysFiles)):
                    sys = sysHighH[k].GetBinContent(j) - sysLowH[k].GetBinContent(j)
		    sys *= sysScale[k]
		    if err2 > 0:
			    syshists[k].SetBinContent(j, sys*sys/err2)
			    syshists[k].SetBinError(j, 0.)
                    err2 += sys*sys
		    err2hist += sys*sys
                if math.fabs(hist.GetXaxis().GetBinCenter(j)) < binmax and err2 > 0:
		    chi2 += res*res/err2
		if i==0:
			base.SetBinError(j, sqrt(err2hist))
		hist.SetBinError(j, sqrt(err2))
	        #if pname=="Truth" and i==1:
		#	hist.SetBinError(j, sqrt(err2hist))
	        '''
	        if sys != 0:
		   chi2 += pow(res/sqrt(err2),2)
		   #chi2 += res*sqrt(base.GetBinContent(j))/sqrt(res*res/base.GetBinContent(j))
		else:
                '''
		if err2 > 0:
			resplot.SetBinContent(j, res/math.sqrt(err2))
			resplot.SetBinError(j, 1)
                else:
			resplot.SetBinContent(j, res/1)
			resplot.SetBinError(j, 0)
                
            if math.fabs(hist.GetXaxis().GetBinCenter(j)) < binmax:
                ndf += 1
            #print i, resplot.GetXaxis().GetBinCenter(j), res, chi2, err2, ndf
        c = TCanvas(self.fname[:-5]+'_'+histvar+'_c1')
        c.SetBottomMargin(0.15)
        c.SetTopMargin(0.075)
        
	if pname == 'rawG4' or pname == 'decoG4' or pname == 'refG4':
            for k in range(len(self.sysFiles)):
                syshists[k].Draw("p")
                c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_'+self.sysFiles[k][3]+'_res.eps')
        #print "Chi-square from residuals is ",chi2, " for " ,ndf," from ", self.histstatedesc[i], binmax
        self.Chi2[histvar][self.histstatenames[i]] = [chi2, ndf]
        
    def MCSPlot(self, pname):

        #print self.fname
	f = TFile(self.fname)
        self.RMS = {}
        self.RMSErr = {}
        self.Chi2 = {}
        self.RMSsysdiff = {}
        self.RMSsyserr  = {}
        
        # create a plot for each histvarname
        for histvar in self.histvarnames:
            self.RMS[histvar]     = {}
            self.RMSErr[histvar]  = {}
            self.RMSsysdiff[histvar] = {}
            self.RMSsyserr[histvar]  = {}
            self.Chi2[histvar]    = {}
            names = [histvar + '_' + x for x in self.histstatenames]
            print names
            hists = [f.Get(histvar+'_'+x) for x in self.histstatenames]
            
            print hists[0]
	    hists[0].SetTitle("")
            # hists[3].Scale(norm)
            # self.formatHist(hist[0], 0)
            resplots = [x.Clone() for x in hists]
            resplots[0].SetTitle('')
            resplots[0].GetYaxis().SetTitle("Normalized Residuals")
        
            # if histvar == 'thetaScatt':
            leg = TLegend(0.55,0.73,0.89,0.92)
            leg.SetLineColor(10)
            # else:
            #    leg = TLegend(0.35,0.2,0.65,0.5)
            for i in range(len(self.histstatedesc)):
                hists[i].Sumw2()
                if histvar=='theta2Scatt':
                    hists[i].Rebin(8)
                    resplots[i].Rebin(8)
                elif histvar=='thetaScatt':
                    hists[i].Rebin(1)
                    resplots[i].Rebin(1)
                else:
                    hists[i].Rebin(1)
                    resplots[i].Rebin(1)
                self.formatHists(hists[i], i)
                self.formatHists(resplots[i], i)
                self.addToRMS(i, hists[i], hists[0], resplots[i], histvar)
                if histvar=='theta2Scatt':
                    hists[i].GetYaxis().SetTitle('Probability per '+str(round(1000*1000*hists[i].GetXaxis().GetBinWidth(4),2))+' mrad^{2}')
                else:
                    hists[i].GetYaxis().SetTitle('Probability per '+str(round(1000*hists[i].GetXaxis().GetBinWidth(4),2))+' mrad')
                leg.AddEntry(hists[i], self.histstatedesc[i], self.histopts[i])
		#print hists[0]
                self.calculateChi2(i, hists[i], hists[0], resplots[i], histvar, pname)

            
            c = TCanvas(self.fname[:-5]+'_'+histvar+'_c1')
            if self.desc[0] == 'XePion':
                t1 = TText(0.18,0.885,"MICE ISIS cycle 2015/03")
                t2 = TText(0.18,0.85,"Xe, "+self.desc[1][2:5]+", MAUS v3.1.2")
            else:
                t1 = TText(0.18,0.885,"MICE ISIS cycle 2015/04")
                t2 = TText(0.18,0.85,"LiH, "+self.desc[1][2:5]+", MAUS v3.1.2")
            t1.SetNDC(1)
            t1.SetTextSize(0.04)
            t1.SetTextFont(42)
            t2.SetNDC(1)
            t2.SetTextSize(0.03)
            t2.SetTextFont(42)
            
            hists[0].GetYaxis().SetRangeUser(4e-4,2.0)
            hists[0].SetTitle(";"+hists[0].GetXaxis().GetTitle()+" (radians);"+hists[0].GetYaxis().GetTitle())
            hists[0].Draw('ep')
            
            c.SetBottomMargin(0.15)
            c.SetTopMargin(0.075)

            for h in hists[1:len(self.histstatedesc)]:
		h.Draw('epsame')
            leg.SetTextSize(0.04)
            leg.Draw('same')
            t1.Draw()
            t2.Draw()
            c.SetLogy()
            c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys.eps')
            c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys.root')
            c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys_pq.jpg')
            c.Clear()
            
            c.SetLogy(0)
            resplots[0].GetYaxis().SetRangeUser(-2,2)
            resplots[0].SetTitle(";"+resplots[0].GetXaxis().GetTitle()+" (radians);"+resplots[0].GetYaxis().GetTitle())
            leg.SetX1NDC(0.5)
            leg.SetX2NDC(0.89)
            leg.SetY1NDC(0.2)
            leg.SetY2NDC(0.4)
            resplots[0].Draw("p")
            for r in resplots:
                r.Draw('psame')
            leg.SetTextSize(0.04)
            leg.Draw('same')
            t1.Draw()
            t2.Draw()
            # pblock.Draw()
            c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys_res_T.eps')
            c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys_res_pq.jpg')

        momhist = f.Get("cor_mom")    
        #mom = [momhist.GetMean() + 19.468, momhist.GetMeanError()]
        #if self.fname.find("LiHMuon_03172")  >= 0:
        #    mom = [momhist.GetMean()*1.107 + 1.05, momhist.GetMeanError()]
        #elif self.fname.find("LiHMuon_03200")  >= 0:
        #    mom = [momhist.GetMean()*1.104 + 1.139, momhist.GetMeanError()]
        #elif self.fname.find("LiHMuon_03240")  >= 0:
        #    mom = [momhist.GetMean()*1.17 - 9.41, momhist.GetMeanError()]
        if self.fname.find("LiHMuon_03172")  >= 0:
            mom = [momhist.GetMean(), momhist.GetMeanError()]
        elif self.fname.find("LiHMuon_03200")  >= 0:
            mom = [momhist.GetMean(), momhist.GetMeanError()]
        elif self.fname.find("LiHMuon_03240")  >= 0:
            mom = [momhist.GetMean(), momhist.GetMeanError()]
        rms = [momhist.GetRMS(), momhist.GetRMSError()]
        summary = []
        syssummary = []
        def sigfig(x):
            if math.fabs(x) > 1e-5:
                return int(math.ceil(math.fabs(math.log(math.fabs(x),10))))
            else:
                return 1
        
        # syssummary.append("p (MeV/c) & "+self.histvarnames[0]+"&"+self.histvarnames[1]+"&"+self.histvarnames[3]+"\\\\")

	if pname!="Truth":
        
         for sys in self.sysFiles:
            # if sys[3] == 'Material': stindx = 1
            # else:
            stindx = 0
            # print sys[3], self.histstatenames[stindx]
            
	    difference0 = self.RMSsysdiff[self.histvarnames[0]][self.histstatenames[stindx]][sys[3]]
            difference1 = self.RMSsysdiff[self.histvarnames[1]][self.histstatenames[stindx]][sys[3]]
            difference3 = self.RMSsysdiff[self.histvarnames[3]][self.histstatenames[stindx]][sys[3]]

            syserr0 = self.RMSsyserr[self.histvarnames[0]][self.histstatenames[stindx]][sys[3]]
            syserr1 = self.RMSsyserr[self.histvarnames[1]][self.histstatenames[stindx]][sys[3]]
            syserr3 = self.RMSsyserr[self.histvarnames[3]][self.histstatenames[stindx]][sys[3]]

            relerr0 = syserr0/self.RMS[self.histvarnames[0]][self.histstatenames[0]]
            relerr1 = syserr1/self.RMS[self.histvarnames[1]][self.histstatenames[0]]
            relerr3 = syserr3/self.RMS[self.histvarnames[3]][self.histstatenames[0]]
            
            syssummary.append(str(round(mom[0],sigfig(mom[1])))+"$\pm$"+str(round(mom[1],sigfig(mom[1])))+\
                              " & "+str(round(difference0,sigfig(difference0)))+\
                              " & "+str(round(syserr0,sigfig(syserr0)))+\
                              " & "+str(round(relerr0,sigfig(relerr0)))+\
                              " & "+str(round(difference1,sigfig(difference1)))+\
                              " & "+str(round(syserr1,sigfig(syserr1)))+\
                              " & "+str(round(relerr1,sigfig(relerr1)))+\
                              " & "+str(round(difference3,sigfig(difference3)))+\
                              " & "+str(round(syserr3,sigfig(syserr3)))+\
                              " & "+str(round(relerr3,sigfig(relerr3)))+"\\\\")
            # print syssummary[-1]
            
         syssummary.append(str(round(mom[0],2))+"$\pm$"+str(round(mom[1],2))+\
                          " & "+ str(round(rms[0],2))+"$\pm$"+str(round(rms[1],2))+\
                          " & "+str(round(self.RMSsysdiff[self.histvarnames[0]][self.histstatenames[0]]['Sum'],2))+\
                          " & "+str(round(self.RMSsyserr[self.histvarnames[0]][self.histstatenames[0]]['Sum'],2))+\
                          " & "+str(round(self.RMSsyserr[self.histvarnames[0]][self.histstatenames[0]]['Sum']/self.RMS[self.histvarnames[0]][self.histstatenames[0]],2))+\
                          " & "+str(round(self.RMSsysdiff[self.histvarnames[1]][self.histstatenames[0]]['Sum'],2))+\
                          " & "+str(round(self.RMSsyserr[self.histvarnames[1]][self.histstatenames[0]]['Sum'],2))+\
                          " & "+str(round(self.RMSsyserr[self.histvarnames[1]][self.histstatenames[0]]['Sum']/self.RMS[self.histvarnames[1]][self.histstatenames[0]],2))+\
                          " & "+str(round(self.RMSsysdiff[self.histvarnames[3]][self.histstatenames[0]]['Sum'],2))+
                          " & "+str(round(self.RMSsyserr[self.histvarnames[3]][self.histstatenames[0]]['Sum'],2))+
                          " & "+str(round(self.RMSsyserr[self.histvarnames[3]][self.histstatenames[0]]['Sum']/self.RMS[self.histvarnames[2]][self.histstatenames[0]],2))+"\\\\")
            
        # summary.append("p (MeV/c) &  &"+str(self.histstatenames[0])+" & "+str(self.histstatenames[1])+" & $\chi^{2}$/ndf & "+
        #                +str(self.histstatenames[2])+" & $\chi^{2}$/ndf \\\\")
        # print mom, self.RMS, self.RMSErr, self.Chi2

         for histvar in self.histvarnames:
            summary.append(str(round(mom[0],2))+"$\pm$"+str(round(mom[1],2))+\
                           "& $\ "+histvar+"$ & "+str(round(self.RMS[histvar][self.histstatenames[0]],2))+ \
                           "$\pm$"+str(round(self.RMSErr[histvar][self.histstatenames[0]],2))+ \
                           "$\pm$"+str(round(self.RMSsyserr[histvar][self.histstatenames[0]]["Sum"],2))+ \
                           " & "+str(round(self.RMS[histvar][self.histstatenames[1]],2))+ \
                           "$\pm$"+str(round(self.RMSErr[histvar][self.histstatenames[1]],2))+ \
                           # "$\pm$"+str(round(self.RMSsyserr[histvar][self.histstatenames[1]]["Sum"],2))+\
                           " & "+str(round(self.Chi2[histvar][self.histstatenames[1]][0],1))+ \
                           " / "+ str(self.Chi2[histvar][self.histstatenames[1]][1])+ \
                           " & "+str(round(self.RMS[histvar][self.histstatenames[2]],2))+ \
                           "$\pm$"+str(round(self.RMSErr[histvar][self.histstatenames[2]],2))+ \
                           #"$\pm$"+str(round(self.RMSsyserr[histvar][self.histstatenames[2]]["Sum"],2))+\
                           " & "+str(round(self.Chi2[histvar][self.histstatenames[2]][0],1))+ \
                           " / "+ str(self.Chi2[histvar][self.histstatenames[2]][1])
                           +"\\\\")
            # print summary[-1]
            
         f.Close()
         return [summary, syssummary]
    #  print "+++++++++++++++ systematics +++++++++++++++++++++++"
    # for i in range(len(RMS)):
    #     print " & ", histvarnames[i]," & ",round(RMSsysdiff[i][0],2)," & ",round(RMSsyserr[i][0],3)," & ",round(RMSsysdiff[i][1],2)," & ",round(RMSsyserr[i][1],3)," & ",round(RMSsysdiff[i][2],2)," & ",round(RMSsyserr[i][2],3),"\\\\"
        
if __name__ == '__main__':
    # baseSet = sys.argv[1]
    raw172 = plotgen("LiHMuon_03172.root")
    raw172.sysFiles = [["./materials/material_0/LiHMu_3172_res_lim0.891.root",
                        "./materials/material_0/LiHMu_3172_res_lim1.109.root",
                        sqrt(0.004/0.218),"Material"],
                       ["./fiducial/LiHMu_3172_R_4.root",
                        "./fiducial/LiHMu_3172_R_6.root",
                       .478/20.,"Fid. Radius"]]
    
    raw200 = plotgen("LiHMuon_03200.root")
    raw200.sysFiles = [["./materials/material_0/LiHMu_3200_res_lim0.891.root",
                        "./materials/material_0/LiHMu_3200_res_lim1.109.root",
                        sqrt(0.004/0.218),"Material"],
                       ["./fiducial/LiHMu_3200_R_4.root",
                        "./fiducial/LiHMu_3200_R_6.root",
                        0.478/20.,"Fid. Radius"]]
    
    raw240 = plotgen("LiHMuon_03240.root")
    raw240.sysFiles = [["./materials/material_0/LiHMu_3240_res_lim0.891.root",
                        "./materials/material_0/LiHMu_3240_res_lim1.109.root",
                        sqrt(0.004/0.218),"Material"],
                       ["./fiducial/LiHMu_3240_R_4.root",
                        "./fiducial/LiHMu_3240_R_6.root",
                        0.478/20.,"Fid. Radius"]]
    
    [out172, sys172] = raw172.MCSPlot("raw")
    [out200, sys200] = raw200.MCSPlot("raw")
    [out240, sys240] = raw240.MCSPlot("raw")
    ''' 
    piraw240 = plotgen("XePion_3240_0.root")
    piraw240.sysFiles = [["TOFsys/XePion_3240_tof_lim26.6_u28.0.root",
                          "TOFsys/XePion_3240_tof_lim27.4_u28.8.root",
                          128./800.,"TOF"],
                         ["material_0/XePion_3240_res_lim0.891.root",
                          "material_0/XePion_3240_res_lim1.109.root",
                          sqrt(0.010/0.218),"Material"],
                       ["../angdef/angdef_0/LiHMu_3172_45.root",
                        "../angdef/angdef_0/LiHMu_3172_135.root",
                        1,"angdef"],
                         ["fiducial_0/XePion_3240_R_5_G_0.root",
                          "fiducial_0/XePion_3240_R_5_G_2.root",
                          0.5/10,"Fid. Pitch"],
                         ["fiducial_0/XePion_3240_R_4_G_1.root",
                          "fiducial_0/XePion_3240_R_6_G_1.root",
                          0.478/20.,"Fid. Radius"]]

    piraw240.histstatenames = ['measdataCobb','refconv_truthsim','refconv_Cobb', 'ref', 'truthsim']
    [piout240, pisys240] = piraw240.MCSPlot("raw")
    '''
    
    print "\n"
    raw172.histstatenames = ['measdataGEANT','recoGEANT','GEANT']
    raw172.histstatedesc = ['Raw Data','Deconvolved, GEANT', 'GEANT Model']
    raw200.histstatenames = ['measdataGEANT','recoGEANT','GEANT']
    raw200.histstatedesc = ['Raw Data','Deconvolved, GEANT', 'GEANT Model']
    raw240.histstatenames = ['measdataGEANT','recoGEANT','GEANT']
    raw240.histstatedesc = ['Raw Data','Deconvolved, GEANT', 'GEANT Model']
    raw172.useeventnorm = True
    raw200.useeventnorm = True
    raw240.useeventnorm = True
    raw172.MCSPlot("rawG4")
    raw200.MCSPlot("rawG4")
    raw240.MCSPlot("rawG4")
    
    print "\n"
    raw172.histstatenames = ['measdataCobb','recoCobb','Cobb']
    raw172.histstatedesc = ['Raw Data','Deconv., Carlisle-Cobb', 'Carlisle-Cobb Model']
    raw200.histstatenames = ['measdataCobb','recoCobb','Cobb']
    raw200.histstatedesc = ['Raw Data','Deconv., Carlisle-Cobb', 'Carlisle-Cobb Model']
    raw240.histstatenames = ['measdataCobb','recoCobb','Cobb']
    raw240.histstatedesc = ['Raw Data','Deconv., Carlisle-Cobb', 'Carlisle-Cobb Model']
    raw172.MCSPlot("rawCC")
    raw200.MCSPlot("rawCC")
    raw240.MCSPlot("rawCC")
    
    print "\n"
    raw172.histstatenames = ['measdataGEANT','ref','GEANT']
    raw172.histstatedesc = ['Raw LiH Data','Empty AFC', 'GEANT Model conv. Empty']
    raw200.histstatenames = ['measdataGEANT','ref','GEANT']
    raw200.histstatedesc = ['Raw LiH Data','Empty AFC', 'GEANT Model conv Empty']
    raw240.histstatenames = ['measdataGEANT','ref','GEANT']
    raw240.histstatedesc = ['Raw LiH Data','Empty AFC', 'GEANT Model conv Empty']
    # piraw240.histstatenames = ['measdatatruthsim','ref','truthsim']
    # piraw240.histstatedesc = ['Raw Data','Empty AFC', 'GEANT Model']
    [ref172, refsys172] = raw172.MCSPlot("refG4")
    [ref200, refsys200] = raw200.MCSPlot("refG4")
    [ref240, refsys240] = raw240.MCSPlot("refG4")
    # [piref240, pirefsys240] = piraw240.MCSPlot("refG4")
    # piraw240.histstatenames = ['measdatatruthsim','ref','Cobb']
    # piraw240.histstatedesc = ['Raw Data','Empty AFC', 'Carlisle-Cobb Model']
    # piraw240.MCSPlot("refCC")
    
    print "\n"
    raw172.histstatenames = ['recoGEANT','GEANT','Cobb']
    raw172.histstatedesc = ['Deconvolved, GEANT', 'GEANT Model', 'Carlisle-Cobb Model']
    raw200.histstatenames = ['recoGEANT','GEANT','Cobb']
    raw200.histstatedesc = ['Deconvolved, GEANT', 'GEANT Model', 'Carlisle-Cobb Model']
    raw240.histstatenames = ['recoGEANT','GEANT','Cobb']
    raw240.histstatedesc = ['Deconvolved, GEANT', 'GEANT Model', 'Carlisle-Cobb Model']
    # piraw240.histstatenames = ['recotruthsim','truthsim','Cobb']
    # piraw240.histstatedesc = ['Deconvolved, GEANT', 'GEANT Model', 'Carlisle-Cobb Model']
    raw172.useeventnorm = False
    raw200.useeventnorm = False
    raw240.useeventnorm = False
    # piraw240.useeventnorm = False
    [G4out172, G4sys172] = raw172.MCSPlot("decoG4")
    [G4out200, G4sys200] = raw200.MCSPlot("decoG4")
    [G4out240, G4sys240] = raw240.MCSPlot("decoG4")
    # [piG4out240, piG4sys240] = piraw240.MCSPlot("decoG4")
    #raw172.histstatenames = ['recoMoliere','Moliere','Cobb']
    #raw172.histstatedesc = ['Deconvolved, Moliere', 'Moliere Model', 'Carlisle-Cobb Model']
    #raw200.histstatenames = ['recoMoliere','Moliere','Cobb']
    #raw200.histstatedesc = ['Deconvolved, Moliere', 'Moliere Model', 'Carlisle-Cobb Model']
    #raw240.histstatenames = ['recoMoliere','Moliere','Cobb']
    #raw240.histstatedesc = ['Deconvolved, Moliere', 'Moliere Model', 'Carlisle-Cobb Model']
    # piraw240.histstatenames = ['recotruthsim','truthsim','Cobb']
    # piraw240.histstatedesc = ['Deconvolved, Moliere', 'Moliere Model', 'Carlisle-Cobb Model']
    raw172.useeventnorm = False
    raw200.useeventnorm = False
    raw240.useeventnorm = False
    # piraw240.useeventnorm = False
    [Moliereout172, Molieresys172] = raw172.MCSPlot("decoG4")
    [Moliereout200, Molieresys200] = raw200.MCSPlot("decoG4")
    [Moliereout240, Molieresys240] = raw240.MCSPlot("decoG4")
    print "\n"
    raw172.histstatenames = ['recoCobb','GEANT','Cobb']
    raw172.histstatedesc = ['Deconv., Carlisle-Cobb','GEANT Model', 'Carlisle-Cobb Model']
    raw200.histstatenames = ['recoCobb','GEANT','Cobb']
    raw200.histstatedesc = ['Deconv., Carlisle-Cobb','GEANT Model', 'Carlisle-Cobb Model']
    raw240.histstatenames = ['recoCobb','GEANT','Cobb']
    raw240.histstatedesc = ['Deconv., Carlisle-Cobb','GEANT Model', 'Carlisle-Cobb Model']
    # piraw240.histstatenames = ['recoCobb','truthsim','Cobb']
    # piraw240.histstatedesc = ['Deconv., Carlisle-Cobb','GEANT Model', 'Carlisle-Cobb Model']
    [CCout172, CCsys172] = raw172.MCSPlot("decoCC")
    [CCout200, CCsys200] = raw200.MCSPlot("decoCC")
    [CCout240, CCsys240] = raw240.MCSPlot("decoCC")
    # [piCCout240, piCCsys240] = piraw240.MCSPlot("decoCC")

    print "Raw deco table" 
    print '\hline'
    print '\hline'
    print out172[0]
    print out172[1]
    print '\hline'
    print out200[0]
    print out200[1]
    print '\hline'
    print out240[0]
    print out240[1]
    # print '\hline'
    # print piout240[0]
    # print piout240[1]

    print '\hline'
    print '\hline'
    print out172[3]
    print out200[3]
    print out240[3]
    # print piout240[3]
    print '\hline'
    print '\hline\n'

    print '\hline'
    print '\hline'
    print ref172[0]
    print ref172[1]
    print '\hline'
    print ref200[0]
    print ref200[1]
    print '\hline'
    print ref240[0]
    print ref240[1]
    # print '\hline'
    # print piref240[0]
    # print piref240[1]

    print '\hline'
    print '\hline'
    print ref172[3]
    print ref200[3]
    print ref240[3]
    # print piref240[3]
    print '\hline'
    print '\hline\n'

    
    print "GEANT4 deco table" 
    print '\hline'
    print '\hline'
    print G4out172[0]
    print G4out172[1]
    print '\hline'
    print G4out200[0]
    print G4out200[1]
    print '\hline'
    print G4out240[0]
    print G4out240[1]
    # print '\hline'
    # print piG4out240[0]
    # print piG4out240[1]

    print '\hline'
    print '\hline'
    print G4out172[3]
    print G4out200[3]
    print G4out240[3]
    # print piG4out240[3]
    print '\hline'
    print '\hline\n'

    print "Moliere deco table" 
    print '\hline'
    print '\hline'
    print Moliereout172[0]
    print Moliereout172[1]
    print '\hline'
    print Moliereout200[0]
    print Moliereout200[1]
    print '\hline'
    print Moliereout240[0]
    print Moliereout240[1]
    # print '\hline'
    # print piMoliereout240[0]
    # print piMoliereout240[1]

    print '\hline'
    print '\hline'
    print Moliereout172[3]
    print Moliereout200[3]
    print Moliereout240[3]
    # print piG4out240[3]
    print '\hline'
    print '\hline\n'
    
    print "Cobb-Carlisle deco table" 
    print '\hline'
    print '\hline'
    print CCout172[0]
    print CCout172[1]
    print '\hline'
    print CCout200[0]
    print CCout200[1]
    print '\hline'
    print CCout240[0]
    print CCout240[1]
    # print '\hline'
    # print piCCout240[0]
    # print piCCout240[1]

    print '\hline'
    print '\hline'
    print CCout172[3]
    print CCout200[3]
    print CCout240[3]
    # print piCCout240[3]
    print '\hline'
    print '\hline'

    
    print '\hline'
    print '\hline'
   
    print 'raw sys'
    for i in range(len(raw172.sysFiles)):
        print raw172.sysFiles[i][3]
        print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
        print sys172[i]
        print sys200[i]
        print sys240[i]
        # print pisys240[i]
        print '\hline\n'
    print '\hline'
    
    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print sys172[-1]
    print sys200[-1]
    print sys240[-1]
    # print pisys240[-1]
    print '\hline\n'
    print '\hline'
    print 'G4 sys'
    for i in range(len(raw172.sysFiles)):
        print raw172.sysFiles[i][3]
        print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
        print G4sys172[i]
        print G4sys200[i]
        print G4sys240[i]
        # print piG4sys240[i]
        print '\hline\n'
    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print sys172[-1]
    print sys200[-1]
    print sys240[-1]
    # print pisys240[-1]
    print '\hline\n'
    print '\hline'
#    print 'Moliere sys'
#    for i in range(len(raw172.sysFiles)):
#        print raw172.sysFiles[i][3]
#        print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
#        print Molieresys172[i]
#        print Molieresys200[i]
#        print Molieresys240[i]
        # print piG4sys240[i]
#        print '\hline\n'
    print '\hline'
    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print G4sys172[-1]
    print G4sys200[-1]
    print G4sys240[-1]
    # print piG4sys240[-1]
    print '\hline\n'
    print '\hline'
    print 'CC sys'
    for i in range(len(raw172.sysFiles)):
        print raw172.sysFiles[i][3]
        print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
        print CCsys172[i]
        print CCsys200[i]
        print CCsys240[i]
        print '\hline\n'
    #print raw172.RMSsyserr
    #print raw200.RMSsyserr
    #print raw240.RMSsyserr
     
    
    
