from load import * 
import numpy as np
import ROOT as r 
import bisect
from setTDRStyle import setTDRStyle
setTDRStyle()
import datetime
import os
import sys
this = sys.modules[__name__]

now = datetime.datetime.now()
dir = "plots_"+now.strftime("%y%m%d_%H%M")
if not os.path.exists(dir+"/bkstll") : os.makedirs(dir+"/bkstll")
if not os.path.exists(dir+"/bkll") : os.makedirs(dir+"/bkll")

import math
def deltaPhi( p1, p2):
    '''Computes delta phi, handling periodic limit conditions.'''
    res = p1 - p2
    while res > math.pi:
        res -= 2*math.pi
    while res < -math.pi:
        res += 2*math.pi
    return res
def deltaR2( e1, p1, e2=None, p2=None):
    """Take either 4 arguments (eta,phi, eta,phi) or two objects that have 'eta', 'phi' methods)"""
    if (e2 == None and p2 == None):
        return deltaR2(e1.eta(),e1.phi(), p1.eta(), p1.phi())
    de = e1 - e2
    dp = deltaPhi(p1, p2)
    return de*de + dp*dp
def deltaR( *args ):
    return math.sqrt( deltaR2(*args) )

################################################################################

def plot_tag_mu_pt() :
    overflow = True
    norm = True 
    logy = True
    setTDRStyle()
    c1 = r.TCanvas()
    xbins,xlow,xhigh = 50,0.,50.
    his = r.TH1F("plot_tag_mu_pt","plot_tag_mu_pt",xbins,xlow,xhigh) 
    for pt,eta in zip(tag_mu_pt,tag_mu_eta) : 
        if abs(eta)<2.5 : 
            if ~overflow or pt<xhigh : his.Fill(pt,1.)
            else : his.Fill(xhigh-1.e-6,1.)
    if norm : 
        his.Scale(1./his.Integral(0,his.GetXaxis().GetNbins()+1))
        his.GetYaxis().SetTitle("Normalised counts")
        if logy : his.GetYaxis().SetRangeUser(0.001,1.0)
        else :    his.GetYaxis().SetRangeUser(0.0,1.0)
    else :
        his.GetYaxis().SetTitle("Events / bin")
    his.SetTitle("")
    his.GetXaxis().SetTitle("p^{tag}_{T} [GeV]")
    his.DrawNormalized("hist")
    his.SetLineWidth(2)
    r.gStyle.SetOptStat(0)
    if logy : c1.SetLogy()
    c1.Update()
    c1.SaveAs("{:s}/{:s}/plot_tag_mu_pt.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

################################################################################

def plot_tag_mu_eta() :
    overflow = True
    norm = True 
    logy = True
    setTDRStyle()
    c1 = r.TCanvas()
    xbins,xhigh = 60,3.
    his = r.TH1F("plot_tag_mu_eta","plot_tag_mu_eta",xbins,-xhigh,xhigh) 
    for pt,eta in zip(tag_mu_pt,tag_mu_eta) : 
        if pt > 7 :
            if ~overflow or abs(eta) < xhigh : his.Fill(eta,1.)
            elif eta > xhigh : his.Fill(xhigh-1.e-6,1.)
            elif eta < -xhigh : his.Fill(-xhigh+1.e-6,1.)
            else : print "error",eta
    if norm : 
        his.Scale(1./his.Integral(0,his.GetXaxis().GetNbins()+1))
        his.GetYaxis().SetTitle("Normalised counts")
        if logy : his.GetYaxis().SetRangeUser(0.001,1.0)
        else :    his.GetYaxis().SetRangeUser(0.0,1.0)
    else :
        his.GetYaxis().SetTitle("Events / bin")
    his.SetTitle("")
    his.GetXaxis().SetTitle("#eta^{tag}")
    his.DrawNormalized("hist")
    his.SetLineWidth(2)
    r.gStyle.SetOptStat(0)
    if logy : c1.SetLogy()
    c1.Update()
    c1.SaveAs("{:s}/{:s}/plot_tag_mu_eta.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

################################################################################

def plot_tag_pt_corr() :
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.15)
    c1 = r.TCanvas()
    his = r.TH2F("plot_tag_pt_corr","plot_tag_pt_corr",40,0.,20.,40,0.,20.) 
    for bd,tag,trig in zip(bd_pt,tag_mu_pt,mu_low) : 
        if trig : his.Fill(tag,bd,1.)
    his.SetTitle(";p^{tag}_{T} [GeV];p^{B}_{T} [GeV]")
    his.GetZaxis().SetNdivisions(505)
    his.DrawNormalized("colz")
    r.gStyle.SetOptStat(0)
    c1.Update()
    c1.SaveAs("{:s}/{:s}/plot_tag_pt_corr.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

################################################################################

def plot_tag_eta_corr() :
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.15)
    c1 = r.TCanvas()
    his = r.TH2F("plot_tag_eta_corr","plot_tag_eta_corr",30,-3.,3.,50,-5.,5.) 
    for bd,tag,trig in zip(bd_eta,tag_mu_eta,mu_low) : 
        if trig : his.Fill(tag,bd,1.)
    his.SetTitle(";#eta^{tag};#eta^{B}")
    his.GetZaxis().SetNdivisions(505)
    his.DrawNormalized("colz")
    r.gStyle.SetOptStat(0)
    c1.Update()
    c1.SaveAs("{:s}/{:s}/plot_tag_eta_corr.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

################################################################################

def plot_lepton_hadron_corr() :

    # utility
    from ROOT import TLorentzVector as P4 
    def add(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2) :
        v1 = P4()
        v1.SetPtEtaPhiM(pt1,eta1,phi1,m1)
        v2 = P4()
        v2.SetPtEtaPhiM(pt2,eta2,phi2,m2)
        v3 = v1+v2
        return v3.Pt(),v3.Eta(),v3.Phi(),v3.M()

    # leptons 
    l_pt = []
    l_eta = []
    l_phi = []
    l_m = []
    for pt1,eta1,phi1,m1,pt2,eta2,phi2,m2 in zip(lp_pt, lp_eta, lp_phi, lp_m, 
                                                 lm_pt, lm_eta, lm_phi, lm_m) :
        pt,eta,phi,m = add(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2)
        l_pt.append(pt)
        l_eta.append(eta)
        l_phi.append(phi)
        l_m.append(m)

    # hadrons 
    h_pt = []
    h_eta = []
    h_phi = []
    h_m = []
    if bkstll :
        for pt1,eta1,phi1,m1,pt2,eta2,phi2,m2 in zip(k_pt, k_eta, k_phi, k_m, 
                                                     pi_pt, pi_eta, pi_phi, pi_m) :
            pt,eta,phi,m = add(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2)
            h_pt.append(pt)
            h_eta.append(eta)
            h_phi.append(phi)
            h_m.append(m)
    else :
        h_pt = k_pt
        h_eta = k_eta
        h_phi = k_phi
        h_m = k_m

    # plot 
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.15)
    c1 = r.TCanvas()
    his = r.TH2F("plot_lepton_hadron_corr","plot_lepton_hadron_corr",40,0.,10.,40,0.,10.) 
    for lpt,hpt,trig in zip(l_pt,h_pt,mu_low) :
        if trig : his.Fill(lpt,hpt,1.)
    his.SetTitle(";p^{leptons}_{T} [GeV];p^{hadrons}_{T} [GeV]")
    his.GetZaxis().SetNdivisions(505)
    his.DrawNormalized("colz")
    r.gStyle.SetOptStat(0)
    c1.Update()
    c1.SaveAs("{:s}/{:s}/plot_lepton_hadron_corr.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

################################################################################

def plot_lepton_corr() :
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.15)
    c1 = r.TCanvas()
    his = r.TH2F("ll","ll",40,0.,10.,40,0.,10.) 
    for lp,lm,trig in zip(lp_pt,lm_pt,mu_low) : 
        if trig : his.Fill(lp,lm,1.)
    his.SetTitle(";p^{l+}_{T} [GeV];p^{l-}_{T} [GeV]")
    his.GetZaxis().SetNdivisions(505)
    his.DrawNormalized("colz")
    r.gStyle.SetOptStat(0)
    c1.Update()
    c1.SaveAs("{:s}/{:s}/plot_lepton_corr.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

    #from scipy.stats.stats import pearsonr
    #print "pearson correlation coeff",pearsonr(lp_pt,lm_pt)
    #from numpy import corrcoef
    #print "correlation matrix",corrcoef(lp_pt,lm_pt)

################################################################################

def plot_hadron_corr() :
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.15)
    c1 = r.TCanvas()
    his = r.TH2F("hadron","hadron",40,0.,10.,40,0.,10.) 
    for lp,lm,trig in zip(pi_pt,k_pt,mu_low) : 
        if trig : his.Fill(lp,lm,1.)
    his.SetTitle(";p^{#pi}_{T} [GeV];p^{K}_{T} [GeV]")
    his.GetZaxis().SetNdivisions(505)
    his.DrawNormalized("colz")
    r.gStyle.SetOptStat(0)
    c1.Update()
    c1.SaveAs("{:s}/{:s}/plot_hadron_corr.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

################################################################################

def plot_object(var) :
    setTDRStyle()
    if "lead" in var : lead = True; var = var.replace("lead","").strip("_")
    else : lead = False
    norm = True
    logy = False
    bins = {
        "pt":{"title":"p_{T} [GeV]","xbins":50,"xlow":0.,"xhigh":10.,"yhigh":0.5},
        "eta":{"title":"#eta","xbins":30,"xlow":-5.,"xhigh":5.,"yhigh":0.1},
        "phi":{"title":"phi","xbins":32,"xlow":-3.2,"xhigh":3.2,"yhigh":0.1},
        }
    particles = \
        [("ll_lead","l^{lead}",r.kBlue),("ll_sub","l^{sub}",r.kRed)] if lead == True else \
        [("lp","l^{#pm}",r.kBlue),("lm","l^{#mp}",r.kRed)]
    particles += [("k","K^{#pm}",r.kOrange)]
    if bkstll == True : particles += [("pi","#pi^{#mp}",r.kGreen)]
    c1 = r.TCanvas()
    hists = odict()
    maximum = 0.
    for (name,title,colour) in particles : 
        hists[name] = r.TH1F(name,name,bins[var]["xbins"],bins[var]["xlow"],bins[var]["xhigh"]) 
        his = hists[name]
        vals = getattr(this,"{:s}_{:s}".format(name,var))
        max = 0. 
        for val,trig in zip(vals,mu_low) : 
            if trig : 
                if val > bins[var]["xlow"] and val < bins[var]["xhigh"] : his.Fill(val)
                elif val < bins[var]["xlow"] : his.Fill(bins[var]["xlow"]+1.e-6)
                elif val > bins[var]["xhigh"] : his.Fill(bins[var]["xhigh"]-1.e-6)
                else : print "Error"
        if his.GetMaximum() > maximum : maximum = his.GetMaximum() 
    legend = r.TLegend(0.7,0.9-0.07*len(particles),0.9,0.9)
    for ip,(name,title,colour) in enumerate(particles) : 
        his = hists[name]
        his.SetTitle("")
        his.GetXaxis().SetTitle(bins[var]["title"])
        his.GetYaxis().SetTitle("Arbitrary units")
        his.GetYaxis().SetNdivisions(505)
        if norm : 
            his.Scale(1./his.Integral(0,his.GetXaxis().GetNbins()+1))
            if logy : his.GetYaxis().SetRangeUser(0.001,bins[var]["yhigh"])
            else : his.GetYaxis().SetRangeUser(0.,bins[var]["yhigh"])
        else : 
            his.GetYaxis().SetRangeUser(0.,maximum*1.1)
        his.SetLineWidth(2)
        his.SetLineColor(colour)
        options = "hist"
        if ip > 0  : options += "same"
        his.Draw(options)
        legend.AddEntry(his,title,"l")
    r.gStyle.SetOptStat(0)
    legend.SetTextSize(0.05)
    legend.Draw("same")
    if logy : c1.SetLogy()
    c1.Update()
    c1.SaveAs("{:s}/{:s}/plot_object_{:s}.pdf".format(dir,"bkstll" if bkstll == True else "bkll",var))

################################################################################

def reco_vs_gen_pt() :
    logz = True
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.15)
    particles = [("lp","l^{#pm}",r.kBlue,20),
                 ("lm","l^{#mp}",r.kRed,24),
                 ("k","K^{#pm}",r.kOrange,21),]
    if bkstll == True : 
        particles += [("pi","#pi^{#mp}",r.kGreen,25)]
    for ip,(name,title,colour,style) in enumerate(particles) : 
        c1 = r.TCanvas()
        his = r.TH2F(name,name,40,0.,5.,40,0.,5.) 
        pt_reco = getattr(this,"{:s}_reco_pt".format(name))
        pt_gen = getattr(this,"{:s}_pt".format(name))
        for reco,gen,trig in zip(pt_reco,pt_gen,mu_low) : 
            if trig : his.Fill(gen,reco if reco > 0. else 0.,1.)
        his.Scale(1./his.Integral(0,his.GetXaxis().GetNbins()+1,
                                  0,his.GetYaxis().GetNbins()+1,))
        his.Draw("colz")
        his.SetTitle(";p^{gen}_{T} [GeV];p^{reco}_{T} [GeV]")
        his.SetMinimum(0.00001)
        his.SetMaximum(0.05)
        r.gStyle.SetOptStat(0)
        if logz : c1.SetLogz()
        c1.Update()
        c1.SaveAs("{:s}/{:s}/plot_reco_vs_gen_pt_{:s}.pdf".format(dir,"bkstll" if bkstll == True else "bkll",name))

################################################################################

def plot_object_eff(lead=False) :
    setTDRStyle()
    particles = \
        [("ll_lead","l^{lead}",r.kBlue,20),("ll_sub","l^{sub}",r.kRed,24)] if lead == True else \
        [("lp","l^{#pm}",r.kBlue,20),("lm","l^{#mp}",r.kRed,24)]
    particles += [("k","K^{#pm}",r.kOrange,21)]
    if bkstll == True : particles += [("pi","#pi^{#mp}",r.kGreen,25)]
    c1 = r.TCanvas()
    legend = r.TLegend(0.3,0.7-0.07*len(particles),0.5,0.7)
    hists = odict()
    for ip,(name,title,colour,style) in enumerate(particles) : 
        hists[name] = r.TEfficiency(name,name,40,0.,5.) 
        his = hists[name]
        gen_pt = getattr(this,"{:s}_pt".format(name))
        reco_pt = getattr(this,"{:s}_reco_pt".format(name))
        for in_acc,reco,gen,trig in zip(all_pass,reco_pt,gen_pt,mu_low) : 
            if trig and in_acc : his.Fill(1 if reco>0.5 else 0,gen)
        options = ""
        if ip > 0  : options += "same"
        his.Draw(options)
        r.gPad.Update()
        his.SetTitle(";p^{gen}_{T} [GeV];Efficiency")
        his.GetPaintedGraph().GetYaxis().SetNdivisions(510)
        his.GetPaintedGraph().GetYaxis().SetRangeUser(0.,1.05)
        his.GetPaintedGraph().GetXaxis().SetRangeUser(0.,5.)
        his.SetMarkerColor(colour)
        his.SetMarkerStyle(style)
        his.SetMarkerSize(1.5)
        his.SetLineColor(colour)
        legend.AddEntry(his,title,"pe")
    r.gStyle.SetOptStat(0)
    legend.SetTextSize(0.05)
    legend.Draw("same")
    c1.Update()
    c1.SaveAs("{:s}/{:s}/plot_object_eff.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

################################################################################

def plot_object_scale() :
    logz = True
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.15)
    particles = [("lp","l^{#pm}",r.kBlue,20),
                 ("lm","l^{#mp}",r.kRed,24),
                 ("k","K^{#pm}",r.kOrange,21),]
    if bkstll == True : 
        particles += [("pi","#pi^{#mp}",r.kGreen,25)]
    for ip,(name,title,colour,style) in enumerate(particles) : 
        c1 = r.TCanvas()
        his = r.TH2F(name,name,40,0.,5.,40,0.,2.) 
        pt_reco = getattr(this,"{:s}_reco_pt".format(name))
        pt_gen = getattr(this,"{:s}_pt".format(name))
        for reco,gen,trig in zip(pt_reco,pt_gen,mu_low) : 
            if trig and gen > 0. : his.Fill(gen,reco/gen if reco > 0. else 0.,1.)
        his.Scale(1./his.Integral(0,his.GetXaxis().GetNbins()+1,
                                  0,his.GetYaxis().GetNbins()+1,))
        his.Draw("colz")
        his.SetTitle(";"+title+" p^{gen}_{T} [GeV];p^{reco}_{T}/p^{gen}_{T}")
        his.SetMinimum(0.00001)
        his.SetMaximum(0.05)
        r.gStyle.SetOptStat(0)
        if logz : c1.SetLogz()
        c1.Update()
        c1.SaveAs("{:s}/{:s}/plot_scale_{:s}.pdf".format(dir,"bkstll" if bkstll == True else "bkll",name))

################################################################################

def eff_vs_hadron_pt_cut() :
    setTDRStyle()
    c1 = r.TCanvas()
    his = r.TEfficiency("","",11,-50.,1050.) 
    effs = {}
    for pt_cut in range(0,1100,100) :
        #k_pass_tmp  = ( k_pt  > pt_cut/1000. ) & ( abs(k_eta)  < eta_cut )
        #pi_pass_tmp = ( pi_pt > pt_cut/1000. ) & ( abs(pi_eta) < eta_cut )
        #all_pass_tmp = lp_pass & lm_pass & k_pass_tmp & ( pi_pass_tmp | ~bkstll )
        k_pass_tmp  = ( k_reco_pt  > pt_cut/1000. ) & ( abs(k_reco_eta)  < eta_reco_cut )
        pi_pass_tmp = ( pi_reco_pt > pt_cut/1000. ) & ( abs(pi_reco_eta) < eta_reco_cut )
        all_pass_tmp = lp_reco_pass & lm_reco_pass & k_pass_tmp & ( pi_pass_tmp | ~bkstll )
        for in_acc,trig in zip(all_pass_tmp,mu_low) :
            if trig : his.Fill(1 if in_acc == True else 0,pt_cut)
    his.SetTitle(";hadron p^{gen}_{T} cut [MeV];Acceptance")
    his.Draw("")
    r.gPad.Update()
    his.GetPaintedGraph().GetYaxis().SetNdivisions(510)
    his.GetPaintedGraph().GetYaxis().SetRangeUser(0.,0.6)
    his.GetPaintedGraph().GetXaxis().SetRangeUser(-50.,1050.)
    his.SetMarkerColor(r.kRed)
    his.SetMarkerStyle(20)
    his.SetMarkerSize(1.)
    his.SetLineColor(r.kRed)
    r.gStyle.SetOptStat(0)
    c1.Update()
    c1.SaveAs("{:s}/{:s}/eff_vs_hadron_pt_cut.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

################################################################################

def eff_vs_lepton_pt_cut() :
    setTDRStyle()
    c1 = r.TCanvas()
    his = r.TEfficiency("","",51,-50.,5050.) 
    effs = {}
    for pt_cut in range(0,5100,100) :
#        lp_pass_tmp = ( lp_pt > pt_cut/1000. ) & ( abs(lp_eta) < eta_cut )
#        lm_pass_tmp = ( lm_pt > pt_cut/1000. ) & ( abs(lm_eta) < eta_cut )
#        all_pass_tmp = lp_pass_tmp & lm_pass_tmp & k_pass & ( pi_pass | ~bkstll )
        lp_pass_tmp = ( lp_reco_pt > pt_cut/1000. ) & ( abs(lp_reco_eta) < eta_reco_cut )
        lm_pass_tmp = ( lm_reco_pt > pt_cut/1000. ) & ( abs(lm_reco_eta) < eta_reco_cut )
        all_pass_tmp = lp_pass_tmp & lm_pass_tmp & k_reco_pass & ( pi_reco_pass | ~bkstll )
        for in_acc,trig in zip(all_pass_tmp,mu_low) :
            if trig : his.Fill(1 if in_acc == True else 0,pt_cut)
    his.SetTitle(";lepton p^{gen}_{T} cut [MeV];Acceptance")
    his.Draw("")
    r.gPad.Update()
    his.GetPaintedGraph().GetYaxis().SetNdivisions(510)
    his.GetPaintedGraph().GetYaxis().SetRangeUser(0.,0.6)
    his.GetPaintedGraph().GetXaxis().SetRangeUser(-50.,5050.)
    his.SetMarkerColor(r.kRed)
    his.SetMarkerStyle(20)
    his.SetMarkerSize(1.)
    his.SetLineColor(r.kRed)
    r.gStyle.SetOptStat(0)
    c1.Update()
    c1.SaveAs("{:s}/{:s}/eff_vs_lepton_pt_cut.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

################################################################################

def deltar_vs_pt() :
    logz = True
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.15)
    particles = [("lp","l^{#pm}",r.kBlue,20),
                 ("lm","l^{#mp}",r.kRed,24),
                 ("k","K^{#pm}",r.kOrange,21),]
    if bkstll == True : 
        particles += [("pi","#pi^{#mp}",r.kGreen,25)]
    for ip,(name,title,colour,style) in enumerate(particles) : 
        c1 = r.TCanvas()
        his = r.TH2F(name,name,40,0.,5.,32,0.,3.2) 
        reco_pt = getattr(this,"{:s}_reco_pt".format(name))
        reco_eta = getattr(this,"{:s}_reco_eta".format(name))
        reco_phi = getattr(this,"{:s}_reco_phi".format(name))
        gen_pt = getattr(this,"{:s}_pt".format(name))
        gen_eta = getattr(this,"{:s}_eta".format(name))
        gen_phi = getattr(this,"{:s}_phi".format(name))
        for     trig,in_acc,\
                gpt,geta,gphi,\
                rpt,reta,rphi in zip(mu_low,all_pass,
                                     gen_pt,gen_eta,gen_phi,
                                     reco_pt,reco_eta,reco_phi) : 
            if trig and in_acc : 
                dr = deltaR(reta,rphi,geta,gphi)
                his.Fill(gpt,dr,1.)
        his.Scale(1./his.Integral(0,his.GetXaxis().GetNbins()+1,
                                  0,his.GetYaxis().GetNbins()+1,))
        his.Draw("colz")
        his.SetTitle(";p^{gen}_{T} [GeV];#DeltaR")
        his.SetMinimum(0.00001)
        his.SetMaximum(0.05)
        r.gStyle.SetOptStat(0)
        if logz : c1.SetLogz()
        c1.Update()
        c1.SaveAs("{:s}/{:s}/plot_deltar_vs_pt_{:s}.pdf".format(dir,"bkstll" if bkstll == True else "bkll",name))

################################################################################

def deltar_vs_pt_binned() :
    logy = True
    setTDRStyle()
    particles = [("lp","l^{#pm}"),
                 ("lm","l^{#mp}"),
                 ("k","K^{#pm}"),]
    if bkstll == True : 
        particles += [("pi","#pi^{#mp}",r.kGreen,25)]
    for ip,(name,title) in enumerate(particles) : 
        c1 = r.TCanvas()
        xbins,xlow,xhigh = 30,0.,0.3
        hists = [r.TH1F(name,name,xbins,xlow,xhigh)] 
        bins = [500,1000,2000,5000]
        for bin in bins :
            hists.append(r.TH1F(name+str(bin),name+str(bin),xbins,xlow,xhigh))
        legend = r.TLegend(0.7,0.6,0.9,0.9)
        reco_pt = getattr(this,"{:s}_reco_pt".format(name))
        reco_eta = getattr(this,"{:s}_reco_eta".format(name))
        reco_phi = getattr(this,"{:s}_reco_phi".format(name))
        gen_pt = getattr(this,"{:s}_pt".format(name))
        gen_eta = getattr(this,"{:s}_eta".format(name))
        gen_phi = getattr(this,"{:s}_phi".format(name))
        for     trig,in_acc,\
                gpt,geta,gphi,\
                rpt,reta,rphi in zip(mu_low,all_pass,
                                     gen_pt,gen_eta,gen_phi,
                                     reco_pt,reco_eta,reco_phi) : 
            if trig and in_acc :
                dr = deltaR(reta,rphi,geta,gphi)
                bin = bisect.bisect_left(bins,gpt*1000.)
                hists[0].Fill(dr if dr < xhigh else xhigh-1.-6,1.)
                hists[bin].Fill(dr if dr < xhigh else xhigh-1.-6,1.)
        for ihis,(his,col) in enumerate(zip(hists,[r.kBlack,r.kBlue,r.kRed,r.kGreen,r.kOrange])) :
            his.SetTitle(";#DeltaR;Arbitrary units")
            his.Scale(1./his.Integral(0,his.GetXaxis().GetNbins()+1))
            his.SetMinimum(0.0001)
            his.SetMaximum(2.0)
            his.SetLineWidth(2)
            his.SetLineColor(col)
            options = "hist"
            if ihis > 0  : options += "same"
            his.Draw(options)
            label = "Inclusive" 
            if ihis > 0 : 
                if ihis < len(bins) : 
                    label = "{:3.1f}-{:3.1f} GeV".format(bins[ihis-1]/1000.,
                                                         bins[ihis]/1000.)
                else : 
                    label = ">{:3.1f} GeV".format(bins[ihis-1]/1000.)
            legend.AddEntry(his,label,"l")
        r.gStyle.SetOptStat(0)
        legend.SetTextSize(0.05)
        legend.Draw("same")
        if logy : c1.SetLogy()
        c1.Update()
        c1.SaveAs("{:s}/{:s}/plot_deltar_vs_pt_binned_{:s}.pdf".format(dir,"bkstll" if bkstll == True else "bkll",name))

################################################################################

def tag_pt_vs_eta(reco=False) :
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.2)
    c1 = r.TCanvas()
    numer = r.TH2F("numer","numer",5,7.,27.,5,0.,2.5) 
    denom = r.TH2F("denom","denom",5,7.,27.,5,0.,2.5) 
    acceptance = all_reco_pass if reco else all_pass
    for in_acc,pt,eta in zip(acceptance,tag_mu_pt,tag_mu_eta) : 
        pt = pt if pt < 27 else 27.-1.e-6
        eta = abs(eta) if abs(eta) < 2.5 else 2.5-1.e-6
        denom.Fill(pt,eta,1.)
        if in_acc: numer.Fill(pt,eta,1.)
    his = r.TH2F(numer)
    his.Divide(denom) 
    title = "Acceptance times efficiency" if reco else "Acceptance" 
    his.SetTitle(";p^{tag}_{T} [GeV];|#eta^{tag}|;"+title)
    his.Draw("colztext")
    his.SetMarkerSize(2.)
    r.gStyle.SetPaintTextFormat("4.2g")
    his.GetZaxis().SetRangeUser(0.,0.2 if reco else 0.7)
    r.gStyle.SetOptStat(0)
    c1.Update()
    c1.SaveAs("{:s}/{:s}/plot_tag_pt_vs_eta.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

#################################################################################
#
#def plot_tag_deltar() :
#    setTDRStyle()
#    r.tdrStyle.SetPadRightMargin(0.15)
#    c1 = r.TCanvas()
#    his = r.TH2F("test","test",40,0.,20.,40,0.,4.,) 
#    for tag,dr,trig in zip(tag_mu_pt,tag_mu_dr_bd_lp,mu_low) : 
#        #if trig and gen > 0. : 
#        his.Fill(tag,dr,1.)
#    his.Draw("colz")
#    his.SetTitle(";tag pt;dr")
#    print "integral",his.Integral(0,his.GetXaxis().GetNbins()+1,
#                                  0,his.GetYaxis().GetNbins()+1)
#    r.gStyle.SetOptStat(0)
#    c1.Update()
#    c1.SaveAs("th2d.pdf")
#    #from scipy.stats.stats import pearsonr
#    #print "pearson correlation coeff",pearsonr(lp_pt,lm_pt)
#    from numpy import corrcoef
#    print "correlation matrix",corrcoef(lp_pt,lm_pt)
#
#################################################################################
#
#def plot_lepton_charge_vs_res() :
#    setTDRStyle()
#    r.tdrStyle.SetPadRightMargin(0.15)
#    c1 = r.TCanvas()
#    his = r.TH2F("ll","ll",30,-1.5,1.5,40,0.,2.,) 
#    for reco,gen,ch,trig in zip(lp_reco_pt,lp_pt,lp_reco_ch,mu_low) : 
#        if trig and gen > 0. : his.Fill(ch,reco/gen,1.)
#    his.Draw("colz")
#    his.SetTitle(";p^{l+}_{T} [GeV];p^{l-}_{T} [GeV]")
#    r.gStyle.SetOptStat(0)
#    c1.Update()
#    c1.SaveAs("th2d.pdf")
#    #from scipy.stats.stats import pearsonr
#    #print "pearson correlation coeff",pearsonr(lp_pt,lm_pt)
#    from numpy import corrcoef
#    print "correlation matrix",corrcoef(lp_pt,lm_pt)
#
#################################################################################
#
#def plot_tag_deltar() :
#    setTDRStyle()
#    r.tdrStyle.SetPadRightMargin(0.15)
#    c1 = r.TCanvas()
#    his = r.TH2F("test","test",40,0.,20.,40,0.,4.,) 
#    for tag,dr,trig in zip(tag_mu_pt,tag_mu_dr_bd_lp,mu_low) : 
#        if trig : his.Fill(tag,dr,1.)
#    his.Draw("colz")
#    his.SetTitle(";tag pt;dr")
#    print "integral",his.Integral(0,his.GetXaxis().GetNbins()+1,
#                                  0,his.GetYaxis().GetNbins()+1)
#    r.gStyle.SetOptStat(0)
#    c1.Update()
#    c1.SaveAs("th2d.pdf")
#
#################################################################################

plot_tag_mu_pt() 
#plot_tag_mu_eta() 
#plot_tag_pt_corr() 
#plot_tag_eta_corr() 
#plot_lepton_hadron_corr()
#plot_lepton_corr() 
#plot_hadron_corr() 
#plot_object("pt")
#plot_object("eta")
#plot_object("phi")
plot_object("pt_lead")
#plot_object("eta_lead")
#plot_object("phi_lead")
#reco_vs_gen_pt() 
#plot_object_eff(True)
#plot_object_scale()
#eff_vs_lepton_pt_cut()
#eff_vs_hadron_pt_cut()
#deltar_vs_pt()
#deltar_vs_pt_binned()
#tag_pt_vs_eta(True)

#plot_lepton_charge_vs_res() 
#plot_tag_deltar()
