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
dir = "plots"#"plots_"+now.strftime("%y%m%d_%H%M")
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
    print "plot_lepton_hadron_corr",his.Integral()
    his.DrawNormalized("colz")
    r.gStyle.SetOptStat(0)
    c1.Update()
    c1.SaveAs("{:s}/{:s}/plot_lepton_hadron_corr.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))

################################################################################

def plot_lepton_corr() :

    from bisect import bisect
    import pickle
    file = open('effs.pkl','rb')
    effs = pickle.load(file)
    file.close()
    option = ["gen-gen","ele-ele","mu-mu","ele-trk","mu-trk"][3]
    
    total = sum(mu_low)
    his1 = r.TH2F("ll","ll",20,0.,5.,20,0.,5.) 
    for lp,lm,trig in zip(ll_lead_pt,ll_sub_pt,mu_low) : 
        if trig : 
            weight = 1.
            if option == "ele-ele" :
                keys,vals = effs["genEles_gsfTrks"] 
                vals.append(vals[-1]) # add overflow 
                #print len(keys),len(vals),lp,bisect(keys,lp),keys
                weight *= vals[bisect(keys,lp)]
                weight *= vals[bisect(keys,lm)]
            elif option == "mu-mu" :
                keys,vals = effs["genMuons_genTrks"]
                vals.append(vals[-1]) # add overflow 
                weight *= vals[bisect(keys,lp)]
                weight *= vals[bisect(keys,lm)]
            elif option == "ele-trk" :
                keys,vals = effs["genEles_gsfTrks"]
                vals.append(vals[-1]) # add overflow 
                weight *= vals[bisect(keys,lp)]
                keys,vals = effs["genEles_genTrks"]
                vals.append(vals[-1]) # add overflow 
                weight *= vals[bisect(keys,lm)]
            elif option == "mu-trk" :
                keys,vals = effs["genMuons_genTrks"]
                vals.append(vals[-1]) # add overflow 
                weight *= vals[bisect(keys,lp)]
                keys,vals = effs["genEles_genTrks"]
                vals.append(vals[-1]) # add overflow 
                weight *= vals[bisect(keys,lm)]
            elif option == "gen-gen" : pass
            else : print "UNKNOWN OPTION!",option
            #print weight 
            his1.Fill(lp,lm,weight)
    xbins = his1.GetXaxis().GetNbins()+1
    ybins = his1.GetYaxis().GetNbins()+1
    his2 = r.TH2F("ll_cumu","ll_cumu",20,0.,5.,20,0.,5.) 
    for xbin in range(xbins) : 
        for ybin in range(ybins) : 
            if xbin < ybin : continue 
            his2.SetBinContent(xbin,ybin,his1.Integral(xbin, xbins, ybin, ybins))
    his1.Scale(1./float(total))
    his2.Scale(1./float(total))
    print "total",total
    print "his1",his1.Integral(1, xbins, 1, ybins)
    print "his2",his2.GetBinContent(1,1)
    
    # plot 
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.15)
    c = r.TCanvas()
    his1.SetTitle(";p^{lead}_{T} [GeV];p^{sub}_{T} [GeV]")
    his1.GetZaxis().SetNdivisions(505)
    his1.SetMaximum(0.01)
    print "plot_lepton_corr",his1.Integral()
    r.gStyle.SetPaintTextFormat("4.2f")
    his1.Draw("colztext")
    r.gStyle.SetOptStat(0)
    c.Update()
    c.SaveAs("{:s}/{:s}/plot_lepton_corr.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))
    del his1
    del c

    # plot 
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.15)
    c = r.TCanvas()
    his2.SetTitle(";p^{lead}_{T} threshold [GeV];p^{sub}_{T} threshold [GeV]")
    his2.SetMaximum(1.0)
    his2.GetZaxis().SetNdivisions(505)
    r.gStyle.SetPaintTextFormat("4.2f")
    his2.Draw("colztext")
    r.gStyle.SetOptStat(0)
    c.Update()
    c.SaveAs("{:s}/{:s}/plot_lepton_corr_cumu.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))
    c.SaveAs("{:s}/{:s}/plot_lepton_corr_cumu.root".format(dir,"bkstll" if bkstll == True else "bkll"))
    del his2
    del c

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
    print "plot_hadron_corr",his.Integral()
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
        "pt":{"title":"p_{T} [GeV]","xbins":50,"xlow":0.,"xhigh":10.,"yhigh":0.16},
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

def plot_object_eff(lead=False,combine_leptons=True) :
    setTDRStyle()
    if combine_leptons == True :
        particles = [("ll","",r.kBlue,20),]
    else : 
        if lead == True : particles = [("ll_lead","e^{lead}",r.kBlue,20),("ll_sub","e^{sub}",r.kRed,24)] 
        else : particles = [("lp","l^{#pm}",r.kBlue,20),("lm","l^{#mp}",r.kRed,24)]
    #particles += [("k","K^{#pm}",r.kOrange,21)]
    #if bkstll == True : particles += [("pi","#pi^{#mp}",r.kGreen,25)]
    c1 = r.TCanvas()
    legend = r.TLegend(0.3,0.7-0.07*len(particles),0.5,0.7)
    hists = odict()
    for ip,(name,title,colour,style) in enumerate(particles) : 
        binning = np.append(np.linspace(0.5,5.,9,endpoint=False),np.linspace(5.,11.,7))
        if ip == 0 : hists[name] = r.TEfficiency(name,name,len(binning)-1,binning) 
        his = hists[name]
        if combine_leptons == True : 
            gen_pt = getattr(this,"lp_pt")
            reco_pt = getattr(this,"lp_reco_pt")
            for in_acc,reco,gen,trig in zip(all_pass,reco_pt,gen_pt,mu_low) : 
                if trig and in_acc : his.Fill(1 if reco>0.1 else 0,min(gen,10.))
            gen_pt = getattr(this,"lm_pt")
            reco_pt = getattr(this,"lm_reco_pt")
            for in_acc,reco,gen,trig in zip(all_pass,reco_pt,gen_pt,mu_low) : 
                if trig and in_acc : his.Fill(1 if reco>0.1 else 0,min(gen,10.))
        else :
            gen_pt = getattr(this,"{:s}_pt".format(name))
            reco_pt = getattr(this,"{:s}_reco_pt".format(name))
            for in_acc,reco,gen,trig in zip(all_pass,reco_pt,gen_pt,mu_low) : 
                if trig and in_acc : his.Fill(1 if reco>0.1 else 0,min(gen,10.))
        options = ""
        if ip > 0  : options += "same"
        his.Draw(options)
        r.gPad.Update()
        his.SetTitle(";p^{gen}_{T} [GeV];Efficiency")
        his.GetPaintedGraph().GetYaxis().SetNdivisions(510)
        his.GetPaintedGraph().GetYaxis().SetRangeUser(0.,1.)
        his.GetPaintedGraph().GetXaxis().SetRangeUser(0.5,11.)
        his.SetMarkerColor(colour)
        his.SetMarkerStyle(style)
        his.SetMarkerSize(1.5)
        his.SetLineColor(colour)
        legend.AddEntry(his,title,"pe")
    r.gStyle.SetOptStat(0)
    legend.SetTextSize(0.05)
    if combine_leptons == False : legend.Draw("same")
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

def eff_vs_hadron_pt_cut(acc = False) :
    setTDRStyle()
    c1 = r.TCanvas()
    his = r.TEfficiency("","",51,-50.,5050.) 
    for pt_cut in range(0,5100,100) :
        if acc == True :
            k_pass_tmp  = ( k_pt  > pt_cut/1000. ) & ( abs(k_eta)  < eta_cut )
            pi_pass_tmp = ( pi_pt > pt_cut/1000. ) & ( abs(pi_eta) < eta_cut )
            all_pass_tmp = lp_pass & lm_pass & k_pass_tmp & ( pi_pass_tmp | ~bkstll )
        else: 
            k_pass_tmp  = ( k_reco_pt  > pt_cut/1000. ) & ( abs(k_reco_eta)  < eta_reco_cut )
            pi_pass_tmp = ( pi_reco_pt > pt_cut/1000. ) & ( abs(pi_reco_eta) < eta_reco_cut )
            all_pass_tmp = lp_reco_pass & lm_reco_pass & k_pass_tmp & ( pi_pass_tmp | ~bkstll )
        for in_acc,trig in zip(all_pass_tmp,mu_low) :
            if trig : his.Fill(1 if in_acc == True else 0,pt_cut)
    his.SetTitle(";hadron p^{gen}_{T} cut [MeV];"+"{:s}".format("Acceptance" if acc == True else "#it{A}#times#epsilon"))
    his.Draw("")
    r.gPad.Update()
    his.GetPaintedGraph().GetYaxis().SetNdivisions(510)
    his.GetPaintedGraph().GetYaxis().SetRangeUser(0.,0.6)
    his.GetPaintedGraph().GetXaxis().SetRangeUser(-50.,5050.)
    his.SetMarkerColor(r.kBlack if acc == True else r.kRed)
    his.SetMarkerStyle(20)
    his.SetMarkerSize(1.)
    his.SetLineColor(r.kBlack if acc == True else r.kRed)
    r.gStyle.SetOptStat(0)
    c1.Update()
    c1.SaveAs("{:s}/{:s}/{:s}_vs_hadron_pt_cut.pdf".format(dir,"bkstll" if bkstll == True else "bkll", "acc" if acc == True else "eff"))

################################################################################

def eff_vs_lepton_pt_cut(only_acc=False,use_effs=False) :
    setTDRStyle()
    c1 = r.TCanvas()
    his = r.TEfficiency("","",51,-50.,5050.) 
    for pt_cut in range(0,5100,100) :
        if only_acc == True :
            lp_pass_tmp = ( ll_lead_pt > pt_cut/1000. ) & ( abs(ll_lead_eta) < eta_cut )
            lm_pass_tmp = ( ll_sub_pt > pt_cut/1000. ) & ( abs(ll_sub_eta) < eta_cut )
            #lp_pass_tmp = ( lp_pt > pt_cut/1000. ) & ( abs(lp_eta) < eta_cut )
            #lm_pass_tmp = ( lm_pt > pt_cut/1000. ) & ( abs(lm_eta) < eta_cut )
            all_pass_tmp = lp_pass_tmp & lm_pass_tmp & k_pass & ( pi_pass | ~bkstll )
        else :
            lp_pass_tmp = ( ll_lead_reco_pt > pt_cut/1000. ) & ( abs(ll_lead_reco_eta) < eta_reco_cut )
            lm_pass_tmp = ( ll_sub_reco_pt > pt_cut/1000. ) & ( abs(ll_sub_reco_eta) < eta_reco_cut )
            #lp_pass_tmp = ( lp_reco_pt > pt_cut/1000. ) & ( abs(lp_reco_eta) < eta_reco_cut )
            #lm_pass_tmp = ( lm_reco_pt > pt_cut/1000. ) & ( abs(lm_reco_eta) < eta_reco_cut )
            all_pass_tmp = lp_pass_tmp & lm_pass_tmp & k_reco_pass & ( pi_reco_pass | ~bkstll )
        # need to weight lp_pass_tmp and lm_pass_tmp
        for in_acc,trig,eff,wei in zip(all_pass_tmp,mu_low,ll_effs,ll_weights) :
            if trig : 
                #weight = wei if use_effs else 1
                #his.Fill(weight if in_acc == True else 0,pt_cut)
                if use_effs : 
                    his.FillWeighted(1 if in_acc == True else 0,eff,pt_cut)
                    his.FillWeighted(0,1.-eff,pt_cut)
                else :
                    his.Fill(1 if in_acc == True else 0,pt_cut)

    his.SetTitle(";lepton p^{gen}_{T} cut [MeV];"+"{:s}".format("Acceptance" if only_acc == True else "#it{A}#times#epsilon"))
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
    c1.SaveAs("{:s}/{:s}/{:s}_vs_lepton_pt_cut.pdf".format(dir,"bkstll" if bkstll == True else "bkll", "acc" if only_acc == True else "eff"))
    print "Acc. x Eff.:",his.GetEfficiency(his.FindFixBin(500.))

################################################################################

def acc_vs_pt_cut() :
    setTDRStyle()
    c1 = r.TCanvas()
    his = r.TEfficiency("","",51,-50.,5050.) 
    for pt_cut in range(0,5100,100) :
        lp_pass_tmp = ( lp_pt > pt_cut/1000. ) & ( abs(lp_eta) < eta_cut )
        lm_pass_tmp = ( lm_pt > pt_cut/1000. ) & ( abs(lm_eta) < eta_cut )
        k_pass_tmp = ( k_pt > pt_cut/1000. ) & ( abs(k_eta) < eta_cut )
        all_pass_tmp = lp_pass_tmp & lm_pass_tmp & k_pass_tmp & ( pi_pass | ~bkstll )
        for in_acc,trig in zip(all_pass_tmp,mu_low) :
            if trig : his.Fill(1 if in_acc == True else 0,pt_cut)
    his.SetTitle(";p^{gen}_{T} cut [MeV];Acceptance")
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
    c1.SaveAs("{:s}/acc_vs_pt_cut.pdf".format(dir))
    c1.SaveAs("{:s}/acc_vs_pt_cut.root".format(dir))

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

def plot_bkll() :

    # utility
    from ROOT import TLorentzVector as P4 
    def add(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2) :
        v1 = P4()
        v1.SetPtEtaPhiM(pt1,eta1,phi1,m1)
        v2 = P4()
        v2.SetPtEtaPhiM(pt2,eta2,phi2,m2)
        v3 = v1+v2
        return v3.Pt(),v3.Eta(),v3.Phi(),v3.M()

    # lp+lm 
    ll_pt = []
    ll_eta = []
    ll_phi = []
    ll_m = []
    for pt1,eta1,phi1,m1,pt2,eta2,phi2,m2 in zip(ll_lead_pt, ll_lead_eta, ll_lead_phi, ll_lead_m, 
                                                 ll_sub_pt, ll_sub_eta, ll_sub_phi, ll_sub_m) :
        pt,eta,phi,m = add(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2)
        ll_pt.append(pt)
        ll_eta.append(eta)
        ll_phi.append(phi)
        ll_m.append(m)

    # k+lp+lm
    kll_pt = []
    kll_eta = []
    kll_phi = []
    kll_m = []
    for pt1,eta1,phi1,m1,pt2,eta2,phi2,m2 in zip(ll_pt, ll_eta, ll_phi, ll_m,
                                                 k_pt, k_eta, k_phi, k_m) :
        pt,eta,phi,m = add(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2)
        kll_pt.append(pt)
        kll_eta.append(eta)
        kll_phi.append(phi)
        kll_m.append(m)
        
    # k+pi+lp+lm
    kpll_pt = []
    kpll_eta = []
    kpll_phi = []
    kpll_m = []
    if bkstll == True :
        for pt1,eta1,phi1,m1,pt2,eta2,phi2,m2 in zip(kll_pt, kll_eta, kll_phi, kll_m,
                                                     pi_pt, pi_eta, pi_phi, pi_m) :
            pt,eta,phi,m = add(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2)
            kpll_pt.append(pt)
            kpll_eta.append(eta)
            kpll_phi.append(phi)
            kpll_m.append(m)

    # fill histos 
    his1 = r.TH2F("plot_kll_ll","plot_kll_ll",60,0.,6.,60,0.,6.) 
    his2 = r.TH1F("plot_kll","plot_kll",60,0.,6.) 
    his3 = r.TH1F("plot_ll","plot_ll",60,0.,6.) 
    for llm,kllm,trig in zip(ll_m,kll_m,mu_low) :
        if trig : 
            his1.Fill(kllm,llm,1.)
            his2.Fill(kllm,1.)
            his3.Fill(llm,1.)
    his4 = r.TH2F("plot_kpll_ll","plot_kpll_ll",60,0.,6.,60,0.,6.) 
    his5 = r.TH1F("plot_kpll","plot_kpll",60,0.,6.) 
    for llm,kpllm,trig in zip(ll_m,kpll_m,mu_low) :
        if trig : 
            print kpllm
            his4.Fill(kpllm,llm,1.)
            his5.Fill(kpllm,1.)

    print len(ll_m),len(kll_m),len(kpll_m)
            
    # plots 
    setTDRStyle()
    r.tdrStyle.SetPadRightMargin(0.15)

    if bkstll == False :
        c = r.TCanvas()
        his1.GetZaxis().SetNdivisions(505)
        his1.SetTitle(";m(Kll) [GeV]; m(ll) [GeV]")
        his1.DrawNormalized("colz")
        r.gStyle.SetOptStat(0)
        c.Update()
        c.SaveAs("{:s}/{:s}/plot_kll_ll.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))
        del his1
        del c

        c = r.TCanvas()
        his2.SetTitle(";m(Kll) [GeV]")
        his2.DrawNormalized()
        r.gStyle.SetOptStat(0)
        c.SetLogy(1)
        c.Update()
        c.SaveAs("{:s}/{:s}/plot_kll.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))
        del his2
        del c

    c = r.TCanvas()
    his3.SetTitle(";m(ll) [GeV]")
    his3.DrawNormalized()
    r.gStyle.SetOptStat(0)
    c.Update()
    c.SaveAs("{:s}/{:s}/plot_ll.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))
    del his3
    del c

    if bkstll == True :
        c = r.TCanvas()
        his4.GetZaxis().SetNdivisions(505)
        his4.SetTitle(";m(K#pill) [GeV]; m(ll) [GeV]")
        his4.DrawNormalized("colz")
        r.gStyle.SetOptStat(0)
        c.Update()
        c.SaveAs("{:s}/{:s}/plot_kpll_ll.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))
        del his4
        del c

        c = r.TCanvas()
        his5.SetTitle(";m(K#pill) [GeV]")
        his5.DrawNormalized()
        r.gStyle.SetOptStat(0)
        c.SetLogy(1)
        c.Update()
        c.SaveAs("{:s}/{:s}/plot_kpll.pdf".format(dir,"bkstll" if bkstll == True else "bkll"))
        del his5
        del c

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

#plot_tag_mu_pt() 
#plot_tag_mu_eta() 
#plot_tag_pt_corr() 
#plot_tag_eta_corr() 
#plot_lepton_hadron_corr()
#plot_lepton_corr() 
#plot_hadron_corr() 
#plot_object("pt")
#plot_object("eta")
#plot_object("phi")
#plot_object("pt_lead")
#plot_object("eta_lead")
#plot_object("phi_lead")
#reco_vs_gen_pt() 
#plot_object_eff(True)
#plot_object_scale()
#eff_vs_hadron_pt_cut()
eff_vs_lepton_pt_cut(only_acc=True,use_effs=True)
#acc_vs_pt_cut()
#deltar_vs_pt()
#deltar_vs_pt_binned()
#tag_pt_vs_eta(True)

#plot_lepton_charge_vs_res() 
#plot_tag_deltar()
#plot_bkll() 
