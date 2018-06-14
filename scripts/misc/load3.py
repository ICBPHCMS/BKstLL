import uproot
import numpy as np
from setTDRStyle import setTDRStyle
import ROOT as r
from math import sqrt
import numpy.ma as ma
import re

################################################################################
# 
class Tree(object) :

    def __init__(self,filename,bkstll) :
        print "Creating tree",filename,"..."
        self.prefix = "branch_"
        self.filename = filename
        self.bkstll = bkstll 
        self.file = uproot.open(filename)
        self.tree = self.file['tree']
        for branch in self.tree.allkeys() : self.setter(branch,self.tree.array(branch))
        self.add() 
        self.branches_withL1 = [ x.replace(self.prefix,"") for x in dir(self) if x.startswith(self.prefix) ]
        self.branches = [ x for x in self.branches_withL1 if "L1" not in x ]
        self.branches_onlyL1 = list(set(self.branches_withL1).difference(set(self.branches)))
        
    def __str__(self) :
        string  = "Filename: {:s}\n".format(self.filename)
        string += "BKstarLL? {:s}\n".format("True" if self.bkstll else "False")
        string += "Branches (noL1):\n"
        for branch in self.branches : string += "  "+branch+"\n"
        return string

    def getter(self,branch) : 
        return getattr(self,self.prefix+branch)

    def setter(self,branch,values) : 
        return setattr(self,self.prefix+branch,values)

    def filter(self,branch,filter) : 
        return ma.array(self.getter(branch),mask=~self.getter(filter)).compressed()

    def add(self) :

        # define mask to identify lead lepton based on gen pT
        lead = [ lp>lm for lp,lm in zip(self.getter("bd_lp_pt"),self.getter("bd_lm_pt"))]

        # add ll_lead and ll_sub variables 
        for var in ["pt","eta","phi","mass","charge"] :
            self.setter("bd_ll_lead_"+var,
                        np.array([ lp if ok else lm for lp,lm,ok in zip(self.getter("bd_lp_"+var),
                                                                        self.getter("bd_lm_"+var),
                                                                        lead) ]))
            self.setter("bd_ll_sub_"+var,
                        np.array([ lp if ~ok else lm for lp,lm,ok in zip(self.getter("bd_lp_"+var),
                                                                         self.getter("bd_lm_"+var),
                                                                         lead) ]))
            self.setter("bd_ll_lead_reco_"+var,
                        np.array([ lp if ok else lm for lp,lm,ok in zip(self.getter("bd_lp_reco_"+var),
                                                                        self.getter("bd_lm_reco_"+var),
                                                                        lead) ]))
            self.setter("bd_ll_sub_reco_"+var,
                        np.array([ lp if ~ok else lm for lp,lm,ok in zip(self.getter("bd_lp_reco_"+var),
                                                                         self.getter("bd_lm_reco_"+var),
                                                                         lead) ]))
            
        # trigger acceptance 
        self.setter("tag_mu_low", ( self.getter("tag_mu_pt") > 7. )  & ( abs(self.getter("tag_mu_eta")) < 2.5 ) )
        self.setter("tag_mu_high", ( self.getter("tag_mu_pt") > 12. ) & ( abs(self.getter("tag_mu_eta")) < 2.5 ) )
        
        # acceptance thresholds
        apply=["standard","open","tight","open_pt"][1]
        pt_cut = {"standard":0.5,"open":0.,"tight":2.,"open_pt":0.}[apply]
        eta_cut = {"standard":2.5,"open":9999.,"tight":0.8,"open_pt":2.5}[apply]
        pt_reco_cut = {"standard":0.5,"open":0.,"tight":2.,"open_pt":0.}[apply]
        eta_reco_cut = {"standard":2.5,"open":9999.,"tight":0.8,"open_pt":2.5}[apply]

        # GEN object acceptance (ignoring trigger)
        self.setter("bd_lp_pass", ( self.getter("bd_lp_pt") > pt_cut ) & ( abs(self.getter("bd_lp_eta")) < eta_cut ) )
        self.setter("bd_lm_pass", ( self.getter("bd_lm_pt") > pt_cut ) & ( abs(self.getter("bd_lm_eta")) < eta_cut ) )
        self.setter("bd_k_pass", ( self.getter("bd_k_pt") > pt_cut ) & ( abs(self.getter("bd_k_eta")) < eta_cut ) )
        self.setter("bd_pi_pass", ( self.getter("bd_pi_pt") > pt_cut ) & ( abs(self.getter("bd_pi_eta")) < eta_cut ) )
        self.setter("bd_all_pass", self.getter("bd_lp_pass") & self.getter("bd_lm_pass") & self.getter("bd_k_pass") & ( self.getter("bd_pi_pass") | ~self.bkstll ) )

        # RECO object acceptance (ignoring trigger and GEN acceptance)
        self.setter("bd_lp_reco_pass", ( self.getter("bd_lp_reco_pt") > pt_reco_cut ) & ( abs(self.getter("bd_lp_reco_eta")) < eta_reco_cut ) )
        self.setter("bd_lm_reco_pass", ( self.getter("bd_lm_reco_pt") > pt_reco_cut ) & ( abs(self.getter("bd_lm_reco_eta")) < eta_reco_cut ) )
        self.setter("bd_k_reco_pass", (  self.getter("bd_k_reco_pt") > pt_reco_cut ) & ( abs( self.getter("bd_k_reco_eta")) < eta_reco_cut ) )
        self.setter("bd_pi_reco_pass", ( self.getter("bd_pi_reco_pt") > pt_reco_cut ) & ( abs(self.getter("bd_pi_reco_eta")) < eta_reco_cut ) )
        self.setter("bd_all_reco_pass", self.getter("bd_lp_reco_pass") & self.getter("bd_lm_reco_pass") & self.getter("bd_k_reco_pass") & ( self.getter("bd_pi_reco_pass") | ~self.bkstll ) )

        # event acceptance given trigger 
        self.setter("evt_low", self.getter("bd_all_pass") & self.getter("tag_mu_low") )
        self.setter("evt_high", self.getter("bd_all_pass") & self.getter("tag_mu_high") )
        self.setter("evt_reco_low", self.getter("bd_all_reco_pass") & self.getter("tag_mu_low") )
        self.setter("evt_reco_high", self.getter("bd_all_reco_pass") & self.getter("tag_mu_high") )

    # print summary of acceptances 
    def acceptance(self) :
        pad = 32
        print "Filename".ljust(pad),self.filename
        print "BKstarLL?".ljust(pad),"True" if self.bkstll else "False"
        print "No requirement on tag muon:"
        print "  N(events)".ljust(pad),self.bd_lp_pt.shape[0]
        print "  N(l+ from probe)".ljust(pad),sum(self.bd_lp_pass)
        print "  N(l- from probe)".ljust(pad),sum(self.bd_lm_pass)
        print "  N(k from probe)".ljust(pad),sum(self.bd_k_pass)
        print "  N(pi from probe)".ljust(pad),sum(self.bd_pi_pass)
        print "  N(all from probe)".ljust(pad),sum(self.bd_all_pass) 
        print "  Acc(all from probe)".ljust(pad),"{:0.3f}".format(sum(self.bd_all_pass)*1. / self.bd_lp_pt.shape[0])
        print "Lower threshold on tag muon:"
        print "  N(tag muon pT > 7 GeV)".ljust(pad),sum(self.tag_mu_low)
        print "  N(l+ from probe)".ljust(pad),sum(self.bd_lp_pass&self.tag_mu_low)
        print "  N(l- from probe)".ljust(pad),sum(self.bd_lm_pass&self.tag_mu_low)
        print "  N(k from probe)".ljust(pad),sum(self.bd_k_pass&self.tag_mu_low)
        print "  N(pi from probe)".ljust(pad),sum(self.bd_pi_pass&self.tag_mu_low)
        print "  N(all from probe)".ljust(pad),sum(self.bd_all_pass&self.tag_mu_low)
        print "  Acc(all from probe | tag muon)".ljust(pad),"{:.3f}".format(sum(self.bd_all_pass&self.tag_mu_low)*1. / sum(self.tag_mu_low)*1.)
        print "Higher threshold on tag muon:"
        print "  N(tag muon pT > 12 GeV)".ljust(pad),sum(self.tag_mu_high)
        print "  N(l+ from probe)".ljust(pad),sum(self.bd_lp_pass&self.tag_mu_high)
        print "  N(l- from probe)".ljust(pad),sum(self.bd_lm_pass&self.tag_mu_high)
        print "  N(k from probe)".ljust(pad),sum(self.bd_k_pass&self.tag_mu_high)
        print "  N(pi from probe)".ljust(pad),sum(self.bd_pi_pass&self.tag_mu_high)
        print "  N(all from probe)".ljust(pad),sum(self.bd_all_pass&self.tag_mu_high)
        print "  Acc(all from probe | tag muon)".ljust(pad),"{:.3f}".format(sum(self.bd_all_pass&self.tag_mu_high)*1. / sum(self.tag_mu_high)*1.)

################################################################################
# 
class Analysis(object) :
    def __init__(self,filenames) :

        self.filenames = [ x[0] for x in filenames ]
        self.bkstll = [ x[1] for x in filenames ]
        self.trees = []
        for filename,bkstll in zip(self.filenames,self.bkstll) :
            self.trees.append( Tree(filename,bkstll) )

    def __str__(self) :
        string  = "Analysis of {:.0f} trees...\n\n".format(len(self.filenames))
        for tree in self.trees : string += tree.__str__()
        return string
    
    def plot(self,
             dir="plots_test",
             norm=True,
             logy=True,
             verbose=False,
             filters=[],#["bd_(.*)(?!reco)(.*)_pt"]
             ) :

        # util
        setTDRStyle()
        binning = {
            "bd_ll_lead_pt":{"xbin":50,"xlow":0.,"xhigh":10.},
            "bd_ll_sub_pt":{"xbin":50,"xlow":0.,"xhigh":10.},
            "bd_ll_lead_reco_pt":{"xbin":50,"xlow":0.,"xhigh":10.},
            "bd_ll_sub_reco_pt":{"xbin":50,"xlow":0.,"xhigh":10.},
            "bd_lp_pt":{"xbin":50,"xlow":0.,"xhigh":10.},
            "bd_lm_pt":{"xbin":50,"xlow":0.,"xhigh":10.},
            "bd_k_pt":{"xbin":50,"xlow":0.,"xhigh":10.},
            "bd_pi_pt":{"xbin":50,"xlow":0.,"xhigh":10.},
            }
        colors = [r.kBlue,r.kRed,r.kGreen,r.kOrange]
        styles = [1,1,9,9]
        widths = [4,2,4,2]

        # branches 
        all_branches = list(set([ branch for tree in self.trees for branch in tree.branches ]))
        if len(filters) > 0 : all_branches = [ b for b in all_branches for f in filters if len(re.findall(f,b)) > 0 ]

        # iterate through branches 
        for ibranch,branch in enumerate(all_branches) : 
            print "Plotting branch",branch,"..."

            # determine tree with most entries in first branch 
            if ibranch == 0 : largest = max([ (index,len(tree.getter(branch))) for index,tree in enumerate(self.trees) ])[0]

            # ignore "event" branch (causes problems)
            if branch == "event" or branch == "lumi" : continue 

            # determine binning, maximum, over using tree with most entries (and remove x > -9. entries!!!)
            values = self.trees[largest].filter(branch,"tag_mu_low")
            values = np.array([ x for x in values if x > -9. ]) #@@ 
            if values.shape[0] > 0 and "bool" in str(type(values[0])) : values = [ 1. if x is True else 0. for x in values ]
            if branch in binning : contents, low_edges = np.histogram(values,bins=np.linspace(binning[branch]["xlow"],
                                                                                              binning[branch]["xhigh"],
                                                                                              binning[branch]["xbin"]+1))
            else : 
                try :
                    _max = max(values)
                    _min = min(values)
                    _diff = np.diff(np.unique(values)).min()
                    _bins = int((_max-_min)/_diff)+1 
                    if _diff < 1. and _bins > 51 : _bins = 51 
                    contents, low_edges = np.histogram(values,bins=np.linspace(_min,_max,_bins))#bins='doane'
                except :
                    print "Error creating np.histogram object!",values[:10]
                    continue 

            underflow = len([ x for x in values if x < low_edges[0] ])
            overflow = len([ x for x in values if x > low_edges[-1] ])
            contents = np.append(contents,[overflow])
            total = sum(contents)
            if norm and total > 0. : contents = [ x*1. / total*1. for x in contents ]
            maximum = max(contents)

            # for each branch, create plot that compares distributions from different trees 
            c1 = r.TCanvas()
            histos = []
            legend = r.TLegend(0.5,0.9-0.05*len(self.trees),0.9,0.9)
            for itree,tree in enumerate(self.trees) :

                # determine histogram content using numpy (quicker!)
                values = tree.filter(branch,"tag_mu_low")
                if values.shape[0] > 0 and "bool" in str(type(values[0])) : values = [ 1. if x is True else 0. for x in values ]
                if branch in binning : contents,_ = np.histogram(values,bins=np.linspace(binning[branch]["xlow"],
                                                                                         binning[branch]["xhigh"],
                                                                                         binning[branch]["xbin"]+1))
                else : contents,_ = np.histogram(values,bins=low_edges)

                name = tree.filename+"_"+branch
                histos.append(r.TH1F(name,branch,len(low_edges)-1,low_edges))
                his = histos[-1]
                for low,val in zip(low_edges,contents) : 
                    his.SetBinContent(his.FindBin(low+1.e-6),val)
                    his.SetBinError(his.FindBin(low+1.e-6),sqrt(val))

                # add overflow to final bin (ignore underflow for now)
                underflow = len([ x for x in values if x < low_edges[0] ])
                overflow = len([ x for x in values if x > low_edges[-1] ])
                his.SetBinContent(his.GetNbinsX(),his.GetBinContent(his.GetNbinsX())+overflow)
                his.SetBinError(his.GetNbinsX(),sqrt(his.GetBinError(his.GetNbinsX())**2.+overflow))
                total = sum(contents)+overflow
                his.SetEntries(total) 

                # normalize?
                if norm and total > 0. : his.Scale(1./total*1.)

                # draw cosmetics
                his.SetTitle(";{:s};a.u.".format(branch))
                his.SetMarkerColor(colors[itree])
                his.SetLineColor(colors[itree])
                his.SetLineStyle(styles[itree])
                his.SetLineWidth(widths[itree])
                if itree == 0 : 
#                    if norm : 
#                        his.SetMaximum(maximum*1./total*1.) if total > 0. else his.SetMaximum(1.0)
#                    else : 
                    his.SetMaximum(maximum*1.1)
                    his.Draw("hist")
                else : 
                    his.Draw("histsame")
                legend.AddEntry(his,tree.filename.replace(".root","")+" ({:.0f})".format(total),"l")

            r.gStyle.SetOptStat(0)
            legend.SetTextSize(0.03)
            legend.Draw("same")
            if logy : c1.SetLogy()
            c1.Update()
            c1.SaveAs("{:s}/{:s}.pdf".format(dir,branch))
            del c1

################################################################################
# print global variables 

if __name__ == '__main__' :

    files = [
#        ("BdKstEEOnlyMuGenFilter.root",True),
        ("BToKEE_decay_first.root",False),
        ("BToKEE_tag_first.root",False),
        ("BToKMuMu_decay_first.root",False),
        ("BToKMuMu_tag_first.root",False),
        ]
    
    a = Analysis(files); #print a 
    a.plot()
#    for tree in a.trees : tree.acceptance(); print 
