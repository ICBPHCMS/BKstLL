import uproot
import numpy as np
from collections import OrderedDict as odict 
import sys
this = sys.modules[__name__]

# ###############################################################################

filename,bkstll = [
    ("../data/BdKstEEOnlyMuGenFilter.root",True),
    ("../data/BToKEE_decay_first.root",False),
    ("../data/BToKEE_tag_first.root",False),
    ("../data/BToKMuMu_decay_first.root",False),
    ("../data/BToKMuMu_tag_first.root",False),
    ][2]

print filename,bkstll
file = uproot.open(filename)
tree = file['tree']
branches = sorted(tree.allkeys(), key=str.lower) 
print ", ".join([ x for x in branches if "L1" not in x ])

# ###############################################################################

# l+
global lp_pt, lp_eta, lp_phi, lp_m, lp_ch
lp_pt, lp_eta, lp_phi, lp_m, lp_ch = tree.arrays(['bd_lp_pt',
                                                  'bd_lp_eta',
                                                  'bd_lp_phi',
                                                  'bd_lp_mass',
                                                  'bd_lp_charge',
                                                  ], outputtype=tuple)

# l-
(lm_pt, lm_eta, lm_phi, lm_m, lm_ch) = tree.arrays(['bd_lm_pt',
                                                    'bd_lm_eta',
                                                    'bd_lm_phi',
                                                    'bd_lm_mass',
                                                    'bd_lm_charge',
                                                    ], outputtype=tuple)

# pion
(pi_pt, pi_eta, pi_phi, pi_m, pi_ch) = tree.arrays(['bd_pi_pt',
                                                    'bd_pi_eta',
                                                    'bd_pi_phi',
                                                    'bd_pi_mass',
                                                    'bd_pi_charge',
                                                    ], outputtype=tuple)

# kaon
(k_pt, k_eta, k_phi, k_m, k_ch) = tree.arrays(['bd_k_pt',
                                               'bd_k_eta',
                                               'bd_k_phi',
                                               'bd_k_mass',
                                               'bd_k_charge',
                                               ], outputtype=tuple)

# l+ RECO
(lp_reco_pt, lp_reco_eta, lp_reco_phi, lp_reco_m, lp_reco_ch) = tree.arrays(['bd_lp_reco_pt',
                                                                             'bd_lp_reco_eta',
                                                                             'bd_lp_reco_phi',
                                                                             'bd_lp_reco_mass',
                                                                             'bd_lp_reco_charge',
                                                                             ], outputtype=tuple)

# l- RECO
(lm_reco_pt, lm_reco_eta, lm_reco_phi, lm_reco_m, lm_reco_ch) = tree.arrays(['bd_lm_reco_pt',
                                                                             'bd_lm_reco_eta',
                                                                             'bd_lm_reco_phi',
                                                                             'bd_lm_reco_mass',
                                                                             'bd_lm_reco_charge',
                                                                             ], outputtype=tuple)


# mask to identify lead lepton based on gen pT
lead = [ lp>lm for lp,lm in zip(lp_pt,lm_pt) ] 

# lead and sublead GEN and RECO leptons, keeping match
#    print "globals",[ x for x,y in globals().items() ]
#    print "locals",[ x for x,y in locals().items() ] 
for var in ["pt","eta","phi","m","ch"] :
    lp_var = getattr(this,"lp_"+var)
    lm_var = getattr(this,"lm_"+var)
    setattr(this,"ll_lead_"+var,np.array([ lp if ok else lm for lp,lm,ok in zip(lp_var,lm_var,lead) ]))
    setattr(this,"ll_sub_"+var,np.array([ lp if ~ok else lm for lp,lm,ok in zip(lp_var,lm_var,lead) ]))
    lp_reco_var = getattr(this,"lp_reco_"+var)
    lm_reco_var = getattr(this,"lm_reco_"+var)
    setattr(this,"ll_lead_reco_"+var,np.array([ lp if ok else lm for lp,lm,ok in zip(lp_reco_var,lm_reco_var,lead) ]))
    setattr(this,"ll_sub_reco_"+var,np.array([ lp if ~ok else lm for lp,lm,ok in zip(lp_reco_var,lm_reco_var,lead) ]))
    
# pion RECO
(pi_reco_pt, pi_reco_eta, pi_reco_phi) = tree.arrays(['bd_pi_reco_pt',
                                                      'bd_pi_reco_eta',
                                                      'bd_pi_reco_phi',
                                                      ], outputtype=tuple)

# kaon RECO
(k_reco_pt, k_reco_eta, k_reco_phi) = tree.arrays(['bd_k_reco_pt',
                                                   'bd_k_reco_eta',
                                                   'bd_k_reco_phi',
                                                   ], outputtype=tuple)

# tag muon 
(tag_mu_pt, tag_mu_eta) = tree.arrays(['tag_mu_pt',
                                       'tag_mu_eta',
                                       ], outputtype=tuple)

# tag muon deltaR
(tag_mu_dr_bd_lp, tag_mu_dr_bd_lm) = tree.arrays(['tag_mu_dr_bd_lp',
                                                  'tag_mu_dr_bd_lm',
                                                  ], outputtype=tuple)

# b -> k*ll
(nbmesons,bd_pt, bd_eta, bd_phi, bd_mass) = tree.arrays(['nbmesons',
                                                         'bd_pt',
                                                         'bd_eta',
                                                         'bd_phi',
                                                         'bd_mass',
                                                         ], outputtype=tuple)

# ###############################################################################

# trigger acceptance 
mu_low  = ( tag_mu_pt > 7. )  & ( abs(tag_mu_eta) < 2.5 )
mu_high = ( tag_mu_pt > 12. ) & ( abs(tag_mu_eta) < 2.5 )

# acceptance thresholds
apply = ["standard","open","tight","open_pt"][0]
pt_cut = {"standard":0.5,"open":0.,"tight":2.,"open_pt":0.}[apply]
eta_cut = {"standard":2.5,"open":9999.,"tight":0.8,"open_pt":2.5}[apply]
pt_reco_cut = {"standard":0.5,"open":0.,"tight":2.,"open_pt":0.}[apply]
eta_reco_cut = {"standard":2.5,"open":9999.,"tight":0.8,"open_pt":2.5}[apply]

# GEN object acceptance (ignoring trigger)
lp_pass = ( lp_pt > pt_cut ) & ( abs(lp_eta) < eta_cut )
lm_pass = ( lm_pt > pt_cut ) & ( abs(lm_eta) < eta_cut )
k_pass  = ( k_pt  > pt_cut ) & ( abs(k_eta)  < eta_cut )
pi_pass = ( pi_pt > pt_cut ) & ( abs(pi_eta) < eta_cut )
all_pass = lp_pass & lm_pass & k_pass & ( pi_pass | ~bkstll )

# RECO object acceptance (ignoring trigger and GEN acceptance)
lp_reco_pass = ( lp_reco_pt > pt_reco_cut ) & ( abs(lp_reco_eta) < eta_reco_cut )
lm_reco_pass = ( lm_reco_pt > pt_reco_cut ) & ( abs(lm_reco_eta) < eta_reco_cut )
k_reco_pass  = (  k_reco_pt > pt_reco_cut ) & ( abs( k_reco_eta) < eta_reco_cut )
pi_reco_pass = ( pi_reco_pt > pt_reco_cut ) & ( abs(pi_reco_eta) < eta_reco_cut )
all_reco_pass = lp_reco_pass & lm_reco_pass & k_reco_pass & ( pi_reco_pass | ~bkstll )

# event acceptance given trigger 
evt_low  = all_pass & mu_low
evt_high = all_pass & mu_high
evt_reco_low  = all_reco_pass & mu_low
evt_reco_high = all_reco_pass & mu_high
