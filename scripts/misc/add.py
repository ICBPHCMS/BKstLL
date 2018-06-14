from load2 import files
from load2 import prefix
from load2 import get
from load2 import load
import numpy as np
import sys

################################################################################
# add new variables to the global scope 

def add(filename) :

    ################################################################################
    # add variables for leading and sub-leading leptons 

#    print __name__
#    print dir(sys.modules['load2'])
#    vars = getattr(sys.modules['load2'],filename)
#    print vars 

    # mask to identify lead lepton based on gen pT
    lead = [ lp>lm for lp,lm in zip(get(filename,"bd_lp_pt"),
                                    get(filename,"bd_lm_pt"))]

    # add ll_lead and ll_sub variables 
    for var in ["pt","eta","phi","mass","charge"] :
        lp_var = get(filename,"bd_lp_"+var)
        lm_var = get(filename,"bd_lm_"+var)
        bd_ll_lead = get(filename,"bd_ll_lead_"+var) 
        bd_ll_lead = np.array([ lp if ok else lm for lp,lm,ok in zip(lp_var,lm_var,lead) ])
        bd_ll_sub = get(filename,"bd_ll_sub_"+var)
        bd_ll_sub = np.array([ lp if ~ok else lm for lp,lm,ok in zip(lp_var,lm_var,lead) ])
        
        lp_reco_var = get(filename,"bd_lp_reco_"+var)
        lm_reco_var = get(filename,"bd_lm_reco_"+var)
        bd_ll_lead_reco = get(filename,"bd_ll_lead_reco_"+var)
        bd_ll_lead_reco = np.array([ lp if ok else lm for lp,lm,ok in zip(lp_reco_var,lm_reco_var,lead) ])
        bd_ll_sub_reco = get(filename,"bd_ll_sub_reco_"+var)
        bd_ll_sub_reco = np.array([ lp if ~ok else lm for lp,lm,ok in zip(lp_reco_var,lm_reco_var,lead) ])

    ################################################################################
    # add variables for leading and sub-leading leptons 

#    # trigger acceptance 
#    tag_mu_pt = var(filename,"tag_mu_pt")
#    #tag_mu_pt = get(filename,"tag_mu_pt"]
#    tag_mu_eta = get(filename,"tag_mu_eta"]
#    var(filename,"mu_low") = ( tag_mu_pt > 7. )  & ( abs(tag_mu_eta) < 2.5 )
#    #get(filename,"mu_low"]  = ( tag_mu_pt > 7. )  & ( abs(tag_mu_eta) < 2.5 )
#    get(filename,"mu_high"] = ( tag_mu_pt > 12. ) & ( abs(tag_mu_eta) < 2.5 )



################################################################################
# add global variables 

if __name__ == '__main__' :
    for filename,bkstll in files[1:2] : 
        print 
        print "Parsing file:",filename,"..."
        load(filename)
        add(filename)
        print 
        print "Globals:\n ","\n  ".join([ name for name,obj in globals().items() if name.startswith(prefix(filename)) ])

    print vars()
