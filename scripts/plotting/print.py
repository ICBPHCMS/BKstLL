from load import * 

print "globals",[ x for x,y in globals().items() ]
print "locals",[ x for x,y in locals().items() ] 

print 
print "lp_pt.shape  ",lp_pt.shape
print "bkstll?      ",bkstll
print 
print "sum(mu_low)  ",sum(mu_low)
print "sum(mu_high) ",sum(mu_high)
print "sum(all_pass)",sum(all_pass) 
print 
print "GEN:"
print "low:"
print "sum(lp_pass) ",sum(lp_pass&mu_low)
print "sum(lm_pass) ",sum(lm_pass&mu_low)
print "sum(k_pass)  ",sum(k_pass &mu_low)
print "sum(pi_pass) ",sum(pi_pass&mu_low)
print "sum(evt_low) ",sum(all_pass&mu_low) 
print "high:"
print "sum(lp_pass) ",sum(lp_pass&mu_high)
print "sum(lm_pass) ",sum(lm_pass&mu_high)
print "sum(k_pass)  ",sum(k_pass &mu_high)
print "sum(pi_pass) ",sum(pi_pass&mu_high)
print "sum(evt_high)",sum(all_pass&mu_high) 
print "7 GeV trig:  ",sum(all_pass&mu_low)*1. / sum(mu_low)*1.
print "12 GeV trig: ",sum(all_pass&mu_high)*1. / sum(mu_high)*1.
print 

print "RECO:"
print "low:"
print "sum(lp_reco_pass) ",sum(lp_reco_pass&mu_low)
print "sum(lm_reco_pass) ",sum(lm_reco_pass&mu_low)
print "sum(k_reco_pass)  ",sum(k_reco_pass &mu_low)
print "sum(pi_reco_pass) ",sum(pi_reco_pass&mu_low)
print "sum(evt_reco_low) ",sum(all_reco_pass&mu_low) 
print "high:"
print "sum(lp_reco_pass) ",sum(lp_reco_pass&mu_high)
print "sum(lm_reco_pass) ",sum(lm_reco_pass&mu_high)
print "sum(k_reco_pass)  ",sum(k_reco_pass &mu_high)
print "sum(pi_reco_pass) ",sum(pi_reco_pass&mu_high)
print "sum(evt_reco_high)",sum(all_reco_pass&mu_high) 
print "7 GeV trig:  ",sum(all_reco_pass&mu_low)*1. / sum(mu_low)*1.
print "12 GeV trig: ",sum(all_reco_pass&mu_high)*1. / sum(mu_high)*1.
print
print "RECO EFFS:"
print "(7 GeV trigger)"
print "lp: ",sum(lp_reco_pass&lp_pass&mu_low)*1. / sum(lp_pass&mu_low)*1.
print "lm: ",sum(lm_reco_pass&lm_pass&mu_low)*1. / sum(lm_pass&mu_low)*1.
print "k:  ",sum(k_reco_pass &k_pass&mu_low)*1.  / sum(k_pass&mu_low)*1.
print "pi: ",sum(pi_reco_pass&pi_pass&mu_low)*1. / sum(pi_pass&mu_low)*1. if bkstll else None
print "eff:",(sum(lp_reco_pass&lp_pass&mu_low)*1. / sum(lp_pass&lp_pass&mu_low)*1.)*\
    (sum(lm_reco_pass&lm_pass&mu_low)*1. / sum(lm_pass&lm_pass&mu_low)*1.)*\
    (sum(k_reco_pass&k_pass&mu_low)*1.  / sum(k_pass&k_pass&mu_low)*1.)
print "(12 GeV trigger)"
print "lp: ",sum(lp_reco_pass&lp_pass&mu_high)*1. / sum(lp_pass&mu_high)*1.
print "lm: ",sum(lm_reco_pass&lm_pass&mu_high)*1. / sum(lm_pass&mu_high)*1.
print "k:  ",sum(k_reco_pass &k_pass&mu_high)*1.  / sum(k_pass&mu_high)*1.
print "pi: ",sum(pi_reco_pass&pi_pass&mu_high)*1. / sum(pi_pass&mu_high)*1. if bkstll else None
print "eff:",(sum(lp_reco_pass&lp_pass&mu_high)*1. / sum(lp_pass&lp_pass&mu_high)*1.)*\
    (sum(lm_reco_pass&lm_pass&mu_high)*1. / sum(lm_pass&lm_pass&mu_high)*1.)*\
    (sum(k_reco_pass&k_pass&mu_high)*1.  / sum(k_pass&k_pass&mu_high)*1.)

#print "lp: ",sum(lp_reco_pass&lp_pass)*1. / sum(lp_pass)*1.
#print "lm: ",sum(lm_reco_pass&lm_pass)*1. / sum(lm_pass)*1.
#print "pi: ",sum(pi_reco_pass&pi_pass)*1. / sum(pi_pass)*1. if sum(pi_pass) > 0 else None
#print "k:  ",sum(k_reco_pass&k_pass)*1.  / sum(k_pass)*1.
#print 
#print "lp: ",sum(lp_reco_pass&lp_pass&mu_low)*1. / sum(lp_pass&mu_low)*1.
#print "lm: ",sum(lm_reco_pass&lm_pass&mu_low)*1. / sum(lm_pass&mu_low)*1.
#print "pi: ",sum(pi_reco_pass&pi_pass&mu_low)*1. / sum(pi_pass&mu_low)*1. if sum(pi_pass) > 0 else None
#print "k:  ",sum(k_reco_pass&k_pass&mu_low)*1.  / sum(k_pass&mu_low)*1.
#print 
#print "l:  ",sum(lp_reco_pass&lp_pass)*1. / sum(lp_pass)*1.
#print "ll: ",sum(lp_reco_pass&lp_pass & lm_reco_pass&lm_pass)*1. / sum(lp_pass & lm_pass)*1.
#print "llk:",sum(lp_reco_pass&lp_pass & lm_reco_pass&lm_pass & k_reco_pass&k_pass)*1. / sum(lp_pass & lm_pass & k_pass)*1.
#print 
#print "l:  ",sum(lp_reco_pass & lp_pass)*1. / sum(lp_pass)*1.
#print "ll: ",sum(lp_reco_pass & lp_pass&lm_pass)*1. / sum(lp_pass&lm_pass)*1.
#print "llk:",sum(lp_reco_pass & lp_pass&lm_pass&k_pass)*1. / sum(lp_pass&lm_pass&k_pass)*1.
#print 
#print "l:  ",sum(lm_reco_pass & lm_pass)*1. / sum(lm_pass)*1.
#print "ll: ",sum(lm_reco_pass & lm_pass&lp_pass)*1. / sum(lm_pass&lp_pass)*1.
#print "llk:",sum(lm_reco_pass & lm_pass&lp_pass&k_pass)*1. / sum(lm_pass&lp_pass&k_pass)*1.
#print 
#print "k:  ",sum(k_reco_pass & k_pass)*1. / sum(k_pass)*1.
#print "kl: ",sum(k_reco_pass & k_pass&lm_pass)*1. / sum(k_pass&lm_pass)*1.
#print "kll:",sum(k_reco_pass & k_pass&lm_pass&lp_pass)*1. / sum(k_pass&lm_pass&lp_pass)*1.

