import uproot

################################################################################
# list of input files 

files = [
    ("BdKstEEOnlyMuGenFilter.root",True),
    ("BToKEE_decay_first.root",False),
    ("BToKEE_tag_first.root",False),
    ("BToKMuMu_decay_first.root",False),
    ("BToKMuMu_tag_first.root",False),
    ]

################################################################################
# utility methods

def prefix(string) : return string.replace('.root','')+'__'

def get(prefix_string,var_string,value=[]) : 
    print __name__
    try :
        return globals()[prefix(prefix_string)+var_string] 
    except :
        globals()[prefix(prefix_string)+var_string] = []
        return globals()[prefix(prefix_string)+var_string]

################################################################################
# method to load branch contents into global variables 

def load(filename,verbose=False) :

    ################################################################################
    # open file 

    file = uproot.open(filename.replace('.root','')+'.root')
    tree = file['tree']
    branches = sorted(tree.allkeys(), key=str.lower) 

    ################################################################################
    # print all branches except L1, then all L1 branches 

    branches_noL1 = [ x for x in branches if "L1" not in x ]
    branches_L1 = list(set(branches).difference(set(branches_noL1)))
    if verbose : 
        print 
        print "Branches (no L1):\n ","\n  ".join(branches_noL1)
        print 
        print "Branches (L1 only):\n ","\n  ".join(branches_L1)

    ################################################################################
    # add global variables containing branch information, with prefix 'filename_'

    for branch in branches_noL1 : globals()[prefix(filename)+branch] = tree.array(branch)
    globals()[filename] = [] # test 

################################################################################
# print global variables 

if __name__ == '__main__' :
    for filename,bkstll in files[:] : 
        print 
        print "Parsing file:",filename,"..."
        load(filename)
        print 
        print "Globals:\n ","\n  ".join([ name for name,obj in globals().items() if name.startswith(prefix(filename)) ])
