import uproot
import sys
this = sys.modules[__name__]

file = uproot.open("BToKEE_decay_first.root")
tree = file['tree']

#branches = tree.allkeys()
#for branch in branches : setattr(this,branch,tree.array(branch))

for n,x in file.allclasses() : setattr(this,repr(n),x)

print dir(this)

