
from parseGENDF_help import   *
import subprocess


#subprocess.run(['cp','u238_tape26','tape26'])
#header, data = readFromTape26('9237')

subprocess.run(['cp','u235_tape26','tape26'])
header, data = readFromTape26('9228')



print("\n")


        
##############################################################################
# Format header
##############################################################################

# CONT
[za,awr,zero,nSig0,minus1,nWordsTitle,lineNum] = header[0]
nSig0,nWordsTitle,lineNum = int(nSig0),int(nWordsTitle),int(lineNum)
assert(zero == 0 and minus1 == -1)

# LIST
[temp,zero1,nGroups,nPhoton,nWordsList,zero2,lineNum] = header[1]
nGroups,nPhoton,nWordsList,lineNum = int(nGroups),int(nPhoton),int(nWordsList),int(lineNum)
assert(zero1 == 0 and zero2 == 0)

restOfFirstLIST = [a for b  in header[2:] for a in b[:-1]]

title   = restOfFirstLIST[0]
Sig0    = restOfFirstLIST[1:nSig0+1]
Ebounds = restOfFirstLIST[nSig0+1:nSig0+nGroups+2]

groups = [group(g,Ebounds[g],Ebounds[g+1],Sig0) for g in range(nGroups)]



##############################################################################
# Split into various reactions
##############################################################################

sigT_lines = [line for line in data if line[-2] == 1    ]
sigF_lines = [line for line in data if line[-2] == 18   ]
sigA_lines = [line for line in data if line[-2] == 3102 ]
nu_lines   = [line for line in data if line[-2] == 3452 ]



##############################################################################
# sigT 
##############################################################################
sigT_groupSplitting = splitTheseGroups(sigT_lines)
reactionHeading = sigT_groupSplitting[0].pop(0)
numLegndr = checkHeading(reactionHeading)

for g in range(nGroups):
    processLIST(sigT_groupSplitting[g].pop(0),g)
    LIST_data = [a for b in sigT_groupSplitting[g] for a in b[:-3]]
    assert(LIST_data.pop() == None)
    groups[g].sigT = LIST_data[int(len(LIST_data)*0.5):][0::numLegndr]

#print(groups[0].sigT)
#print(groups[1].sigT)
#print("\n")




##############################################################################
# sigF 
##############################################################################
sigF_groupSplitting = splitTheseGroups(sigF_lines)
reactionHeading = sigF_groupSplitting[0].pop(0)
numLegndr = checkHeading(reactionHeading)

for g in range(nGroups):
    processLIST(sigF_groupSplitting[g].pop(0),g)
    LIST_data = [a for b in sigF_groupSplitting[g] for a in b[:-3]]
    assert(LIST_data.pop() == None)
    groups[g].sigF = LIST_data[int(len(LIST_data)*0.5):][0::numLegndr]

#print(groups[0].sigF)
#print(groups[1].sigF)
#print("\n")




##############################################################################
# sigA 
##############################################################################
sigA_groupSplitting = splitTheseGroups(sigA_lines)
reactionHeading = sigA_groupSplitting[0].pop(0)
numLegndr = checkHeading(reactionHeading)

for g in range(nGroups):
    processLIST(sigA_groupSplitting[g].pop(0),g)
    LIST_data = [a for b in sigA_groupSplitting[g] for a in b[:-2]]
    assert(LIST_data.pop() == None)
    groups[g].sigA = LIST_data[int(len(LIST_data)*0.5):][0::numLegndr]

#print(groups[0].sigA)
#print(groups[1].sigA)
#print("\n")




##############################################################################
# nu-bar
##############################################################################
nu_groupSplitting = splitTheseGroups(nu_lines)
reactionHeading = nu_groupSplitting[0].pop(0)
numLegndr = checkHeading(reactionHeading)

for g in range(nGroups):
    processLIST(nu_groupSplitting[g].pop(0),g)
    LIST_data = [a for b in nu_groupSplitting[g] for a in b[:-2]]
    assert(LIST_data.pop() == None and len(LIST_data)==3)
    groups[g].nuBar = LIST_data[int(len(LIST_data)*0.5):][0]

#print(groups[0].nuBar)
#print(groups[1].nuBar)



verbose = True
if verbose:
    for group in groups: print("SIGT",["%.4e"%g for g in group.sigT])
    print()
    for group in groups: print("SIGF",["%.4e"%g for g in group.sigF])
    print()
    for group in groups: print("SIGA",["%.4e"%g for g in group.sigA])
    print()
    for group in groups: print("NU  ",["%.4e"%group.nuBar])




print("\n")



subprocess.run(['rm','tape26'])














