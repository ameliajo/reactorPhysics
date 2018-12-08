

class group:
    def __init__(self,ID,E_low,E_high):
        self.ID     = ID
        self.E_low  = E_low
        self.E_high = E_high

def convertToExpPos(x):
    split = x.split('+')
    return float(split[0])*10.0**float(split[1])

def convertToExpNeg(x):
    split = x.split('-')
    return float(split[0])*10.0**(-1*float(split[1]))


def number(x):
    if not x: return x
    if '+' in x: return convertToExpPos(x)
    if '-' in x[1:]: return convertToExpNeg(x)
    try: return float(x)
    except: return x


def splitTheseGroups(lines):
    groupSplitting, dilution = [], []
    for line in lines:
        dilution.append(line)
        if None in line:
            groupSplitting.append(dilution)
            dilution = []
    return groupSplitting

def checkHeading(reactionHeading):
    if len(reactionHeading) == 9:
        [zaNew,awrNew,numLegndr,nSig0,breakupFlag,nGroups,MF,MT,lineNum] = reactionHeading
    else:
        [zaNew,awrNew,numLegndr,nSig0,breakupFlag,nGroups,MF_MT,lineNum] = reactionHeading
        MF, MT = float(str(MF_MT)[0]), float(str(MF_MT)[1:])

    nSig0,numLegndr,nGroups,lineNum,MF,MT = \
    int(nSig0),int(numLegndr),int(nGroups),int(lineNum),int(MF),int(float(MT))
    if MT != 452:
        assert(zaNew == za);        assert(lineNum == 1); 
        assert(nSig0 == len(Sig0)); assert(nGroups == len(Ebounds)-1)
    return numLegndr


def processLIST(LIST_for_this_group_for_this_reaction,g):
    if len(LIST_for_this_group_for_this_reaction) == 9:
        [temp,zero,numSecPos,IG2LO,nWordsList,groupIndex,MF,MT,lineNum] = \
        LIST_for_this_group_for_this_reaction
    else:
        [temp,zero,numSecPos,IG2LO,nWordsList,groupIndex,MF_MT,lineNum] = \
        LIST_for_this_group_for_this_reaction 
        MF, MT = str(MF_MT)[0], str(MF_MT)[1:]

    numSecPos,IG2LO,nWordsList,groupIndex,lineNum = \
    int(numSecPos),int(IG2LO),int(nWordsList),int(groupIndex),int(lineNum) 
    assert(g+1 == groupIndex and int(zero) == 0)





print("\n")

##############################################################################
# Read in the values from tape26
##############################################################################
header = []
data = []
with open("tape26", "r") as f:
    for line in f:
        split = line.split()
        if split[-2] == '1451':
            split[-3] = split[-3][:-4] # Get rid of the 1451 ID tag in header
            split = split[:-2]+[split[-1]]
            header.append([number(x) for x in split])

        elif len(split) > 4 and split[-4][-4:] == '1395':
            split[-4] = split[-4][:-4]
            split = [None if val == '' else val for val in split]
            data.append([number(x) for x in split])

        # For instances where MF is more than two characters, and accidentally
        # hits our friendly MT
        elif len(split) > 4 and split[-3][-4:] == '1395':
            split[-3] = split[-3][:-4]
            split = [None if val == '' else val for val in split]
            data.append([number(x) for x in split])



        
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

groups = [group(g,Ebounds[g],Ebounds[g+1]) for g in range(nGroups)]



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

    for group in groups: print("SIGT","%.4e"%group.sigT[0],"%.4e"%group.sigT[1],"%.4e"%group.sigT[2],"%.4e"%group.sigT[3],"%.4e"%group.sigT[4],"%.4e"%group.sigT[5],"%.4e"%group.sigT[6])
    print()

    for group in groups:
        print("SIGF","%.4e"%group.sigF[0],"%.4e"%group.sigF[1],"%.4e"%group.sigF[2],"%.4e"%group.sigF[3],"%.4e"%group.sigF[4],"%.4e"%group.sigF[5],"%.4e"%group.sigF[6])
    print()

    for group in groups:
        print("SIGA","%.4e"%group.sigA[0],"%.4e"%group.sigA[1],"%.4e"%group.sigA[2],"%.4e"%group.sigA[3],"%.4e"%group.sigA[4],"%.4e"%group.sigA[5],"%.4e"%group.sigA[6])
    print()

    for group in groups:
        print("NU  ","%.4e"%group.nuBar)




print("\n")

















