

class group:
    def __init__(self,ID,E_low,E_high):
        self.ID     = ID
        self.E_low  = E_low
        self.E_high = E_high
        self.sigT   = False

    def __str__(self):
        return  \
        "ID:  "+str(self.ID)+\
        "    E range:   "+str('%.2E'%self.E_low)+" - "+str('%.2E'%self.E_high)+\
        "    sigT(dil0) = "+str('%.2E'%self.sigT[0])+ \
        "    sigF(dil0) = "+str('%.2E'%self.sigF[0]) \
        if self.sigT and self.sigF else \
        "ID:  "+str(self.ID)+\
        "    E range:   "+str('%.2E'%self.E_low)+" - "+str('%.2E'%self.E_high)+\
        "    sigT(dil0) = "+str('%.2E'%self.sigT[0])+ \
        "    sigF(dil0) = "+str(None)\
        if self.sigT and not self.sigF else \
        "ID:  "+str(self.ID)+\
        "    E range:   "+str('%.2E'%self.E_low)+" - "+str('%.2E'%self.E_high)+\
        "    sigT(dil0) = "+str(None)+ \
        "    sigF(dil0) = "+str('%.2E'%self.sigF[0]) \
        if not self.sigT and self.sigF else \
        "    sigT(dil0) = "+str(None)+ \
        "    sigF(dil0) = "+str(None)




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



##############################################################################
# Read in the values from tape26
##############################################################################
header = []
array = []
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
            array.append([number(x) for x in split])

        elif len(split) > 4 and split[-3][-4:] == '1395':
            split[-3] = split[-3][:-4]
            array.append([number(x) for x in split])
        
##############################################################################
# Format header
##############################################################################
[za,awr,zero,nSig0,minus1,nWordsTitle,lineNum] = header[0]
nSig0,nWordsTitle,lineNum = int(nSig0),int(nWordsTitle),int(lineNum)
assert(zero == 0 and minus1 == -1)
#print(header[0])
print("ZA",za,"AWR",awr,"#Sig0",nSig0,"#Words",nWordsTitle)
print()


[temp,zero1,nGroups,nPhoton,nWordsList,zero2,lineNum] = header[1]
nGroups,nPhoton,nWordsList,lineNum = int(nGroups),int(nPhoton),int(nWordsList),int(lineNum)
assert(zero1 == 0 and zero2 == 0)
#print(header[1])
print("TEMP",temp,"#Groups",nGroups,"#PhotonGroups",nPhoton,"#Words",nWordsList)
print()

restOfTitle = []
for entry in header[2:]:
    restOfTitle += entry[:-1]
title = restOfTitle[0]
Sig0 = restOfTitle[1:nSig0+1]
Ebounds = restOfTitle[nSig0+1:nSig0+nGroups+2]
print("Dilution Vals",Sig0)
print("energy bounds",Ebounds)
print()


groups = []
for g in range(nGroups):
    groups.append(group(g,Ebounds[g],Ebounds[g+1]))





##############################################################################
# Split into various reactions
##############################################################################

sigT_lines = [line for line in array if line[-2] == 1   ]
sigF_lines = [line for line in array if line[-2] == 18  ]
sigG_lines = [line for line in array if line[-2] == 102 ]


##############################################################################
# Split sigT into various dilutions
##############################################################################
sigT_groupSplitting = []
dilution = []
for line in sigT_lines:
    dilution.append(line)
    if None in line:
        sigT_groupSplitting.append(dilution)
        dilution = []

reactionHeading = sigT_groupSplitting[0].pop(0)
[za,awr,numLegndr,nSig0,breakupFlag,nGroups,MF,MT,lineNum] = reactionHeading
nSig0,numLegndr,nGroups,lineNum = int(nSig0),int(numLegndr),int(nGroups),int(lineNum)
print("ZA",za,"AWR",awr,"#Legndr",numLegndr,"#Sig0",nSig0,"#Groups",nGroups)
print()
assert(lineNum == 1)


##############################################################################
# loop through all groups
##############################################################################
for g in range(nGroups):
    [temp,zero,numSecPos,IG2LO,nWordsList,groupIndex,MF,MT,lineNum] = sigT_groupSplitting[g][0]
    numSecPos,IG2LO,nWordsList,groupIndex,lineNum = int(numSecPos),int(IG2LO),int(nWordsList),int(groupIndex),int(lineNum) 
    assert(g+1 == groupIndex)

    fluxSigma = []
    for line in sigT_groupSplitting[g][1:]:
        fluxSigma += line[:-3]
    assert(fluxSigma.pop() == None)

    flux   = fluxSigma[:int(len(fluxSigma)*0.5)]
    sigmaT = fluxSigma[int(len(fluxSigma)*0.5):]
    sigmaT = [sigmaT[i]/flux[i] for i in range(nSig0)]

    groups[g].sigT = sigmaT





##############################################################################
# Split sigF into various dilutions
##############################################################################

sigF_groupSplitting = []
dilution = []
for line in sigF_lines:
    dilution.append(line)
    if None in line:
        sigF_groupSplitting.append(dilution)
        dilution = []


reactionHeading = sigF_groupSplitting[0].pop(0)
[za,awr,numLegndr,nSig0,breakupFlag,nGroups,MF,MT,lineNum] = reactionHeading
nSig0,numLegndr,nGroups,lineNum = int(nSig0),int(numLegndr),int(nGroups),int(lineNum)

##############################################################################
# loop through all groups
##############################################################################
for g in range(nGroups):
    [temp,zero,numSecPos,IG2LO,nWordsList,groupIndex,MF,MT,lineNum] = sigF_groupSplitting[g][0]
    numSecPos,IG2LO,nWordsList,groupIndex,lineNum = int(numSecPos),int(IG2LO),int(nWordsList),int(groupIndex),int(lineNum) 

    assert(g+1 == groupIndex)

    fluxSigma = []
    for line in sigF_groupSplitting[g][1:]:
        fluxSigma += line[:-3]
    assert(fluxSigma.pop() == None)

    flux   = fluxSigma[:int(len(fluxSigma)*0.5)]
    sigmaF = fluxSigma[int(len(fluxSigma)*0.5):]
    sigmaF = [sigmaT[i]/flux[i] for i in range(nSig0)]

    groups[g].sigF = sigmaF





























