

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
[za,awr,zero,nSig0,minus1,nWordsTitle,lineNum] = header[0]
nSig0,nWordsTitle,lineNum = int(nSig0),int(nWordsTitle),int(lineNum)
assert(zero == 0 and minus1 == -1)
#print(header[0])
#print("ZA",za,"AWR",awr,"#Sig0",nSig0,"#Words",nWordsTitle)
#print()


[temp,zero1,nGroups,nPhoton,nWordsList,zero2,lineNum] = header[1]
nGroups,nPhoton,nWordsList,lineNum = int(nGroups),int(nPhoton),int(nWordsList),int(lineNum)
assert(zero1 == 0 and zero2 == 0)
#print(header[1])
#print("TEMP",temp,"#Groups",nGroups,"#PhotonGroups",nPhoton,"#Words",nWordsList)
#print()

restOfTitle = []
for entry in header[2:]:
    restOfTitle += entry[:-1]
title = restOfTitle[0]
Sig0 = restOfTitle[1:nSig0+1]
Ebounds = restOfTitle[nSig0+1:nSig0+nGroups+2]
#print("Dilution Vals",Sig0)
#print("energy bounds",Ebounds)
#print()


groups = []
for g in range(nGroups):
    groups.append(group(g,Ebounds[g],Ebounds[g+1]))





##############################################################################
# Split into various reactions
##############################################################################

sigT_lines = [line for line in data if line[-2] == 1    ]
sigF_lines = [line for line in data if line[-2] == 18   ]
sigA_lines = [line for line in data if line[-2] == 3102 ]
nu_lines   = [line for line in data if line[-2] == 3452 ]



def checkHeading(reactionHeading):
    [zaNew,awrNew,numLegndr,nSig0,breakupFlag,nGroups,MF,MT,lineNum] = reactionHeading
    nSig0,numLegndr,nGroups,lineNum,MF,MT = int(nSig0),int(numLegndr),int(nGroups),int(lineNum),int(MF),int(MT)
    #print("ZA",za,"AWR",awr,"#Legndr",numLegndr,"#Sig0",nSig0,"#Groups",nGroups,"MF",MF,"MT",MT)
    #print()
    assert(zaNew == za)
    assert(lineNum == 1); 
    assert(nSig0 == len(Sig0))
    assert(nGroups == len(Ebounds)-1)
    return numLegndr



##############################################################################
# Split sigT into various dilutions
##############################################################################
sigT_groupSplitting = splitTheseGroups(sigT_lines)

reactionHeading = sigT_groupSplitting[0].pop(0)
numLegndr = checkHeading(reactionHeading)


##############################################################################

for g in range(nGroups):
    #[temp,zero,numSecPos,IG2LO,nWordsList,groupIndex,MF,MT,lineNum] = sigT_groupSplitting[g][0]
    [temp,zero,numSecPos,IG2LO,nWordsList,groupIndex,MF,MT,lineNum] = sigT_groupSplitting[g].pop(0)
    numSecPos,IG2LO,nWordsList,groupIndex,lineNum = int(numSecPos),int(IG2LO),int(nWordsList),int(groupIndex),int(lineNum) 
    #print("TEMP",temp,"#Vals",nWordsList,"Group#",groupIndex,"MF",MF,"MT",MT,"line#",lineNum)
    assert(g+1 == groupIndex and int(zero) == 0)
    LIST_data = []
    for line in sigT_groupSplitting[g]:
        LIST_data += line[:-3]
    assert(LIST_data.pop() == None)

    sigmaT = LIST_data[int(len(LIST_data)*0.5):][0::numLegndr]

    groups[g].sigT = sigmaT

#print(groups[0].sigT)
#print(groups[1].sigT)



#print("\n")




##############################################################################
# Split sigF into various dilutions
##############################################################################

sigF_groupSplitting = splitTheseGroups(sigF_lines)

reactionHeading = sigF_groupSplitting[0].pop(0)
[za,awr,numLegndr,nSig0,breakupFlag,nGroups,MF,MT,lineNum] = reactionHeading
nSig0,numLegndr,nGroups,lineNum = int(nSig0),int(numLegndr),int(nGroups),int(lineNum)
#print("ZA",za,"AWR",awr,"#Legndr",numLegndr,"#Sig0",nSig0,"#Groups",nGroups,"MF",MF,"MT",MT)
#print()

##############################################################################

for g in range(nGroups):
    [temp,zero,numSecPos,IG2LO,nWordsList,groupIndex,MF,MT,lineNum] = sigF_groupSplitting[g][0]
    numSecPos,IG2LO,nWordsList,groupIndex,lineNum = int(numSecPos),int(IG2LO),int(nWordsList),int(groupIndex),int(lineNum) 
    #print("TEMP",temp,"#Vals",nWordsList,"Group#",groupIndex,"MF",MF,"MT",MT,"line#",lineNum)

    assert(g+1 == groupIndex and int(zero) == 0)

    fluxSigma = []
    for line in sigF_groupSplitting[g][1:]:
        fluxSigma += line[:-3]
    assert(fluxSigma.pop() == None)

    flux   = fluxSigma[:int(len(fluxSigma)*0.5)]
    sigmaF = fluxSigma[int(len(fluxSigma)*0.5):][0::numLegndr]

    groups[g].sigF = sigmaF



#print(groups[0].sigF)
#print(groups[1].sigF)



#print("\n")




##############################################################################
# Split sigA into various dilutions
##############################################################################

sigA_groupSplitting = splitTheseGroups(sigA_lines)

reactionHeading = sigA_groupSplitting[0].pop(0)
[za,awr,numLegndr,nSig0,breakupFlag,nGroups,MF_MT,lineNum] = reactionHeading
MF, MT = str(MF_MT)[0], str(MF_MT)[1:]
nSig0,numLegndr,nGroups,lineNum = int(nSig0),int(numLegndr),int(nGroups),int(lineNum)
#print("ZA",za,"AWR",awr,"#Legndr",numLegndr,"#Sig0",nSig0,"#Groups",nGroups,"MF",MF,"MT",MT)
#print()

##############################################################################

for g in range(nGroups):
    [temp,zero,numSecPos,IG2LO,nWordsList,groupIndex,MF_MT,lineNum] = sigA_groupSplitting[g][0]
    MF, MT = str(MF_MT)[0], str(MF_MT)[1:]
    numSecPos,IG2LO,nWordsList,groupIndex,lineNum = int(numSecPos),int(IG2LO),int(nWordsList),int(groupIndex),int(lineNum) 
    #print("TEMP",temp,"#Vals",nWordsList,"Group#",groupIndex,"MF",MF,"MT",MT,"line#",lineNum)

    assert(g+1 == groupIndex and int(zero) == 0)

    fluxSigma = []
    for line in sigA_groupSplitting[g][1:]:
        fluxSigma += line[:-2]
    assert(fluxSigma.pop() == None)

    flux   = fluxSigma[:int(len(fluxSigma)*0.5)]
    sigmaA = fluxSigma[int(len(fluxSigma)*0.5):][0::numLegndr]

    groups[g].sigA = sigmaA



#print(groups[0].sigA)
#print(groups[1].sigA)




#print("\n")





##############################################################################
# Split nu-bar
##############################################################################


nu_groupSplitting = splitTheseGroups(nu_lines)


reactionHeading = nu_groupSplitting[0].pop(0)
[za,awr,numLegndr,nSig0,breakupFlag,nGroups,MF_MT,lineNum] = reactionHeading
MF, MT = str(MF_MT)[0], str(MF_MT)[1:]
nSig0,numLegndr,nGroups,lineNum = int(nSig0),int(numLegndr),int(nGroups),int(lineNum)
#print("ZA",za,"AWR",awr,"#Legndr",numLegndr,"#Sig0",nSig0,"#Groups",nGroups,"MF",MF,"MT",MT)
#print()

##############################################################################

for g in range(nGroups):
    [temp,zero,numSecPos,IG2LO,nWordsList,groupIndex,MF_MT,lineNum] = nu_groupSplitting[g][0]
    MF, MT = str(MF_MT)[0], str(MF_MT)[1:]
    numSecPos,IG2LO,nWordsList,groupIndex,lineNum = int(numSecPos),int(IG2LO),int(nWordsList),int(groupIndex),int(lineNum) 
    #print("TEMP",temp,"#Vals",nWordsList,"Group#",groupIndex,"MF",MF,"MT",MT,"line#",lineNum)

    assert(g+1 == groupIndex and int(zero) == 0)

    fluxSigma = []
    for line in nu_groupSplitting[g][1:]:
        fluxSigma += line[:-2]
    assert(fluxSigma.pop() == None)

    flux  = fluxSigma[:int(len(fluxSigma)*0.5)]
    nuBar = fluxSigma[int(len(fluxSigma)*0.5):][1]
    assert(len(fluxSigma)==3)

    groups[g].nuBar = nuBar



print(groups[0].nuBar)
print(groups[1].nuBar)





for group in groups:
    print(group.sigT[0])






print("\n")

















