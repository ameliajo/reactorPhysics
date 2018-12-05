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
#print(za,awr,nSig0,nWordsTitle)
#print()


[temp,zero1,nGroups,nPhoton,nWordsList,zero2,lineNum] = header[1]
nGroups,nPhoton,nWordsList,lineNum = int(nGroups),int(nPhoton),int(nWordsList),int(lineNum)
assert(zero1 == 0 and zero2 == 0)
#print(header[1])
#print(temp,nGroups,nPhoton,nWordsList)
#print()

restOfTitle = []
for entry in header[2:]:
    restOfTitle += entry[:-1]
title = restOfTitle[0]
Sig0 = restOfTitle[1:nSig0+1]
Ebounds = restOfTitle[nSig0+1:nSig0+nGroups+2]
#print(Sig0)
#print(Ebounds)
#print()




##############################################################################
# Split into various reactions
##############################################################################

sigT_lines = [line for line in array if line[-2] == 1   ]
sigE_lines = [line for line in array if line[-2] == 2   ]
sigF_lines = [line for line in array if line[-2] == 18  ]
sigG_lines = [line for line in array if line[-2] == 102 ]



sigT = {}

##############################################################################
# Split sigT into various dilutions
##############################################################################

sigT_dilutionVecs = []
dilution = []
for line in sigT_lines:
    dilution.append(line)
    if None in line:
        sigT_dilutionVecs.append(dilution)
        dilution = []


[za,awr,numLegndr,nSig0,breakupFlag,nGroups,MF,MT,lineNum] = sigT_dilutionVecs[0][0]
numLegndr,nGroups,lineNum = int(numLegndr),int(nGroups),int(lineNum)
assert(lineNum == 1)
#print(sigT_dilutionVecs[0][0])
#print(za,awr,numLegndr,nSig0,breakupFlag,nGroups,MF,MT,lineNum)
#print()



##############################################################################
# look at sigT dilution #1
##############################################################################

[temp,zero,numSecPos,indexToLowestNonzeroGroup,nWordsList,groupIndex,MF,MT,lineNum] = sigT_dilutionVecs[0][1]
#print(sigT_dilutionVecs[0][1])
#print(temp,zero,numSecPos,indexToLowestNonzeroGroup,nWordsList,groupIndex,MF,MT,lineNum)
#print()

fluxSigma = []
for line in sigT_dilutionVecs[0][2:]:
    fluxSigma += line[:-3]
assert(fluxSigma.pop() == None)
flux  = fluxSigma[:int(len(fluxSigma)*0.5)]
sigma = fluxSigma[int(len(fluxSigma)*0.5):]
print(flux)
print(sigma)








































