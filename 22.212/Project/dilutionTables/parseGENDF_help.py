
class Group:
    def __init__(self,ID,E_low,E_high,Sig0):
        self.ID     = ID
        self.E_low  = E_low
        self.E_high = E_high
        self.dilutionVals = Sig0

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

    assert(lineNum == 1); 
    return numLegndr


def readFromTape26(matNum):
    header = []
    data = []
    with open("tape26", "r") as f:
        for line in f:
            split = line.split()
            if split[-2] == '1451':
                split[-3] = split[-3][:-4] # Get rid of the 1451 ID tag in header
                split = split[:-2]+[split[-1]]
                header.append([number(x) for x in split])

            elif len(split) > 4 and split[-4][-4:] == matNum:
                split[-4] = split[-4][:-4]
                split = [None if val == '' else val for val in split]
                data.append([number(x) for x in split])

            # For instances where MF is more than two characters, and accidentally
            # hits our friendly MT
            elif len(split) > 4 and split[-3][-4:] == matNum:
                split[-3] = split[-3][:-4]
                split = [None if val == '' else val for val in split]
                data.append([number(x) for x in split])
    return header, data




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




