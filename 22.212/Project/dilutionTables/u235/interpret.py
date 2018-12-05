def convertToExpPos(x):
    split = x.split('+')
    return float(split[0])*10.0**float(split[1])

def convertToExpNeg(x):
    split = x.split('-')
    return float(split[0])*10.0**(-1*float(split[1]))


def number(x):
    if '+' in x: return convertToExpPos(x)
    if '-' in x[1:]: return convertToExpNeg(x)
    try: return float(x)
    except: return x



header = []
array = []
with open("tape26", "r") as f:
    for line in f:
        split = line.split()
        if split[-2] == '1451':
            split[-3] = split[-3][:-4]
            if split[-3] == '':
                del split[-3]
            split = split[:-2]+[split[-1]]
            header.append([number(x) for x in split])
        elif len(split) > 4 and split[-4][-4:] == '1395':
            split[-4] = split[-4][:-4]
            if split[-4] == '':
                del split[-3]
            array.append([number(x) for x in split])
        elif len(split) > 4 and split[-3][-4:] == '1395':
            split[-3] = split[-3][:-4]
            if split[-3] == '':
                del split[-3]
            array.append([number(x) for x in split])

a0 = header[0]
za = a0[0]
awr = a0[1]
assert(a0[2]==0)
nSig0 = int(a0[3])
assert(a0[4]==-1)
nWordsTitle = a0[5]
print(a0)
print(za,awr,nSig0,nWordsTitle)
print()






a1 = header[1]
temp = a1[0]
assert(a1[1] == 0)
nGroups = int(a1[2])
nPhoton = int(a1[3])
nWordsList = a1[4]
assert(a1[5] == 0)
print(a1)
print(temp,nGroups,nPhoton,nWordsList)
print()

restOfTitle = []
for entry in header[2:]:
    restOfTitle += entry[:-1]
title = restOfTitle[0]
Sig0 = restOfTitle[1:nSig0+1]
Ebounds = restOfTitle[nSig0+1:nSig0+nGroups+2]
print(restOfTitle)
print(Sig0)
print(Ebounds)




print()


