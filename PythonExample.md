  * Source code of a example program, seqext.py.
```
#!/usr/bin/python
#filename: seqext.py

import sys
import getopt
import gbparsy

def help(iCode=0):
    print """Extract specific feature from a GenBank flat file.

Usage: seqext [-h] [-i Genbank_file] [-f Feature to extract (default: gene)] [-q Qualifier to print | -s Word to search]

-h  print help and exit"""
    sys.exit(iCode)

sNorBase = "ACGTRYMKWSBDHVNacgtrymkwsbdhvn";
sComBase = "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn";

def getRevCom(sSequence):
    lSequence = map(lambda x: sComBase[sNorBase.find(x)], sSequence)
    lSequence.reverse()

    return "".join(lSequence)

def getQualValue(sQualifier, dFeature):
    if not sQualifier: return ""

    for tQualifier in dFeature["qualifiers"]:
        if tQualifier[0] == sQualifier: return tQualifier[1]

    return ""

def searchQualValue(sWord, dFeature):
    if not sWord: return ""

    for tQualifier in dFeature["qualifiers"]:
        if sWord in tQualifier[1]: return tQualifier[1]

    return ""

def getSequence(sWholeSequence, dFeature):
    sSequence = sWholeSequence[dFeature["start"] - 1: dFeature["end"]]
    if dFeature["direction"] == "C":
        sSequence = getRevCom(sSequence)

    return sSequence
    
try:
    lOpts, lArgs = getopt.getopt(sys.argv[1:], "f:i:q:s:h")
except:
    help()

sFileName = None
sFeature = "gene"
sQualifier = None
sWord = None
    
for sOpt, sVal in lOpts:
    if sOpt == "-h":
        help(1)
    elif sOpt == "-i":
        sFileName = sVal
    elif sOpt == "-f":
        sFeature = sVal
    elif sOpt == "-q":
        sQualifier = sVal
    elif sOpt == "-s":
        sWord = sVal
    else:
        help()

oGBParsy = gbparsy.gbparsy()
oGBParsy.parse(sFileName) # parse a GBF file which contains more than one GBF sequence data

# If you want to get a list of SeqRecord instances, use getBioData method instead of getRawData.
lSeqData = oGBParsy.getRawData()
   
for dSeqData in lSeqData: # dSeqData has a parsed data of a GBF sequence data
    # start of user process
    for dFeature in dSeqData["features"]:
        if dFeature["feature"] == sFeature:
            if sWord:
                sValue = searchQualValue(sWord, dFeature)
                if sValue:
                    sys.stdout.write(">%s_%i-%i_%s %s\n%s\n" % (
                        dFeature["feature"],
                        dFeature["start"],
                        dFeature["end"],
                        dFeature["direction"],
                        sValue,
                        getSequence(dSeqData["sequence"], dFeature)))
            else:
                sys.stdout.write(">%s_%i-%i_%s %s\n%s\n" % (
                    dFeature["feature"],
                    dFeature["start"],
                    dFeature["end"],
                    dFeature["direction"],
                    getQualValue(sQualifier, dFeature),
                    getSequence(dSeqData["sequence"], dFeature)))
    # end of user process
```