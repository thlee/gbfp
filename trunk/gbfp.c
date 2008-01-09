#define __EXTENSIONS__

#include <stdio.h>
#include <stdlib.h>
#include <regex.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <sys/types.h>
 
#include "gbfp.h"
 
const char sVer[] = "0.3.0";
const char sNorBase[] = "ACGTRYMKWSBDHVNacgtrymkwsbdhvn";
const char sComBase[] = "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn";
const unsigned int iBaseLen = 30;

static void rtrim(string sLine) {
    register int i;
    
    for (i = (strlen(sLine) - 1); i >= 0; i--) {
        if (! isspace(sLine[i])) {
            sLine[++i] = '\0';
            break;
        }
    }
}

static void rtrim2(string sLine, char cRemove) {
    register int i;
    
    for (i = (strlen(sLine) - 1); i >= 0; i--) {
        if (sLine[i] == cRemove) {
            sLine[i] = '\0';
            break;
        }
    }
}

static int Pos2Num(string sPositions, int aiPositions[]) {
    register int i;
    int iNum = 0;

    for (i = strlen(sPositions); i >= 0; i--) {
        if (isdigit(*(sPositions + i))) aiPositions[(aiPositions[iNum] - 1 == i) ? iNum : ++iNum] = i;
        else *(sPositions + i) = '\0';
    }

    return iNum;
}

static int Positions2Numbers(string sPositions, unsigned long *lStart, unsigned long *lEnd) {
    int aiPositions[16] = {-2,};
    int iNum;

    iNum = Pos2Num(sPositions, aiPositions);
   
    if (iNum == 2) {
        *lStart = atol(sPositions + aiPositions[2]);
        *lEnd = atol(sPositions + aiPositions[1]);

        return 1;
    } else if (iNum == 1) {
        *lStart = *lEnd = atol(sPositions + aiPositions[1]);
       
        return 1;
    } else {
        fprintf(stderr, "Warning: cannot parse '%s'\n", sPositions);

        return 0;
    }

}

static void LocusParser(string sLocusStr, gbdata *ptGBData) {
    /*    
    01-05      'LOCUS'
    06-12      spaces
    13-28      Locus name
    29-29      space
    30-40      Length of sequence, right-justified
    41-41      space
    42-43      bp
    44-44      space
    45-47      spaces, ss- (single-stranded), ds- (double-stranded), or
               ms- (mixed-stranded)
    48-53      NA, DNA, RNA, tRNA (transfer RNA), rRNA (ribosomal RNA),
               mRNA (messenger RNA), uRNA (small nuclear RNA), snRNA,
               snoRNA. Left justified.
    54-55      space
    56-63      'linear' followed by two spaces, or 'circular'
    64-64      space
    65-67      The division code (see Section 3.3)
    68-68      space
    69-79      Date, in the form dd-MMM-yyyy (e.g., 15-MAR-1991)
    */
    
    regex_t ptRegExLocus;
    regmatch_t ptRegMatch[7];
    
    const char sLocus[] = "LOCUS +([a-z|A-Z|0-9|_]+) +([0-9]+) bp +([a-z|A-Z|-]+) +([a-z]+) +([A-Z]{3}) (.+)";
    char sTemp[LINELEN];
    unsigned int i, iLen;
    
    
    struct tData {
        char cType;
        void *Pointer;
    };
    
    struct tData tDatas[] = {
        {STRING, NULL},
        {LONG, NULL},
        {STRING, NULL},
        {STRING, NULL},
        {STRING, NULL},
        {STRING, NULL}};

    tDatas[0].Pointer = ptGBData->sLocusName;
    tDatas[1].Pointer = &(ptGBData->lLength);
    tDatas[2].Pointer = ptGBData->sType;
    tDatas[3].Pointer = ptGBData->sTopology;
    tDatas[4].Pointer = ptGBData->sDivisionCode;
    tDatas[5].Pointer = ptGBData->sDate;

    i = regcomp(&ptRegExLocus, sLocus, REG_EXTENDED | REG_ICASE);
    if (i != 0) {
        regerror(i, &ptRegExLocus, sTemp, LINELEN);
        fprintf(stderr, "%s\n", sTemp);
        exit(1);
    }
    
    rtrim(sLocusStr);
        
    if (regexec(&ptRegExLocus, sLocusStr, 7, ptRegMatch, 0) == 0) {
        for (i = 0; i < 6; i++) {
            iLen = ptRegMatch[i + 1].rm_eo - ptRegMatch[i + 1].rm_so;
            switch (tDatas[i].cType) {
            case STRING:
                memcpy(tDatas[i].Pointer, (sLocusStr + ptRegMatch[i + 1].rm_so), iLen);
                *((string ) tDatas[i].Pointer + iLen) = '\0';
                /*
                printf("%i, %s\n", iLen, (string ) tDatas[i].Pointer);
                */
                break;
            case LONG:
                memcpy(sTemp, (sLocusStr + ptRegMatch[i + 1].rm_so), iLen);
                sTemp[iLen] = '\0';
                *((unsigned long *) tDatas[i].Pointer) = atol(sTemp);
                break;
            default:
                perror("Unknown Data Type!!!");
            }
        }
    } else {
        regerror(i, &ptRegExLocus, sTemp, LINELEN);
        fprintf(stderr, "%s\n", sTemp);
        exit(1);
    }
}

static string checkComplement(string sLocation) {
    string sPosition;

    for (;isspace(*sLocation);sLocation++);

    for (sPosition = sLocation; *sPosition; sPosition++) {
        /* Check the 1st and the 2nd characters of 'complement' */
        if (*sPosition == 'c' && *(sPosition + 1) == 'o') {
            rtrim2(sLocation, ')');
            return sPosition + 11;
        }
    }

    return sLocation;
}

static string checkJoin(string sLocation) {
    string sPosition;

    for (;isspace(*sLocation);sLocation++);

    for (sPosition = sLocation; *sPosition; sPosition++) {
        /* Check the 1st and the 2nd characters of 'complement' */
        if (*sPosition == 'j' && *(sPosition + 1) == 'o') {
            rtrim2(sLocation, ')');
            return sPosition + 5;
        }
    }

    return sLocation;
}

/* Parsing a string which contains location information */

static void LocationParser(string sLocation, feature *pFeature) {
    string sTemp;
    string sString = NULL;
    
    unsigned int iLocationNum = 1;
    
    /* Evalue sequence direction
    sString has location and join informations
    */

    sString = checkComplement(sLocation);

    if (sLocation == sString) pFeature->cDirection = NORMAL;
    else pFeature->cDirection = REVCOM;

    /* Remove 'join' string
    sString has location informations
    */

    sString = checkJoin(sString);

    sTemp = sString - 1;
    while((sTemp = strchr((sTemp + 1), ','))) iLocationNum++;
    pFeature->ptLocation = malloc(iLocationNum * sizeof(*(pFeature->ptLocation)));

    iLocationNum = 0;
    sLocation = strtok_r(sString, ",", &sTemp);
    if (Positions2Numbers(sLocation,
        &(((pFeature->ptLocation)+iLocationNum)->lStart),
        &(((pFeature->ptLocation)+iLocationNum)->lEnd)) == 1) iLocationNum++;

    while((sLocation = strtok_r(NULL, ",", &sTemp))) {
        if (Positions2Numbers(sLocation,
            &(((pFeature->ptLocation)+iLocationNum)->lStart),
            &(((pFeature->ptLocation)+iLocationNum)->lEnd))) iLocationNum++;
    }
    
    pFeature->lStart = (pFeature->ptLocation)->lStart;
    pFeature->lEnd = ((pFeature->ptLocation)+(iLocationNum - 1))->lEnd;
    pFeature->iLocationNum = iLocationNum;
}

static string parseQualifier(string sQualifier, string *psValue) {
    string sPosition;

    for (;isspace(*sQualifier); sQualifier++);

    if ((sPosition = strchr(sQualifier, '=')) == NULL) {
        *psValue = sQualifier + strlen(sQualifier);
        return sQualifier;
    }

    *sPosition = '\0';

    sPosition++;

    for (;isspace(*sQualifier); sQualifier++);

    if (*sPosition == '"') {
        sPosition++;
        rtrim2(sPosition, '"');
    }

    *psValue = sPosition;

    return sQualifier;
}

static void QualifierParser(string sQualifier, feature *pFeature) {
    string sValue;
    string sTemp = NULL;
    string sString = NULL;
    qualifier *ptQualifier;
    
    pFeature->ptQualifier = malloc(INITQUALIFIERNUM * sizeof(qualifier));
    ptQualifier = pFeature->ptQualifier;

    /* Parse the 1st qualifier string */
    sString = strtok_r(sQualifier, "\n", &sTemp);

    sQualifier = parseQualifier(sString, &sValue);

    ptQualifier->sQualifier = sQualifier;
    ptQualifier->sValue = sValue;
    ptQualifier++;
    
    /* Parse the rest qualifier string */
    while((sString = strtok_r(NULL, "\n", &sTemp)) != NULL) {
        sQualifier = parseQualifier(sString, &sValue);
        ptQualifier->sQualifier = sQualifier;
        ptQualifier->sValue = sValue;
        ptQualifier++;
    }

    pFeature->iQualifierNum = ptQualifier - pFeature->ptQualifier;
    pFeature->ptQualifier =  realloc(pFeature->ptQualifier, pFeature->iQualifierNum * sizeof(qualifier));
}

static unsigned int SequenceParser(string sSequence, string sSequence2) {
    register unsigned int i = 0;
    register unsigned int j = 0;
    register char c;
    
    while((c = *(sSequence + i++)) != '\0')
        if (isalpha(c) != 0) *(sSequence2 + j++) = c;

    *(sSequence2 + j) = '\0';
        
    return j;
}

static void RevCom(string sSequence) {
    char c;
    unsigned int k;
    unsigned long i, j;
    
    for (i = 0, j = strlen(sSequence) - 1; i < j; i++, j--) {
	c = *(sSequence + i);
	*(sSequence + i) = 'X';
	for (k = 0; k < iBaseLen; k++)
	    if (*(sNorBase + k) == *(sSequence + j)) {
		*(sSequence + i) = *(sComBase + k);
		break;
	    }
	*(sSequence + j) = 'X';
	for (k = 0; k < iBaseLen; k++)
	    if (*(sNorBase + k) == c) {
		*(sSequence + j) = *(sComBase + k);
		break;
	}
    }
}

string getSequence(string sSequence, feature *ptFeature) {
    unsigned long lSeqLen = 1; /* For the '\0' characher */
    unsigned long lStart, lEnd;
    unsigned int i;
    string sSequenceTemp;
    
    for (i = 0; i < ptFeature->iLocationNum; i++)
        lSeqLen += (((ptFeature->ptLocation) + i)->lEnd - ((ptFeature->ptLocation) + i)->lStart + 1);
    
    sSequenceTemp = malloc(lSeqLen * sizeof(char));
    
    lSeqLen = 0;
    
    for (i = 0; i < ptFeature->iLocationNum; i++) {
        lStart = ((ptFeature->ptLocation) + i)->lStart;
        lEnd = ((ptFeature->ptLocation) + i)->lEnd;
        strncpy(sSequenceTemp + lSeqLen, sSequence + lStart - 1, lEnd - lStart + 1);
        lSeqLen += (lEnd - lStart + 1);
    }
    
    *(sSequenceTemp + lSeqLen) = '\0';
    
    if (ptFeature->cDirection == REVCOM) RevCom(sSequenceTemp);

    return sSequenceTemp;
}

void freeGBData(gbdata **pptGBData) {
    gbdata *ptGBData = NULL;
    feature *tpFeatures = NULL;
    unsigned int iFeatureNumber = 0;
    unsigned int iFeatureNum = 0;
    unsigned int iSeqPos = 0;
  
    for (iSeqPos = 0; *(pptGBData + iSeqPos) != NULL; iSeqPos++) {
        ptGBData = *(pptGBData + iSeqPos);

        tpFeatures = ptGBData->ptFeatures;
        iFeatureNumber = ptGBData->iFeatureNumber;
    
        for (iFeatureNum = 0; iFeatureNum < iFeatureNumber; iFeatureNum++) {
            /*
            printf("%i, %i\n", iFeatureNum, (tpFeatures+iFeatureNum)->iQualifierNum);
            */
            free((tpFeatures+iFeatureNum)->ptLocation);
            free(((tpFeatures+iFeatureNum)->ptQualifier)->sQualifier);
        }
    }

    free(tpFeatures);
}

static gbdata *_GBFF_Parser(FILE *FSeqFile) {
    char sLine[LINELEN] = {'\0',};
    char sLocation[LINELEN] = {'\0',};
    string sQualifier = NULL;
    string sQualifierTemp = NULL;
   
    unsigned int iReadPos = INELSE;
    unsigned int iFeatureNum = 0;
    unsigned int iFeatureMem = INITFEATURENUM;
    unsigned int i = 0;
    unsigned long lSeqLen;
    
    feature *pFeatures = NULL;
    feature *pFeature = NULL;
    gbdata *ptGBData = NULL;
    
    pFeatures = (feature *) malloc(iFeatureMem * sizeof(feature));

    /* Confirming GBFF File with LOCUS line */
    while(fgets(sLine, LINELEN, FSeqFile)) {
        if (memcmp(sLine, "LOCUS", 5) == 0) {
            ptGBData = malloc(sizeof(gbdata));
            break;
        }
    }

    /* If there is a no LOCUS line, next statement return NULL value to end parsing */
    if (ptGBData == NULL) return NULL;
   
    /* Parse LOCUS line */ 
    LocusParser(sLine, ptGBData);

    while(fgets(sLine, LINELEN, FSeqFile))
        if (memcmp(sLine, "FEATURES", 8) == 0) break;

    /* Parse FEATURES */
    while(fgets(sLine, LINELEN, FSeqFile)) {
        rtrim(sLine);

        if (memcmp(sLine, "BASE COUNT", 10) == 0 || memcmp(sLine, "ORIGIN", 6) == 0) {
            if (iFeatureNum == iFeatureMem) {
                iFeatureMem += INITFEATURENUM;
                pFeatures = realloc(pFeatures, sizeof(feature) * iFeatureMem);
            }

            if (strlen(sLocation) != 0) LocationParser(sLocation, (pFeatures + iFeatureNum - 1));
            if (sQualifier < sQualifierTemp) {
                *sQualifierTemp++ = '\n';
                *sQualifierTemp = '\0';
                sQualifier = realloc(sQualifier, (sQualifierTemp - sQualifier + 1) * sizeof(*sQualifier));
                QualifierParser(sQualifier, (pFeatures + iFeatureNum - 1));
            }
            
            break;
        }
        
        if (memcmp(sLine + 5, "               ", 15) != 0) {
            if (iFeatureNum == iFeatureMem) {
                iFeatureMem += INITFEATURENUM;
                pFeatures = realloc(pFeatures, sizeof(feature) * iFeatureMem);
            }

            if (strlen(sLocation) != 0) LocationParser(sLocation, (pFeatures + iFeatureNum - 1));
            if (sQualifier < sQualifierTemp) {
                *sQualifierTemp++ = '\n';
                *sQualifierTemp = '\0';
                /*
                printf("=====\n%s=====", sQualifier);
                */
                sQualifier = realloc(sQualifier, (sQualifierTemp - sQualifier + 1) * sizeof(*sQualifier));
                QualifierParser(sQualifier, (pFeatures + iFeatureNum - 1));
            }

            *sLocation = '\0';
            sQualifier = malloc(sizeof(*sQualifier) * LINELEN);
            sQualifierTemp = sQualifier;
            
            iReadPos = INFEATURE;
            
            memcpy((pFeatures + iFeatureNum)->sFeature, (sLine + 5), 15);
            *(((pFeatures + iFeatureNum)->sFeature) + 15) = '\0';
            rtrim((pFeatures + iFeatureNum)->sFeature);
            strcpy(sLocation, (sLine + 21));

            /* Feature Initalize */
            pFeature = pFeatures + iFeatureNum;
            pFeature->iNumber = iFeatureNum;
            pFeature->cDirection = NORMAL;
            pFeature->iLocationNum = 0;
            pFeature->lStart = 0;
            pFeature->lEnd = 0;
            pFeature->iQualifierNum = 0;
            pFeature->ptLocation = NULL;
            pFeature->ptQualifier = NULL;

            iFeatureNum++;
        } else if (*(sLine + QUALIFIERSTART) == '/') {
            iReadPos = INQUALIFIER;
            if (sQualifier < sQualifierTemp) *sQualifierTemp++ = '\n';
            i = strlen(sLine) - (QUALIFIERSTART + 1);
            memcpy(sQualifierTemp, sLine + (QUALIFIERSTART + 1), i);
            sQualifierTemp += i;
        } else {
            if (iReadPos == INFEATURE) {
                strcpy((sLocation + strlen(sLocation)), (sLine + QUALIFIERSTART));
            } else if (iReadPos == INQUALIFIER) {
                i = strlen(sLine) - QUALIFIERSTART;
                memcpy(sQualifierTemp, sLine + QUALIFIERSTART, i);
                sQualifierTemp += i;
            }
        }
    }
    
    if (memcmp(sLine, "ORIGIN", 6) != 0)
        while(fgets(sLine, LINELEN, FSeqFile))
            if (memcmp(sLine, "ORIGIN", 6) == 0) break;

    lSeqLen = 0;
    ptGBData->sSequence = malloc((ptGBData->lLength + 1) * sizeof(char));
    
    while(fgets(sLine, LINELEN, FSeqFile)) {
        if (memcmp(sLine, "//", 2) == 0) break;

        lSeqLen += SequenceParser(sLine, ptGBData->sSequence + lSeqLen);
    }

    ptGBData->iFeatureNumber = iFeatureNum;
    ptGBData->ptFeatures = pFeatures;
    
    return ptGBData;
}

gbdata **parseGBFF(string spFileName) {
    int iGBFSeqPos = 0;
    unsigned int iGBFSeqNum = INITGBFSEQNUM;
    gbdata **pptGBDatas;
    FILE *FSeqFile;
 
    if (spFileName == NULL) {
        FSeqFile = stdin;
    } else {
        if (access(spFileName, F_OK) != 0) {
            /* perror(spFileName); */
            return NULL;
        } else {
            FSeqFile = fopen(spFileName, "r");
        }
    }
    
    pptGBDatas = malloc(iGBFSeqNum * sizeof(gbdata *));

    do {
        if (iGBFSeqNum == iGBFSeqPos) {
            iGBFSeqNum += INITGBFSEQNUM;
            pptGBDatas = realloc(pptGBDatas, iGBFSeqNum * sizeof(gbdata *));
        }
        *(pptGBDatas + iGBFSeqPos) = _GBFF_Parser(FSeqFile);
    } while (*(pptGBDatas + iGBFSeqPos++) != NULL);
    
    if (spFileName) fclose(FSeqFile);

    return pptGBDatas;
}
