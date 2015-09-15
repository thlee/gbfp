  * Source code of a example program, seqext.
```
/* filename: seqext.c */

#define __EXTENSIONS__

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <regex.h>
#include <unistd.h>
#include <string.h>
#include <gbfp.h>
 
void help(void) {
        printf("Extract specific feature from a GenBank flat file.\n"
        "\n"
        "Usage: seqext [-h] [-i Genbank_file] [-f Feature to extract (default: gene)] [-q Qualifier to print | -s Word to search]\n"
        "\n"
        "-h  print help and exit\n");
}

static char *getQualValue(char *sQualifier, gb_feature *ptFeature) {
    static char *null = "";
    gb_qualifier *i;

    for (i = ptFeature->ptQualifier; (i - ptFeature->ptQualifier) < ptFeature->iQualifierNum; i++) 
        if (strcmp(sQualifier, i->sQualifier) == 0)
            return i->sValue;
    return null;
}

static char *searchQualValue(char *sWord, gb_feature *ptFeature) {
    gb_qualifier *i;

    for (i = ptFeature->ptQualifier; (i - ptFeature->ptQualifier) < ptFeature->iQualifierNum; i++) 
        if (strstr(sWord, i->sValue) != NULL)
            return i->sValue;
    return NULL;
}

int main(int argc, char *argv[]) {
    char *sFileName = NULL;
    char *sFeature = "gene";
    char *sQualifier = NULL;
    char *sSequence = NULL;
    char *sWord = NULL;
    char *sValue = NULL;
    
    gb_data **pptSeqData, *ptSeqData;
    gb_feature *ptFeature;

    unsigned int iOpt;
    unsigned int i, j;
    
    while((iOpt = getopt(argc, argv, "f:i:q:s:h")) != -1) {
        switch(iOpt) {
        case 'h':
            help();
            exit(0);
            break;
        case 'i':
            sFileName = optarg;
            break;
        case 'f':
            sFeature = optarg;
            break;
        case 'q':
            sQualifier = optarg;
            break;
        case 's':
            sWord = optarg;
            break;
        default:
            help();
            exit(0);
        }
    }
 
    pptSeqData = parseGBFF(sFileName); /* parse a GBF file which contains more than one GBF sequence data */
    for (i = 0; (ptSeqData = *(pptSeqData + i)) != NULL; i++) { /* ptSeqData points a parsed data of a GBF sequence data */
        /* start of user process */
        for (j = 0; j < ptSeqData->iFeatureNum; j++) {
            ptFeature = (ptSeqData->ptFeatures + j);
            if (strcmp(sFeature, ptFeature->sFeature) == 0) {
                if (sWord != NULL) {
                    if ((sValue = searchQualValue(sWord, ptFeature)) != NULL) {
                        printf(">%s_%li-%li_%c %s\n%s\n", \
                            ptFeature->sFeature, \
                            ptFeature->lStart, \
                            ptFeature->lEnd, \
                            ptFeature->cDirection, \
                            sValue == NULL ? "" : sValue, \
                            sSequence = getSequence(ptSeqData->sSequence, ptFeature));
                        free(sSequence);
                    }
                } else {
                    printf(">%s_%li-%li_%c %s\n%s\n", \
                        ptFeature->sFeature, \
                        ptFeature->lStart, \
                        ptFeature->lEnd, \
                        ptFeature->cDirection, \
                        sQualifier == NULL ? "" : getQualValue(sQualifier, ptFeature), \
                        sSequence = getSequence(ptSeqData->sSequence, ptFeature));
                    free(sSequence);
                }
            }
        }
        /* end of user process */
    }
    freeGBData(pptSeqData); /* release memory space */

    return 0;
}

```

  * Makefile to compile 'seqext.c'
```
CC = gcc
CFLAGS = -O2 -Wall -ansi -pedantic
DEBUG = -ggdb -pg

prefix = /usr/local
bindir = $(prefix)/bin

all: seqext

install: seqext
        if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
        cp -f seqext $(bindir)
        chmod 755 $(bindir)/seqext

uninstall:
        rm -f $(bindir)/seqext

seqext: seqext.c
        $(CC) $(CFLAGS) -D_GNU_SOURCE -o seqext seqext.c -lgbfp

debug: seqext.c
        $(CC) $(CFLAGS) $(DEBUG) -D_GNU_SOURCE -o seqext seqext.c -lgbfp

clean:
        rm -f seqext
```