#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <regex.h>
#include <unistd.h>
#include <string.h>
#include <zlib.h>
#include <getopt.h>
#include "gbfp.h"
 
void help(void) {
        printf("Convert GenBank sequence into FASTA sequence.\n\nUsage: gbk2fna [-h] [-V] [-f Feature_to_extract] [-q Qualifier_to_print] [-i Genbank_file]\n\n-V  print version information and exit\n");
}

int main(int argc, char *argv[]) {
    char *psFileName = NULL;
    char *psFeature = "gene";
    char *psQualifier = NULL;
    
    struct tGBFFData **pptSeqData;
    struct tGBFFData *ptSeqData;
    struct tFeature *ptFeature;
    struct tQualifier *ptQualifier;

    unsigned int iOpt;
    unsigned int i, j, k;
    
    while((iOpt = getopt(argc, argv, "f:i:q:h")) != -1) {
        switch(iOpt) {
        case 'h':
            help();
            exit(0);
            break;
        case 'i':
            psFileName = optarg;
            break;
        case 'f':
            psFeature = optarg;
            break;
        case 'q':
            psQualifier = optarg;
            break;
        default:
            help();
            exit(0);
        }
    }
 
    
    pptSeqData = GBFF_Parser(psFileName);
    if (pptSeqData != NULL) {
        for (i = 0; (ptSeqData = *(pptSeqData + i)) != NULL; i++) {
            for (j = 0; j < ptSeqData->iFeatureNumber; j++) {
                ptFeature = (ptSeqData->ptFeatures + j);
                if (strcmp(psFeature, ptFeature->sFeature) == 0) {
                    if (psQualifier == NULL) {
                        printf(">%s_%li_%li\n", \
                            ptFeature->sFeature, \
                            ptFeature->lStart, \
                            ptFeature->lEnd);
                        printf("%s\n", GBFF_Get_Sequence(ptSeqData->psSequence, ptSeqData->ptFeatures + j));
                    } else {
                        for (k = 0; k < ptFeature->iQualifierNum; k++) {
                            ptQualifier = (ptFeature->ptQualifier + k);
                            if (strcmp(psQualifier, ptQualifier->sQualifier) == 0) {
                                printf(">%s_%li_%li %s\n", \
                                    ptFeature->sFeature, \
                                    ptFeature->lStart, \
                                    ptFeature->lEnd, \
                                    ptQualifier->psValue);
                                printf("%s\n", GBFF_Get_Sequence(ptSeqData->psSequence, ptSeqData->ptFeatures + j));
                                break;
                            }
                        }
                    }
                }
            }
        }
        GBFF_Free(pptSeqData);
    }

    return 0;
}
