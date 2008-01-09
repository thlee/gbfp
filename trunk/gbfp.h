#define LINELEN             65536
#define MEGA                1048576
#define INITGBFSEQNUM       5
#define INITFEATURENUM      64
#define INITQUALIFIERNUM    128
#define FEATURELEN          16
#define QUALIFIERLEN        16
#define LOCUSLEN            16
#define TYPELEN             9
#define TOPOLOGYSTRLEN      8
#define DIVISIONCODELEN     3
#define DATESTRLEN          11
#define QUALIFIERSTART      21

#define INELSE              0
#define INFEATURE           1
#define INQUALIFIER         2

#define NORMAL              'N'
#define REVCOM              'C'
#define LINEAR              'L'
#define CIRCULAR            'C'

#define CHARACTER           'C'
#define LONG                'L'
#define STRING              'S'

typedef char *string;

typedef struct tLocation {
    unsigned long lStart;
    unsigned long lEnd;
} location;

typedef struct tQualifier {
    string sQualifier;
    string sValue;
} qualifier;

typedef struct tFeature {
    char sFeature[FEATURELEN + 1];
    char cDirection;
    unsigned long lStart;
    unsigned long lEnd;
    unsigned int iNumber;
    unsigned int iLocationNum;
    unsigned int iQualifierNum;
    location *ptLocation;
    qualifier *ptQualifier;
} feature;

typedef struct tVersion {
    string sVersion;
    string sGI;
} version;

typedef struct tGBFFData {
    char sLocusName[LOCUSLEN + 1];
    unsigned long lLength;
    char sType[TYPELEN + 1];
    char sTopology[TOPOLOGYSTRLEN + 1];
    char sDivisionCode[DIVISIONCODELEN + 1];
    char sDate[DATESTRLEN + 1];
    feature *ptFeatures;
    unsigned int iFeatureNumber;
    string sDef;
    string sAccession;
    version tVersion;
    string sKeywords;
    string sSequence;
} gbdata;

gbdata **parseGBFF(string spFileName);
void freeGBData(gbdata **pptGBFFData);
string getSequence(string sSequence, feature *ptFeature);
