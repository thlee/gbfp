#include "Python.h"
#include "gbfp.h"

static PyObject *ErrorObject;

PyObject* MakeLocationList(struct tLocation *ptLocation, unsigned int iLocationNum) {
    unsigned int i;
    struct tLocation *ptLocData;
    PyObject *LocationList;
    PyObject *LocationTuple;

    LocationList = PyList_New(0);

    for (i = 0; i < iLocationNum; i++) {
        ptLocData = ptLocation + i;
        
        LocationTuple = PyTuple_New(2);
        PyTuple_SetItem(LocationTuple, 0, PyInt_FromLong((long) ptLocData->lStart));
        PyTuple_SetItem(LocationTuple, 1, PyInt_FromLong((long) ptLocData->lEnd));
        PyList_Append(LocationList, LocationTuple);
    }

    return LocationList;
}

PyObject* MakeQualifierList(struct tQualifier *ptQualifier, unsigned int iQualifierNum) {
    unsigned int i;
    struct tQualifier *ptQualData;
    PyObject *QualifierList;
    PyObject *QualifierTuple;

    QualifierList = PyList_New(0);

    for (i = 0; i < iQualifierNum; i++) {
        ptQualData = ptQualifier + i;
        
        QualifierTuple = PyTuple_New(2);
        PyTuple_SetItem(QualifierTuple, 0, PyString_FromString(ptQualData->psQualifier));
        PyTuple_SetItem(QualifierTuple, 1, PyString_FromString(ptQualData->psValue));
        PyList_Append(QualifierList, QualifierTuple);
    }

    return QualifierList;
}

PyObject* MakeFeatureDict(struct tFeature *ptFeature) {
    PyObject *FeatureDict;

    FeatureDict =  PyDict_New();

    PyDict_SetItemString(FeatureDict, "feature", PyString_FromString((char *) &(ptFeature->sFeature)));
    PyDict_SetItemString(FeatureDict, "direction", PyString_FromStringAndSize((char *) &(ptFeature->cDirection), 1));
    PyDict_SetItemString(FeatureDict, "start", PyInt_FromLong((long) ptFeature->lStart));
    PyDict_SetItemString(FeatureDict, "end", PyInt_FromLong((long) ptFeature->lEnd));
    PyDict_SetItemString(FeatureDict, "number", PyInt_FromLong((long) ptFeature->iNumber));
    PyDict_SetItemString(FeatureDict, "location_num", PyInt_FromLong((long) ptFeature->iLocationNum));
    PyDict_SetItemString(FeatureDict, "qualifier_num", PyInt_FromLong((long) ptFeature->iQualifierNum));
    PyDict_SetItemString(FeatureDict, "location", MakeLocationList(ptFeature->ptLocation, ptFeature->iLocationNum));
    PyDict_SetItemString(FeatureDict, "qualifier", MakeQualifierList(ptFeature->ptQualifier, ptFeature->iQualifierNum));

   return FeatureDict; 
}

PyObject* MakeGBFFDataDict(struct tGBFFData *ptGBFFData) {
    int i;

    PyObject *GBFFDataDict;
    PyObject *FeatureList;

    GBFFDataDict =  PyDict_New();

    PyDict_SetItemString(GBFFDataDict, "locus_name", PyString_FromString((char *) &(ptGBFFData->sLocusName)));
    PyDict_SetItemString(GBFFDataDict, "length", PyInt_FromLong((long) ptGBFFData->lLength));
    PyDict_SetItemString(GBFFDataDict, "type", PyString_FromString((char *) &(ptGBFFData->sType)));
    PyDict_SetItemString(GBFFDataDict, "topology", PyString_FromString((char *) &(ptGBFFData->sTopology)));
    PyDict_SetItemString(GBFFDataDict, "division_code", PyString_FromString((char *) &(ptGBFFData->sDivisionCode)));
    PyDict_SetItemString(GBFFDataDict, "date", PyString_FromString((char *) &(ptGBFFData->sDate)));
    PyDict_SetItemString(GBFFDataDict, "feature_num", PyInt_FromLong((long) ptGBFFData->iFeatureNumber));
    PyDict_SetItemString(GBFFDataDict, "sequence", PyString_FromString(ptGBFFData->psSequence));

    FeatureList = PyList_New(0);

    for(i = 0; i < ptGBFFData->iFeatureNumber; i++) {
        PyList_Append(FeatureList, MakeFeatureDict((ptGBFFData->ptFeatures) + i));
    }

    PyDict_SetItemString(GBFFDataDict, "features", FeatureList);

    return GBFFDataDict;
}

static PyObject* parse(PyObject *self, PyObject *args) {
    int i;
    char *psFileName;

    struct tGBFFData **pptGBFFData;
    struct tGBFFData *ptGBFFData;
    struct tFeature *ptFeature;

    PyObject *GBFFDataList;

    /* Parsing arguments */
    if (!PyArg_ParseTuple(args, "s", &psFileName)) return NULL;
    
    /* Parsing with C function */
    pptGBFFData = GBFF_Parser(psFileName);

    if (pptGBFFData != NULL) {
        /* Convert datas from C to Python */
        GBFFDataList = PyList_New(0);

        for (i = 0; *(pptGBFFData + i) != NULL; i++) \
            PyList_Append(GBFFDataList, MakeGBFFDataDict(*(pptGBFFData + i)));

        GBFF_Free(pptGBFFData);
        
        return GBFFDataList;
    } else {
        PyErr_SetString(PyExc_IOError, "File not found !!!");

        return NULL;
    }
}         
 
static struct PyMethodDef gbfp_methods[] = { 
    {"parse", parse, METH_VARARGS},
    {NULL, NULL}
};

void initgbfp() {
    PyObject *parser; 
    parser = Py_InitModule("gbfp", gbfp_methods);
    ErrorObject = Py_BuildValue("s", "GBFF parser error !!!");
} 
