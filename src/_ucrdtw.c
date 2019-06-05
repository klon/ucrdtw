#include <Python.h>
#include <numpy/arrayobject.h>
#include "ucrdtw.h"
#include <math.h>

/* Docstrings */
static char module_docstring[] = "This module implements fast nearest-neighbor retrieval under the dynamic time warping (DTW).";
static char ucrdtw_docstring[] = "Calculate the nearest neighbor of a times series in a larger time series expressed as location and distance, "
                                 "using the UCR suite optimizations.";
static char ucrdtw_all_docstring[] = "Calculate the top n nearest neighbor of a times series in a larger time series expressed as location and distance, "
                                 "using the UCR suite optimizations.";

/* Available functions */
static PyObject* ucrdtw_ucrdtw(PyObject *self, PyObject *args);

static PyObject* ucrdtw_all(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = { { "ucrdtw", ucrdtw_ucrdtw, METH_VARARGS, ucrdtw_docstring }, 
                                        { "ucrdtw_all", ucrdtw_all, METH_VARARGS, ucrdtw_all_docstring },
                                        { NULL, NULL, 0, NULL } };

/* Initialize the module */
#if PY_MAJOR_VERSION >= 3
  #define MOD_ERROR_VAL NULL
  #define MOD_SUCCESS_VAL(val) val
  #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
          ob = PyModule_Create(&moduledef);
#else
  #define MOD_ERROR_VAL
  #define MOD_SUCCESS_VAL(val)
  #define MOD_INIT(name) void init##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          ob = Py_InitModule3(name, methods, doc);
#endif

MOD_INIT(_ucrdtw)
{
    PyObject *m;

    MOD_DEF(m, "_ucrdtw", module_docstring,
            module_methods)

    if (m == NULL)
        return MOD_ERROR_VAL;

    import_array();

    return MOD_SUCCESS_VAL(m);
}

/*
Calculate the nearest n neighbors of a times series in a larger time series expressed as location and distance,
using the UCR suite optimizations.

Arguments:
data - a list of floats or an ndarray for the time series to process
query - a list of floats or an ndarray for the time series to search for
warp_width - Allowed warp width as a fraction of query size
n - number of top candidatews
verbose - Optional boolean to print stats
*/
static PyObject* ucrdtw_all(PyObject* self, PyObject* args) {
    PyObject* data_obj = NULL;
    PyObject* query_obj = NULL;
    double warp_width = -1;
    int n = -1;
    PyObject* verbose_obj = NULL;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOdi|O", &data_obj, &query_obj, &warp_width, &n, &verbose_obj)) {
        return NULL;
    }

    /* Interpret the input objects as numpy arrays. */
    if (!PyList_Check(data_obj) && !PyArray_Check(data_obj)) {
        PyErr_SetString(PyExc_TypeError, "Data argument must be a list or ndarray");
        return NULL;
    }
    PyObject* data_array = PyArray_FROM_OTF(data_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (data_array == NULL) {
        Py_XDECREF(data_array);
        PyErr_SetString(PyExc_TypeError, "Data argument must be a list or ndarray");
        return NULL;
    }

    if (!PyList_Check(query_obj) && !PyArray_Check(query_obj)) {
        Py_XDECREF(data_array);
        PyErr_SetString(PyExc_TypeError, "Query argument must be a list or ndarray");
        return NULL;
    }
    PyObject* query_array = PyArray_FROM_OTF(query_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (query_array == NULL) {
        Py_XDECREF(data_array);
        Py_XDECREF(query_array);
        PyErr_SetString(PyExc_TypeError, "Query argument must be a list or ndarray");
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double* data = (double*) PyArray_DATA(data_array);
    int data_size = (long) PyArray_DIM(data_array, 0);

    double* query = (double*) PyArray_DATA(query_array);
    int query_size = (int) PyArray_DIM(query_array, 0);

    int verbose = verbose_obj != NULL ? PyObject_IsTrue(verbose_obj) : 0;

    /* Call the external C function to compute the best DTW location and distance. */
    struct candidate *candidates = malloc(n * sizeof(candidate));
    int p_size = -1;
    int status = ucrdtwa(data, data_size, query, query_size, warp_width, verbose, n, candidates, &p_size);

    /* Clean up. */
    Py_XDECREF(data_array);
    Py_XDECREF(query_array);

    if (status) {
        PyErr_SetString(PyExc_RuntimeError, "ucrdtw could not allocate memory");
        return NULL;
    }

   for (int wtf = 0; wtf < p_size; wtf++) {
         printf("RETURNED result : %.12f \n", candidates[wtf].distance);
        }

   /* Build the output list up. */
    PyObject *lst = PyList_New(p_size);
    if (!lst)
        return NULL;
    for (int i = 0; i < p_size; i++) {
        PyObject* num = Py_BuildValue("Ld", candidates[i].index, sqrt(candidates[i].distance));
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
    PyList_SET_ITEM(lst, i, num);   // reference to num stolen
    }
    return lst;
}

/*
Calculate the nearest neighbor of a times series in a larger time series expressed as location and distance,
using the UCR suite optimizations.

Arguments:
data - a list of floats or an ndarray for the time series to process
query - a list of floats or an ndarray for the time series to search for
warp_width - Allowed warp width as a fraction of query size
verbose - Optional boolean to print stats
*/
static PyObject* ucrdtw_ucrdtw(PyObject* self, PyObject* args) {
    PyObject* data_obj = NULL;
    PyObject* query_obj = NULL;
    double warp_width = -1;
    PyObject* verbose_obj = NULL;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOd|O", &data_obj, &query_obj, &warp_width, &verbose_obj)) {
        return NULL;
    }

    /* Interpret the input objects as numpy arrays. */
    if (!PyList_Check(data_obj) && !PyArray_Check(data_obj)) {
        PyErr_SetString(PyExc_TypeError, "Data argument must be a list or ndarray");
        return NULL;
    }
    PyObject* data_array = PyArray_FROM_OTF(data_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (data_array == NULL) {
        Py_XDECREF(data_array);
        PyErr_SetString(PyExc_TypeError, "Data argument must be a list or ndarray");
        return NULL;
    }

    if (!PyList_Check(query_obj) && !PyArray_Check(query_obj)) {
        Py_XDECREF(data_array);
        PyErr_SetString(PyExc_TypeError, "Query argument must be a list or ndarray");
        return NULL;
    }
    PyObject* query_array = PyArray_FROM_OTF(query_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (query_array == NULL) {
        Py_XDECREF(data_array);
        Py_XDECREF(query_array);
        PyErr_SetString(PyExc_TypeError, "Query argument must be a list or ndarray");
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double* data = (double*) PyArray_DATA(data_array);
    int data_size = (long) PyArray_DIM(data_array, 0);

    double* query = (double*) PyArray_DATA(query_array);
    int query_size = (int) PyArray_DIM(query_array, 0);

    int verbose = verbose_obj != NULL ? PyObject_IsTrue(verbose_obj) : 0;

    /* Call the external C function to compute the best DTW location and distance. */
    long long location = -1;
    double distance = -1;
    int status = ucrdtw(data, data_size, query, query_size, warp_width, verbose, &location, &distance);

    /* Clean up. */
    Py_XDECREF(data_array);
    Py_XDECREF(query_array);

    if (status) {
        PyErr_SetString(PyExc_RuntimeError, "ucrdtw could not allocate memory");
        return NULL;
    }

    /* Build the output tuple */
    PyObject* ret = Py_BuildValue("Ld", location, distance);
    return ret;
}
