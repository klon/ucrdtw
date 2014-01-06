#include <Python.h>
#include <numpy/arrayobject.h>
#include "ucrdtw.h"

/* Docstrings */
static char module_docstring[] = "This module implements fast nearest-neighbor retrieval under the dynamic time warping (DTW).";
static char ucrdtw_docstring[] = "Calculate the nearest neighbor of a times series in a larger time series expressed as location and distance, "
                                 "using the UCR suite optimizations.";

/* Available functions */
static PyObject* ucrdtw_ucrdtw(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = { { "ucrdtw", ucrdtw_ucrdtw, METH_VARARGS, ucrdtw_docstring }, { NULL, NULL, 0, NULL } };

/* Initialize the module */
PyMODINIT_FUNC init_ucrdtw(void) {
    PyObject* m = Py_InitModule3("_ucrdtw", module_methods, module_docstring);
    if (m == NULL) {
        return;
    }

    /* Load numpy functionality. */
    import_array();
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
        Py_XDECREF(verbose_obj);
        PyErr_SetString(PyExc_TypeError, "Data argument must be a list or ndarray");
        return NULL;
    }

    if (!PyList_Check(query_obj) && !PyArray_Check(query_obj)) {
        PyErr_SetString(PyExc_TypeError, "Query argument must be a list or ndarray");
        return NULL;
    }
    PyObject* query_array = PyArray_FROM_OTF(query_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (query_array == NULL) {
        Py_XDECREF(data_array);
        Py_XDECREF(query_array);
        Py_XDECREF(verbose_obj);
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
    Py_DECREF(data_array);
    Py_DECREF(query_array);
    Py_XDECREF(verbose_obj);

    if (status) {
        PyErr_SetString(PyExc_RuntimeError, "ucrdtw could not allocate memory");
        return NULL;
    }

    /* Build the output tuple */
    PyObject* ret = Py_BuildValue("Ld", location, distance);
    return ret;
}
