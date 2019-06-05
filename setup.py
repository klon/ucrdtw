from distutils.core import setup, Extension
import numpy.distutils.misc_util

c_ext = Extension('_ucrdtw', ['src/_ucrdtw.c', 'src/ucrdtw.c'])

setup(
    name='ucrdtw',
    version='1.0',
    ext_modules=[c_ext],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)
