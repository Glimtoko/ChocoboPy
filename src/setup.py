# python setup.py build_ext --inplace

from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        "lagrangian_hydro",
        ["lagrangian_hydro.pyx"],
        include_dirs=[numpy.get_include()],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
    )
]

setup(
    name='lagrangian_hydro',
    ext_modules=cythonize(ext_modules),
)

# setup(
#     name='Hydrocode',
#     ext_modules=cythonize("lagrangian_hydro.pyx"),
#     include_dirs=[numpy.get_include()],
#     zip_safe=False,
# )
