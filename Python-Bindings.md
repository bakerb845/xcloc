# Python

Python3 bindings are being crafted to interface to the xcloc library.  Unless the library path is explicitly specified then it is important to understand that Python will scrape the system's load library path, LD\_LIBRARY\_PATH in attempt to find the xcloc shared library.  After the library is obtained a [ctypes](https://docs.python.org/3/library/ctypes.html) interface is used to access xcloc.  

In the future the bindings to the MPI variant of xcloc will be added.  In this instance the user will need to remember to call the script with mpirun and have the [MPI for Python](http://mpi4py.scipy.org/docs/) library.


