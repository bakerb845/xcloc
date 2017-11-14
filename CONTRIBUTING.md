# Ways to Help

This code is being actively developed.  There exists enough to handle event location serially but the full parallel code is not yet there.  If you want to get involved a few sub-tasks that come to mind are:

1. Try things, break things, and make an [issue](https://github.com/bakerb845/xcloc/issues).
2. If you'd like to write some code it would be worth reworking the RMS filter so that it works with OpenMP.
3. An envelope computation would be useful.  It's sketched out but not lashed in.
4. Be brave and set the floats and float complexes in the xcfft.h to void * then try out the mixed precision code.
5. Tightly integrate with the travel-time computations in https://github.com/bakerb845/fteik.  
6. Extend the HDF5 IO to spherical coordinates and/or unstructured grids.
7. Ask.  There's always stuff to do.

