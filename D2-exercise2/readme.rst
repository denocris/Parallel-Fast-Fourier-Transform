============
D2-exercise2
============

1. Write a fftw parallel wrapper 
================================

Suppose you are on a machine for which the available fftw implementation does provide the fftw_mpi interface.
In such a case, you need to use the fftw serial interface and manage the communication with mpi yourself.

Take a look to the template provided_code/C/fft_wrapper.c. 
There you can see that you need to call 1d fftw for each dimension, but you need to be careful in managing the fftw in the distributed
direction. The All2all communication needs to be done properly, such that the sent data is positioned properly.

Remember that the data are distributed dividing the orginal 3D domain in a 1D decomposition among the outer dimension (see the following picture).

|1D-decomp|



2. Coding style
================
You can use c preprocessor macro to choose, at compilation time, whether to use the fftw_mpi calls, or the homemade one.
To do this, you can use 

.. code:: c

	#ifdef __HOMEMADE
	.....
	#else
	....
	#endif

.. |1D-decomp| image:: ./1D-decomposition.png
   :alt: Frequencies
   :scale:  100%
   :align: middle

