# infCS

Code to reproduce numerical experiments from [Non-uniform recovery guarantees for binary measurements and infinite-dimensional compressed sensing](https://arxiv.org/abs/1909.01143) . 

The [cww](https://github.com/vegarant/cww) library from Antun is extended to allow the numerical implementation of infinite wavelet coefficient bandwidth. See [Generalized sampling and infinite-dimensionalcompressed sensing](https://www.repository.cam.ac.uk/bitstream/handle/1810/284101/BAACHGSCS_Rev10.pdf?sequence=1) for details on the discretization of the infinite problem setting.

### Dependencies

- [cww](https://github.com/vegarant/cww) Continuous Walsh sampling and wavelet reconstruction
- [fastwht](https://bitbucket.org/vegarant/fastwht/src/master/) Fast Walsh transform
- [SPGL-1](https://github.com/vegarant/spgl1) A solver for large-scale sparse reconstruction 

### Getting started

The script “multiple_examples” chooses the hyperparameter which are altered for the different experiments in the paper. It is possible to change other hyperparameter in “Example_handle_1D” and “Example_handle_1D_flip”. The two scripts are closely related to "Ex_CS_1d" in cww. 
