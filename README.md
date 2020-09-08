# infCS

Code to reproduce numerical experiments from [Non-uniform recovery guarantees for binary measurements and infinite-dimensional compressed sensing](https://arxiv.org/abs/1909.01143) . 

The [cww](https://github.com/vegarant/cww) library from Antun and the [GS-wavelets](https://github.com/milana-gataric/gs-wavelets) library from Gataric and Poon are extended to allow the numerical implementation of infinite wavelet coefficient bandwidth. See [A Practical Guide to the Recovery of Wavelet Coefficients from Fourier Measurements](https://epubs.siam.org/doi/abs/10.1137/15M1018630) for information about the implementation of the recovery method in the Fourier setting and [Generalized sampling and infinite-dimensional compressed sensing](https://www.repository.cam.ac.uk/bitstream/handle/1810/284101/BAACHGSCS_Rev10.pdf?sequence=1) for details on the discretization of the infinite problem setting.

### Dependencies

- [cww](https://github.com/vegarant/cww) Continuous Walsh sampling and wavelet reconstruction
- [fastwht](https://bitbucket.org/vegarant/fastwht/src/master/) Fast Walsh transform
- [SPGL-1](https://github.com/vegarant/spgl1) A solver for large-scale sparse reconstruction 
- [GS-wavelets](https://github.com/milana-gataric/gs-wavelets) Generalized sampling for the reconstruction of wavelet coefficients from Fourier measurements (with the related [Publication](https://epubs.siam.org/doi/abs/10.1137/15M1018630))

### Getting started

The script “multiple_examples” chooses the hyperparameter which are altered for the different experiments in the paper. It is possible to change other hyperparameter in “Example_handle_1D” and “Example_handle_1D_flip”. The two scripts are an adaptation of "Ex_CS_1d" in cww. 

For the Fourier examples use "CS_Fourier" in the file the wavelet as well as the number of samples can be chosen accordingly. The script is an adaptation of "GSCS_1D" in GS-wavelets.
