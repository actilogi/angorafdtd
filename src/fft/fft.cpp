//Definitions of the FFT wrappers

//Whatever FFT is used, it should have the following:
// 1) 2D FFT support (for imaging)
// 2) 3D FFT support (for random-medium generation)
// 3) complex<double> type support at the input
// 4) Different FFT size than the input size

#include "headers.h"

//#include "fft.h"

// #include <gsl/gsl_fft_complex.h>
// #include <fftw3.h>

//KissFFT library header
#include "kissfft/kiss_fftnd.h"
// Although KissFFT does not *natively* support arbitrary FFT sizes, this can be partially alleviated by resizing the input to the desired FFT size and applying zero padding.


void perform_1D_FFT(const Array<complex<kiss_fft_scalar>,1>& input, Array<complex<kiss_fft_scalar>,1>& output, const int& N_FFT, const int& inverse_fft)
{//computes the 1D Fourier transform
	kiss_fft_cfg mycfg=kiss_fft_alloc(N_FFT,inverse_fft,NULL,NULL); //computes inverse FFT if inverse_fft is not 0  NOTE:inverse FFT is NOT scaled by 1/(N_FFT)
	//Here, we assume that kiss_fft_cpx is bit-compatible to the C++ complex<T> datatype
	//kiss_fft_cpx is just a struct with two elements of type T, so it is sufficiently simple for this assumption
	const kiss_fft_cpx* input_kissfft = reinterpret_cast<const kiss_fft_cpx*>(input.data());
	kiss_fft_cpx* output_kissfft = reinterpret_cast<kiss_fft_cpx*>(output.data());
	kiss_fft(mycfg,input_kissfft,output_kissfft);
	free(mycfg);
}

void perform_2D_FFT(const Array<complex<kiss_fft_scalar>,2>& input, Array<complex<kiss_fft_scalar>,2>& output, const int& N1, const int& N2, const int& inverse_fft)
{//computes the 2D Fourier transform
	const int ndims = 2;
	const int dims[ndims] = {N1,N2};
	kiss_fftnd_cfg mycfg=kiss_fftnd_alloc(dims,ndims,inverse_fft,NULL,NULL); //computes inverse FFT if inverse_fft is not 0   NOTE:inverse FFT is NOT scaled by 1/(N1*N2)
	//Here, we assume that kiss_fft_cpx is bit-compatible to the C++ complex<T> datatype
	//kiss_fft_cpx is just a struct with two elements of type T, so it is sufficiently simple for this assumption
	const kiss_fft_cpx* input_kissfft = reinterpret_cast<const kiss_fft_cpx*>(input.data());
	kiss_fft_cpx* output_kissfft = reinterpret_cast<kiss_fft_cpx*>(output.data());
	kiss_fftnd(mycfg,input_kissfft,output_kissfft);
	free(mycfg); //mycfg is allocated dynamically, so free it
}

void perform_3D_FFT(const Array<complex<kiss_fft_scalar>,3>& input, Array<complex<kiss_fft_scalar>,3>& output, const int& N1, const int& N2, const int& N3, const int& inverse_fft)
{//computes the 3D Fourier transform
	const int ndims = 3;
	const int dims[ndims] = {N1,N2,N3};
	kiss_fftnd_cfg mycfg=kiss_fftnd_alloc(dims,ndims,inverse_fft,NULL,NULL); //computes inverse FFT if inverse_fft is not 0   NOTE:inverse FFT is NOT scaled by 1/(N1*N2)
	//Here, we assume that kiss_fft_cpx is bit-compatible to the C++ complex<T> datatype
	//kiss_fft_cpx is just a struct with two elements of type T, so it is sufficiently simple for this assumption
	const kiss_fft_cpx* input_kissfft = reinterpret_cast<const kiss_fft_cpx*>(input.data());
	kiss_fft_cpx* output_kissfft = reinterpret_cast<kiss_fft_cpx*>(output.data());
	kiss_fftnd(mycfg,input_kissfft,output_kissfft);
	free(mycfg); //mycfg is allocated dynamically, so free it
}
