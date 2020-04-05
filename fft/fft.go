package fft

import (
	"math"
	"math/bits"
	"math/cmplx"
)

//FFT computes the fast fourier transform of the data x.
//The fast fourier transform is computed by using an out of place radix-2 fft for lengths which are powers of 2 and using the Chirp-z transform followed by radix-2 for lengths which are not a power of 2.
func FFT(x []complex128) []complex128 {
	n := len(x)
	if isp2(uint64(n)) {
		return fft2(x, n, 1)
	}
	return chirpz(x)
}

//IFFT computes the inverse fast fourier transform of the data x.
func IFFT(x []complex128) []complex128 {
	n := len(x)
	for i := 0; i < n; i++ {
		x[i] = cmplx.Conj(x[i])
	}
	out := FFT(x)
	for i := 0; i < n; i++ {
		x[i] = cmplx.Conj(x[i])
		out[i] = cmplx.Conj(out[i])
	}
	for i := range out {
		out[i] /= complex(float64(n), 0)
	}
	return out
}

func chirpz(x []complex128) []complex128 {
	n := len(x)
	N := int(clp2(2*uint64(n) - 1))
	y := make([]complex128, N)
	factors := make([]complex128, n)
	for i := range factors {
		sin, cos := math.Sincos(math.Pi / float64(n) * float64(i*i))
		factors[i] = complex(cos, sin)
	}
	for i, v := range x {
		y[i] = v * cmplx.Conj(factors[i])
	}
	z := make([]complex128, N)
	for i := range x {
		z[i] = factors[i]
		if i != 0 {
			z[N-i] = factors[i]
		}
	}

	y = fft2(y, N, 1)
	z = fft2(z, N, 1)

	for j := range y {
		y[j] *= z[j]
	}

	y = ifft2(y)

	for i := range factors {
		y[i] *= cmplx.Conj(factors[i])
	}

	return y[:n]
}

func clp2(x uint64) uint64 {
	nlz := bits.LeadingZeros64(x - 1)
	return 1 << uint(64-nlz)
}

//Returns true if
func isp2(x uint64) bool {
	return (x & (x - 1)) == 0
}

func fft2(x []complex128, N int, s int) []complex128 {
	out := make([]complex128, N)
	if N <= 1 {
		copy(out, x)
		return out
	}
	out1 := fft2(x, N/2, 2*s)
	out2 := fft2(x[s:], N/2, 2*s)

	for k := 0; k < N/2; k++ {
		sin, cos := math.Sincos(-2 * math.Pi * float64(k) / float64(N))
		twiddle := complex(cos, sin)
		out[k] = out1[k] + twiddle*out2[k]
		out[N/2+k] = out1[k] - twiddle*out2[k]
	}
	return out
}

func ifft2(x []complex128) []complex128 {
	n := len(x)
	for i := 0; i < n; i++ {
		x[i] = cmplx.Conj(x[i])
	}
	out := fft2(x, n, 1)
	for i := 0; i < n; i++ {
		x[i] = cmplx.Conj(x[i])
		out[i] = cmplx.Conj(out[i])
	}
	for i := range out {
		out[i] /= complex(float64(n), 0)
	}
	return out
}
