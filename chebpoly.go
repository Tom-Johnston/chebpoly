package chebpoly

import "math"
import "github.com/mjibson/go-dsp/fft"
import "fmt"

//Chebpts(n, dl, du) gives the n CHebyshev points of the second kind suitably scaled to the domain [dl,du]
func Chebpts(n uint, domainLower float64, domainUpper float64) []float64 {
	if domainLower > domainUpper {
		fmt.Println("Warning: domainLower is greater than domainUpper. Switching them round...")
		domainLower, domainUpper = domainUpper, domainLower
	}
	points := make([]float64, n)
	N := int(n)
	for i, _ := range points {
		points[i] = (1+math.Cos(math.Pi*float64(N-1-i)/float64(N-1)))*(domainUpper-domainLower)/2 + domainLower
	}

	return points
}

type Chebpoly struct {
	DomainLower float64
	DomainUpper float64

	Coeffs []float64 //Coeffs[k] is the coefficient of kth Chebyshev polynomial of the fist kind
}

func (poly Chebpoly) Length() int {
	return len(poly.Coeffs)
}

//Interp(points) returns the Chebyshev interpolate taking value points[i] at the ith Chebyshev point of the second kind.
func Interp(values []float64, domainLower float64, domainUpper float64) *Chebpoly {
	n := len(values)
	poly := new(Chebpoly)

	if domainLower > domainUpper {
		fmt.Println("Warning: domainLower is greater than domainUpper. Switching them round...")
		domainLower, domainUpper = domainUpper, domainLower
	}
	poly.DomainLower = domainLower
	poly.DomainUpper = domainUpper

	expandedData := make([]float64, 2*(n-1))
	for i := 0; i < n-1; i++ {
		expandedData[n-1-i] = values[i]
		expandedData[n-1+i] = values[i]
	}
	expandedData[0] = values[n-1]
	dctValues := fft.IFFTReal(expandedData)
	expandedData[0] = real(dctValues[0])
	expandedData[n-1] = real(dctValues[n-1])
	for i := 1; i < n-1; i++ {
		expandedData[i] = 2 * real(dctValues[i])
	}

	poly.Coeffs = expandedData[0:n]

	return poly
}

//Evaluate(x) returns the value of the Chebypoly poly at the point x. It uses Clenshaw's algorithm.
func (poly Chebpoly) Evaluate(x float64) float64 {
	scaledX := (x-poly.DomainLower)/(poly.DomainUpper-poly.DomainLower)*2 - 1
	bk2 := 0.0
	bk1 := 0.0
	for i := poly.Length() - 1; i > 1; i -= 2 {
		bk2 = poly.Coeffs[i] + 2*scaledX*bk1 - bk2
		bk1 = poly.Coeffs[i-1] + 2*scaledX*bk2 - bk1
	}
	if ((poly.Length() - 1) & 1) != 0 {
		//We are of odd degree so bk2 is b3 and bk1 is b2
		bk2 = poly.Coeffs[1] + 2*scaledX*bk1 - bk2
		return poly.Coeffs[0] + scaledX*bk2 - bk1
	} else {
		//We are of even degree so bk2 is b2 and bk1 is b1
		return poly.Coeffs[0] + scaledX*bk1 - bk2
	}
}

//Cumsum returns the Chebpoly of the indefinite integral.
//The constant of integration is taken so that Cumsum(DomainLower) = 0
func (poly Chebpoly) Cumsum() Chebpoly {
	n := poly.Length()
	integral := Chebpoly{DomainLower: poly.DomainLower, DomainUpper: poly.DomainUpper}
	integral.Coeffs = make([]float64, poly.Length()+1)
	scaleFactor := (poly.DomainUpper - poly.DomainLower) / 2 //Multiply all the coefficients by this factor to account for the change of variables in the integral.
	b0 := 0.0                                                //This will be the zero coefficient of the integrand
	integral.Coeffs[1] = scaleFactor * (poly.Coeffs[0] - poly.Coeffs[2]/2)
	b0 += integral.Coeffs[1]
	for i := 2; i < n-1; i++ {
		integral.Coeffs[i] = scaleFactor * (poly.Coeffs[i-1] - poly.Coeffs[i+1]) / (2 * float64(i))
		if (i % 2) == 0 {
			b0 -= integral.Coeffs[i]
		} else {
			b0 += integral.Coeffs[i]
		}
	}
	integral.Coeffs[n-1] = scaleFactor * poly.Coeffs[n-2] / (2 * float64(n-1))
	integral.Coeffs[n] = scaleFactor * poly.Coeffs[n-1] / (2 * float64(n))
	if (n % 2) == 1 {
		b0 -= integral.Coeffs[n-1]
		b0 += integral.Coeffs[n]
	} else {
		b0 += integral.Coeffs[n-1]
		b0 -= integral.Coeffs[n]
	}
	integral.Coeffs[0] = b0
	return integral
}
