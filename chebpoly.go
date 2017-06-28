package chebpoly

import "math"
import "github.com/mjibson/go-dsp/fft"
import "github.com/gonum/matrix/mat64"
import "fmt"
import "sort"

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

func (poly Chebpoly) Sum() float64 {
	sum := 0.0
	for i := 0; i < poly.Length(); i += 2 {
		sum += -poly.Coeffs[i] * 2 / float64(i*i-1)
	}
	return sum * (poly.DomainUpper - poly.DomainLower) / 2
}

func (poly Chebpoly) Diff() Chebpoly {
	n := poly.Length()
	derivative := Chebpoly{DomainLower: poly.DomainLower, DomainUpper: poly.DomainUpper}
	derivative.Coeffs = make([]float64, n-1)
	scaleFactor := 2 / (poly.DomainUpper - poly.DomainLower)
	derivative.Coeffs[n-2] = scaleFactor * 2 * float64(n-1) * poly.Coeffs[n-1]
	derivative.Coeffs[n-3] = scaleFactor * 2 * float64(n-2) * poly.Coeffs[n-2]
	for i := n - 4; i > 0; i-- {
		derivative.Coeffs[i] = derivative.Coeffs[i+2] + scaleFactor*2*float64(i+1)*poly.Coeffs[i+1]
	}
	derivative.Coeffs[0] = derivative.Coeffs[2]/2 + scaleFactor*poly.Coeffs[1]
	return derivative
}

func (poly Chebpoly) Roots() []float64 {

	roots := poly.getApproximateRoots()
	roots = poly.refineRoots(roots)
	sort.Float64s(roots)
	return roots
}

func (poly Chebpoly) getApproximateRoots() []float64 {
	n := poly.Length()

	//Truncate the series. We are only after the roots approximately anyway.

	//Calculate the 1 norm
	norm := 0.0
	for _, v := range poly.Coeffs {
		norm += math.Abs(v)
	}
	cutoffValue := 1e-13*norm + 1e-15
	cutOffPoint := 1
	//Iterate back through the coefficients to find the first point that we
	for i := 0; i < n-1; i++ {
		if math.Abs(poly.Coeffs[n-1-i]) >= cutoffValue {
			cutOffPoint = n - i
			break
		}
	}

	coeffs := poly.Coeffs[:cutOffPoint]
	n = len(coeffs)
	//We will find the roots by solving the Colleague matrix eigenvalue problem if it is at most maxEigProbSize else we will split the interval in 2 and try again.

	const maxEigProbSize = 50 //Note that making this too small can cause problems as the splitting can reach intervals of lenth around machine precision.

	const complexTolerance = 1e-14
	const realTolerance = 1e-14

	if n == 1 && coeffs[0] != 0 {
		return []float64{}
	}
	if n == 1 {
		// We will return the middle of the domain like chebfun...
		fmt.Println("Warning: Finding roots of the 0 polynomial")
		return []float64{(poly.DomainLower + poly.DomainUpper) / 2}
	}
	if n == 2 {
		root := -coeffs[0] / coeffs[1]
		if math.Abs(root) <= 1+realTolerance {
			return []float64{(1+root)*(poly.DomainUpper-poly.DomainLower)/2 + poly.DomainLower}
		} else {
			return []float64{}
		}
	}
	if n <= maxEigProbSize {
		roots := make([]float64, 0, n-1)
		a := make([]float64, (n-1)*(n-1))
		a[1] = 1
		for i := 1; i < n-2; i++ {
			a[i*(n-1)+i-1] = 0.5
			a[i*(n-1)+i+1] = 0.5
		}
		for i := 0; i < n-1; i++ {
			a[(n-1)*(n-2)+i] = -coeffs[i] / (2 * coeffs[n-1])
		}
		a[(n-1)*(n-1)-2] = -coeffs[n-3]/(2*coeffs[n-1]) + 0.5
		A := mat64.NewDense(n-1, n-1, a)
		var eigen mat64.Eigen
		eigen.Factorize(A, false, false)
		eigenValues := make([]complex128, n-1)
		eigen.Values(eigenValues)

		for _, v := range eigenValues {
			if math.Abs(imag(v)) < complexTolerance && math.Abs(real(v)) <= 1+realTolerance {
				scaledV := (1+real(v))*(poly.DomainUpper-poly.DomainLower)/2 + poly.DomainLower
				roots = append(roots, scaledV)
			}
		}
		return roots
	} else {
		//Split
		const splitPoint = -0.004849834917525 // This is an arbitrary number copied from Chefun.

		scaledSplitPoint := (1+splitPoint)*(poly.DomainUpper-poly.DomainLower)/2 + poly.DomainLower

		valuesLeft := make([]float64, n)
		for i, v := range Chebpts(uint(n), poly.DomainLower, scaledSplitPoint) {
			valuesLeft[i] = poly.Evaluate(v)
		}
		valuesRight := make([]float64, n)
		for i, v := range Chebpts(uint(n), scaledSplitPoint, poly.DomainUpper) {
			valuesRight[i] = poly.Evaluate(v)
		}

		polyLeft := Interp(valuesLeft, poly.DomainLower, scaledSplitPoint)
		polyRight := Interp(valuesRight, scaledSplitPoint, poly.DomainUpper)

		return append(polyLeft.getApproximateRoots(), polyRight.getApproximateRoots()...)

	}
}

func (poly Chebpoly) refineRoots(roots []float64) []float64 {
	return roots
}
