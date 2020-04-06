package chebpoly

import (
	"fmt"
	"math"
	"sort"

	"github.com/Tom-Johnston/chebpoly/fft"
)

//Chebpts gives the n Chebyshev points of the second kind suitably scaled to the domain [dl,du].
func Chebpts(n uint, domainLower float64, domainUpper float64) []float64 {
	if domainLower > domainUpper {
		domainLower, domainUpper = domainUpper, domainLower
	}
	points := make([]float64, n)
	if n == 1 {
		points[0] = 0
		return points
	}

	N := int(n)
	for i := range points {
		points[i] = (1+math.Cos(math.Pi*float64(N-1-i)/float64(N-1)))*(domainUpper-domainLower)/2 + domainLower
	}

	return points
}

//Chebpoly contains a polynomial approximation to a real valued function on a closed bounded interval represented as a sum of Chebyshev polynomials of the first kind.
type Chebpoly struct {
	DomainLower float64
	DomainUpper float64

	Coeffs []float64 //coeffs[k] is the coefficient of kth Chebyshev polynomial of the fist kind
}

//Length returns the number of Chebyshev polynomials in the approximation including those with coefficient 0.
//The degree of the polynomial is one less than the length.
func (poly Chebpoly) Length() int {
	return len(poly.Coeffs)
}

//Copy returns a deep copy of the Chebpoly.
func (poly Chebpoly) Copy() Chebpoly {
	copyCoeffs := make([]float64, len(poly.Coeffs))
	copy(copyCoeffs, poly.Coeffs)
	return Chebpoly{DomainLower: poly.DomainLower, DomainUpper: poly.DomainUpper, Coeffs: copyCoeffs}
}

//String prints basic information about the Chebpoly.
func (poly Chebpoly) String() string {
	if poly.Length() < 10 {
		return fmt.Sprintf("Interval: [%f %f]  Degree: %d  coeffs: %v", poly.DomainLower, poly.DomainUpper, poly.Length()-1, poly.Coeffs)
	}
	return fmt.Sprintf("Interval: [%f %f]  Degree: %d  coeffs: [%f %f %f %f %f %f %f %f %f %f ...]", poly.DomainLower, poly.DomainUpper, poly.Length()-1, poly.Coeffs[0], poly.Coeffs[1], poly.Coeffs[2], poly.Coeffs[3], poly.Coeffs[4], poly.Coeffs[5], poly.Coeffs[6], poly.Coeffs[7], poly.Coeffs[8], poly.Coeffs[9])
}

//New returns the Chebpoly on the given domain with the specified coefficients.
func New(coeffs []float64, domainLower float64, domainUpper float64) Chebpoly {
	if domainLower > domainUpper {
		domainLower, domainUpper = domainUpper, domainLower
	}

	if len(coeffs) == 0 {
		return Chebpoly{DomainLower: domainLower, DomainUpper: domainUpper, Coeffs: []float64{0}}
	}
	return Chebpoly{DomainLower: domainLower, DomainUpper: domainUpper, Coeffs: coeffs}
}

//Interp returns the Chebyshev approximation  on [domainLower, domainUpper] which takes the value values[i] at the ith Chebyshev point of the second kind scaled to [domainLower, domainUpper].
func Interp(values []float64, domainLower float64, domainUpper float64) Chebpoly {
	n := len(values)
	poly := new(Chebpoly)

	if domainLower > domainUpper {
		domainLower, domainUpper = domainUpper, domainLower
	}
	poly.DomainLower = domainLower
	poly.DomainUpper = domainUpper

	if n == 1 {
		poly.Coeffs = []float64{values[0]}
		return *poly
	}

	expandedData := make([]complex128, 2*(n-1))
	for i := 0; i < n-1; i++ {
		expandedData[n-1-i] = complex(values[i], 0)
		expandedData[n-1+i] = complex(values[i], 0)
	}
	expandedData[0] = complex(values[n-1], 0)
	dctValues := fft.IFFT(expandedData)
	coeffs := make([]float64, n)
	coeffs[0] = real(dctValues[0])
	coeffs[n-1] = real(dctValues[n-1])
	for i := 1; i < n-1; i++ {
		coeffs[i] = 2 * real(dctValues[i])
	}

	poly.Coeffs = coeffs

	return *poly
}

//Adaptive creates a Chebyshev approximation to the function f on [domainLower, domainUpper].
func Adaptive(f func(float64) float64, domainLower float64, domainUpper float64) Chebpoly {
	const precision = 1e-15
	const maxSize = 65536

	if domainLower > domainUpper {
		domainLower, domainUpper = domainUpper, domainLower
	}

	const startingDegree = 16
	currentValues := make([]float64, startingDegree+1)
	points := Chebpts(startingDegree+1, domainLower, domainUpper)
	maxValue := 0.0
	for i, v := range points {
		w := f(v)
		currentValues[i] = w
		if math.Abs(w) > maxValue {
			maxValue = math.Abs(w)
		}
	}
	currentDegree := startingDegree
	currentPoly := Interp(currentValues, domainLower, domainUpper)
	for true {
		points := Chebpts(uint(2*currentDegree+1), domainLower, domainUpper)
		currentValues = append(currentValues, make([]float64, currentDegree)...)
		for i := currentDegree; i >= 1; i-- {
			currentValues[2*i] = currentValues[i]
		}

		imprecise := false
		for i := 0; i < currentDegree; i++ {
			w := f(points[2*i+1])
			currentValues[2*i+1] = w
			if math.Abs(w) > maxValue {
				maxValue = math.Abs(w)
			}
		}
		for i := 0; i < currentDegree; i++ {
			if !imprecise && math.Abs(Evaluate(currentPoly, points[2*i+1])-currentValues[2*i+1]) > precision*maxValue+precision+1e-16*float64(currentDegree) {
				imprecise = true
			}
		}

		if imprecise {
			//We aren't accurate on the new points so improve
			currentPoly = Interp(currentValues, domainLower, domainUpper)
			currentDegree = 2 * currentDegree
			if currentDegree > maxSize {
				panic("Unable to approximate with the maximum number of points.")
			}
		} else {
			break
		}
	}

	for i := len(currentPoly.Coeffs) - 1; i >= 0; i-- {
		if currentPoly.Coeffs[i] < precision {
			currentPoly.Coeffs = currentPoly.Coeffs[:i]
		} else {
			break
		}
	}

	return currentPoly
}

//Add returns the sum of two Chebpolys on the same domain.
func Add(a, b Chebpoly) Chebpoly {
	if a.DomainLower != b.DomainLower || a.DomainUpper != b.DomainUpper {
		panic("The chebpolys have different domains")
	}
	if a.Length() < b.Length() {
		return Add(b, a)
	}
	sum := a.Copy()
	for i := 0; i < b.Length(); i++ {
		sum.Coeffs[i] += b.Coeffs[i]
	}
	return sum
}

//Multiply returns the product of two Chebpolys on the same domain.
func Multiply(a, b Chebpoly) Chebpoly {
	if a.DomainLower != b.DomainLower || a.DomainUpper != b.DomainUpper {
		panic("The chebpolys have different domains")
	}
	coeffs := make([]float64, a.Length()+b.Length()-1)
	for i := 0; i < a.Length(); i++ {
		for j := 0; j < b.Length(); j++ {
			coeffs[i+j] += 0.5 * a.Coeffs[i] * b.Coeffs[j]
			if i > j {
				coeffs[i-j] += 0.5 * a.Coeffs[i] * b.Coeffs[j]
			} else {
				coeffs[j-i] += 0.5 * a.Coeffs[i] * b.Coeffs[j]
			}
		}
	}
	return Chebpoly{DomainLower: a.DomainLower, DomainUpper: b.DomainUpper, Coeffs: coeffs}
}

//Values returns the values of the approximation at the n + 1 Chebyshev points of the second kind where n is the degree of the approximation.
func (poly Chebpoly) Values() []float64 {
	n := len(poly.Coeffs)
	vals := make([]float64, n)
	expandedData := make([]complex128, 2*(n-1))
	expandedData[0] = complex(poly.Coeffs[0], 0)
	for i := 1; i < n-1; i++ {
		expandedData[i] = complex(poly.Coeffs[i]/2, 0)
		expandedData[n-1+i] = complex(poly.Coeffs[n-1-i]/2, 0)
	}
	expandedData[n-1] = complex(poly.Coeffs[n-1], 0)
	fftOut := fft.FFT(expandedData)
	for i := range vals {
		vals[i] = real(fftOut[n-1-i])
	}
	return vals
}

//Evaluate returns the value of the Chebypoly poly at the point x. It uses Clenshaw's algorithm.
func Evaluate(poly Chebpoly, x float64) float64 {
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
	}
	//We are of even degree so bk2 is b2 and bk1 is b1
	return poly.Coeffs[0] + scaledX*bk1 - bk2
}

//Cumsum returns the chebpoly of the indefinite integral.
//The constant of integration is taken so that Cumsum(domainLower) = 0
func Cumsum(poly Chebpoly) Chebpoly {
	n := poly.Length()
	integral := Chebpoly{DomainLower: poly.DomainLower, DomainUpper: poly.DomainUpper}
	integral.Coeffs = make([]float64, poly.Length()+1)
	scaleFactor := (poly.DomainUpper - poly.DomainLower) / 2 //Multiply all the coefficients by this factor to account for the change of variables in the integral.
	b0 := 0.0                                                //This will be the zero coefficient of the integrand
	if n > 2 {
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
	} else if n == 2 {
		integral.Coeffs[1] = scaleFactor * poly.Coeffs[0]
		integral.Coeffs[2] = scaleFactor * poly.Coeffs[1] / 4
	} else if n == 1 {
		integral.Coeffs[1] = scaleFactor * poly.Coeffs[0]
	}
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

//Sum returns the integral of the chebpoly over the domain.
func Sum(poly Chebpoly) float64 {
	sum := 0.0
	for i := 0; i < poly.Length(); i += 2 {
		sum += -poly.Coeffs[i] * 2 / float64(i*i-1)
	}
	return sum * (poly.DomainUpper - poly.DomainLower) / 2
}

//Diff returns the chebpoly representing the derivative of the given chebpoly
func Diff(poly Chebpoly) Chebpoly {
	n := poly.Length()
	derivative := Chebpoly{DomainLower: poly.DomainLower, DomainUpper: poly.DomainUpper}
	if n == 1 {
		//poly is a constant so derivative is 0.
		derivative.Coeffs = make([]float64, 1)
		return derivative
	}
	coeffs := make([]float64, n+1) //We will add a couple of extra zero entries here to make the code simpler. Note that these entries will probably always be in memory but this is not really a problem.
	scaleFactor := 2 / (poly.DomainUpper - poly.DomainLower)
	for i := n - 2; i > 0; i-- {
		coeffs[i] = coeffs[i+2] + scaleFactor*2*float64(i+1)*poly.Coeffs[i+1]
	}
	coeffs[0] = coeffs[2]/2 + scaleFactor*poly.Coeffs[1]
	derivative.Coeffs = coeffs[:n-1]
	return derivative
}

//Split returns the Chebpolys from doing separate Chebyshev approximations for each interval created by splitting [domainLower, domainUpper] at the splitPoints.
func Split(poly Chebpoly, splitPoints ...float64) []Chebpoly {
	n := poly.Length()
	sections := make([]Chebpoly, len(splitPoints)+1)
	values := make([]float64, n)
	points := Chebpts(uint(n), 0, 1) //Note that we manipulate the same set of points instead of calling cos multiple times.
	scaleFactor := splitPoints[0] - poly.DomainLower
	for i, v := range points {
		values[i] = Evaluate(poly, v*scaleFactor+poly.DomainLower)
	}
	sections[0] = Interp(values, poly.DomainLower, splitPoints[0])

	for i := 1; i < len(splitPoints); i++ {
		scaleFactor = (splitPoints[i] - splitPoints[i-1])
		for j, v := range points {
			values[j] = Evaluate(poly, v*scaleFactor+splitPoints[i-1])
		}
		sections[i] = Interp(values, splitPoints[i-1], splitPoints[i])
	}

	scaleFactor = poly.DomainUpper - splitPoints[len(splitPoints)-1]
	for i, v := range points {
		values[i] = Evaluate(poly, v*scaleFactor+splitPoints[len(splitPoints)-1])
	}
	sections[len(splitPoints)] = Interp(values, splitPoints[len(splitPoints)-1], poly.DomainUpper)
	return sections
}

//Roots returns all the real roots of the chebpoly in [-1,1]
//This could be extended to provide all the roots of the chebpoly but it isn't.
func Roots(poly Chebpoly) []float64 {

	roots := poly.getApproximateRoots(1e-15)
	roots = poly.refineRoots(roots)
	sort.Float64s(roots)
	return roots
}

//This finds the roots in the [-1,1] but they may not be perfectly accurate.
func (poly Chebpoly) getApproximateRoots(tol float64) []float64 {
	//We will be quite generous here and then we will double check that a root is at least close to 0 in the refinement.The danger is repeating a root too many times.

	n := poly.Length()
	//Truncate the series. We are only after the roots approximately anyway.

	//Calculate the 1 norm
	norm := 0.0
	for _, v := range poly.Coeffs {
		norm += math.Abs(v)
	}
	cutoffValue := 1e-13 * norm
	cutOffPoint := 1
	//Iterate back through the coefficients to find the first point that we
	for i := 0; i < n-1; i++ {
		if math.Abs(poly.Coeffs[n-1-i]) > cutoffValue {
			cutOffPoint = n - i
			break
		}
	}

	coeffs := poly.Coeffs[:cutOffPoint]
	n = len(coeffs)
	//We will find the roots by solving the Colleague matrix eigenvalue problem if it is at most maxEigProbSize else we will split the interval in 2 and try again.
	const maxEigProbSize = 50 //Note that making this too small can cause problems as the splitting can reach intervals of lenth around machine precision.

	if n == 1 && coeffs[0] != 0 {
		return []float64{}
	}
	if n == 1 {
		panic("Finding roots of the 0 polynomial")
	}
	if n == 2 {
		root := -coeffs[0] / coeffs[1]
		if math.Abs(root) <= 1+tol {
			return []float64{(1+root)*(poly.DomainUpper-poly.DomainLower)/2 + poly.DomainLower}
		}
		return []float64{}
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
		hess(a, n-1)
		eigenValues := eig(a, n-1)
		for _, v := range eigenValues {
			if math.Abs(imag(v)) < tol && math.Abs(real(v)) <= 1+tol {
				scaledV := (1+real(v))*(poly.DomainUpper-poly.DomainLower)/2 + poly.DomainLower
				roots = append(roots, scaledV)
			}
		}
		return roots
	}
	//Split
	const splitPoint = -0.004849834917525 // This is an arbitrary number copied from Chefun.

	scaledSplitPoint := (1+splitPoint)*(poly.DomainUpper-poly.DomainLower)/2 + poly.DomainLower

	sections := Split(poly, scaledSplitPoint)

	for i := range sections {
		maxCoeff := 0.0
		for _, v := range sections[i].Coeffs {
			if math.Abs(v) > maxCoeff {
				maxCoeff = math.Abs(v)
			}
		}
		scale := 1 / maxCoeff
		for j := range sections[i].Coeffs {
			sections[i].Coeffs[j] *= scale
		}
	}

	return append(sections[0].getApproximateRoots(10*tol), sections[1].getApproximateRoots(10*tol)...)
}

func (poly Chebpoly) refineRoots(roots []float64) []float64 {
	const maxIterations = 20

	derivative := Diff(poly)

	for i := len(roots) - 1; i >= 0; i-- {
		v := roots[i]
		for i := 0; i < maxIterations; i++ {
			if Evaluate(derivative, v) != 0 {
				change := Evaluate(poly, v) / Evaluate(derivative, v)
				if math.Abs(change) > 1e-15 {
					v = v - change
				} else {
					break
				}
			}
		}
		if math.Abs(Evaluate(poly, v)) < 1e-12+1e-16*float64(poly.Length()) {
			roots[i] = v
		} else {
			roots[i] = roots[len(roots)-1]
			roots = roots[:len(roots)-1]
		}

	}
	return roots
}

func (poly Chebpoly) evaluateCmplx(x complex128) complex128 {
	var bk2 complex128
	var bk1 complex128
	complexCoeffs := make([]complex128, len(poly.Coeffs))
	for i := range complexCoeffs {
		complexCoeffs[i] = complex(poly.Coeffs[i], 0)
	}
	for i := poly.Length() - 1; i > 1; i -= 2 {
		bk2 = complexCoeffs[i] + 2*x*bk1 - bk2
		bk1 = complexCoeffs[i-1] + 2*x*bk2 - bk1
	}
	if ((poly.Length() - 1) & 1) != 0 {
		//We are of odd degree so bk2 is b3 and bk1 is b2
		bk2 = complexCoeffs[1] + 2*x*bk1 - bk2
		return complexCoeffs[0] + x*bk2 - bk1
	}
	//We are of even degree so bk2 is b2 and bk1 is b1
	return complexCoeffs[0] + x*bk1 - bk2
}

//Extremum is a local maximum or minimum of the Chebyshev approximation. The end points are considered to be extrema.
//Point contains the x value of the extremum.
//Value contains the value of the Chebyshev approximation at the extremum.
//Maximum is true if the extremum is a local maximum.
type Extremum struct {
	Point float64
	Value float64

	Maximum bool
}

type byPoint []Extremum

func (a byPoint) Len() int           { return len(a) }
func (a byPoint) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a byPoint) Less(i, j int) bool { return a[i].Point < a[j].Point }

//Extrema returns the local maxmima and minima of the Chebyshev approximation. The end points are considered to be extrema.
func Extrema(poly Chebpoly) []Extremum {
	const realTolerance = 1e-13

	n := poly.Length()
	if n == 1 {
		//TODO What should I actually do here?
		//I'm just going to error because handling these stupid cases is annoying.
		panic("Uncountably many extrema as function is constant")
	}
	diff := Diff(poly)
	criticalPoints := Roots(diff)
	extrema := make([]Extremum, 0, len(criticalPoints)+2)

	//End points
	if w := Evaluate(diff, poly.DomainLower); w > 0 { //Note that the criticalPoints are sorted
		e := Extremum{Point: poly.DomainLower, Value: Evaluate(poly, poly.DomainLower), Maximum: false}
		extrema = append(extrema, e)
	} else if w < 0 {
		e := Extremum{Point: poly.DomainLower, Value: Evaluate(poly, poly.DomainLower), Maximum: true}
		extrema = append(extrema, e)
	}

	if w := Evaluate(diff, poly.DomainUpper); w > 0 { //Note that the criticalPoints are sorted
		e := Extremum{Point: poly.DomainUpper, Value: Evaluate(poly, poly.DomainUpper), Maximum: true}
		extrema = append(extrema, e)
	} else if w < 0 {
		e := Extremum{Point: poly.DomainUpper, Value: Evaluate(poly, poly.DomainUpper), Maximum: false}
		extrema = append(extrema, e)
	}

	//Remove duplicates
	for i := len(criticalPoints) - 1; i > 0; i-- {
		if math.Abs(criticalPoints[i]-criticalPoints[i-1]) < realTolerance {
			criticalPoints[i] = criticalPoints[len(criticalPoints)-1]
			criticalPoints = criticalPoints[:len(criticalPoints)-1]
		}
	}

	//This will classify the criticalPoints using the Higher Order Derivative Test. Note that we don't care about points of inflection so we don't test for these. They will be added to extrema as they are classified and sorted at the end.
	d := diff
	derivativeNumber := 1
	for len(criticalPoints) > 0 {
		d = Diff(d) //Compute the next deritive
		derivativeNumber++
		for i := len(criticalPoints) - 1; i >= 0; i-- {
			if w := Evaluate(d, criticalPoints[i]); w > 0 && (derivativeNumber%2) == 0 {
				//Local Minimum
				e := Extremum{Point: criticalPoints[i], Value: Evaluate(poly, criticalPoints[i]), Maximum: false}
				extrema = append(extrema, e)
			} else if w < 0 && (derivativeNumber%2) == 0 {
				//Local Maximum
				e := Extremum{Point: criticalPoints[i], Value: Evaluate(poly, criticalPoints[i]), Maximum: true}
				extrema = append(extrema, e)
			} else if w != 0 {
				//We have classified the point as either a local minimum, a local maximum or a point of inflection and dealt with it so we can remove it.
				criticalPoints[i] = criticalPoints[len(criticalPoints)-1]
				criticalPoints = criticalPoints[:len(criticalPoints)-1]
			}
		}
	}

	sort.Sort(byPoint(extrema))

	return extrema
}

//MaxAndMin returns the maximum and the minimum of the poly.
func MaxAndMin(poly Chebpoly) (float64, float64) {
	max := Evaluate(poly, poly.DomainLower)
	min := max
	if poly.Length() == 1 {
		return max, min
	}
	if v := Evaluate(poly, poly.DomainUpper); v > max {
		max = v
	} else if v < min {
		min = v
	}
	criticalPoints := Roots(Diff(poly))
	for _, w := range criticalPoints {
		if v := Evaluate(poly, w); v > max {
			max = v
		} else if v < min {
			min = v
		}
	}
	return max, min
}
