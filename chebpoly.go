package chebpoly

import "math"
import "github.com/mjibson/go-dsp/fft"
import "github.com/gonum/matrix/mat64"
import "fmt"
import "sort"

//Chebpts(n, dl, du) gives the n CHebyshev points of the second kind suitably scaled to the domain [dl,du]
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
	for i, _ := range points {
		points[i] = (1+math.Cos(math.Pi*float64(N-1-i)/float64(N-1)))*(domainUpper-domainLower)/2 + domainLower
	}

	return points
}

type chebpoly struct {
	domainLower float64
	domainUpper float64

	coeffs []float64 //coeffs[k] is the coefficient of kth Chebyshev polynomial of the fist kind
}

func AddChebpoly(a, b chebpoly) chebpoly {
	if a.domainLower != b.domainLower || a.domainUpper != b.domainUpper {
		panic("The chebpolys have different domains")
	}
	if a.Length() < b.Length() {
		return AddChebpoly(b, a)
	}
	sum := a.Copy()
	for i := 0; i < b.Length(); i++ {
		sum.coeffs[i] += b.coeffs[i]
	}
	return sum
}

func MultiplyChebpoly(a, b chebpoly) chebpoly {
	if a.domainLower != b.domainLower || a.domainUpper != b.domainUpper {
		panic("The chebpolys have different domains")
	}
	coeffs := make([]float64, a.Length()+b.Length()-1)
	for i := 0; i < a.Length(); i++ {
		for j := 0; j < b.Length(); j++ {
			coeffs[i+j] += 0.5 * a.coeffs[i] * b.coeffs[j]
			if i > j {
				coeffs[i-j] += 0.5 * a.coeffs[i] * b.coeffs[j]
			} else {
				coeffs[j-i] += 0.5 * a.coeffs[i] * b.coeffs[j]
			}
		}
	}
	return chebpoly{domainLower: a.domainLower, domainUpper: b.domainUpper, coeffs: coeffs}
}

func (poly chebpoly) String() string {
	if poly.Length() < 10 {
		return fmt.Sprintf("Interval: [%f %f]  Degree: %d  coeffs: %v", poly.domainLower, poly.domainUpper, poly.Length()-1, poly.coeffs)
	} else {
		return fmt.Sprintf("Interval: [%f %f]  Degree: %d  coeffs: [%f %f %f %f %f %f %f %f %f %f ...]", poly.domainLower, poly.domainUpper, poly.Length()-1, poly.coeffs[0], poly.coeffs[1], poly.coeffs[2], poly.coeffs[3], poly.coeffs[4], poly.coeffs[5], poly.coeffs[6], poly.coeffs[7], poly.coeffs[8], poly.coeffs[9])
	}
}

func (poly chebpoly) Copy() chebpoly {
	copyCoeffs := make([]float64, len(poly.coeffs))
	copy(copyCoeffs, poly.coeffs)
	return chebpoly{domainLower: poly.domainLower, domainUpper: poly.domainUpper, coeffs: copyCoeffs}
}

func (poly chebpoly) Length() int {
	return len(poly.coeffs)
}

func New(coeffs []float64, domainLower float64, domainUpper float64) chebpoly {
	if domainLower > domainUpper {
		domainLower, domainUpper = domainUpper, domainLower
	}

	if len(coeffs) == 0 {
		return chebpoly{domainLower: domainLower, domainUpper: domainUpper, coeffs: []float64{0}}
	} else {
		return chebpoly{domainLower: domainLower, domainUpper: domainUpper, coeffs: coeffs}
	}
}

//Interp(points) returns the Chebyshev interpolate taking value points[i] at the ith Chebyshev point of the second kind.
func Interp(values []float64, domainLower float64, domainUpper float64) chebpoly {
	n := len(values)
	poly := new(chebpoly)

	if domainLower > domainUpper {
		domainLower, domainUpper = domainUpper, domainLower
	}
	poly.domainLower = domainLower
	poly.domainUpper = domainUpper

	if n == 1 {
		poly.coeffs = []float64{values[0]}
		return *poly
	}

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

	poly.coeffs = expandedData[0:n]

	return *poly
}

//Evaluate(x) returns the value of the Chebypoly poly at the point x. It uses Clenshaw's algorithm.
func (poly chebpoly) Evaluate(x float64) float64 {
	scaledX := (x-poly.domainLower)/(poly.domainUpper-poly.domainLower)*2 - 1
	bk2 := 0.0
	bk1 := 0.0
	for i := poly.Length() - 1; i > 1; i -= 2 {
		bk2 = poly.coeffs[i] + 2*scaledX*bk1 - bk2
		bk1 = poly.coeffs[i-1] + 2*scaledX*bk2 - bk1
	}
	if ((poly.Length() - 1) & 1) != 0 {
		//We are of odd degree so bk2 is b3 and bk1 is b2
		bk2 = poly.coeffs[1] + 2*scaledX*bk1 - bk2
		return poly.coeffs[0] + scaledX*bk2 - bk1
	} else {
		//We are of even degree so bk2 is b2 and bk1 is b1
		return poly.coeffs[0] + scaledX*bk1 - bk2
	}
}

//Cumsum returns the chebpoly of the indefinite integral.
//The constant of integration is taken so that Cumsum(domainLower) = 0
func (poly chebpoly) Cumsum() chebpoly {
	n := poly.Length()
	integral := chebpoly{domainLower: poly.domainLower, domainUpper: poly.domainUpper}
	integral.coeffs = make([]float64, poly.Length()+1)
	scaleFactor := (poly.domainUpper - poly.domainLower) / 2 //Multiply all the coefficients by this factor to account for the change of variables in the integral.
	b0 := 0.0                                                //This will be the zero coefficient of the integrand
	if n > 2 {
		integral.coeffs[1] = scaleFactor * (poly.coeffs[0] - poly.coeffs[2]/2)
		b0 += integral.coeffs[1]
		for i := 2; i < n-1; i++ {
			integral.coeffs[i] = scaleFactor * (poly.coeffs[i-1] - poly.coeffs[i+1]) / (2 * float64(i))
			if (i % 2) == 0 {
				b0 -= integral.coeffs[i]
			} else {
				b0 += integral.coeffs[i]
			}
		}
		integral.coeffs[n-1] = scaleFactor * poly.coeffs[n-2] / (2 * float64(n-1))
		integral.coeffs[n] = scaleFactor * poly.coeffs[n-1] / (2 * float64(n))
	} else if n == 2 {
		integral.coeffs[1] = scaleFactor * poly.coeffs[0]
		integral.coeffs[2] = scaleFactor * poly.coeffs[1] / 4
	} else if n == 1 {
		integral.coeffs[1] = scaleFactor * poly.coeffs[0]
	}
	if (n % 2) == 1 {
		b0 -= integral.coeffs[n-1]
		b0 += integral.coeffs[n]
	} else {
		b0 += integral.coeffs[n-1]
		b0 -= integral.coeffs[n]
	}
	integral.coeffs[0] = b0
	return integral
}

//Sum returns the integral of the chebpoly over the domain.
func (poly chebpoly) Sum() float64 {
	sum := 0.0
	for i := 0; i < poly.Length(); i += 2 {
		sum += -poly.coeffs[i] * 2 / float64(i*i-1)
	}
	return sum * (poly.domainUpper - poly.domainLower) / 2
}

//Diff returns the chebpoly representing the derivative of the given chebpoly
func (poly chebpoly) Diff() chebpoly {
	n := poly.Length()
	derivative := chebpoly{domainLower: poly.domainLower, domainUpper: poly.domainUpper}
	if n == 1 {
		//poly is a constant so derivative is 0.
		derivative.coeffs = make([]float64, 1)
		return derivative
	}
	coeffs := make([]float64, n+1) //We will add a couple of extra zero entries here to make the code simpler. Note that these entries will probably always be in memory but this is not really a problem.
	scaleFactor := 2 / (poly.domainUpper - poly.domainLower)
	for i := n - 2; i > 0; i-- {
		coeffs[i] = coeffs[i+2] + scaleFactor*2*float64(i+1)*poly.coeffs[i+1]
	}
	coeffs[0] = coeffs[2]/2 + scaleFactor*poly.coeffs[1]
	derivative.coeffs = coeffs[:n-1]
	return derivative
}

const complexTolerance = 1e-14
const realTolerance = 1e-14

//Roots returns all the real roots of the chebpoly in [-1,1]
//This could be extended to provide all the roots of the chebpoly but it isn't.
func (poly chebpoly) Roots() []float64 {

	roots := poly.getApproximateRoots()
	roots = poly.refineRoots(roots)
	sort.Float64s(roots)
	return roots
}

//This finds the roots in the [-1,1] but they may not be perfectly accurate.
func (poly chebpoly) getApproximateRoots() []float64 {
	n := poly.Length()
	//Truncate the series. We are only after the roots approximately anyway.

	//Calculate the 1 norm
	norm := 0.0
	for _, v := range poly.coeffs {
		norm += math.Abs(v)
	}
	cutoffValue := 1e-13 * norm
	cutOffPoint := 1
	//Iterate back through the coefficients to find the first point that we
	for i := 0; i < n-1; i++ {
		if math.Abs(poly.coeffs[n-1-i]) > cutoffValue {
			cutOffPoint = n - i
			break
		}
	}

	coeffs := poly.coeffs[:cutOffPoint]
	n = len(coeffs)
	//We will find the roots by solving the Colleague matrix eigenvalue problem if it is at most maxEigProbSize else we will split the interval in 2 and try again.

	const maxEigProbSize = 50 //Note that making this too small can cause problems as the splitting can reach intervals of lenth around machine precision.

	if n == 1 && coeffs[0] != 0 {
		return []float64{}
	}
	if n == 1 {
		// We will return the middle of the domain like chebfun...
		panic("Finding roots of the 0 polynomial")
	}
	if n == 2 {
		root := -coeffs[0] / coeffs[1]
		if math.Abs(root) <= 1+realTolerance {
			return []float64{(1+root)*(poly.domainUpper-poly.domainLower)/2 + poly.domainLower}
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
				scaledV := (1+real(v))*(poly.domainUpper-poly.domainLower)/2 + poly.domainLower
				roots = append(roots, scaledV)
			}
		}
		return roots
	} else {
		//Split
		const splitPoint = -0.004849834917525 // This is an arbitrary number copied from Chefun.

		scaledSplitPoint := (1+splitPoint)*(poly.domainUpper-poly.domainLower)/2 + poly.domainLower

		sections := poly.Split(scaledSplitPoint)

		return append(sections[0].getApproximateRoots(), sections[1].getApproximateRoots()...)
	}
}

func (poly chebpoly) refineRoots(roots []float64) []float64 {
	const maxIterations = 10

	derivative := poly.Diff()

	for i, v := range roots {
		for i := 0; i < maxIterations; i++ {
			if derivative.Evaluate(v) != 0 {
				change := poly.Evaluate(v) / derivative.Evaluate(v)
				if math.Abs(change) > 1e-15 {
					v = v - change
				} else {
					break
				}
			}
		}
		roots[i] = v
	}
	return roots
}

type extremum struct {
	Point float64
	Value float64

	Maximum bool
}

type byPoint []extremum

func (a byPoint) Len() int           { return len(a) }
func (a byPoint) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a byPoint) Less(i, j int) bool { return a[i].Point < a[j].Point }

func (poly chebpoly) Extrema() []extremum {
	n := poly.Length()
	if n == 1 {
		//TODO What should I actually do here?
		//I'm just going to error because handling these stupid cases is annoying.
		panic("Uncountably many extrema as function is constant")
	}
	diff := poly.Diff()
	criticalPoints := diff.Roots()
	extrema := make([]extremum, 0, len(criticalPoints)+2)

	//End points
	if w := diff.Evaluate(poly.domainLower); w > 0 { //Note that the criticalPoints are sorted
		e := extremum{Point: poly.domainLower, Value: poly.Evaluate(poly.domainLower), Maximum: false}
		extrema = append(extrema, e)
	} else if w < 0 {
		e := extremum{Point: poly.domainLower, Value: poly.Evaluate(poly.domainLower), Maximum: true}
		extrema = append(extrema, e)
	}

	if w := diff.Evaluate(poly.domainUpper); w > 0 { //Note that the criticalPoints are sorted
		e := extremum{Point: poly.domainUpper, Value: poly.Evaluate(poly.domainUpper), Maximum: true}
		extrema = append(extrema, e)
	} else if w < 0 {
		e := extremum{Point: poly.domainUpper, Value: poly.Evaluate(poly.domainUpper), Maximum: false}
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
		d = d.Diff() //Compute the next deritive
		derivativeNumber++
		for i := len(criticalPoints) - 1; i >= 0; i-- {
			if w := d.Evaluate(criticalPoints[i]); w > 0 && (derivativeNumber%2) == 0 {
				//Local Minimum
				e := extremum{Point: criticalPoints[i], Value: poly.Evaluate(criticalPoints[i]), Maximum: false}
				extrema = append(extrema, e)
			} else if w < 0 && (derivativeNumber%2) == 0 {
				//Local Maximum
				e := extremum{Point: criticalPoints[i], Value: poly.Evaluate(criticalPoints[i]), Maximum: true}
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

func (poly chebpoly) MaxAndMin() (float64, float64) {
	max := poly.Evaluate(poly.domainLower)
	min := max
	if poly.Length() == 1 {
		return max, min
	}
	if v := poly.Evaluate(poly.domainUpper); v > max {
		max = v
	} else if v < min {
		min = v
	}
	criticalPoints := poly.Diff().Roots()
	for _, w := range criticalPoints {
		if v := poly.Evaluate(w); v > max {
			max = v
		} else if v < min {
			min = v
		}
	}
	return max, min
}

func Adaptive(f func(float64) float64, domainLower float64, domainUpper float64) chebpoly {
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
			if !imprecise && math.Abs(currentPoly.Evaluate(points[2*i+1])-currentValues[2*i+1]) > 1e-13*maxValue {
				imprecise = true
			}
		}

		if imprecise {
			//We aren't accurate on the new points so improve
			currentPoly = Interp(currentValues, domainLower, domainUpper)
			currentDegree = 2 * currentDegree
		} else {
			break
		}
	}

	return currentPoly
}

func (poly chebpoly) Split(splitPoints ...float64) []chebpoly {
	n := poly.Length()
	sections := make([]chebpoly, len(splitPoints)+1)
	values := make([]float64, n)
	points := Chebpts(uint(n), 0, 1) //Note that we manipulate the same set of points instead of calling cos multiple times.
	scaleFactor := splitPoints[0] - poly.domainLower
	for i, v := range points {
		values[i] = poly.Evaluate(v*scaleFactor + poly.domainLower)
	}
	sections[0] = Interp(values, poly.domainLower, splitPoints[0])

	for i := 1; i < len(splitPoints); i++ {
		scaleFactor = (splitPoints[i] - splitPoints[i-1])
		for j, v := range points {
			values[j] = poly.Evaluate(v*scaleFactor + splitPoints[i-1])
		}
		sections[i] = Interp(values, splitPoints[i-1], splitPoints[i])
	}

	scaleFactor = poly.domainUpper - splitPoints[len(splitPoints)-1]
	for i, v := range points {
		values[i] = poly.Evaluate(v*scaleFactor + splitPoints[len(splitPoints)-1])
	}
	sections[len(splitPoints)] = Interp(values, splitPoints[len(splitPoints)-1], poly.domainUpper)
	return sections
}

func (poly chebpoly) Abs() chebfun {
	roots := poly.Roots()
	//Remove duplicates
	for i := len(roots) - 1; i > 0; i-- {
		if math.Abs(roots[i]-roots[i-1]) < realTolerance {
			roots = append(roots[:i], roots[i+1:]...)
		}
	}
	fun := poly.Split(roots...)
	for i, _ := range fun {
		v := fun[i]
		if v.Evaluate((v.domainLower+v.domainUpper)/2) < 0 {
			for j, u := range v.coeffs {
				v.coeffs[j] = -u
			}
		}
	}
	return fun
}

//Note that chebfuns are right-continuous and most be sorted with no overlap.
type chebfun []chebpoly

func (fun chebfun) Copy() chebfun {
	copy := make([]chebpoly, len(fun))
	for i := 0; i < len(fun); i++ {
		copy[i] = fun[i].Copy()
	}
	return copy
}

func (a chebfun) Len() int           { return len(a) }
func (a chebfun) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a chebfun) Less(i, j int) bool { return a[i].domainLower < a[j].domainLower }

func NewChebfun(parts []chebpoly) (chebfun, bool) {
	fun := chebfun(parts).Copy()
	sort.Sort(fun)
	hscale := fun[len(fun)-1].domainUpper - fun[0].domainLower
	for i := 0; i < len(parts)-1; i++ {
		if fun[i].domainUpper < fun[i+1].domainLower-realTolerance*hscale {
			fun = append(fun, chebpoly{domainLower: fun[i].domainUpper, domainUpper: fun[i+1].domainLower, coeffs: []float64{0}})
		} else if fun[i].domainUpper > fun[i+1].domainLower+realTolerance*hscale {
			return nil, false
		} else {
			v := (fun[i+1].domainLower + fun[i].domainUpper) / 2
			fun[i+1].domainLower = v
			fun[i].domainUpper = v
		}
	}
	sort.Sort(fun)
	return fun, true
}

func alignChebfuns(funA, funB chebfun) (chebfun, chebfun) {
	a := funA.Copy()
	b := funB.Copy()

	hscale := math.Max(a[len(a)-1].domainUpper-a[0].domainLower, b[len(b)-1].domainUpper-b[0].domainLower)

	aSplitPoints := make([]float64, len(a)+1)
	aSplitPoints[0] = a[0].domainLower
	bSplitPoints := make([]float64, len(b)+1)
	bSplitPoints[0] = b[0].domainLower

	for i := 0; i < len(a); i++ {
		aSplitPoints[i+1] = a[i].domainUpper
	}
	for i := 0; i < len(b); i++ {
		bSplitPoints[i+1] = b[i].domainUpper
	}

	bIndex := 0
	for i := 0; i < len(aSplitPoints); i++ {
		if (i == len(aSplitPoints)-1 || aSplitPoints[i+1]-aSplitPoints[i] > 2*hscale*realTolerance) && (i == 0 || aSplitPoints[i]-aSplitPoints[i-1] > 2*hscale*realTolerance) {
			//We can move this splitPoint
			for bIndex < len(bSplitPoints) {
				if bSplitPoints[bIndex] < aSplitPoints[i]-hscale*realTolerance {
					bIndex++
				} else if bSplitPoints[bIndex] > aSplitPoints[i]+hscale*realTolerance {
					break
				} else if (bIndex == len(bSplitPoints)-1 || aSplitPoints[i+1]-aSplitPoints[i] > 2*hscale*realTolerance) && (i == 0 || aSplitPoints[i]-aSplitPoints[i-1] > 2*hscale*realTolerance) {
					v := (aSplitPoints[i] + bSplitPoints[bIndex]) / 2
					aSplitPoints[i], bSplitPoints[bIndex] = v, v
					bIndex++
				}
			}
		}
	}

	a[0].domainLower = aSplitPoints[0]
	b[0].domainLower = bSplitPoints[0]

	for i := 0; i < len(a); i++ {
		a[i].domainLower = aSplitPoints[i]
		a[i].domainUpper = aSplitPoints[i+1]
	}
	for i := 0; i < len(b); i++ {
		b[i].domainLower = bSplitPoints[i]
		b[i].domainUpper = bSplitPoints[i+1]
	}

	//Extend a
	if b[0].domainLower < a[0].domainLower {
		a = append(a, chebpoly{})
		copy(a[1:], a[0:])
		a[0] = chebpoly{domainLower: b[0].domainLower, domainUpper: a[0].domainLower, coeffs: []float64{0}}
	}
	if b[len(b)-1].domainUpper > a[len(a)-1].domainUpper {
		a = append(a, chebpoly{domainLower: a[len(a)-1].domainUpper, domainUpper: b[len(b)-1].domainUpper, coeffs: []float64{0}})
	}
	//Extend b
	if a[0].domainLower < b[0].domainLower {
		b = append(b, chebpoly{})
		copy(b[1:], b[0:])
		b[0] = chebpoly{domainLower: a[0].domainLower, domainUpper: b[0].domainLower, coeffs: []float64{0}}
	}
	if a[len(a)-1].domainUpper > b[len(b)-1].domainUpper {
		b = append(b, chebpoly{domainLower: b[len(b)-1].domainUpper, domainUpper: a[len(a)-1].domainUpper, coeffs: []float64{0}})
	}
	//Split a
	var splitPoints []float64
	if len(a) > len(b) {
		splitPoints = make([]float64, 0, len(a)+1)
	} else {
		splitPoints = make([]float64, 0, len(b)+1)
	}

	pointIndex := 0
	i := 0
	for i < len(a) {
		for pointIndex < len(bSplitPoints) {
			if bSplitPoints[pointIndex] < a[i].domainUpper && bSplitPoints[pointIndex] > a[i].domainLower {
				splitPoints = append(splitPoints, bSplitPoints[pointIndex])
				pointIndex++
			} else if bSplitPoints[pointIndex] <= a[i].domainLower {
				pointIndex++
			} else {
				break
			}
		}
		if len(splitPoints) > 0 {
			a = append(a[0:i], append(a[i].Split(splitPoints...), a[i+1:]...)...)
		}
		i += len(splitPoints) + 1
		splitPoints = splitPoints[:0]
	}
	//Split b
	pointIndex = 0
	i = 0
	for i < len(b) {
		for pointIndex < len(aSplitPoints) {
			if aSplitPoints[pointIndex] < b[i].domainUpper && aSplitPoints[pointIndex] > b[i].domainLower {
				splitPoints = append(splitPoints, aSplitPoints[pointIndex])
				pointIndex++
			} else if aSplitPoints[pointIndex] <= b[i].domainLower {
				pointIndex++
			} else {
				break
			}
		}
		if len(splitPoints) > 0 {
			b = append(b[0:i], append(b[i].Split(splitPoints...), b[i+1:]...)...)
		}
		i += len(splitPoints) + 1
		splitPoints = splitPoints[:0]
	}
	return a, b

}

func Add(a, b chebfun) chebfun {
	a2, b2 := alignChebfuns(a, b)
	for i := 0; i < len(a2); i++ {
		a2[i] = AddChebpoly(a2[i], b2[i])
	}

	return a2
}

func Multiply(a, b chebfun) chebfun {
	a2, b2 := alignChebfuns(a, b)
	for i := 0; i < len(a2); i++ {
		a2[i] = MultiplyChebpoly(a2[i], b2[i])
	}

	return a2
}

func (fun chebfun) TotalLength() int {
	length := 0
	for _, v := range fun {
		length += v.Length()
	}
	return length
}

func (fun chebfun) Evaluate(x float64) float64 {
	f := func(i int) bool { return fun[i].domainUpper > x }
	index := sort.Search(len(fun), f)
	if index == 0 && fun[0].domainLower > x {
		return 0
	} else if index == len(fun) {
		return 0
	}
	return fun[index].Evaluate(x)
}

func (fun chebfun) Cumsum() chebfun {
	newPieces := make([]chebpoly, len(fun))

	sum := 0.0 //This will be added to
	for i, v := range fun {
		newPieces[i] = v.Cumsum()
		newPieces[i].coeffs[0] += sum
		sum += v.Sum()
	}
	cumsum := newPieces
	return cumsum
}

func (fun chebfun) Sum() float64 {
	sum := 0.0
	for _, v := range fun {
		sum += v.Sum()
	}
	return sum
}

func (fun chebfun) MaxAndMin() (float64, float64) {
	max, min := fun[0].MaxAndMin()
	tmpmax := 0.0
	tmpmin := 0.0
	for i := 1; i < len(fun); i++ {
		tmpmax, tmpmin = fun[i].MaxAndMin()
		if tmpmax > max {
			max = tmpmax
		}
		if tmpmin < min {
			min = tmpmin
		}
	}
	return max, min
}

//chebfun says there is a root when there a discontinuity is across 0. Do I want to follow suit?
// func (fun chebfun) Roots() []float64{
// 	roots := make([]float64, 0, fun.TotalLength());
// 	roots = append(roots, fun)
// 	for i := 1 ; i < len(fun); i++{
//
// 	}
// }
