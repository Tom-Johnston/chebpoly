package chebpoly

import (
	"math"
	"testing"
)

const tol = 1e-13

//The function f(x) = x on [-0.5, 2].
var x Chebpoly = Chebpoly{coeffs: []float64{0.75, 1.25}, domainLower: -0.5, domainUpper: 2}

//The polynomial T_4(x) = 8x^4 - 8x^2 + 1 on the standard domain of [-1,1].
var t4 Chebpoly = Chebpoly{coeffs: []float64{0, 0, 0, 0, 1}, domainLower: -1, domainUpper: 1}

//The indefinite integral of the polynomial T_4(x) with the constant chosen so that it is 0 at -1.
var t4Cumsum Chebpoly = Chebpoly{coeffs: []float64{-1.0 / 15, 0, 0, -1.0 / 6, 0, 1.0 / 10}, domainLower: -1, domainUpper: 1}

//The polynomial
var t4Scaled Chebpoly = Chebpoly{coeffs: []float64{0, 0, 0, 0, 1}, domainLower: 0.5, domainUpper: 2}

var t5 Chebpoly = Chebpoly{coeffs: []float64{0, 0, 0, 0, 0, 1}, domainLower: -1, domainUpper: 1}

var t5Scaled Chebpoly = Chebpoly{coeffs: []float64{0, 0, 0, 0, 0, 1}, domainLower: 0.5, domainUpper: 2}
var t5ScaledCumsum Chebpoly = Chebpoly{coeffs: []float64{1.0 / 32, 0, 0, 0, -3.0 / 32, 0, 1.0 / 16}, domainLower: 0.5, domainUpper: 2}

var t3plust5 Chebpoly = Chebpoly{coeffs: []float64{0, 0, 0, 1, 0, 1}, domainLower: -1, domainUpper: 1}

var extendedt3plust5 Chebpoly = Chebpoly{coeffs: []float64{0, 0, 0, 1, 0, 1, 0, 0, 0, 0}, domainLower: -1, domainUpper: 1}

func floatsNearlyEqual(found, expected float64, t *testing.T) {
	t.Helper()
	difference := math.Abs(found - expected)
	if difference > tol {
		t.Errorf("Fail - Found: %f Expected: %f Diff: %f", found, expected, difference)
	}
}

func floatArraysNearlyEqual(found, expected []float64, t *testing.T) {
	t.Helper()
	if len(found) != len(expected) {
		t.Logf("Length doesn't match - Found: %d Expected: %d", len(found), len(expected))
		t.FailNow()
	}
	maxDifference := 0.0
	maxPosition := 0
	for i, v := range found {
		if math.Abs(v-expected[i]) > maxDifference {
			maxDifference = math.Abs(v - expected[i])
			maxPosition = i
		}
	}
	if maxDifference > tol {
		t.Errorf("Values don't match - Max Diff: %g Position: %d", maxDifference, maxPosition)
	}
}

func extremaNearlyEqual(found, expected []Extremum, t *testing.T) {
	t.Helper()
	if len(found) != len(expected) {
		t.Logf("Length doesn't match - Found: %d Expected: %d", len(found), len(expected))
		t.FailNow()
	}

	maxPointDifference := 0.0
	maxPointPosition := 0

	maxValueDifference := 0.0
	maxValuePosition := 0
	for i, v := range found {
		if math.Abs(v.Point-expected[i].Point) > maxPointDifference {
			maxPointDifference = math.Abs(v.Point - expected[i].Point)
			maxPointPosition = i
		}
		if math.Abs(v.Value-expected[i].Value) > maxValueDifference {
			maxValueDifference = math.Abs(v.Value - expected[i].Value)
			maxValuePosition = i
		}
		if v.Maximum != expected[i].Maximum {
			t.Errorf("Type of extremum disagreement - Found: %v Expected: %v", v.Maximum, expected[i].Maximum)
		}
	}
	if maxPointDifference > tol || maxValueDifference > tol {
		t.Errorf("Values don't match - Max Point Diff: %g at %d Max Value Diff: %g at %d", maxPointDifference, maxPointPosition, maxValueDifference, maxValuePosition)
	}
}

func chebpolysNearlyEqual(found, expected Chebpoly, t *testing.T) {
	t.Helper()
	if math.Abs(found.domainLower-expected.domainLower) > tol || math.Abs(found.domainUpper-expected.domainUpper) > tol {
		t.Errorf("Domains don't match - Found: [%f %f] Expected: [%f %f]", found.domainLower, found.domainUpper, expected.domainLower, expected.domainUpper)
	}
	if len(found.coeffs) != len(expected.coeffs) {
		t.Errorf("Degrees don't match - Found: %d Expected: %d", len(found.coeffs)-1, len(expected.coeffs)-1)
	}
	maxDifference := 0.0
	maxPosition := 0
	for i, v := range found.coeffs {
		if math.Abs(v-expected.coeffs[i]) > maxDifference {
			maxDifference = math.Abs(v - expected.coeffs[i])
			maxPosition = i
		}
	}
	if maxDifference > tol {
		t.Errorf("coeffs don't match - Max Diff: %g Position: %d", maxDifference, maxPosition)
	}
}

func TestChebpts(t *testing.T) {

	output := Chebpts(20, -1.0, 1.0)
	correctResult := []float64{-1.000000000000000, -0.986361303402722, -0.945817241700635, -0.879473751206489, -0.789140509396394, -0.677281571625741, -0.546948158122427, -0.401695424652969, -0.245485487140799, -0.082579345472332, 0.082579345472332, 0.245485487140799, 0.401695424652969, 0.546948158122427, 0.677281571625741, 0.789140509396394, 0.879473751206489, 0.945817241700635, 0.986361303402722, 1.000000000000000}
	floatArraysNearlyEqual(output, correctResult, t)

	output = Chebpts(20, 0.5, 2)
	correctResult = []float64{0.500000000000000, 0.510229022447958, 0.540637068724524, 0.590394686595133, 0.658144617952705, 0.742038821280694, 0.839788881408180, 0.948728431510273, 1.065885884644401, 1.188065490895751, 1.311934509104249, 1.434114115355599, 1.551271568489727, 1.660211118591820, 1.757961178719306, 1.841855382047295, 1.909605313404867, 1.959362931275476, 1.989770977552042, 2.000000000000000}
	floatArraysNearlyEqual(output, correctResult, t)

	output = Chebpts(20, 2, 0.5)
	correctResult = []float64{0.500000000000000, 0.510229022447958, 0.540637068724524, 0.590394686595133, 0.658144617952705, 0.742038821280694, 0.839788881408180, 0.948728431510273, 1.065885884644401, 1.188065490895751, 1.311934509104249, 1.434114115355599, 1.551271568489727, 1.660211118591820, 1.757961178719306, 1.841855382047295, 1.909605313404867, 1.959362931275476, 1.989770977552042, 2.000000000000000}
	floatArraysNearlyEqual(output, correctResult, t)
}

func TestInterp(t *testing.T) {
	output := Interp([]float64{1, -1, 1, -1, 1}, -1, 1)
	chebpolysNearlyEqual(output, t4, t)

	output = Interp([]float64{-2.000000000000000, -0.326351822333070, 1.439692620785908, 0.500000000000000, -0.266044443118978, 0.266044443118978, -0.500000000000000, -1.439692620785908, 0.326351822333070, 2.000000000000000}, -1, 1)
	chebpolysNearlyEqual(output, extendedt3plust5, t)

}

func TestEvaluate(t *testing.T) {
	floatsNearlyEqual(t5.Evaluate(Chebpts(6, -1, 1)[2]), -1.0, t)

	floatsNearlyEqual(t4.Evaluate(Chebpts(5, -1, 1)[2]), 1.0, t)

	floatsNearlyEqual(t5Scaled.Evaluate(Chebpts(6, 0.5, 2)[2]), -1.0, t)
}

func TestCumsum(t *testing.T) {

	chebpolysNearlyEqual(t4.Cumsum(), t4Cumsum, t)

	chebpolysNearlyEqual(t5Scaled.Cumsum(), t5ScaledCumsum, t)
}

func TestSum(t *testing.T) {

	floatsNearlyEqual(t4.Sum(), -2.0/15, t)

	floatsNearlyEqual(t5.Sum(), 0.0, t)

	floatsNearlyEqual(t4Scaled.Sum(), -1.0/10, t)

}

func TestDiff(t *testing.T) {

	chebpolysNearlyEqual(t4Cumsum.Diff(), t4, t)

	chebpolysNearlyEqual(t5ScaledCumsum.Diff(), t5Scaled, t)

}

func TestRoots(t *testing.T) {

	poly := Chebpoly{domainLower: -1, domainUpper: 1, coeffs: []float64{2.4}}
	correctResult := []float64{}
	floatArraysNearlyEqual(poly.Roots(), correctResult, t)

	poly = Chebpoly{domainLower: -1, domainUpper: 1, coeffs: []float64{0, 1}}
	correctResult = []float64{0}
	floatArraysNearlyEqual(poly.Roots(), correctResult, t)

	poly = Chebpoly{domainLower: -1, domainUpper: 1, coeffs: []float64{2, 1}}
	correctResult = []float64{}
	floatArraysNearlyEqual(poly.Roots(), correctResult, t)

	correctResult = []float64{math.Cos(4*math.Pi/5 + math.Pi/10), math.Cos(3*math.Pi/5 + math.Pi/10), 0, math.Cos(1*math.Pi/5 + math.Pi/10), math.Cos(math.Pi / 10)}
	floatArraysNearlyEqual(t5.Roots(), correctResult, t)

	correctResult = []float64{0.536707612778635, 0.809161060780646, 1.25, 1.690838939219355, 1.963292387221366}
	floatArraysNearlyEqual(t5Scaled.Roots(), correctResult, t)

	correctResult = make([]float64, 201)
	for i := range correctResult {
		correctResult[i] = -1.0 + float64(i)*(0.01)
	}
	s100 := Adaptive(sin100, -1, 1)
	floatArraysNearlyEqual(s100.Roots(), correctResult, t)

}

func TestExtrema(t *testing.T) {
	correctResult := []Extremum{
		Extremum{Point: -1, Value: -1, Maximum: false},
		Extremum{Point: math.Cos(math.Pi * 4 / 5), Value: 1, Maximum: true},
		Extremum{Point: math.Cos(math.Pi * 3 / 5), Value: -1, Maximum: false},
		Extremum{Point: math.Cos(math.Pi * 2 / 5), Value: 1, Maximum: true},
		Extremum{Point: math.Cos(math.Pi * 1 / 5), Value: -1, Maximum: false},
		Extremum{Point: 1, Value: 1, Maximum: true},
	}

	extremaNearlyEqual(t5.Extrema(), correctResult, t)
}

func TestMaxAndMin(t *testing.T) {
	max, min := t5.MaxAndMin()
	floatsNearlyEqual(max, 1, t)
	floatsNearlyEqual(min, -1, t)

	max, min = t5Scaled.MaxAndMin()
	floatsNearlyEqual(max, 1, t)
	floatsNearlyEqual(min, -1, t)
}

func cos5(x float64) float64 {
	return math.Cos(5 * x)
}

func sin100(x float64) float64 {
	return math.Sin(100 * x * math.Pi)
}

func TestAdaptive(t *testing.T) {
	poly := Adaptive(cos5, -1.0, 1.0)
	correctValue := math.Sin(5) / 2.5
	floatsNearlyEqual(poly.Sum(), correctValue, t)

	correctExtrema := []Extremum{
		Extremum{Point: -1, Value: cos5(-1), Maximum: true},
		Extremum{Point: -math.Pi / 5, Value: -1, Maximum: false},
		Extremum{Point: 0, Value: 1, Maximum: true},
		Extremum{Point: math.Pi / 5, Value: -1, Maximum: false},
		Extremum{Point: 1, Value: cos5(1), Maximum: true},
	}
	extremaNearlyEqual(poly.Extrema(), correctExtrema, t)

	max, min := poly.MaxAndMin()
	floatsNearlyEqual(max, 1, t)
	floatsNearlyEqual(min, -1, t)

	correctRoots := []float64{
		-3 * math.Pi / 10,
		-math.Pi / 10,
		math.Pi / 10,
		3 * math.Pi / 10,
	}
	t.Log(poly.Roots())
	floatArraysNearlyEqual(poly.Roots(), correctRoots, t)
}

//Benchmarks

var result float64

func BenchmarkChebpts(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Chebpts(100, -1, 1)
	}
}

func BenchmarkInterp(b *testing.B) {
	points := Chebpts(100, -1, 1)
	for i, v := range points {
		points[i] = sin100(v)
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Interp(points, -1, 1)
	}
}

func BenchmarkEvaluate(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		s100.Evaluate(0.6)
	}
}

func BenchmarkCumsum(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		s100.Cumsum()
	}
}

func BenchmarkSum(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		s100.Sum()
	}
}

func BenchmarkDiff(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		s100.Diff()
	}
}

func BenchmarkRoots(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		s100.Roots()
	}
}

func BenchmarkExtrema(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		s100.Extrema()
	}
}

func BenchmarkMaxAndMin(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		s100.MaxAndMin()
	}
}
