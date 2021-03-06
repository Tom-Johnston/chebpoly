package chebpoly

import (
	"math"
	"testing"
)

const tol = 1e-13

//The function f(x) = x on [-0.5, 2].
var x Chebpoly = Chebpoly{Coeffs: []float64{0.75, 1.25}, DomainLower: -0.5, DomainUpper: 2}

//The polynomial T_4(x) = 8x^4 - 8x^2 + 1 on the standard domain of [-1,1].
var t4 Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 0, 1}, DomainLower: -1, DomainUpper: 1}

//The indefinite integral of the polynomial T_4(x) with the constant chosen so that it is 0 at -1.
var t4Cumsum Chebpoly = Chebpoly{Coeffs: []float64{-1.0 / 15, 0, 0, -1.0 / 6, 0, 1.0 / 10}, DomainLower: -1, DomainUpper: 1}

//The polynomial
var t4Scaled Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 0, 1}, DomainLower: 0.5, DomainUpper: 2}

var t5 Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 0, 0, 1}, DomainLower: -1, DomainUpper: 1}

var t5Scaled Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 0, 0, 1}, DomainLower: 0.5, DomainUpper: 2}
var t5ScaledCumsum Chebpoly = Chebpoly{Coeffs: []float64{1.0 / 32, 0, 0, 0, -3.0 / 32, 0, 1.0 / 16}, DomainLower: 0.5, DomainUpper: 2}

var t3plust5 Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 1, 0, 1}, DomainLower: -1, DomainUpper: 1}

var extendedt3plust5 Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 1, 0, 1, 0, 0, 0, 0}, DomainLower: -1, DomainUpper: 1}

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
	if math.Abs(found.DomainLower-expected.DomainLower) > tol || math.Abs(found.DomainUpper-expected.DomainUpper) > tol {
		t.Errorf("Domains don't match - Found: [%f %f] Expected: [%f %f]", found.DomainLower, found.DomainUpper, expected.DomainLower, expected.DomainUpper)
	}
	if len(found.Coeffs) != len(expected.Coeffs) {
		t.Errorf("Degrees don't match - Found: %d Expected: %d", len(found.Coeffs)-1, len(expected.Coeffs)-1)
	}
	maxDifference := 0.0
	maxPosition := 0
	for i, v := range found.Coeffs {
		if math.Abs(v-expected.Coeffs[i]) > maxDifference {
			maxDifference = math.Abs(v - expected.Coeffs[i])
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
	floatsNearlyEqual(Evaluate(t5, Chebpts(6, -1, 1)[2]), -1.0, t)

	floatsNearlyEqual(Evaluate(t4, Chebpts(5, -1, 1)[2]), 1.0, t)

	floatsNearlyEqual(Evaluate(t5Scaled, Chebpts(6, 0.5, 2)[2]), -1.0, t)
}

func TestCumsum(t *testing.T) {

	chebpolysNearlyEqual(Cumsum(t4), t4Cumsum, t)

	chebpolysNearlyEqual(Cumsum(t5Scaled), t5ScaledCumsum, t)
}

func TestSum(t *testing.T) {

	floatsNearlyEqual(Sum(t4), -2.0/15, t)

	floatsNearlyEqual(Sum(t5), 0.0, t)

	floatsNearlyEqual(Sum(t4Scaled), -1.0/10, t)

}

func TestDiff(t *testing.T) {

	chebpolysNearlyEqual(Diff(t4Cumsum), t4, t)

	chebpolysNearlyEqual(Diff(t5ScaledCumsum), t5Scaled, t)

}

func TestRoots(t *testing.T) {

	poly := Chebpoly{DomainLower: -1, DomainUpper: 1, Coeffs: []float64{2.4}}
	correctResult := []float64{}
	floatArraysNearlyEqual(Roots(poly), correctResult, t)

	poly = Chebpoly{DomainLower: -1, DomainUpper: 1, Coeffs: []float64{0, 1}}
	correctResult = []float64{0}
	floatArraysNearlyEqual(Roots(poly), correctResult, t)

	poly = Chebpoly{DomainLower: -1, DomainUpper: 1, Coeffs: []float64{2, 1}}
	correctResult = []float64{}
	floatArraysNearlyEqual(Roots(poly), correctResult, t)

	correctResult = []float64{math.Cos(4*math.Pi/5 + math.Pi/10), math.Cos(3*math.Pi/5 + math.Pi/10), 0, math.Cos(1*math.Pi/5 + math.Pi/10), math.Cos(math.Pi / 10)}
	floatArraysNearlyEqual(Roots(t5), correctResult, t)

	correctResult = []float64{0.536707612778635, 0.809161060780646, 1.25, 1.690838939219355, 1.963292387221366}
	floatArraysNearlyEqual(Roots(t5Scaled), correctResult, t)

	correctResult = make([]float64, 201)
	for i := range correctResult {
		correctResult[i] = -1.0 + float64(i)*(0.01)
	}
	s100 := Adaptive(sin100, -1, 1)
	floatArraysNearlyEqual(Roots(s100), correctResult, t)

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

	extremaNearlyEqual(Extrema(t5), correctResult, t)
}

func TestMaxAndMin(t *testing.T) {
	max, min := MaxAndMin(t5)
	floatsNearlyEqual(max, 1, t)
	floatsNearlyEqual(min, -1, t)

	max, min = MaxAndMin(t5Scaled)
	floatsNearlyEqual(max, 1, t)
	floatsNearlyEqual(min, -1, t)

	s := Adaptive(xsinx, -1, 1)
	max, min = MaxAndMin(s)
	floatsNearlyEqual(max, 0.57923032735651960910655, t)
	floatsNearlyEqual(min, 0, t)
}

func cos5(x float64) float64 {
	return math.Cos(5 * x)
}

func sin100(x float64) float64 {
	return math.Sin(100 * x * math.Pi)
}

func xsinx(x float64) float64 {
	return x * math.Sin(x*math.Pi)
}

func TestAdaptive(t *testing.T) {
	poly := Adaptive(cos5, -1.0, 1.0)
	correctValue := math.Sin(5) / 2.5
	floatsNearlyEqual(Sum(poly), correctValue, t)

	correctExtrema := []Extremum{
		Extremum{Point: -1, Value: cos5(-1), Maximum: true},
		Extremum{Point: -math.Pi / 5, Value: -1, Maximum: false},
		Extremum{Point: 0, Value: 1, Maximum: true},
		Extremum{Point: math.Pi / 5, Value: -1, Maximum: false},
		Extremum{Point: 1, Value: cos5(1), Maximum: true},
	}
	extremaNearlyEqual(Extrema(poly), correctExtrema, t)

	max, min := MaxAndMin(poly)
	floatsNearlyEqual(max, 1, t)
	floatsNearlyEqual(min, -1, t)

	correctRoots := []float64{
		-3 * math.Pi / 10,
		-math.Pi / 10,
		math.Pi / 10,
		3 * math.Pi / 10,
	}
	t.Log(Roots(poly))
	floatArraysNearlyEqual(Roots(poly), correctRoots, t)
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
		Evaluate(s100, 0.6)
	}
}

func BenchmarkCumsum(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Cumsum(s100)
	}
}

func BenchmarkSum(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Sum(s100)
	}
}

func BenchmarkDiff(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Diff(s100)
	}
}

func BenchmarkRoots(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Roots(s100)
	}
}

func BenchmarkExtrema(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Extrema(s100)
	}
}

func BenchmarkMaxAndMin(b *testing.B) {
	s100 := Adaptive(sin100, -1, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		MaxAndMin(s100)
	}
}
