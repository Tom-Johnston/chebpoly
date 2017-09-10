package chebpoly

import "testing"
import "math"

const tol = 1e-13

var x chebpoly = chebpoly{coeffs: []float64{0.75, 1.25}, domainLower: -0.5, domainUpper: 2}

var s100 chebpoly = chebpoly{coeffs: []float64{0, -0.154290704028225, 0, -0.152568403440665, 0, -0.148391473929028, 0, -0.140345381974426, 0, -0.126473522812062, 0, -0.104580652037872, 0, -0.072787348681248, 0, -0.030396242447856, 0, 0.020968775379588, 0, 0.076187842328998, 0, 0.125961809127669, 0, 0.157435430279415, 0, 0.157008546711987, 0, 0.115723867368825, 0, 0.036268545968470, 0, -0.061483609528938, 0, -0.141355879131755, 0, -0.162779647030418, 0, -0.103422694466848, 0, 0.017396789940632, 0, 0.133719597657129, 0, 0.163752486893221, 0, 0.071287631631195, 0, -0.084406330587917, 0, -0.170701322789670, 0, -0.093304643332630, 0, 0.086165468950018, 0, 0.173896049004493, 0, 0.050635977430853, 0, -0.143987256125872, 0, -0.141435680910038, 0, 0.075164844458731, 0, 0.177526624584611, 0, -0.021548487864416, 0, -0.187386167733437, 0, 0.003928443986477, 0, 0.192677150775782, 0, -0.029671200913063, 0, -0.190366542693884, 0, 0.102004112571603, 0, 0.144140250398605, 0, -0.196030233129984, 0, 0.003334537327747, 0, 0.194749402093528, 0, -0.198857364307938, 0, 0.021068993112999, 0, 0.202711360106877, 0, -0.364232422013750, 0, 0.441210883358808, 0, -0.230487850646075}, domainLower: -1, domainUpper: 1}

var s100Cumsum chebpoly = chebpoly{coeffs: []float64{0.008423812166895, 0, -0.000430575146890, 0, -0.000522116188955, 0, -0.000670507662883, 0, -0.000866991197648, 0, -0.001094643538709, 0, -0.001324720973193, 0, -0.001513968079764, 0, -0.001605156807108, 0, -0.001533862970817, 0, -0.001244349169967, 0, -0.000715309571631, 0, 0.000008893407655, 0, 0.000793936141215, 0, 0.001418845025006, 0, 0.001629202591623, 0, 0.001248004212544, 0, 0.000315055410274, 0, -0.000824402118938, 0, -0.001589730057993, 0, -0.001454035096456, 0, -0.000357534395668, 0, 0.001050736991614, 0, 0.001692325676295, 0, 0.000898906168768, 0, -0.000773966794570, 0, -0.001725674156564, 0, -0.000812320185690, 0, 0.001100536353336, 0, 0.001677786496179, 0, -0.000021263126799, 0, -0.001746778430393, 0, -0.000799701407233, 0, 0.001508144791281, 0, 0.001219394704919, 0, -0.001366532940857, 0, -0.001310754908259, 0, 0.001502353727627, 0, 0.001057206195926, 0, -0.001874170867086, 0, -0.000263350861419, 0, 0.002074210265418, 0, -0.001186695062248, 0, -0.001112877120731, 0, 0.002236402081827, 0, -0.001221813096783, 0, -0.000987186777141, 0, 0.003015658415535, 0, -0.004195017215482, 0, 0.003427034357168, 0, -0.001152439253230}, domainLower: -1, domainUpper: 1}

var t4 chebpoly = chebpoly{coeffs: []float64{0, 0, 0, 0, 1}, domainLower: -1, domainUpper: 1}
var t4Cumsum chebpoly = chebpoly{coeffs: []float64{-1.0 / 15, 0, 0, -1.0 / 6, 0, 1.0 / 10}, domainLower: -1, domainUpper: 1}

var t4Scaled chebpoly = chebpoly{coeffs: []float64{0, 0, 0, 0, 1}, domainLower: 0.5, domainUpper: 2}

var t5 chebpoly = chebpoly{coeffs: []float64{0, 0, 0, 0, 0, 1}, domainLower: -1, domainUpper: 1}

var t5Scaled chebpoly = chebpoly{coeffs: []float64{0, 0, 0, 0, 0, 1}, domainLower: 0.5, domainUpper: 2}
var t5ScaledCumsum chebpoly = chebpoly{coeffs: []float64{1.0 / 32, 0, 0, 0, -3.0 / 32, 0, 1.0 / 16}, domainLower: 0.5, domainUpper: 2}

var t3t5 chebpoly = chebpoly{coeffs: []float64{0, 0, 0, 1, 0, 1}, domainLower: -1, domainUpper: 1}

var extendedt3t5 chebpoly = chebpoly{coeffs: []float64{0, 0, 0, 1, 0, 1, 0, 0, 0, 0}, domainLower: -1, domainUpper: 1}

func floatsNearlyEqual(found, expected float64, t *testing.T) {
	t.Helper()
	difference := math.Abs(found - expected)
	if difference > tol {
		t.Logf("Fail - Found: %f Expected: %f Diff: %f", found, expected, difference)
		t.Fail()
	} else {
		t.Logf("Pass - Found: %f Expected: %f Diff: %f", found, expected, difference)
	}
}

func floatArraysNearlyEqual(found, expected []float64, t *testing.T) {
	t.Helper()
	if len(found) != len(expected) {
		t.Logf("Length doesn't match - Found: %d Expected: %d", len(found), len(expected))
		t.Fail()
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
		t.Logf("Values don't match - Max Diff: %g Position: %d", maxDifference, maxPosition)
		t.Fail()
	} else {
		t.Logf("Values match - Max Diff: %g Position: %d", maxDifference, maxPosition)
	}
}

func extremaNearlyEqual(found, expected []extremum, t *testing.T) {
	t.Helper()
	if len(found) != len(expected) {
		t.Logf("Length doesn't match - Found: %d Expected: %d", len(found), len(expected))
		t.Fail()
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
			t.Logf("Type of extremum disagreement - Found: %v Expected: %v", v.Maximum, expected[i].Maximum)
			t.Fail()
		}
	}
	if maxPointDifference > tol || maxValueDifference > tol {
		t.Logf("Values don't match - Max Point Diff: %g at %d Max Value Diff: %g at %d", maxPointDifference, maxPointPosition, maxValueDifference, maxValuePosition)
		t.Fail()
	} else {
		t.Logf("Values match - Max Point Diff: %g at %d Max Value Diff: %g at %d", maxPointDifference, maxPointPosition, maxValueDifference, maxValuePosition)
	}
}

func chebpolysNearlyEqual(found, expected chebpoly, t *testing.T) {
	t.Helper()
	if math.Abs(found.domainLower-expected.domainLower) > tol || math.Abs(found.domainUpper-expected.domainUpper) > tol {
		t.Logf("Domains don't match - Found: [%f %f] Expected: [%f %f]", found.domainLower, found.domainUpper, expected.domainLower, expected.domainUpper)
		t.Fail()
	}
	if len(found.coeffs) != len(expected.coeffs) {
		t.Logf("Degrees don't match - Found: %d Expected: %d", len(found.coeffs)-1, len(expected.coeffs)-1)
		t.Fail()
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
		t.Logf("coeffs don't match - Max Diff: %g Position: %d", maxDifference, maxPosition)
		t.Fail()
	} else {
		t.Logf("coeffs match - Max Diff: %g Position: %d", maxDifference, maxPosition)
	}
}

func chebfunsNearlyEqual(found, expected chebfun, t *testing.T) {
	t.Helper()
	if len(found) != len(expected) {
		t.Logf("Contain a different number of chebpoly - Found: %d Expected: %d", len(found), len(expected))
	}
	for i, v := range found {
		u := expected[i]
		if math.Abs(v.domainLower-u.domainLower) > tol || math.Abs(v.domainUpper-u.domainUpper) > tol {
			t.Logf("Fail in %d - Domains don't match - Found: [%f %f] Expected: [%f %f]", i, v.domainLower, v.domainUpper, u.domainLower, u.domainUpper)
			t.Fail()
		}
		if len(v.coeffs) != len(u.coeffs) {
			t.Logf("Fail in %d - Degrees don't match - Found: %d Expected: %d", i, len(v.coeffs)-1, len(u.coeffs)-1)
			t.Fail()
		}
		maxDifference := 0.0
		maxPosition := 0
		for j, w := range v.coeffs {
			if math.Abs(w-u.coeffs[j]) > maxDifference {
				maxDifference = math.Abs(w - u.coeffs[j])
				maxPosition = j
			}
		}
		if maxDifference > tol {
			t.Logf("Fail in %d - coeffs don't match - Max Diff: %g Position: %d", i, maxDifference, maxPosition)
			t.Fail()
		} else {
			t.Logf("Pass %d - coeffs match - Max Diff: %g Position: %d", i, maxDifference, maxPosition)
		}
	}
}

func TestChebpts(t *testing.T) {

	t.Log("Chebpts")
	output := Chebpts(20, -1.0, 1.0)
	correctResult := []float64{-1.000000000000000, -0.986361303402722, -0.945817241700635, -0.879473751206489, -0.789140509396394, -0.677281571625741, -0.546948158122427, -0.401695424652969, -0.245485487140799, -0.082579345472332, 0.082579345472332, 0.245485487140799, 0.401695424652969, 0.546948158122427, 0.677281571625741, 0.789140509396394, 0.879473751206489, 0.945817241700635, 0.986361303402722, 1.000000000000000}

	floatArraysNearlyEqual(output, correctResult, t)

	t.Log("Chebpts - scaled")
	output = Chebpts(20, 0.5, 2)
	correctResult = []float64{0.500000000000000, 0.510229022447958, 0.540637068724524, 0.590394686595133, 0.658144617952705, 0.742038821280694, 0.839788881408180, 0.948728431510273, 1.065885884644401, 1.188065490895751, 1.311934509104249, 1.434114115355599, 1.551271568489727, 1.660211118591820, 1.757961178719306, 1.841855382047295, 1.909605313404867, 1.959362931275476, 1.989770977552042, 2.000000000000000}

	floatArraysNearlyEqual(output, correctResult, t)

	t.Log("Chebpts - domain switched")
	output = Chebpts(20, 2, 0.5)
	correctResult = []float64{0.500000000000000, 0.510229022447958, 0.540637068724524, 0.590394686595133, 0.658144617952705, 0.742038821280694, 0.839788881408180, 0.948728431510273, 1.065885884644401, 1.188065490895751, 1.311934509104249, 1.434114115355599, 1.551271568489727, 1.660211118591820, 1.757961178719306, 1.841855382047295, 1.909605313404867, 1.959362931275476, 1.989770977552042, 2.000000000000000}

	floatArraysNearlyEqual(output, correctResult, t)
}

func TestInterp(t *testing.T) {

	t.Log("TestInterp - t4")
	output := Interp([]float64{1, -1, 1, -1, 1}, -1, 1)
	chebpolysNearlyEqual(output, t4, t)

	t.Log("TestInterp - extendedt3t5")
	output = Interp([]float64{-2.000000000000000, -0.326351822333070, 1.439692620785908, 0.500000000000000, -0.266044443118978, 0.266044443118978, -0.500000000000000, -1.439692620785908, 0.326351822333070, 2.000000000000000}, -1, 1)
	chebpolysNearlyEqual(output, extendedt3t5, t)

	t.Log("TestInterp - sin(100x)")

	points := Chebpts(100, -1, 1)
	for i, v := range points {
		points[i] = math.Sin(100 * v)
	}

	output = Interp(points, -1, 1)

	chebpolysNearlyEqual(output, s100, t)

	t.Log("TestInterp - sin(100x) with bounds switched")
	output = Interp(points, 1, -1)
	chebpolysNearlyEqual(output, s100, t)
}

func TestEvaluate(t *testing.T) {

	floatsNearlyEqual(t5.Evaluate(Chebpts(6, -1, 1)[2]), -1.0, t)

	floatsNearlyEqual(t4.Evaluate(Chebpts(5, -1, 1)[2]), 1.0, t)

	floatsNearlyEqual(t5Scaled.Evaluate(Chebpts(6, 0.5, 2)[2]), -1.0, t)

	floatsNearlyEqual(s100.Evaluate(0.6), -0.499401732854600, t)

}

func TestCumsum(t *testing.T) {
	chebpolysNearlyEqual(t4.Cumsum(), t4Cumsum, t)

	chebpolysNearlyEqual(t5Scaled.Cumsum(), t5ScaledCumsum, t)

	chebpolysNearlyEqual(s100.Cumsum(), s100Cumsum, t)

	//TODO Tests for the split cases
}

func TestSum(t *testing.T) {

	floatsNearlyEqual(t4.Sum(), -2.0/15, t)

	floatsNearlyEqual(t5.Sum(), 0.0, t)

	floatsNearlyEqual(t4Scaled.Sum(), -1.0/10, t)

}

func TestDiff(t *testing.T) {

	chebpolysNearlyEqual(t4Cumsum.Diff(), t4, t)

	chebpolysNearlyEqual(t5ScaledCumsum.Diff(), t5Scaled, t)

	//Note that differentiating s100Cumsum doesn't give s100 but this is likely due to errors in the terms of s100Cumsum. Indeed, pasting these values into chebfun the difference in the derivative's coefficients is about 2.2e-14.
	chebpolysNearlyEqual(s100.Cumsum().Diff(), s100, t)

}

func TestRoots(t *testing.T) {

	poly := chebpoly{domainLower: -1, domainUpper: 1, coeffs: []float64{2.4}}
	correctResult := []float64{}
	floatArraysNearlyEqual(poly.Roots(), correctResult, t)

	poly = chebpoly{domainLower: -1, domainUpper: 1, coeffs: []float64{0, 1}}
	correctResult = []float64{0}
	floatArraysNearlyEqual(poly.Roots(), correctResult, t)

	poly = chebpoly{domainLower: -1, domainUpper: 1, coeffs: []float64{2, 1}}
	correctResult = []float64{}
	floatArraysNearlyEqual(poly.Roots(), correctResult, t)

	correctResult = []float64{math.Cos(4*math.Pi/5 + math.Pi/10), math.Cos(3*math.Pi/5 + math.Pi/10), 0, math.Cos(1*math.Pi/5 + math.Pi/10), math.Cos(math.Pi / 10)}
	floatArraysNearlyEqual(t5.Roots(), correctResult, t)

	correctResult = []float64{0.536707612778635, 0.809161060780646, 1.25, 1.690838939219355, 1.963292387221366}
	floatArraysNearlyEqual(t5Scaled.Roots(), correctResult, t)

	t.Log("Roots - s100")
	correctResult = []float64{-0.973476537473530, -0.941963922846026, -0.909997630679600, -0.880898165476357, -0.846951011003893, -0.818161004985501, -0.785125636034162, -0.752213060303187, -0.722848191200046, -0.693342479264069, -0.661825115146115, -0.627104608999527, -0.593792976919849, -0.562609279390930, -0.532149578262066, -0.501865177209623, -0.471554903526460, -0.441138543138549, -0.410583514910297, -0.379877544362281, -0.349016870625636, -0.318000644447448, -0.286828379638309, -0.255499328351622, -0.224013236948676, -0.192371973822785, -0.160581414207517, -0.128652968862129, -0.096604334095511, -0.064459334959752, -0.032246984407967, 0, 0.032246984407966, 0.064459334959752, 0.096604334095511, 0.128652968862129, 0.160581414207517, 0.192371973822785, 0.224013236948676, 0.255499328351621, 0.286828379638309, 0.318000644447447, 0.349016870625636, 0.379877544362281, 0.410583514910296, 0.441138543138549, 0.471554903526460, 0.501865177209623, 0.532149578262066, 0.562609279390930, 0.593792976919849, 0.627104608999527, 0.661825115146115, 0.693342479264069, 0.722848191200045, 0.752213060303187, 0.785125636034162, 0.818161004985501, 0.846951011003893, 0.880898165476357, 0.909997630679600, 0.941963922846026, 0.973476537473530}

	floatArraysNearlyEqual(s100.Roots(), correctResult, t)

	//TODO Test a scaled one that splits.
}

func TestExtrema(t *testing.T) {
	correctResult := []extremum{
		extremum{Point: -1, Value: -1, Maximum: false},
		extremum{Point: math.Cos(math.Pi * 4 / 5), Value: 1, Maximum: true},
		extremum{Point: math.Cos(math.Pi * 3 / 5), Value: -1, Maximum: false},
		extremum{Point: math.Cos(math.Pi * 2 / 5), Value: 1, Maximum: true},
		extremum{Point: math.Cos(math.Pi * 1 / 5), Value: -1, Maximum: false},
		extremum{Point: 1, Value: 1, Maximum: true},
	}
	extremaNearlyEqual(t5.Extrema(), correctResult, t)

	//TODO Scaled Extrema

}

func TestMaxAndMin(t *testing.T) {
	t.Log("MaxAndMin - t5")
	max, min := t5.MaxAndMin()
	floatsNearlyEqual(max, 1, t)
	floatsNearlyEqual(min, -1, t)

	t.Log("MaxAndMin - t5Scaled")
	max, min = t5Scaled.MaxAndMin()
	floatsNearlyEqual(max, 1, t)
	floatsNearlyEqual(min, -1, t)

	t.Log("MaxAndMin - s100")
	max, min = s100.MaxAndMin()
	floatsNearlyEqual(max, 1.391205486041127, t)
	floatsNearlyEqual(min, -1.391205486041127, t)
}

func f(x float64) float64 {
	return math.Cos(5 * x)
}

func TestAdaptive(t *testing.T) {
	poly := Adaptive(f, -1.0, 1.0)
	correctValue := math.Sin(5) / 2.5
	floatsNearlyEqual(poly.Sum(), correctValue, t)

	correctExtrema := []extremum{
		extremum{Point: -1, Value: f(-1), Maximum: true},
		extremum{Point: -math.Pi / 5, Value: -1, Maximum: false},
		extremum{Point: 0, Value: 1, Maximum: true},
		extremum{Point: math.Pi / 5, Value: -1, Maximum: false},
		extremum{Point: 1, Value: f(1), Maximum: true},
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
	floatArraysNearlyEqual(poly.Roots(), correctRoots, t)
}

//|x| on [-0.5,2]
var absX chebfun = []chebpoly{chebpoly{domainLower: -0.5, domainUpper: 0, coeffs: []float64{0.25, -0.25}}, chebpoly{domainLower: 0, domainUpper: 2, coeffs: []float64{1, 1}}}

var sgnX chebfun = []chebpoly{chebpoly{domainLower: -1, domainUpper: 0, coeffs: []float64{-1}}, chebpoly{domainLower: 0, domainUpper: 1, coeffs: []float64{1}}}

//max(0,|x| - 1) on [-2,2]
var max0absXminus1 chebfun = []chebpoly{chebpoly{domainLower: -2, domainUpper: -1, coeffs: []float64{0.5, -0.5}}, chebpoly{domainLower: -1, domainUpper: 1, coeffs: []float64{0}}, chebpoly{domainLower: 1, domainUpper: 2, coeffs: []float64{0.5, 0.5}}}

func TestAbs(t *testing.T) {
	chebfunsNearlyEqual(x.Abs(), absX, t)
}

func TestChebfunSum(t *testing.T) {
	floatsNearlyEqual(absX.Sum(), 2.125, t)
	floatsNearlyEqual(sgnX.Sum(), 0.0, t)

}

func TestChebfunCumSum(t *testing.T) {
	correctResult := []chebpoly{chebpoly{domainLower: -0.5, domainUpper: 0, coeffs: []float64{5.0 / 64, 1.0 / 16, -1.0 / 64}}, chebpoly{domainLower: 0, domainUpper: 2, coeffs: []float64{7.0 / 8, 1, 1.0 / 4}}}
	chebfunsNearlyEqual(absX.Cumsum(), correctResult, t)

	correctResult = []chebpoly{chebpoly{domainLower: -1, domainUpper: 0, coeffs: []float64{-0.5, -0.5}}, chebpoly{domainLower: 0, domainUpper: 1, coeffs: []float64{-0.5, 0.5}}}
	chebfunsNearlyEqual(sgnX.Cumsum(), correctResult, t)
}

func TestChebfunMaxAndMin(t *testing.T) {
	max, min := absX.MaxAndMin()
	floatsNearlyEqual(max, 2, t)
	floatsNearlyEqual(min, 0, t)

	max, min = sgnX.MaxAndMin()
	floatsNearlyEqual(max, 1, t)
	floatsNearlyEqual(min, -1, t)
}

func TestNewChebfun(t *testing.T){
	parts := []chebpoly{chebpoly{domainLower: -1, domainUpper: 1, coeffs: []float64{0}}, chebpoly{domainLower: 0.5, domainUpper: 2, coeffs: []float64{0}}}
	fun, pass := NewChebfun(parts);
	if fun != nil || pass != false{
		t.Fail()
	}

	parts = []chebpoly{sgnX[1], sgnX[0]}
	fun, pass = NewChebfun(parts)
	if pass == false{
		t.Fail()
	}
	chebfunsNearlyEqual(fun, sgnX, t)

	parts = []chebpoly{chebpoly{domainLower: -2, domainUpper: -1, coeffs: []float64{0.5, -0.5}}, chebpoly{domainLower: 1, domainUpper: 2, coeffs: []float64{0.5, 0.5}}}
	fun, pass = NewChebfun(parts)
	if pass == false{
		t.Fail()
	}
	chebfunsNearlyEqual(fun, max0absXminus1,t )
}

func TestChebfunEvaluate(t *testing.T){
	floatsNearlyEqual(sgnX.Evaluate(-1), -1, t)
	floatsNearlyEqual(sgnX.Evaluate(0), 1, t)
	floatsNearlyEqual(sgnX.Evaluate(1), 0, t)

	floatsNearlyEqual(absX.Evaluate(-0.3), 0.3, t)
	floatsNearlyEqual(absX.Evaluate(0.3), 0.3, t)
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
		points[i] = math.Sin(100 * v)
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Interp(points, -1, 1)
	}
}

func BenchmarkEvaluate(b *testing.B) {
	for i := 0; i < b.N; i++ {
		s100.Evaluate(0.6)
	}
}

func BenchmarkCumsum(b *testing.B) {
	for i := 0; i < b.N; i++ {
		s100.Cumsum()
	}
}

func BenchmarkSum(b *testing.B) {
	for i := 0; i < b.N; i++ {
		s100.Sum()
	}
}

func BenchmarkDiff(b *testing.B) {
	for i := 0; i < b.N; i++ {
		s100.Diff()
	}
}

func BenchmarkRoots(b *testing.B) {
	for i := 0; i < b.N; i++ {
		s100.Roots()
	}
}

func BenchmarkExtrema(b *testing.B) {
	for i := 0; i < b.N; i++ {
		s100.Extrema()
	}
}

func BenchmarkMaxAndMin(b *testing.B) {
	for i := 0; i < b.N; i++ {
		s100.MaxAndMin()
	}
}

func BenchmarkAbs(b *testing.B) {
	for i := 0; i < b.N; i++ {
		s100.Abs()
	}
}
