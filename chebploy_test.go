package chebpoly

import "testing"
import "math"
import "fmt"

const tol = 1e-14

var s100 Chebpoly = Chebpoly{Coeffs: []float64{0, -0.154290704028225, 0, -0.152568403440665, 0, -0.148391473929028, 0, -0.140345381974426, 0, -0.126473522812062, 0, -0.104580652037872, 0, -0.072787348681248, 0, -0.030396242447856, 0, 0.020968775379588, 0, 0.076187842328998, 0, 0.125961809127669, 0, 0.157435430279415, 0, 0.157008546711987, 0, 0.115723867368825, 0, 0.036268545968470, 0, -0.061483609528938, 0, -0.141355879131755, 0, -0.162779647030418, 0, -0.103422694466848, 0, 0.017396789940632, 0, 0.133719597657129, 0, 0.163752486893221, 0, 0.071287631631195, 0, -0.084406330587917, 0, -0.170701322789670, 0, -0.093304643332630, 0, 0.086165468950018, 0, 0.173896049004493, 0, 0.050635977430853, 0, -0.143987256125872, 0, -0.141435680910038, 0, 0.075164844458731, 0, 0.177526624584611, 0, -0.021548487864416, 0, -0.187386167733437, 0, 0.003928443986477, 0, 0.192677150775782, 0, -0.029671200913063, 0, -0.190366542693884, 0, 0.102004112571603, 0, 0.144140250398605, 0, -0.196030233129984, 0, 0.003334537327747, 0, 0.194749402093528, 0, -0.198857364307938, 0, 0.021068993112999, 0, 0.202711360106877, 0, -0.364232422013750, 0, 0.441210883358808, 0, -0.230487850646075}, DomainLower: -1, DomainUpper: 1}

var s100Cumsum Chebpoly = Chebpoly{Coeffs: []float64{0.008423812166895, 0, -0.000430575146890, 0, -0.000522116188955, 0, -0.000670507662883, 0, -0.000866991197648, 0, -0.001094643538709, 0, -0.001324720973193, 0, -0.001513968079764, 0, -0.001605156807108, 0, -0.001533862970817, 0, -0.001244349169967, 0, -0.000715309571631, 0, 0.000008893407655, 0, 0.000793936141215, 0, 0.001418845025006, 0, 0.001629202591623, 0, 0.001248004212544, 0, 0.000315055410274, 0, -0.000824402118938, 0, -0.001589730057993, 0, -0.001454035096456, 0, -0.000357534395668, 0, 0.001050736991614, 0, 0.001692325676295, 0, 0.000898906168768, 0, -0.000773966794570, 0, -0.001725674156564, 0, -0.000812320185690, 0, 0.001100536353336, 0, 0.001677786496179, 0, -0.000021263126799, 0, -0.001746778430393, 0, -0.000799701407233, 0, 0.001508144791281, 0, 0.001219394704919, 0, -0.001366532940857, 0, -0.001310754908259, 0, 0.001502353727627, 0, 0.001057206195926, 0, -0.001874170867086, 0, -0.000263350861419, 0, 0.002074210265418, 0, -0.001186695062248, 0, -0.001112877120731, 0, 0.002236402081827, 0, -0.001221813096783, 0, -0.000987186777141, 0, 0.003015658415535, 0, -0.004195017215482, 0, 0.003427034357168, 0, -0.001152439253230}, DomainLower: -1, DomainUpper: 1}

var t4 Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 0, 1}, DomainLower: -1, DomainUpper: 1}
var t4Cumsum Chebpoly = Chebpoly{Coeffs: []float64{-1.0 / 15, 0, 0, -1.0 / 6, 0, 1.0 / 10}, DomainLower: -1, DomainUpper: 1}

var t4Scaled Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 0, 1}, DomainLower: 0.5, DomainUpper: 2}

var t5 Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 0, 0, 1}, DomainLower: -1, DomainUpper: 1}

var t5Scaled Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 0, 0, 1}, DomainLower: 0.5, DomainUpper: 2}
var t5ScaledCumsum Chebpoly = Chebpoly{Coeffs: []float64{1.0 / 32, 0, 0, 0, -3.0 / 32, 0, 1.0 / 16}, DomainLower: 0.5, DomainUpper: 2}

var t3t5 Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 1, 0, 1}, DomainLower: -1, DomainUpper: 1}

var extendedt3t5 Chebpoly = Chebpoly{Coeffs: []float64{0, 0, 0, 1, 0, 1, 0, 0, 0, 0}, DomainLower: -1, DomainUpper: 1}

func floatArraysNearlyEqual(array1, array2 []float64) (bool, string) {
	if len(array1) != len(array2) {
		return false, fmt.Sprintf("len doesn't match - Array1: %d Array2: %d", len(array1), len(array2))
	}
	maxDifference := 0.0
	maxPosition := 0
	for i, v := range array1 {
		if math.Abs(v-array2[i]) > maxDifference {
			maxDifference = math.Abs(v - array2[i])
			maxPosition = i
		}
	}
	if maxDifference > tol {
		return false, fmt.Sprintf("Values don't match - Max Diff: %g Max Diff Position: %d", maxDifference, maxPosition)
	} else {
		return true, fmt.Sprintf("Max Diff: %g Max Diff Position: %d", maxDifference, maxPosition)
	}
}

func chebpolysNearlyEqual(f, g Chebpoly) (bool, string) {
	if math.Abs(f.DomainLower-g.DomainLower) > tol || math.Abs(f.DomainUpper-g.DomainUpper) > tol {
		return false, "Domains don't match"
	}
	b, err := floatArraysNearlyEqual(f.Coeffs, g.Coeffs)
	return b, err
}

func TestChebpts(t *testing.T) {

	t.Log("Chebpts")
	output := Chebpts(20, -1.0, 1.0)
	correctResult := []float64{-1.000000000000000, -0.986361303402722, -0.945817241700635, -0.879473751206489, -0.789140509396394, -0.677281571625741, -0.546948158122427, -0.401695424652969, -0.245485487140799, -0.082579345472332, 0.082579345472332, 0.245485487140799, 0.401695424652969, 0.546948158122427, 0.677281571625741, 0.789140509396394, 0.879473751206489, 0.945817241700635, 0.986361303402722, 1.000000000000000}

	b, err := floatArraysNearlyEqual(output, correctResult)
	t.Log(err)
	if !b {
		t.Fail()
	}

	t.Log("Chebpts - scaled")
	output = Chebpts(20, 0.5, 2)
	correctResult = []float64{0.500000000000000, 0.510229022447958, 0.540637068724524, 0.590394686595133, 0.658144617952705, 0.742038821280694, 0.839788881408180, 0.948728431510273, 1.065885884644401, 1.188065490895751, 1.311934509104249, 1.434114115355599, 1.551271568489727, 1.660211118591820, 1.757961178719306, 1.841855382047295, 1.909605313404867, 1.959362931275476, 1.989770977552042, 2.000000000000000}

	b, err = floatArraysNearlyEqual(output, correctResult)
	t.Log(err)
	if !b {
		t.Fail()
	}

	t.Log("Chebpts - domain switched")
	output = Chebpts(20, 2, 0.5)
	correctResult = []float64{0.500000000000000, 0.510229022447958, 0.540637068724524, 0.590394686595133, 0.658144617952705, 0.742038821280694, 0.839788881408180, 0.948728431510273, 1.065885884644401, 1.188065490895751, 1.311934509104249, 1.434114115355599, 1.551271568489727, 1.660211118591820, 1.757961178719306, 1.841855382047295, 1.909605313404867, 1.959362931275476, 1.989770977552042, 2.000000000000000}

	b, err = floatArraysNearlyEqual(output, correctResult)
	t.Log(err)
	if !b {
		t.Fail()
	}
}

func TestInterp(t *testing.T) {

	output := Interp([]float64{1, -1, 1, -1, 1}, -1, 1)
	b, err := chebpolysNearlyEqual(*output, t4)
	t.Log(err)
	if !b {
		t.Fail()
	}

	output = Interp([]float64{-2.000000000000000, -0.326351822333070, 1.439692620785908, 0.500000000000000, -0.266044443118978, 0.266044443118978, -0.500000000000000, -1.439692620785908, 0.326351822333070, 2.000000000000000}, -1, 1)
	b, err = chebpolysNearlyEqual(*output, extendedt3t5)
	t.Log(err)
	if !b {
		t.Fail()
	}

	//Test interpolating sin(100x) in 100 points on [-1,1]

	t.Log("TestInterp - sin(100x)")

	points := Chebpts(100, -1, 1)
	for i, v := range points {
		points[i] = math.Sin(100 * v)
	}

	output = Interp(points, -1, 1)

	b, err = chebpolysNearlyEqual(*output, s100)
	t.Log(err)
	if !b {
		t.Fail()
	}

	t.Log("Starting sin(100x) test with the domain bounds switched")
	output = Interp(points, 1, -1)
	b, err = chebpolysNearlyEqual(*output, s100)
	t.Log(err)
	if !b {
		t.Fail()
	}
}

func TestEvaluate(t *testing.T) {

	if math.Abs(t5.Evaluate(Chebpts(6, -1, 1)[2])+1) > tol {
		t.FailNow()
	}

	if math.Abs(t4.Evaluate(Chebpts(5, -1, 1)[2])-1) > tol {
		t.Log(t4.Evaluate(Chebpts(5, -1, 1)[2]))
		t.FailNow()
	}

	if math.Abs(t5Scaled.Evaluate(Chebpts(6, 0.5, 2)[2])+1) > tol {
		t.FailNow()
	}

	if math.Abs(s100.Evaluate(0.6)+0.499401732854600) > tol {
		t.FailNow()
	}
}

func TestCumsum(t *testing.T) {
	b, err := chebpolysNearlyEqual(t4.Cumsum(), t4Cumsum)
	t.Log(err)
	if !b {
		t.Fail()
	}

	b, err = chebpolysNearlyEqual(t5Scaled.Cumsum(), t5ScaledCumsum)
	t.Log(err)
	if !b {
		t.Fail()
	}

	b, err = chebpolysNearlyEqual(s100.Cumsum(), s100Cumsum)
	t.Log(err)
	if !b {
		t.Fail()
	}
}

func TestSum(t *testing.T) {
	if math.Abs(t4.Sum()+2.0/15) > tol {
		t.FailNow()
	}

	if math.Abs(t5.Sum()) > tol {
		t.FailNow()
	}
	if math.Abs(t4Scaled.Sum()+1.0/10) > tol {
		t.FailNow()
	}
}

func TestDiff(t *testing.T) {
	b, err := chebpolysNearlyEqual(t4Cumsum.Diff(), t4)
	t.Log(err)
	if !b {
		t.Fail()
	}

	b, err = chebpolysNearlyEqual(t5ScaledCumsum.Diff(), t5Scaled)
	t.Log(err)
	if !b {
		t.Fail()
	}

	//Note that differentiating s100Cumsum doesn't give s100 but this is likely due to errors in the terms of s100Cumsum. Indeed, pasting these values into chebfun the difference in the derivative's coefficients is about 2.2e-14.
	b, err = chebpolysNearlyEqual(s100.Cumsum().Diff(), s100)
	t.Log(err)
	if !b {
		t.Fail()
	}
}

func TestRoots(t *testing.T) {
	poly := Chebpoly{DomainLower: -1, DomainUpper: 1, Coeffs: []float64{0}}
	correctResult := []float64{0}
	b, err := floatArraysNearlyEqual(poly.Roots(), correctResult)
	t.Log(err)
	if !b {
		t.Fail()
	}

	poly = Chebpoly{DomainLower: -1, DomainUpper: 1, Coeffs: []float64{2.4}}
	correctResult = []float64{}
	b, err = floatArraysNearlyEqual(poly.Roots(), correctResult)
	t.Log(err)
	if !b {
		t.Fail()
	}

	poly = Chebpoly{DomainLower: -1, DomainUpper: 1, Coeffs: []float64{0, 1}}
	correctResult = []float64{0}
	b, err = floatArraysNearlyEqual(poly.Roots(), correctResult)
	t.Log(err)
	if !b {
		t.Fail()
	}

	poly = Chebpoly{DomainLower: -1, DomainUpper: 1, Coeffs: []float64{2, 1}}
	correctResult = []float64{}
	b, err = floatArraysNearlyEqual(poly.Roots(), correctResult)
	t.Log(err)
	if !b {
		t.Fail()
	}

	correctResult = []float64{math.Cos(4*math.Pi/5 + math.Pi/10), math.Cos(3*math.Pi/5 + math.Pi/10), 0, math.Cos(1*math.Pi/5 + math.Pi/10), math.Cos(math.Pi / 10)}
	b, err = floatArraysNearlyEqual(t5.Roots(), correctResult)
	t.Log(err)
	if !b {
		t.Fail()
	}

	correctResult = []float64{0.536707612778635, 0.809161060780646, 1.25, 1.690838939219355, 1.963292387221366}
	b, err = floatArraysNearlyEqual(t5Scaled.Roots(), correctResult)
	t.Log(err)
	if !b {
		t.Fail()
	}

	t.Log("s100")
	correctResult = []float64{-0.973476537473530, -0.941963922846026, -0.909997630679600, -0.880898165476357, -0.846951011003893, -0.818161004985501, -0.785125636034162, -0.752213060303187, -0.722848191200046, -0.693342479264069, -0.661825115146115, -0.627104608999527, -0.593792976919849, -0.562609279390930, -0.532149578262066, -0.501865177209623, -0.471554903526460, -0.441138543138549, -0.410583514910297, -0.379877544362281, -0.349016870625636, -0.318000644447448, -0.286828379638309, -0.255499328351622, -0.224013236948676, -0.192371973822785, -0.160581414207517, -0.128652968862129, -0.096604334095511, -0.064459334959752, -0.032246984407967, 0, 0.032246984407966, 0.064459334959752, 0.096604334095511, 0.128652968862129, 0.160581414207517, 0.192371973822785, 0.224013236948676, 0.255499328351621, 0.286828379638309, 0.318000644447447, 0.349016870625636, 0.379877544362281, 0.410583514910296, 0.441138543138549, 0.471554903526460, 0.501865177209623, 0.532149578262066, 0.562609279390930, 0.593792976919849, 0.627104608999527, 0.661825115146115, 0.693342479264069, 0.722848191200045, 0.752213060303187, 0.785125636034162, 0.818161004985501, 0.846951011003893, 0.880898165476357, 0.909997630679600, 0.941963922846026, 0.973476537473530}

	b, err = floatArraysNearlyEqual(s100.Roots(), correctResult)
	t.Log(err)
	if !b {
		t.Fail()
	}

	//TODO Test a scaled one that splits.
}

func extremaNearlyEqual(array1, array2 []extremum) (bool, string) {
	if len(array1) != len(array2) {
		return false, fmt.Sprintf("len doesn't match - Extrema1: %d Extrema2: %d", len(array1), len(array2))
	}
	maxPointDifference := 0.0
	maxPointPosition := 0

	maxValueDifference := 0.0
	maxValuePosition := 0
	for i, v := range array1 {
		if math.Abs(v.Point-array2[i].Point) > maxPointDifference {
			maxPointDifference = math.Abs(v.Point - array2[i].Point)
			maxPointPosition = i
		}
		if math.Abs(v.Value-array2[i].Value) > maxValueDifference {
			maxValueDifference = math.Abs(v.Value - array2[i].Value)
			maxValuePosition = i
		}
		if v.Maximum != array2[i].Maximum {
			return false, fmt.Sprintf("Type of extremum disagreement - Extrema 1: %v Extrema 2: %v", v.Maximum, array2[i].Maximum)
		}
	}
	if maxPointDifference > tol || maxValueDifference > tol {
		return false, fmt.Sprintf("Values don't match - Max Point Diff: %g at %d Max Value Diff: %g at %d", maxPointDifference, maxPointPosition, maxValueDifference, maxValuePosition)
	} else {
		return true, fmt.Sprintf("Max Point Diff: %g at %d Max Value Diff: %g at %d", maxPointDifference, maxPointPosition, maxValueDifference, maxValuePosition)
	}
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
	b, err := extremaNearlyEqual(t5.Extrema(), correctResult)
	t.Log(err)
	if !b {
		t.Fail()
	}
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
