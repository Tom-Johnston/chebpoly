package chebpoly

import (
	"math/cmplx"
	"testing"
)

func cmplxArraysNearlyEqual(found, expected []complex128, t *testing.T) {
	t.Helper()
	if len(found) != len(expected) {
		t.Logf("Length doesn't match - Found: %d Expected: %d", len(found), len(expected))
		t.Log("Found: ", found)
		t.Log("Expected: ", expected)
		t.FailNow()
	}
	maxDifference := 0.0
	maxPosition := 0
	for i, v := range found {
		if diff := cmplx.Abs(v - expected[i]); diff > maxDifference {
			maxDifference = diff
			maxPosition = i
		}
	}
	if maxDifference > tol {
		t.Logf("Values don't match - Max Diff: %g Position: %d", maxDifference, maxPosition)
		t.Log("Found: ", found)
		t.Log("Expected: ", expected)
		t.Fail()
	}
}

func TestHessEig(t *testing.T) {
	input := []float64{7, 3, 4, -11, -9, -2, -6, 4, -5, 7, 1, 12, -1, -9, 2, 2, 9, 1, -8, 0, -1, 5, 0, 8, -4, 3, -5, 7, 2, 10, 6, 1, 4, -11, -7, -1}
	hess(input, 6)
	output := eig(input, 6)
	expected := []complex128{1 - 2i, 1 + 2i, 3, 4, 5 - 6i, 5 + 6i}
	cmplxArraysNearlyEqual(output, expected, t)
}
