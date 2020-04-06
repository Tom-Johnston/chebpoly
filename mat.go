package chebpoly

import (
	"fmt"
	"math"
	"math/cmplx"
	"sort"
	"strings"
)

func formatMatrix(mat []float64, n int) string {
	var sb strings.Builder
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			sb.WriteString(fmt.Sprintf("%16f \t", mat[i*n+j]))
		}
		sb.WriteString("\n")
	}
	sb.WriteString("\n")
	return sb.String()
}

func hess(mat []float64, n int) {
	u := make([]float64, n)
	for k := 0; k < n-2; k++ {
		norm := 0.0
		for i := k + 2; i < n; i++ {
			norm += mat[i*n+k] * mat[i*n+k]
		}
		if norm == 0 {
			continue
		}
		alpha := norm + mat[(k+1)*n+k]*mat[(k+1)*n+k]
		if mat[(k+1)*n+k] < 0 {
			alpha = math.Sqrt(alpha)
		} else {
			alpha = -math.Sqrt(alpha)
		}
		u = u[:n-k-1]
		norm += (mat[(k+1)*n+k] - alpha) * (mat[(k+1)*n+k] - alpha)
		norm = math.Sqrt(norm)
		u[0] = (mat[(k+1)*n+k] - alpha) / norm
		for i := 1; i < n-k-1; i++ {
			u[i] = mat[(k+i+1)*n+k] / norm
		}
		for i := k; i < n; i++ {
			scale := 0.0
			for j := 0; j < n-k-1; j++ {
				scale += mat[(k+j+1)*n+i] * u[j]
			}
			scale *= 2
			for j := 0; j < n-k-1; j++ {
				mat[(k+j+1)*n+i] -= u[j] * scale
			}
		}

		for i := 0; i < n; i++ {
			scale := 0.0
			for j := 0; j < n-k-1; j++ {
				scale += mat[i*n+(k+j+1)] * u[j]
			}
			scale *= 2
			for j := 0; j < n-k-1; j++ {
				mat[i*n+(k+j+1)] -= u[j] * scale
			}
		}
	}
}

//mat contains a matrix in upper Hessenberg form where the entires are in the order
func eig(mat []float64, n int) []complex128 {
	// printMatrix(mat, n)

	const eps = 1e-15
	maxIter := 100

	eigs := make([]complex128, 0, n)

	u := make([]float64, 3)
	p := n

	iter := 0
	for p > 2 && iter < maxIter {
		iter++
		q := p - 1
		s := mat[(q-1)*n+(q-1)] + mat[(p-1)*n+(p-1)]
		t := mat[(q-1)*n+(q-1)]*mat[(p-1)*n+(p-1)] - mat[(q-1)*n+p-1]*mat[(p-1)*n+q-1]
		x := mat[0]*mat[0] + mat[1]*mat[n] - s*mat[0] + t
		y := mat[n] * (mat[0] + mat[n+1] - s)
		z := mat[n] * mat[2*n+1]
		for k := 0; k < p-2; k++ {
			var alpha float64
			if x > 0 {
				alpha = -math.Sqrt(x*x + y*y + z*z)
			} else {
				alpha = math.Sqrt(x*x + y*y + z*z)
			}
			u[0] = x - alpha
			u[1] = y
			u[2] = z

			//Normalise u.
			norm := math.Sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2])
			if norm > 0 {
				u[0] /= norm
				u[1] /= norm
				u[2] /= norm
			}

			r := k - 1
			if r == -1 {
				r = 0
			}
			for i := r; i < n; i++ {
				scale := 2 * (u[0]*mat[k*n+i] + u[1]*mat[(k+1)*n+i] + u[2]*mat[(k+2)*n+i])
				mat[k*n+i] -= u[0] * scale
				mat[(k+1)*n+i] -= u[1] * scale
				mat[(k+2)*n+i] -= u[2] * scale
			}

			r = k + 3
			if r == p {
				r = p - 1
			}
			for i := 0; i <= r; i++ {
				scale := 2 * (u[0]*mat[i*n+k] + u[1]*mat[i*n+k+1] + u[2]*mat[i*n+k+2])
				mat[i*n+k] -= u[0] * scale
				mat[i*n+k+1] -= u[1] * scale
				mat[i*n+k+2] -= u[2] * scale
			}
			x = mat[(k+1)*n+k]
			y = mat[(k+2)*n+k]
			if k < p-3 {
				z = mat[(k+3)*n+k]
			}
		}
		norm := math.Sqrt(x*x + y*y)
		if norm == 0 {
			norm = 1
		}
		c := x / norm
		s = -y / norm
		for k := p - 3; k < n; k++ {
			mat[(p-2)*n+k], mat[(p-1)*n+k] = c*mat[(p-2)*n+k]-s*mat[(p-1)*n+k], s*mat[(p-2)*n+k]+c*mat[(p-1)*n+k]
		}
		for k := 0; k < p; k++ {
			mat[k*n+p-2], mat[k*n+p-1] = c*mat[k*n+p-2]-s*mat[k*n+p-1], s*mat[k*n+p-2]+c*mat[k*n+p-1]
		}
		start := -1
		subdivide := false
		for i := 0; i < p-3; i++ {
			if math.Abs(mat[(i+1)*n+i]) < eps*(math.Max(math.Abs(mat[i*n+i]), math.Abs(mat[(i+1)*n+i+1]))) {
				//Subdivide
				subdivide = true
				dist := i - start
				tmp := make([]float64, dist*dist)
				for j := range tmp {
					tmp[j] = mat[(start+1+(j/dist))*n+start+1+(j%dist)]
				}
				e := eig(tmp, dist)
				eigs = append(eigs, e...)
				start = i
			}
		}
		if subdivide {
			dist := p - 1 - start
			tmp := make([]float64, dist*dist)
			for j := range tmp {
				tmp[j] = mat[(start+1+(j/dist))*n+start+1+(j%dist)]
			}
			e := eig(tmp, dist)
			eigs = append(eigs, e...)
			sort.Slice(eigs, func(i, j int) bool {
				if math.Abs(cmplx.Abs(eigs[i])-cmplx.Abs(eigs[j])) > tol {
					return cmplx.Abs(eigs[i]) < cmplx.Abs(eigs[j])
				}
				if math.Abs(real(eigs[i])-real(eigs[j])) > tol {
					return real(eigs[i]) < real(eigs[j])
				}

				return imag(eigs[i]) < imag(eigs[j])
			})

			return eigs
		}
		if math.Abs(mat[(p-1)*n+q-1]) < eps*(math.Abs(mat[(q-1)*n+(q-1)])+math.Abs(mat[(p-1)*n+(p-1)]))+eps {

			mat[(p-1)*n+q-1] = 0
			eigs = append(eigs, complex(mat[(p-1)*n+p-1], 0))
			p--
			iter = 0
		} else if math.Abs(mat[(p-2)*n+q-2]) < eps*(math.Abs(mat[(q-2)*n+(q-2)])+math.Abs(mat[(q-1)*n+(q-1)]))+eps {
			mat[(p-2)*n+q-2] = 0
			eigs = append(eigs, eig2x2(mat[(q-1)*n+(q-1)], mat[(q-1)*n+q], mat[q*n+q-1], mat[q*n+q])...)
			p -= 2
			iter = 0
		}
	}

	if p > 2 {
		panic("Failed to converge")
	}

	if p == 1 {
		eigs = append(eigs, complex(mat[0], 0))
	} else {
		eigs = append(eigs, eig2x2(mat[0], mat[1], mat[n], mat[n+1])...)
	}
	sort.Slice(eigs, func(i, j int) bool {
		if math.Abs(cmplx.Abs(eigs[i])-cmplx.Abs(eigs[j])) > tol {
			return cmplx.Abs(eigs[i]) < cmplx.Abs(eigs[j])
		}
		if math.Abs(real(eigs[i])-real(eigs[j])) > tol {
			return real(eigs[i]) < real(eigs[j])
		}

		return imag(eigs[i]) < imag(eigs[j])
	})
	return eigs
}

func eig2x2(a, b, c, d float64) []complex128 {
	tr := a + d
	det := a*d - b*c
	if dis := tr*tr - 4*det; dis < 0 {
		sqrt := math.Sqrt(-dis)
		return []complex128{complex(tr/2, sqrt/2), complex(tr/2, -sqrt/2)}
	}
	sqrt := math.Sqrt(tr*tr - 4*det)
	return []complex128{complex((tr+sqrt)/2, 0), complex((tr-sqrt)/2, 0)}
}
