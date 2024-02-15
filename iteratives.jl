using LinearAlgebra
using Random

function Jacobi(A, b, x0)
	D = Diagonal(A)
	E = -(LowerTriangular(A) - D)
	F = -(UpperTriangular(A) - D)
	M = D
	N = E + F
	Minv = inv(M)
	
	ancien = x0
	res = Minv * (N*ancien + b)
	norme = norm(res - ancien)

	while(norme >= 10^(-15))
		ancien = res
		res = Minv * (N*ancien + b)
		norme = norm(res - ancien)
	end

	return res
end

function GaussSeidel(A, b, x0)
	D = Diagonal(A)
	E = -(LowerTriangular(A) - D)
	F = -(UpperTriangular(A) - D)
	M = D - E
	N = F
	Minv = inv(M)
	
	ancien = x0
	res = Minv * (N*ancien + b)
	norme = norm(res - ancien)

	while(norme >= 10^(-15))
		ancien = res
		res = Minv * (N*ancien + b)
		norme = norm(res - ancien)
	end

	return res
end

function SOR(A, b, w, x0)
	D = Diagonal(A)
	E = -(LowerTriangular(A) - D)
	F = -(UpperTriangular(A) - D)
	M = D/w - E
	N = (1/w - 1)D + F
	Minv = inv(M)
	
	ancien = x0
	res = Minv * (N*ancien + b)
	norme = norm(res - ancien)

	while(norme >= 10^(-15))
		ancien = res
		res = Minv * (N*ancien + b)
		norme = norm(res - ancien)
	end

	return res
end

A = kron(4I(3), Tridiagonal(-ones((3, 3))) + 5I) + kron(Tridiagonal(-ones((3, 3))) + I, I(3))
b = ones((9, 1))
x0 = rand(1:9, 9)
resJ = Jacobi(A, b, x0)
resGS = GaussSeidel(A, b, x0)
w = 0.5
resS = SOR(A, b, w, x0)

abs(norm(A\b - resJ))
abs(norm(A\b - resGS))
abs(norm(A\b - resS))