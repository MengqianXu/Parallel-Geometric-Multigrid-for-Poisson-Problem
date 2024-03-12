using LinearAlgebra
using Random

# function Jacobi(A, b, x0)
# 	D = Diagonal(A)
# 	E = -(LowerTriangular(A) - D)
# 	F = -(UpperTriangular(A) - D)
# 	M = D
# 	N = E + F
# 	Minv = inv(M)
	
# 	ancien = x0
# 	res = Minv * (N*ancien + b)
# 	norme = norm(res - ancien)

# 	while(norme >= 10^(-15))
# 		ancien = res
# 		res = Minv * (N*ancien + b)
# 		norme = norm(res - ancien)
# 	end

# 	return res
# end

# function GaussSeidel(A, b, x0)
# 	D = Diagonal(A)
# 	E = -(LowerTriangular(A) - D)
# 	F = -(UpperTriangular(A) - D)
# 	M = D - E
# 	N = F
# 	Minv = inv(M)
	
# 	ancien = x0
# 	res = Minv * (N*ancien + b)
# 	norme = norm(res - ancien)

# 	while(norme >= 10^(-15))
# 		ancien = res
# 		res = Minv * (N*ancien + b)
# 		norme = norm(res - ancien)
# 	end

# 	return res
# end

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

function Jacobi(A, b, x0, tol = 10^(-15))
	ancien = copy(x0)
	nouveau = copy(b)
	n = length(b)
	
	for i = 1:n
		for j = 1:n
			if j != i 
				nouveau[i] -= A[i, j]*ancien[j]
			end
		end
		nouveau[i] /= A[i, i]
	end
	
	norme = norm(nouveau - ancien)
	while norme >= tol
		ancien = copy(nouveau)
		nouveau = copy(b)
		
		for i = 1:n
			for j = 1:n
				if j != i 
					nouveau[i] -= A[i, j]*ancien[j]
				end
			end
			nouveau[i] /= A[i, i]
		end
		
		norme = norm(nouveau - ancien)
	end

	return nouveau
end

function GaussSeidel(A, b, x0, tol = 10^(-15))
	ancien = copy(x0)
	nouveau = copy(b)
	n = length(b)
	
	for i = 1:n
		for j = 1:n
			if j != i
				if j < i 
					nouveau[i] -= A[i, j]*nouveau[j]
				else
					nouveau[i] -= A[i, j]*ancien[j]
				end
			end
		end
		nouveau[i] /= A[i, i]
	end
	
	norme = norm(nouveau - ancien)
	while norme >= 10^(-15)
		ancien = copy(nouveau)
		nouveau = copy(b)
		
		for i = 1:n
			for j = 1:n
				if j != i
					if j < i 
						nouveau[i] -= A[i, j]*nouveau[j]
					else
						nouveau[i] -= A[i, j]*ancien[j]
					end
				end
			end
			nouveau[i] /= A[i, i]
		end
		
		norme = norm(nouveau - ancien)
	end

	return nouveau
end

function Multigrid(A, F, N::Int64, x0, pre, post)
	if N == 1
		return (1/4)*F
	else
		U = copy(x0)

		for i = 1:pre
			U = Jacobi(A, F, U)
		end

		newN = floor(Int64, N/2)
		newA = Creer_A(newN)
		newF = zeros(newN*newN)
		newU = zeros(newN*newN)

		for i = 1:newN
			for j = 1:newN
				newF[j + N*(i - 1)] = F[2j + N*(2i - 1)]
				newU[j + N*(i - 1)] = U[2j + N*(2i - 1)]
			end
		end
		
		newU = Multigrid(newA, newF, newN, newU, pre, post)

		for i = 1:newN
			for j = 1:newN
				U[2j + N*(2i - 1)] = newU[j + N*(i - 1)]
			end
		end

		for i = 1:pre
			U = Jacobi(A, F, U)
		end

		return U
	end
end