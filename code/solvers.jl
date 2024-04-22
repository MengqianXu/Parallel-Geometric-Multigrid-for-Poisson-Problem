using LinearAlgebra
using Random
using SparseArrays

# function Jacobi(A, b, x0, tol = 10^(-15), MaxIter = 10000)
# 	D = Diagonal(A)
# 	E = -(LowerTriangular(A) - D)
# 	F = -(UpperTriangular(A) - D)
# 	M = D
# 	N = E + F
# 	Minv = inv(M)
	
# 	ancien = x0
# 	res = Minv * (N*ancien + b)
# 	norme = norm(res - ancien)
# 	iteration = 1

# 	while (norme >= tol) && (iteration < MaxIter)
# 		ancien = res
# 		res = Minv * (N*ancien + b)
# 		norme = norm(res - ancien)
# 		iteration += 1
# 	end

# 	return res
# end

# function GaussSeidel(A, b, x0, tol = 10^(-15), MaxIter = 10000)
# 	D = Diagonal(A)
# 	E = -(LowerTriangular(A) - D)
# 	F = -(UpperTriangular(A) - D)
# 	M = D - E
# 	N = F
# 	Minv = inv(M)
	
# 	ancien = x0
# 	res = Minv * (N*ancien + b)
# 	norme = norm(res - ancien)
# 	iteration = 1

# 	while (norme >= tol) && (iteration < MaxIter)
# 		ancien = res
# 		res = Minv * (N*ancien + b)
# 		norme = norm(res - ancien)
# 		iteration += 1
# 	end

# 	return res
# end

# function SOR(A, b, w, x0, tol = 10^(-15), MaxIter = 10000)
# 	D = Diagonal(A)
# 	E = -(LowerTriangular(A) - D)
# 	F = -(UpperTriangular(A) - D)
# 	M = D/w - E
# 	N = (1/w - 1)D + F
# 	Minv = inv(M)
	
# 	ancien = x0
# 	res = Minv * (N*ancien + b)
# 	norme = norm(res - ancien)
# 	iteration = 1

# 	while (norme >= tol) && (iteration < MaxIter)
# 		ancien = res
# 		res = Minv * (N*ancien + b)
# 		norme = norm(res - ancien)
# 		iteration += 1
# 	end

# 	return res
# end

# function Jacobi(A, b, x0, tol = 10^(-15), MaxIter = 10000)
# 	ancien = copy(x0)
# 	nouveau = copy(b)
# 	n = length(b)
	
# 	for i = 1:n
# 		for j = 1:n
# 			if j != i
# 				nouveau[i] -= A[i, j]*ancien[j]
# 			end
# 		end
# 		nouveau[i] /= A[i, i]
# 	end
# 	norme = norm(nouveau - ancien)
# 	iteration = 1

# 	while (norme >= tol) && (iteration < MaxIter)
# 		ancien = copy(nouveau)
# 		nouveau = copy(b)
		
# 		for i = 1:n
# 			for j = 1:n
# 				if j != i
# 					nouveau[i] -= A[i, j]*ancien[j]
# 				end
# 			end
# 			nouveau[i] /= A[i, i]
# 		end
		
# 		norme = norm(nouveau - ancien)
# 		iteration += 1
# 	end

# 	return nouveau
# end

# function GaussSeidel(A, b, x0, tol = 10^(-15), MaxIter = 10000)
# 	ancien = copy(x0)
# 	nouveau = copy(b)
# 	n = length(b)
	
# 	for i = 1:n
# 		for j = 1:n
# 			if j != i
# 				if j < i 
# 					nouveau[i] -= A[i, j]*nouveau[j]
# 				else
# 					nouveau[i] -= A[i, j]*ancien[j]
# 				end
# 			end
# 		end
# 		nouveau[i] /= A[i, i]
# 	end
	
# 	norme = norm(nouveau - ancien)
# 	iteration = 1

# 	while (norme >= tol) && (iteration < MaxIter)
# 		ancien = copy(nouveau)
# 		nouveau = copy(b)
		
# 		for i = 1:n
# 			for j = 1:n
# 				if j != i
# 					if j < i 
# 						nouveau[i] -= A[i, j]*nouveau[j]
# 					else
# 						nouveau[i] -= A[i, j]*ancien[j]
# 					end
# 				end
# 			end
# 			nouveau[i] /= A[i, i]
# 		end
		
# 		norme = norm(nouveau - ancien)
# 		iteration += 1
# 	end

# 	return nouveau
# end

# function SOR(A, b, x0, tol = 10^(-15), MaxIter = 10000)
# 	ancien = copy(x0)
# 	nouveau = w*copy(b)
# 	n = length(b)
	
# 	for i = 1:n
# 		for j = 1:n
# 			if j != i
# 				if j < i 
# 					nouveau[i] -= w*A[i, j]*nouveau[j]
# 				else
# 					nouveau[i] -= w*A[i, j]*ancien[j]
# 				end
# 			end
# 		end
# 		nouveau[i] /= A[i, i]
# 		nouveau[i] += (1 - w)*ancien[i]
# 	end
	
# 	norme = norm(nouveau - ancien)
# 	iteration = 1

# 	while (norme >= tol) && (iteration < MaxIter)
# 		ancien = copy(nouveau)
# 		nouveau = w*copy(b)
		
# 		for i = 1:n
# 			for j = 1:n
# 				if j != i
# 					if j < i 
# 						nouveau[i] -= w*A[i, j]*nouveau[j]
# 					else
# 						nouveau[i] -= w*A[i, j]*ancien[j]
# 					end
# 				end
# 			end
# 			nouveau[i] /= A[i, i]
# 			nouveau[i] += (1 - w)*ancien[i]
# 		end
		
# 		norme = norm(nouveau - ancien)
# 		iteration += 1
# 	end

# 	return nouveau
# end

function Jacobi(A, b, x0, tol = 10^(-15), MaxIter = 10000)
	ancien = copy(x0)
	nouveau = copy(b)
	C, L, V = findnz(A)
	n = length(C)
	
	prec = L[1]
	for k = 1:n
		i = L[k]
		j = C[k]
		v = V[k]
		
		if i != prec
			nouveau[prec] /= A[prec, prec]
			prec = i
		end
		if j != i
			nouveau[i] -= v*ancien[j]
		end
	end
	nouveau[prec] /= A[prec, prec]
	iteration = 1
	
	norme = norm(nouveau - ancien)
	
	while (norme >= tol) && (iteration < MaxIter)
		ancien = copy(nouveau)
		nouveau = copy(b)

		prec = L[1]
		for k = 1:n
			i = L[k]
			j = C[k]
			v = V[k]
			if i != prec
				nouveau[prec] /= A[prec, prec]
				prec = i
			end
			if j != i
				nouveau[i] -= v*ancien[j]
			end
		end
		nouveau[prec] /= A[prec, prec]

		norme = norm(nouveau - ancien)
		iteration += 1
	end

	return nouveau
end

function GaussSeidel(A, b, x0, tol = 10^(-15), MaxIter = 10000)
	ancien = copy(x0)
	nouveau = copy(b)
	C, L, V = findnz(A)
	n = length(C)
	
	prec = L[1]
	for k = 1:n
		i = L[k]
		j = C[k]
		v = V[k]
		
		if i != prec
			nouveau[prec] /= A[prec, prec]
			prec = i
		end
		if j != i
			if j < i 
				nouveau[i] -= A[i, j]*nouveau[j]
			else
				nouveau[i] -= A[i, j]*ancien[j]
			end
		end
	end
	nouveau[prec] /= A[prec, prec]
	
	norme = norm(nouveau - ancien)
	iteration = 1
	while (norme >= tol) && (iteration < MaxIter)
		ancien = copy(nouveau)
		nouveau = copy(b)
		
		prec = L[1]
		for k = 1:n
			i = L[k]
			j = C[k]
			v = V[k]
			
			if i != prec
				nouveau[prec] /= A[prec, prec]
				prec = i
			end
			if j != i
				if j < i 
					nouveau[i] -= A[i, j]*nouveau[j]
				else
					nouveau[i] -= A[i, j]*ancien[j]
				end
			end
		end
		nouveau[prec] /= A[prec, prec]
		
		norme = norm(nouveau - ancien)
		iteration += 1
	end

	return nouveau
end

function SOR(A, b, w, x0, tol = 10^(-15), MaxIter = 10000)
	ancien = x0
	nouveau = w*b
	C, L, V = findnz(A)
	n = length(C)

	prec = L[1]
	for k = 1:n
		i = L[k]
		j = C[k]
		v = V[k]
		
		if i != prec
			nouveau[prec] /= A[prec, prec]
			nouveau[prec] += (1 - w)*ancien[prec]
			prec = i
		end
		if j != i
			if j < i 
				nouveau[i] -= w*A[i, j]*nouveau[j]
			else
				nouveau[i] -= w*A[i, j]*ancien[j]
			end
		end
	end
	nouveau[prec] /= A[prec, prec]
	nouveau[prec] += (1 - w)*ancien[prec]
	
	norme = norm(nouveau - ancien)
	iteration = 1

	while (norme >= tol) && (iteration < MaxIter)
		ancien = nouveau
		nouveau = w*b;

		prec = L[1]
		for k = 1:n
			i = L[k]
			j = C[k]
			v = V[k]
			
			if i != prec
				nouveau[prec] /= A[prec, prec]
				nouveau[prec] += (1 - w)*ancien[prec]
				prec = i
			end
			if j != i
				if j < i 
					nouveau[i] -= w*A[i, j]*nouveau[j]
				else
					nouveau[i] -= w*A[i, j]*ancien[j]
				end
			end
		end
		nouveau[prec] /= A[prec, prec]
		nouveau[prec] += (1 - w)*ancien[prec]

		norme = norm(nouveau - ancien)
		iteration += 1
	end

	return nouveau
end

function MultigridJ(A, F, N::Int64, x0, pre, post)
	# println("début multigrid, N = ", N)
	if N == 1
		return (1/4)*F
	else
		U = copy(x0)

		# println("Pre smoothing")
		for i = 1:pre
			# @show i
			U = Jacobi(A, F, U)
		end

		newN = floor(Int64, N/2)
		newA = Creer_A(newN)
		newF = zeros(newN*newN)
		newU = zeros(newN*newN)

		# println("Fine to coarse")
		for i = 1:newN
			# @show i
			for j = 1:newN
				# @show j
				newF[j + newN*(i - 1)] = F[2j + N*(2i - 1)]
				newU[j + newN*(i - 1)] = U[2j + N*(2i - 1)]
			end
		end
		
		newU = MultigridJ(newA, newF, newN, newU, pre, post)

		# println("Coarse to fine")
		for i = 1:newN
			# @show i
			for j = 1:newN
				# @show j
				U[2j + N*(2i - 1)] = newU[j + newN*(i - 1)]
			end
		end

		# println("Post smoothing")
		for i = 1:post
			# @show i
			U = Jacobi(A, F, U)
		end

		# println("fin multigrid, N = ", N, "\n")
		return U
	end
end

function MultigridGS(A, F, N::Int64, x0, pre, post)
	# println("début multigrid, N = ", N)
	if N == 1
		return (1/4)*F
	else
		U = copy(x0)

		# println("Pre smoothing")
		for i = 1:pre
			# @show i
			U = GaussSeidel(A, F, U)
		end

		newN = floor(Int64, N/2)
		newA = Creer_A(newN)
		newF = zeros(newN*newN)
		newU = zeros(newN*newN)

		# println("Fine to coarse")
		for i = 1:newN
			# @show i
			for j = 1:newN
				# @show j
				newF[j + newN*(i - 1)] = F[2j + N*(2i - 1)]
				newU[j + newN*(i - 1)] = U[2j + N*(2i - 1)]
			end
		end
		
		newU = MultigridGS(newA, newF, newN, newU, pre, post)

		# println("Coarse to fine")
		for i = 1:newN
			# @show i
			for j = 1:newN
				# @show j
				U[2j + N*(2i - 1)] = newU[j + newN*(i - 1)]
			end
		end

		# println("Post smoothing")
		for i = 1:post
			# @show i
			U = GaussSeidel(A, F, U)
		end

		# println("fin multigrid, N = ", N, "\n")
		return U
	end
end