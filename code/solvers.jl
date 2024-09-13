using LinearAlgebra
using Random
using SparseArrays

include("problem.jl")

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

function Jacobi(A, b, x0, MaxIter = 1000000, tol = 10^(-15))
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
			# @show i 
			# @show j 
			# @show v 
			# println("")
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

	return nouveau, iteration
end

function GaussSeidel(A, b, x0, MaxIter = 1000000, tol = 10^(-15))
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

	return nouveau, iteration
end

function SOR(A, b, w, x0, MaxIter = 1000000, tol = 10^(-15))
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

	return nouveau, iteration
end

function Prolongation(U, N)
	newN = 2N + 1
	temp = reshape(U, (N, N))
	res = zeros(newN, newN)
	for i = 1:N
		for j = 1:N
			res[2i, 2j] += temp[i, j]
			res[2i - 1, 2j - 1] += temp[i, j]/4
			res[2i - 1, 2j + 1] += temp[i, j]/4
			res[2i + 1, 2j - 1] += temp[i, j]/4
			res[2i + 1, 2j + 1] += temp[i, j]/4
			res[2i, 2j - 1] += temp[i, j]/2
			res[2i, 2j + 1] += temp[i, j]/2
			res[2i - 1, 2j] += temp[i, j]/2
			res[2i + 1, 2j] += temp[i, j]/2
		end
	end
	return res[:], newN
end

function Restriction(U, N)
	newN = floor(Int64, N/2)
	temp = reshape(U, (N, N))
	res = zeros(newN, newN)
	for i = 1:newN
		for j = 1:newN
			res[i, j] += temp[2i, 2j]/4
			res[i, j] += temp[2i, 2j - 1]/8
			res[i, j] += temp[2i, 2j + 1]/8
			res[i, j] += temp[2i - 1, 2j]/8
			res[i, j] += temp[2i + 1, 2j]/8
			res[i, j] += temp[2i - 1, 2j - 1]/16
			res[i, j] += temp[2i - 1, 2j + 1]/16
			res[i, j] += temp[2i + 1, 2j - 1]/16
			res[i, j] += temp[2i + 1, 2j + 1]/16
		end
	end
	return res[:], newN
end

function CycleJ(A, F, x0, N::Int64, pre, post, nb = 1)
	U, _ = Jacobi(A, F, x0, pre)
	r, newN = Restriction(F - A*U, N)
	if newN == 1
		d = (r[1]/A[1])*ones(1)
	else
		tempA = Creer_A(newN)
		tempX = zeros(newN*newN)
		for i = 1:nb
			d = CycleJ(tempA, r, tempX, newN, pre, post, nb)
		end
	end
	d, _ = Prolongation(d, newN)
	U += d
	U, _ = Jacobi(A, F, U, post)
	return U
end

function CycleGS(A, F, x0, N::Int64, pre, post, nb = 1)
	U, _ = GaussSeidel(A, F, x0, pre)
	r, newN = Restriction(F - A*U, N)
	if newN == 1
		d = (r[1]/A[1])*ones(1)
	else
		tempA = Creer_A(newN)
		tempX = zeros(newN*newN)
		for i = 1:nb
			d = CycleGS(tempA, r, tempX, newN, pre, post, nb)
		end
	end
	d, _ = Prolongation(d, newN)
	U += d
	U, _ = GaussSeidel(A, F, U, post)
	return U
end

function MultigridJ(F, N::Int64, pre, post, nbc = 1, nbm = 2)
	tempN = N + 1
	k = -1
	while tempN != 1
		k += 1
		tempN /= 2
	end
	tempN = floor(Int64, tempN)
	tempF = F[N*(2^k - 1) + 2^k]*ones(tempN)
	A = Creer_A(tempN)
	U = (F[1]/A[1])*ones(tempN)
	while tempN != N
		U, tempN = Prolongation(U, tempN)
		k -= 1
		tempF = zeros(tempN*tempN)
		for i = 1:tempN
			for j = 1:tempN
				tempF[j + (i - 1)tempN] = F[N*((2^k)*i - 1) + (2^k)*j]
			end
		end
		A = Creer_A(tempN)
		for i = 1:nbm
			U = CycleJ(A, tempF, U,  tempN, pre, post, nbc)
		end
	end
	return U
end

function MultigridGS(F, N::Int64, pre, post, nbc = 1, nbm = 2)
	tempN = N + 1
	k = -1
	while tempN != 1
		k += 1
		tempN /= 2
	end
	tempN = floor(Int64, tempN)
	tempF = F[N*(2^k - 1) + 2^k]*ones(tempN)
	A = Creer_A(tempN)
	U = (F[1]/A[1])*ones(tempN)
	while tempN != N
		U, tempN = Prolongation(U, tempN)
		k -= 1
		tempF = zeros(tempN*tempN)
		for i = 1:tempN
			for j = 1:tempN
				tempF[j + (i - 1)tempN] = F[N*((2^k)*i - 1) + (2^k)*j]
			end
		end
		A = Creer_A(tempN)
		for i = 1:nbm
			U = CycleGS(A, tempF, U, tempN, pre, post, nbc)
		end
	end
	return U
end