using LinearAlgebra
using Random
using SparseArrays
using MPI

include("problem.jl")

function decomposition_2_facteurs(n)
	i = floor(Int64, sqrt(n))
	while n%i != 0
		i -= 1
	end
	j = floor(Int64, n/i)
	if i > j
		return i, j
	else
		return j, i
	end
end

function Jacobi_MPI(F, U0, N, tol = 10^(-15), MaxIter = 10000)
	# Initialisation des paramètres
	# MPI.Init()
	comm = MPI.COMM_WORLD
	rank = MPI.Comm_rank(comm)
	size = MPI.Comm_size(comm)
	a, b = decomposition_2_facteurs(size)
	x = ceil(Int64, N/a)
	y = ceil(Int64, N/b)
	d = rank%a
	c = floor(Int64, (rank - d)/a)
	h = 1/(N + 1)


	# Initialisation des variables utilisées
	tempU = zeros(1)
	Ulocal = zeros(x*y)
	Flocal = zeros(x*y)

	if rank == 0
		tempU0 = reshape(U0, (N, N))
		tempF0 = reshape(F, (N, N))
		tempU = NaN64*ones(a*x, b*y)
		tempF = NaN64*ones(a*x, b*y)
		Ures = zeros(N, N)
		pos = 1

		for m = 0:(a - 1)
			for n = 0:(b - 1)
				for i = 1:x
					colonne = m*x + i
					for j = 1:y
						ligne = n*y + j
						if (ligne <= N) && (colonne <= N)
							tempU[pos] = tempU0[ligne, colonne]
							tempF[pos] = tempF0[ligne, colonne]
						end
						pos += 1
					end
				end
			end
		end
		
		tempU = tempU[:]
		tempF = tempF[:]
		MPI.Scatter!(tempU, Ulocal, 0, comm)
		MPI.Scatter!(tempF, Flocal, 0, comm)
	end

	Ulocal = reshape(Ulocal, (y, x))
	Flocal = reshape(Flocal, (y, x))
	@show Flocal, rank
	Ulocal_new = zeros(y, x)
	requeteS = MPI.MultiRequest(4)
	requeteR = MPI.MultiRequest(4)
	compteur = 1

	if c > 0
		tabgaucheS = Ulocal[:, 1]
		tabgaucheR = zeros(y)
		MPI.Isend(tabgaucheS, rank - a, 0, comm, requeteS[compteur])
		MPI.Irecv!(tabgaucheR, rank - a, 0, comm, requeteR[compteur])
		compteur += 1
	end

	if c < (a - 1)
		tabdroiteS = Ulocal[:, x]
		tabdroiteR = zeros(y)
		MPI.Isend(tabdroiteS, rank + a, 0, comm, requeteS[compteur])
		MPI.Irecv!(tabdroiteR, rank + a, 0, comm, requeteR[compteur])
		compteur += 1
	end

	if d > 0
		tabhautS = Ulocal[1, :]
		tabhautR = zeros(x)
		MPI.Isend(tabhautS, rank - 1, 0, comm, requeteS[compteur])
		MPI.Irecv!(tabhautR, rank - 1, 0, comm, requeteR[compteur])
		compteur += 1
	end

	if d < (b - 1)
		tabbasS = Ulocal[y, :]
		tabbasR = zeros(x)
		MPI.Isend(tabbasS, rank + 1, 0, comm, requeteS[compteur])
		MPI.Irecv!(tabbasR, rank + 1, 0, comm, requeteR[compteur])
		compteur += 1
	end

	MPI.Waitall(requeteS)
	MPI.Waitall(requeteR)	


	# Première itération de Jacobi
	norme_carre_local = 0
	for i = 1:x
		colonne = c*x + i 
		for j = 1:y
			ligne = d*y + j
			if (ligne <= N) && (colonne <= N)
				res = (h^2)Flocal[j, i]
				if i > 1
					res += Ulocal[j, i - 1]
				else
					if c > 0
						res += tabgaucheR[j]
					end
				end

				if (colonne + 1) <= N
					if i < x
						res += Ulocal[j, i + 1]
					else
						if c < (a - 1)
							res += tabdroiteR[j]
						end
					end
				end

				if j > 1
					res += Ulocal[j - 1, i]
				else
					if d > 0
						res += tabhautR[i]
					end
				end

				if (ligne + 1) <= N
					if j < y
						res += Ulocal[j + 1, i]
					else
						if d < (b - 1)
							res += tabbasR[i]
						end
					end
				end
				Ulocal_new[j, i] = res/4
				norme_carre_local += (Ulocal_new[j, i] - Ulocal[j, i])^2
			end
		end
	end
	
	Ulocal = Ulocal[:]
	Ulocal_new = Ulocal_new[:]
	norme_carre = MPI.Allreduce(norme_carre_local, +, comm)
	norme = sqrt(norme_carre)
	iteration = 1
	@show norme

	# Iterations de Jacobi
	while (norme >= tol) && (iteration <= MaxIter)
		@show rank, iteration
		Ulocal = copy(reshape(Ulocal_new, (y, x)))
		Ulocal_new = reshape(Ulocal_new, (y, x))
		compteur = 1
		
		if c > 0
			tabgaucheS = Ulocal[:, 1]
			MPI.Isend(tabgaucheS, rank - a, 0, comm, requeteS[compteur])
			MPI.Irecv!(tabgaucheR, rank - a, 0, comm, requeteR[compteur])
			compteur += 1
		end
	
		if c < (a - 1)
			tabdroiteS = Ulocal[:, x]
			MPI.Isend(tabdroiteS, rank + a, 0, comm, requeteS[compteur])
			MPI.Irecv!(tabdroiteR, rank + a, 0, comm, requeteR[compteur])
			compteur += 1
		end
	
		if d > 0
			tabhautS = Ulocal[1, :]
			MPI.Isend(tabhautS, rank - 1, 0, comm, requeteS[compteur])
			MPI.Irecv!(tabhautR, rank - 1, 0, comm, requeteR[compteur])
			compteur += 1
		end
	
		if d < (b - 1)
			tabbasS = Ulocal[y, :]
			MPI.Isend(tabbasS, rank + 1, 0, comm, requeteS[compteur])
			MPI.Irecv!(tabbasR, rank + 1, 0, comm, requeteR[compteur])
			compteur += 1
		end
	
		MPI.Waitall(requeteS)
		MPI.Waitall(requeteR)

		norme_carre_local = 0
		for i = 1:x
			colonne = c*x + i 
			for j = 1:y
				ligne = d*y + j
				if (ligne <= N) && (colonne <= N)
					res = (h^2)Flocal[j, i]
					if i > 1
						res += Ulocal[j, i - 1]
					else
						if c > 0
							res += tabgaucheR[j]
						end
					end
	
					if (colonne + 1) <= N
						if i < x
							res += Ulocal[j, i + 1]
						else
							if c < (a - 1)
								res += tabdroiteR[j]
							end
						end
					end
	
					if j > 1
						res += Ulocal[j - 1, i]
					else
						if d > 0
							res += tabhautR[i]
						end
					end
	
					if (ligne + 1) <= N
						if j < y
							res += Ulocal[j + 1, i]
						else
							if d < (b - 1)
								res += tabbasR[i]
							end
						end
					end
					Ulocal_new[j, i] = res/4
					norme_carre_local += (Ulocal_new[j, i] - Ulocal[j, i])^2
				end
			end
		end
		
		Ulocal = Ulocal[:]
		Ulocal_new = Ulocal_new[:]
		norme_carre = MPI.Allreduce(norme_carre_local, +, comm)
		norme = sqrt(norme_carre)
		iteration += 1
	end
	@show norme

	MPI.Gather!(Ulocal_new, tempU, comm)
	# MPI.Finalize()

	if rank == 0
		pos = 1
		tempU = reshape(tempU, (a*x, b*y))
		for m = 0:(a - 1)
			for n = 0:(b - 1)
				for i = 1:x
					colonne = m*x + i 
					for j = 1:y
						ligne = n*y + j
						if (ligne <= N) && (colonne <= N)
							Ures[ligne, colonne] = tempU[pos]
						end
						pos += 1
					end
				end
			end
		end
		return Ures[:]
	end
end