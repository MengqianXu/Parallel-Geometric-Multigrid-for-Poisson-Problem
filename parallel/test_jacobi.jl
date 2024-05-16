using LinearAlgebra
using SparseArrays
using Random
using MPI

function u(p, q, x, y)
    return sin(p * π * x) * sin(q * π * y)
end

function Creer_A(N)
	return sparse(kron(I(N), Tridiagonal(-ones(N, N)) + 5I(N)) + kron(Tridiagonal(-ones(N, N)) + I(N), I(N)))
end

function Creer_F(p, q, N)
    f = zeros(N*N)
    dx = 1.0 / (N + 1)
    dy = 1.0 / (N + 1)
    for i = 1:N
        x = i * dx
        for j = 1:N
            y = j * dy
            f[N*(i - 1) + j] = u(p, q, x, y) * ((π^2) * (p^2 + q^2))
        end
    end
    return ((1 / (N + 1))^2) * f
end
function Creer_F(p, q, N)
    f = zeros(N*N)
    dx = 1.0 / (N + 1)
    dy = 1.0 / (N + 1)
    for i = 1:N
        x = i * dx
        for j = 1:N
            y = j * dy
            f[N*(i - 1) + j] = u(p, q, x, y) * ((π^2) * (p^2 + q^2))
        end
    end
    return ((1 / (N + 1))^2) * f
end

function Jacobi(A, b, x0, MaxIter, tol = 10^(-15))
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
			#@show i 
			#@show j 
			#@show v 
			#println("")
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

function decomposition(n)
    i::Int64 = floor(sqrt(n))
    while n % i != 0
        i -= 1
    end
    j::Int64 = floor(n / i)
    if i > j
        a = i
        b = j
    else
        a = j
        b = i
    end

    return a,b
end


function parallel_jacobi_solver(F,U0,N::Int64,MaxIter,tol=1e-15)
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    
    a, b = decomposition(size)
    x = ceil(Int64, N / a)
    y = ceil(Int64, N / b)
    c = rand() % a
    d = floor(Int64, (rank - c) / a)

    U_local = zeros(x*y)
    F_local = zeros(x*y)
    if rank == 0
        temp_U0 = reshape(U0,(N,N))
        temp_F0 = reshape(F,(N,N))
        temp_U = NaN64*ones(a*x, b*y)
		temp_F = NaN64*ones(a*x, b*y)
        pos = 1
        for m = 0:(a - 1)
			for n = 0:(b - 1)
				for i = 1:x
					colonne = m*x + i
					for j = 1:y
						ligne = n*y + j
						if (ligne <= N) && (colonne <= N)
							temp_U[pos] = temp_U0[ligne, colonne]
							temp_F[pos] = temp_F0[ligne, colonne]
						end
						pos += 1
					end
				end
			end
		end

        temp_U = temp_U[:]
        temp_F = tempF[:]
        Scatterv!(temp_U,U_local,comm)
        Scatterv!(temp_F,F_local,comm)
    end 
    U_local = reshape(U_local, (y, x))
	F_local = reshape(F_local, (y, x))
    U_new = zeros(y, x)
    iteration = 0

    requeteS = MPI.MultiRequest(4)
	requeteR = MPI.MultiRequest(4)


    nb = 1
    while iteration < MaxIter
        iteration += 1
        if d < a - 1
            rank_droit = rank + a
            if rank_droit < size
                TabDroit_R = zeros(y)
                MPI.Isend(U_local[:, x] , comm, dest = rank_droit,requeteS[nb])
                MPI.Irecv!(TabDroit_R, comm, source=rank_droit, requeteS[nb])
                U_local[:,x]  = TabDroit_R 
                nb +=1
            end
            
        end
        
        if d > 0
            rank_gauch = rank - a
            if rank_gauch >= 0
                TabGauch_R = zeros(y)
                MPI.Isend(U_local[:,1] , comm, dest = rank_gauch,requeteS[nb])
                MPI.Irecv!(TabGauch_R, comm, source=rank_gauch, requeteS[nb])
                U_local[:,1]  = TabGauch_R 
                nb +=1
            end
        end
        
        if c < b - 1
            rank_bas = rank + 1
            if rank_bas < size
                TabBas_R = zeros(x)
                MPI.Isend(U_local[y,:], comm, dest = rank_bas,requeteS[nb])
                MPI.Irecv!(TabBas_R, comm, source=rank_bas, requeteS[nb])
                U_local[y,:]  = TabBas_R 
                nb +=1
            end
        end
        
        if c > 0
            rank_haut = rank - 1
            if rank_haut >= 0
                TabHaut_R = zeros(x)
                MPI.Isend(U_local[1,:], comm, dest = rank_haut,requeteS[nb])
                MPI.Irecv!(TabHaut_R,comm, source = rank_haut, requeteS[nb])
                U_local[1,:] = TabHaut_R
                nb +=1
            end
        end
    
        #for i = 1:x
            #for j = 1:y
                #idx = (i - 1) * y + j
                #U_new[idx] = (F_local[idx] - sum_neighbours(U_old, i, j, x, y)) * 1/4
            #end
        #end

        norme_local = 0.0
        for i = 1:x
            colonne = d*x + i 
            for j = 1:y
                ligne = c*y + j
                if (ligne <= N) && (colonne <= N)
                    sum = (1/(N + 1)^2) * F_local[j, i]
                    if i > 1
                        sum += U_local[j, i - 1]
                    else
                        if d > 0
                            sum += TabGauch_R[j]
                        end
                    end
    
                    if (colonne + 1) <= N
                        if i < x
                            sum += U_local[j, i + 1]
                        else
                            if d < (a - 1)
                                sum += TabDroit_R[j]
                            end
                        end
                    end
    
                    if j > 1
                        sum += U_local[j - 1, i]
                    else
                        if c > 0
                            sum += TabHaut_R[i]
                        end
                    end
    
                    if (ligne + 1) <= N
                        if j < y
                            sum += U_local[j + 1, i]
                        else
                            if c < (b - 1)
                                sum += TabBas_R[i]
                            end
                        end
                    end
                    U_new[j, i] = sum / 4
                    norme_local += norm(U_new[j, i] - U_local[j, i])
                end
            end
        end
        norm = MPI.Reduce(norme_local, MPI.SUM, comm)
        if norm < tol
            break
        end
        
	    MPI.Waitall(requeteS)
	    MPI.Waitall(requeteR)
        U_local .= U_new
    end
    
    # Gather U_local from all processes to U
    U_local = U_local[:]
    temp_U = zeros(a*x*b*y)
    MPI.Gather!(U_local, temp_U, 0, comm)
    
    U = zeros(N, N)
    if rank == 0
		pos = 1
		temp_U = reshape(temp_U, (a*x, b*y))
		for m = 0:(a - 1)
			for n = 0:(b - 1)
				for i = 1:x
					colonne = m*x + i 
					for j = 1:y
						ligne = n*y + j
						if (ligne <= N) && (colonne <= N)
							U[ligne, colonne] = temp_U[pos]
						end
						pos += 1
					end
				end
			end
		end
		return U[:]
	end
    
    MPI.Finalize()
    return U
end

# Calculate the sum of neighbours of U(i, j)
function sum_neighbours(U, i, j, x, y)
    sum = 0.0
    if i > 1
        sum += U[(i - 2) * y + j]
    end
    if i < x
        sum += U[i * y + j]
    end
    if j > 1
        sum += U[(i - 1) * y + j - 1]
    end
    if j < y
        sum += U[(i - 1) * y + j + 1]
    end
    return sum
end
    





function main(N)
    p = 1
    q = 1
    A = Creer_A(N)
    F = Creer_F(p, q, N)
    x0 = rand(1:9, (N*N, 1))
    x0 = ones(N * N)
    MaxIter = 1
    # Solve with parallel Jacobi method
    println("Solving with parallel Jacobi...")
    @time x_parallel_jacobi = parallel_jacobi_solver(F, x0, N, MaxIter)
    x_parallel = reshape(x_parallel_jacobi, (N*N, 1))

    println("Solving with non-parallel Jacobi...")
    @time x_jacobi = Jacobi(A, F, x0,MaxIter)
    #println("Final result vector size: $(length(x_jacobi))")

    println("non-parallel Jacobi",x_jacobi)
    println("parallel Jacobi",x_parallel)
    
    # Compare solutions
    println("difference between solutions:", maximum(abs.(x_jacobi - x_parallel)))
end

N = 7
main(N)