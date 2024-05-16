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

function Jacobi(A, b, x0, tol = 10^(-15), MaxIter )
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

"""
    split_count(N::Integer, n::Integer)

Return a vector of `n` integers which are approximately equally sized and sum to `N`.
"""
function split_count(N::Integer, n::Integer)
    q,r = divrem(N, n)
    return [i <= r ? q+1 : q for i = 1:n]
end

function parallel_jacobi_solver(p, q, N,MaxIter)
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    f = Creer_F(p, q, N)
    f_local = f[rank  * N + 1 : rank * N + N ] 
    #@show length(f)
    #@show length(f_local)
    x_local = ones(N,N)
    x_new = ones(N,1)
    if rank == 0
        output = similar(x_local)
        M_counts = [N for i = 1:size]
        N_counts = split_count(N, size)
        # store sizes in 2 * comm_size Array
        sizes = vcat(M_counts', N_counts')
        # store number of values to send to each rank in comm_size length Vector
        counts = vec(prod(sizes, dims=1))
        output_vbuf = VBuffer(output, counts) # VBuffer for gather
    else
        output_vbuf = VBuffer(nothing)
    end

    for iter = 1:MaxIter
        #@show iter
        x_local_old = copy(x_local)
        #@show length(x_local_old)
        for i = 2:N
            for j = 1:N
                left = j == 1 ? x_local_old[i, 1] : x_local_old[i, j-1]
                right = j == N ? x_local_old[i, N] : x_local_old[i, j+1]
                up = i == N ? x_local_old[N, j] : x_local_old[i+1, j]
                down = i == 2 ? x_local_old[1, j] : x_local_old[i-1, j]
                x_local[i, j] = 0.25 * (f_local[i-1] + left + right + up + down)
                x_new = [x_local[i,i] for i in 1:N]
            end
        end
    end
    for i = 0:size-1
        if rank == i
            #@show i
            #@show rank x_local x_new
            #@show length(x_local)
        end
        MPI.Barrier(comm)
    end

    # Collecter les vecteurs de solution locaux de toutes les unités de traitement en vecteurs de solution globaux
    MPI.Gatherv!(x_new,output_vbuf, 0, comm) 


    if rank == 0 
        println("Final result")
        println("================")
        #@show output_vbuf
        global_solution = output_vbuf.data
        #@show global_solution
        return global_solution
    end
    MPI.Finalize()
    
end





function main(N)
    p = 1
    q = 1
    A = Creer_A(N)
    F = Creer_F(p, q, N)
    x0 = rand(1:9, (N*N, 1))
    x0 = zeros(N * N)
    MaxIter = 1
    # Solve with parallel Jacobi method
    println("Solving with parallel Jacobi...")
    @time x_parallel_jacobi = parallel_jacobi_solver(p, q, N,MaxIter)
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

