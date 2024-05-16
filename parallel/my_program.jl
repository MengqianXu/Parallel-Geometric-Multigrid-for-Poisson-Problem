using LinearAlgebra
using SparseArrays
using Random
using MPI



function u(p, q, x, y)
	function g(x, y)
		return sin(p*pi*x)*sin(q*pi*y)
	end
	return g(x, y)
end
p = rand(-5:5) 
q = rand(-5:5)

function Creer_A(N)
	return sparse(kron(I(N), Tridiagonal(-ones(N, N)) + 5I(N)) + kron(Tridiagonal(-ones(N, N)) + I(N), I(N)))
end

function Creer_F(p, q, N)
	F = zeros(N*N)
	for i = 1:N
		x = i/(N + 1)
		for j = 1:N
			y = j/(N + 1)
			F[N*(i - 1) + j] = u(p, q, x, y)*((pi^2)*((p^2) + (q^2)))
		end
	end
	return ((1/(N + 1))^2)F
end


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
			#@show i 
			#@show j 
			#@show v 
			println("")
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

function find_factor(n)
    for i = 2:n
        if n % i == 0
            return i
        end
    end
    return n  
end

function isprime(n)
    if n <= 1
        return false
    end
    if n == 2
        return true
    end
    if n % 2 == 0
        return false
    end
    for i = 3:2:isqrt(n)
        if n % i == 0
            return false
        end
    end
    return true
end

function calculate_optimal_process(N)
    is_prime = isprime(N)
    if !is_prime
        factor = find_factor(N)
        if factor == 0
            error("Error")
        end
        npx = factor
        npy = N ÷ factor
    else
        npx = N
        npy = 1  
    end
    return npx, npy
end


function calculate_subdomain(rank, npx, npy, nx, ny)
    myidx = rank % npx
    myidy = div(rank, npx)
    
    rx = (nx - 1) % npx
    ry = (ny - 1) % npy
    
    nlx = ifelse(myidx < rx, div(nx - 1, npx) + 2, div(nx - 1, npx) + 1)
    nly = ifelse(myidy < ry, div(ny - 1, npy) + 2, div(ny - 1, npy) + 1)
    
    x0 = myidx * div(nx - 1, npx) + min(myidx, rx)+1
    y0 = myidy * div(ny - 1, npy) + min(myidy, ry)+1
    
    println("Rank $(rank): x0 = $(x0), y0 = $(y0), nlx = $(nlx), nly = $(nly)")
    return x0, y0, nlx, nly
end





function parallel_jacobi(A, b, N, npx, npy, tol=1e-15, max_iter=10000)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    #println("Size",size)
    x0, y0, nlx, nly = calculate_subdomain(rank, npx, npy, N, N)
    
    @assert 1 <= x0 <= N*N && 1 <= y0 <= N*N && (x0 + nlx - 1) <= N*N && (y0 + nly - 1) <= N*N
    
    local_A = A[x0:x0+nlx-1, y0:y0+nly-1]
    local_b = b[x0:x0+nlx-1]
    
    local_u = zeros(nlx, nly)  # Local solution vector, including boundaries
    new_local_u = similar(local_u)  # Used to store the updated local solution vector
    println("Initial local_u:")
    println(local_u)
    iteration = 0
    norm_diff = Inf


    while norm_diff > tol && iteration < max_iter
        new_local_u .= local_u
        # Perform Jacobi iteration over interior points
        for i = 2:nlx-1
            for j = 2:nly-1
                sum = local_A[i, j-1] * local_u[i, j-1] +
                      local_A[i, j+1] * local_u[i, j+1] +
                      local_A[i-1, j] * local_u[i-1, j] +
                      local_A[i+1, j] * local_u[i+1, j]
                new_local_u[i, j] = (local_b[i] - sum) / local_A[i, i]
            end
        end
        
        # Update local solution vector
        local_u .= new_local_u
        
        # Exchange boundary data
        # left
        if x0 > 1 
            dest_rank = mod(rank - 1, size)
            left_boundary = local_u[2:end-1, 1]
            println("Sending left boundary from rank $rank to rank $dest_rank: $left_boundary")
            MPI.Send(left_boundary, dest_rank, 0, comm)
            MPI.Recv!(local_u[1, :], comm, source=dest_rank, tag=0)
            println("Received left boundary at rank $rank: $(local_u[1, :])")
        end
        
        # right boundary 
        if x0 + nlx - 1 < N
            dest_rank = mod(rank + 1, size)
            right_boundary = local_u[2:end-1, end]
            println("Sending right boundary from rank $rank to rank $dest_rank: $right_boundary")
            MPI.Send(right_boundary, dest_rank, 0, comm)
            MPI.Recv!(local_u[end, :], comm, source=dest_rank, tag=0)
            println("Received right boundary at rank $rank: $(local_u[end, :])")
        end
        
        # top boundary exchange
        if y0 > 1
            dest_rank = mod(rank - npx, size)
            top_boundary = local_u[1, 2:end-1]
            println("Sending top boundary from rank $rank to rank $dest_rank: $top_boundary")
            MPI.Send(top_boundary, dest_rank, 0, comm)
            MPI.Recv!(local_u[:, 1], comm, source=dest_rank, tag=0)
            println("Received top boundary at rank $rank: $(local_u[:, 1])")
        end
        
        # bottom boundary exchange
        if y0 + nly - 1 < N
            dest_rank = mod(rank + npx, size)
            bottom_boundary = local_u[end, 2:end-1]
            println("Sending bottom boundary from rank $rank to rank $dest_rank: $bottom_boundary")
            MPI.Send(bottom_boundary, dest_rank, 0, comm)
            MPI.Recv!(local_u[:, end], comm, source=dest_rank, tag=0)
            println("Received bottom boundary at rank $rank: $(local_u[:, end])")
        end

        local_u .= new_local_u
        println("Local solution at iteration $iteration:")
        println(local_u)

        norm_diff = norm(new_local_u - local_u)
        sendbuf = [norm_diff]
        recvbuf = similar(sendbuf)
        MPI.Allreduce!(sendbuf, recvbuf, MPI.MIN, comm)
        iteration += 1
    end
    
    # Gather local solutions to the root process
    gathered_local_u = MPI.Gather(local_u[2:end-1, 2:end-1], 0, comm)
    
    if rank == 0
        result = vcat(gathered_local_u...)
        return result
    else
        return nothing
    end
end


function benchmark_jacobi_methods(N)
    MPI.Init()
    A = Creer_A(N)
    F = Creer_F(p, q, N)
    x0 = rand(1:9, (N*N, 1))
    x0 = zeros(N * N)
    #println("Solving with non-parallel Jacobi...")
    #@time x_jacobi = Jacobi(A, F, x0)
    #println("Final result vector size: $(length(x_jacobi))")
    #println("non-parallel Jacobi",x_jacobi)
    
    # Solve with parallel Jacobi method
    # npx =  Number of processes in x-direction
    # npy =  Number of processes in y-direction
    #npx,npy = calculate_optimal_process(N)
    npx = 3
    npy = 2
    println("Solving with parallel Jacobi...")
    @time x_parallel_jacobi = parallel_jacobi(A, F, N, npx, npy)
    println("Final result vector size: $(length(x_parallel_jacobi))")
    #println("parallel Jacobi",x_parallel_jacobi)
    # Compare solutions
    #println("difference between solutions:", maximum(abs.(x_jacobi - x_parallel_jacobi)))
    MPI.Finalize()
end

N = 12
benchmark_jacobi_methods(N)



