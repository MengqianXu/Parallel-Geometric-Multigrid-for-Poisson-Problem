using MPI


function calculate_subdomain(rank, npx, npy, N)
    myidx = rank % npx
    myidy = div(rank, npx)
    
    rx = (N - 1) % npx
    ry = (N - 1) % npy
    
    nlx = ifelse(myidx < rx, (N - 1) ÷ npx + 2, (N - 1) ÷ npx + 1)
    nly = ifelse(myidy < ry, (N - 1) ÷ npy + 2, (N - 1) ÷ npy + 1)
    
    x0 = myidx * ((N - 1) ÷ npx) + min(myidx, rx) + 1
    y0 = myidy * ((N - 1) ÷ npy) + min(myidy, ry) + 1
    
    return x0, y0, nlx, nly
end


function parallel_jacobi(A, b, N, npx, npy, tol=1e-15, max_iter=10000)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    
 
    x0, y0, nlx, nly = calculate_subdomain(rank, npx, npy, N)
    @assert 1 <= x0 <= N && 1 <= y0 <= N && (x0 + nlx - 1) <= N && (y0 + nly - 1) <= N
    
    # Extract local matrix A and vector b
    local_A = A[x0:x0+nlx-1, y0:y0+nly-1]
    local_b = b[x0:x0+nlx-1]
    
    local_u = zeros(nlx, nly)  # local solution vector, including boundaries
    new_local_u = similar(local_u)  # Used to store the updated local solution vector
    
    iteration = 0
    norm_diff = Inf
    
    while norm_diff > tol && iteration < max_iter
       
        new_local_u .= local_u
        for i = 2:nlx-1
            for j = 2:nly-1
                sum = local_A[i, j-1] * local_u[i, j-1] +
                      local_A[i, j+1] * local_u[i, j+1] +
                      local_A[i-1, j] * local_u[i-1, j] +
                      local_A[i+1, j] * local_u[i+1, j]
                new_local_u[i, j] = (local_b[i] - sum) / local_A[i, i]
            end
        end
        
        #Update the current local solution vector
        local_u .= new_local_u
        
        # Send and receive boundary data, blocking communication
        if x0 > 1
            MPI.Sendrecv(local_u[2, :], rank-1, 0,
                         local_u[1, :], rank-1, 0, comm)
        end
        if x0 + nlx - 1 < N
            MPI.Sendrecv(local_u[nlx-1, :], rank+1, 0,
                         local_u[nlx, :], rank+1, 0, comm)
        end
        if y0 > 1
            left_boundary = local_u[:, 2]
            right_boundary = similar(left_boundary)
            MPI.Sendrecv(left_boundary, rank-npx, 0,
                         right_boundary, rank-npx, 0, comm)
            local_u[:, 1] .= right_boundary
        end
        if y0 + nly - 1 < N
            right_boundary = local_u[:, nly-1]
            left_boundary = similar(right_boundary)
            MPI.Sendrecv(right_boundary, rank+npx, 0,
                         left_boundary, rank+npx, 0, comm)
            local_u[:, nly] .= left_boundary
        end
        
        
        iteration += 1
    end
    
    # Collect local solution vectors of all processes to the main process
    gathered_local_u = MPI.Gather(local_u[2:end-1, 2:end-1], 0, comm)
    
    if rank == 0
        # The main process splices all local solution vectors to obtain the final result vector
        result = vcat(gathered_local_u...)
        return result
    else
        return nothing
    end
end



function benchmark_jacobi_methods(N)
 
    A = Creer_A(N)
    F = Creer_F(p, q, N)
    x0 = rand(1:9, (N*N, 1))
    x0 = zeros(N * N)
    println("Solving with non-parallel Jacobi...")
    @time x_jacobi = Jacobi(A, F, x0)
    println("Final result vector size: $(length(x_jacobi))")
    #println("non-parallel Jacobi",x_jacobi)
    
    # Solve with parallel Jacobi method
    npx = 4  # Number of processes in x-direction
    npy = 2  # Number of processes in y-direction
  
    
    println("Solving with parallel Jacobi...")
    @time x_parallel_jacobi = parallel_jacobi(A, b, N, npx, npy)
    println("Final result vector size: $(length(x_parallel_jacobi))")
    println("parallel Jacobi",x_parallel_jacobi)
    # Compare solutions
    println("difference between solutions:", maximum(abs.(x_jacobi - x_parallel_jacobi)))
end


N = 10
benchmark_jacobi_methods(N)

