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
    

    local_A = A[x0:x0+nlx-1, y0:y0+nly-1]
    local_b = b[x0:x0+nlx-1]
    
    local_u = zeros(nlx * nly)  # local solution vector, including boundaries
    new_local_u = similar(local_u)  # Used to store the updated local solution vector
    
    iteration = 0
    norm_diff = Inf
    
    while norm_diff > tol && iteration < max_iter
        # Update internal points
        for idx = 1:nlx * nly
            i = div(idx - 1, nly) + 1
            j = mod(idx - 1, nly) + 1
            
            if i > 1 && i < nlx && j > 1 && j < nly
                #internal point update
                local_sum = local_A[i, j-1] * local_u[idx-1] +
                            local_A[i, j+1] * local_u[idx+1] +
                            local_A[i-1, j] * local_u[idx-nly] +
                            local_A[i+1, j] * local_u[idx+nly]
                new_local_u[idx] = (local_b[i] - local_sum) / local_A[i, j]
            else
                # Boundary points retain their current values first
                new_local_u[idx] = local_u[idx]
            end
        end
        
        # Update boundary points (outermost circle)
        for j = 1:nly
            if j > 1 && j < nly
                # lower
                idx = 1 + (j - 1)
                new_local_u[idx] = local_b[1] / local_A[1, j]
                
                # upper
                idx = nlx * nly - (nly - j)
                new_local_u[idx] = local_b[nlx] / local_A[nlx, j]
            end
        end
        
        for i = 2:nlx-1
            # left
            idx = (i - 1) * nly + 1
            new_local_u[idx] = local_b[i] / local_A[i, 1]
            
            # right
            idx = i * nly
            new_local_u[idx] = local_b[i] / local_A[i, nly]
        end
        
        iteration += 1
    end
    

    gathered_local_u = MPI.Gather(local_u, 0, comm)
    
    if rank == 0
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

