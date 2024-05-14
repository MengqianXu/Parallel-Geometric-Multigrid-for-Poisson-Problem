using MPI

MPI.Init()

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
    println("Size",size)
    x0, y0, nlx, nly = calculate_subdomain(rank, npx, npy, N)
    
    @assert 1 <= x0 <= N && 1 <= y0 <= N && (x0 + nlx - 1) <= N && (y0 + nly - 1) <= N
    
    local_A = A[x0:x0+nlx-1, y0:y0+nly-1]
    local_b = b[x0:x0+nlx-1]
    
    local_u = zeros(nlx, nly)  # Local solution vector, including boundaries
    new_local_u = similar(local_u)  # Used to store the updated local solution vector
    
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
                new_local_u[i, j] = (local_b[i] - sum) / local_A[i, j]
            end
        end
        
        # Update local solution vector
        local_u .= new_local_u
        
        # Exchange boundary data
        if x0 > 1
            dest_rank = rank - 1
            left_boundary = local_u[2:end-1, 1]
            MPI.send(left_boundary, dest_rank, 0, comm)
            MPI.recv!(local_u[1, 1:end], dest_rank, 0, comm)
        end
        
        if x0 + nlx - 1 < N
            dest_rank = rank + 1
            right_boundary = local_u[2:end-1, end]
            MPI.send(right_boundary, dest_rank, 0, comm)
            MPI.recv!(local_u[end, 1:end], dest_rank, 0, comm)
        end
        
        if y0 > 1
            dest_rank = rank - npx
            top_boundary = local_u[1, 2:end-1]
            MPI.send(top_boundary, dest_rank, 0, comm)
            MPI.recv!(local_u[1:end, 1], dest_rank, 0, comm)
        end
        
        if y0 + nly - 1 < N
            dest_rank = rank + npx
            bottom_boundary = local_u[end, 2:end-1]
            MPI.send(bottom_boundary, dest_rank, 0, comm)
            MPI.recv!(local_u[1:end, end], dest_rank, 0, comm)
        end
        
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
 
    A = Creer_A(N)
    F = Creer_F(p, q, N)
    x0 = rand(1:9, (N*N, 1))
    x0 = zeros(N * N)
    println("Solving with non-parallel Jacobi...")
    @time x_jacobi = Jacobi(A, F, x0)
    println("Final result vector size: $(length(x_jacobi))")
    #println("non-parallel Jacobi",x_jacobi)
    
    # Solve with parallel Jacobi method
    # npx =  Number of processes in x-direction
    # npy =  Number of processes in y-direction
    npx,npy = calculate_optimal_process(N)
    
    println("Solving with parallel Jacobi...")
    @time x_parallel_jacobi = parallel_jacobi(A, F, N, npx, npy)
    println("Final result vector size: $(length(x_parallel_jacobi))")
    println("parallel Jacobi",x_parallel_jacobi)
    # Compare solutions
    #println("difference between solutions:", maximum(abs.(x_jacobi - x_parallel_jacobi)))
end

N = 7
benchmark_jacobi_methods(N)
MPI.Finalize()


