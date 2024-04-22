using Distributed
using MPI
using LinearAlgebra
using SparseArrays
using Random

MPI.Init()

function parallel_jacobi(A, b, x0, tol=1e-15, max_iter=10000)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    
    n = size  # 总进程数与矩阵大小相同
    local_n = Int(floor(sqrt(n)))  # 每行或每列的进程数
    @assert local_n^2 == n "The number of processes must be a perfect square."
    
    myidx = rank % local_n  # 当前进程在行方向的索引
    myidy = rank ÷ local_n  # 当前进程在列方向的索引
    
    # 计算每个进程负责的局部矩阵大小
    global_n = size
    local_rows = Int(floor(global_n / local_n))
    local_cols = Int(floor(global_n / local_n))
    
    # 创建局部矩阵和局部向量
    local_A = zeros(local_rows, local_cols)
    local_b = zeros(local_rows)
    local_x0 = zeros(local_rows)
    
    # 使用 MPI 分发全局矩阵 A 和向量 b 给各个进程的局部矩阵和向量
    MPI.Scatter!(comm, A, local_A, root=0)
    MPI.Scatter!(comm, b, local_b, root=0)
    MPI.Bcast!(comm, x0, root=0)
    
    iteration = 0
    norm_diff = Inf
    local_x = copy(local_x0)
    new_local_x = similar(local_x)
    
    while norm_diff > tol && iteration < max_iter
        for i = 1:local_rows
            local_sum = 0.0
            for j = 1:local_cols
                if i != j
                    local_sum += local_A[i, j] * local_x[j]
                end
            end
            new_local_x[i] = (local_b[i] - local_sum) / local_A[i, i]
        end
        
        norm_diff = norm(new_local_x - local_x)
        local_x .= new_local_x
        
        # 全局归约以计算所有进程的差值范数
        global_norm_diff = MPI.Allreduce(norm_diff, MPI.SUM, comm)
        
        if rank == 0
            println("Iteration: $iteration, Norm difference: $global_norm_diff")
        end
        
        iteration += 1
    end
    
    return local_x
end


N = 3
A = Creer_A(N)
F = Creer_F(p, q, N)
x0 = rand(1:9, (N*N, 1))
x0 = zeros(N * N)
result = parallel_jacobi(A, F, x0)
println("Parallel Jacobi result:")
println(result)







function ParallelMultigrid(A, F, N::Int64, x0, pre, post)
    println("start Multigrid,N = ", N)

    if N == 1
        return (1/4) * F
    else
        U = copy(x0)

        println("preprocessing")
        for i = 1:pre
            println("Preprocessing iteration ", i)
            U = parallel_jacobi(A, F, U)
        end

        newN = div(N, 2)
        local_n = div(newN, size)
        local_start = (rank - 1) * local_n + 1
        local_end = rank * local_n

        newA = Creer_A(newN)
        newF = zeros(newN * newN)
        newU = zeros(newN * newN)

        println("fine mesh to coarse mesh")
        for i = local_start:local_end
            for j = 1:newN
                newF[j + newN * (i - 1)] = F[2*j + N * (2*i - 1)]
                newU[j + newN * (i - 1)] = U[2*j + N * (2*i - 1)]
            end
        end

        newU = ParallelMultigrid(newA, newF, newN, newU, pre, post)

        MPI.Allgather(newU, U, comm)

        println("coarse mesh to fine mesh")
        for i = local_start:local_end
            for j = 1:newN
                U[2*j + N * (2*i - 1)] = newU[j + newN * (i - 1)]
            end
        end

        println("Post-processing")
        for i = 1:post
            println("Post-processing iteration ", i)
            U = parallel_jacobi(A, F, U)
        end

        println("end Multigrid N = ", N, "\n")
        return U
    end
end




time1 = @elapsed begin
    N = 3
    A = Creer_A(N)
    F = Creer_F(p, q, N)
    U = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
    x0 = rand(1:9, (N*N, 1))
    x0 = zeros(N * N)
    result = ParallelMultigrid(A, F1, N, x0, 1, 1)
    println("Multigrid:", result)
end
println("Elapsed time: ", time1, " seconds")



time2 = @elapsed begin
    N = 3
    A = Creer_A(N)
    F1 = Creer_F(p, q, N)
    U = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
    x0 = rand(1:9, (N*N, 1))
    x0 = zeros(N * N)
    resMJ = MultigridJ(A, F1, N, x0, 1, 1)
    println("Multigrid original：", resMJ)
end
println("Elapsed time(origine): ", time2, " seconds")

MPI.Barrier(comm)
MPI.Finalize()
