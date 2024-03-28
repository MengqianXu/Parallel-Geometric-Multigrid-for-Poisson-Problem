using IterativeSolvers
using LinearAlgebra
using MPI
using SparseArrays
using Random
using Plots
using Distributed
using Base.Threads

function Multigrid(A, F, N::Int64, x0, pre, post)
	if N == 1
		return (1/4)*F
	else
		U = copy(x0)

		for i = 1:pre
			U = Jacobi(A, F, U)
		end

		newN = floor(Int64, N/2)
		newA = Creer_A(newN)
		newF = zeros(newN*newN)
		newU = zeros(newN*newN)

		for i = 1:newN
			for j = 1:newN
				newF[j + N*(i - 1)] = F[2j + N*(2i - 1)]
				newU[j + N*(i - 1)] = U[2j + N*(2i - 1)]
			end
		end
		
		newU = Multigrid(newA, newF, newN, newU, pre, post)

		for i = 1:newN
			for j = 1:newN
				U[2j + N*(2i - 1)] = newU[j + N*(i - 1)]
			end
		end

		for i = 1:pre
			U = Jacobi(A, F, U)
		end

		return U
	end
end


function Multigrid_distributed(A, F, N::Int64, x0, pre, post)
    if N == 1
        return (1/4)*F
    else
        U = copy(x0)
        @sync begin
            @distributed for i = 1:pre
                U = Jacobi(A, F, U)
            end
        end
        newN = floor(Int64, N/2)
		newA = Creer_A(newN)
		newF = zeros(newN*newN)
		newU = zeros(newN*newN)

		for i = 1:newN
			for j = 1:newN
				newF[j + N*(i - 1)] = F[2j + N*(2i - 1)]
				newU[j + N*(i - 1)] = U[2j + N*(2i - 1)]
			end
		end
		
		newU = Multigrid_distributed(newA, newF, newN, newU, pre, post)

		for i = 1:newN
			for j = 1:newN
				U[2j + N*(2i - 1)] = newU[j + N*(i - 1)]
			end
		end
        @sync begin
            @distributed for i = 1:pre
                U = Jacobi(A, F, U)
            end
        end
        return U
    end
end




function Multigrid_threads(A, F, N::Int64, x0, pre, post)
    if N == 1
        return (1/4)*F
    else
        U = copy(x0)
        Threads.@threads for i = 1:pre
            U = Jacobi(A, F, U)
        end
        newN = floor(Int64, N/2)
        newA = Creer_A(newN)
        newF = zeros(newN*newN)
        newU = zeros(newN*newN)

        for i = 1:newN
            for j = 1:newN
                newF[j + N*(i - 1)] = F[2j + N*(2i - 1)]
                newU[j + N*(i - 1)] = U[2j + N*(2i - 1)]
            end
        end

        newU = Multigrid_threads(newA, newF, newN, newU, pre, post)

        for i = 1:newN
            for j = 1:newN
                U[2j + N*(2i - 1)] = newU[j + N*(i - 1)]
            end
        end

        Threads.@threads for i = 1:pre
            U = Jacobi(A, F, U)
        end

        return U
    end
end




time1 = @elapsed begin
    N = 3
    A = Creer_A(N)
    F1 = Creer_F(p, q, N)
    U1 = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
    x0 = rand(1:9, (N*N, 1))
    res_origine = Multigrid(A, F1, N, x0, 1, 1)
    #println(res_origine)
end
println("Elapsed time(origine): ", time1, " seconds")


time2 = @elapsed begin
    N = 3
    A = Creer_A(N)
    F1 = Creer_F(p, q, N)
    U1 = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
    x0 = rand(1:9, (N*N, 1))
    res = Multigrid_distributed(A, F1, N, x0, 1, 1)
    #println(res)
end
println("Elapsed time(distributed): ", time2, " seconds")


time_threads = @elapsed begin
    N = 3
    A = Creer_A(N)
    F1 = Creer_F(p, q, N)
    U1 = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
    x0 = rand(1:9, (N*N, 1))
    res_threads = Multigrid_threads(A, F1, N, x0, 1, 1)
    #println(res_threads)
end

println("Elapsed time (threads): ", time_threads, " seconds")


MPI.Init()

function Multigrid_MPI(A, F, N::Int64, x0, pre, post, comm=MPI.COMM_WORLD)
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    
    if N == 1
        return (1/4)*F
    else
        U = copy(x0)
        for i = 1:pre
            U = Jacobi(A, F, U)
        end

        newN = floor(Int64, N/2)
        newA = Creer_A(newN)
        newF = zeros(newN*newN)
        newU = zeros(newN*newN)

        for i = 1:newN
            for j = 1:newN
                newF[j + newN*(i - 1)] = F[2j + N*(2i - 1)]
                newU[j + newN*(i - 1)] = U[2j + N*(2i - 1)]
            end
        end
       
        
        newU_global = zeros(Float64, newN, newN, size)
        newU_global = MPI.Allgather(newU, comm)
        local_newU = newU_global[:, :, rank]

        for i = 1:post
            local_newU = Jacobi(newA, newF, local_newU)
        end

        MPI.Gather(local_newU, newU_global, comm)

        if rank == 0
            newU = sum(newU_global, dims=3)[:, :, 1]
        end

        for i = 1:newN
            for j = 1:newN
                U[2j + N*(2i - 1)] = newU[j + newN*(i - 1)]
            end
        end

        for i = 1:pre
            U = Jacobi(A, F, U)
        end

        return U
    end
end


time_mpi = @elapsed begin
    N = 3
    A = Creer_A(N)
    F1 = Creer_F(p, q, N)
    x0 = rand(1:9, (N*N, 1))
    res_mpi = Multigrid_MPI(A, F1, N, x0, 1, 1)
    #println(res_mpi)
end

println("Elapsed time (MPI): ", time_mpi, " seconds")

# 终止 MPI
MPI.Finalize()
