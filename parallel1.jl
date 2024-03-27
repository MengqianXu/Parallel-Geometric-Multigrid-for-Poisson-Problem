using Distributed
using MPI
using LinearAlgebra



comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

function Multigrid(A, F, N::Int64, x0, pre, post)
    if N == 1
        return (1/4)*F
    else
        U = copy(x0)

        @sync begin
            @distributed for i = 1:pre
                U = parallel_jacobi(A, F, U)
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

        newU = Multigrid(newA, newF, newN, newU, pre, post)

        for i = 1:newN
            for j = 1:newN
                U[2j + N*(2i - 1)] = newU[j + N*(i - 1)]
            end
        end

        @sync begin
            @distributed for i = 1:pre
                U = parallel_jacobi(A, F, U)
            end
        end

        return U
    end
end



function parallel_jacobi(A, b, x0, tol = 1e-15, max_iters = 1000)
    ancien = copy(x0)
    nouveau = copy(b)
    n = length(b)
    iter = 0
    while iter < max_iters
        iter += 1
        for i = 1:n
            sum = 0.0
            for j = 1:n
                if j != i
                    sum += A[i, j] * ancien[j]
                end
            end
            nouveau[i] = (b[i] - sum) / A[i, i]
        end
        
        MPI.Allreduce(nouveau, ancien, MPI.SUM, comm)
        
        if norm(nouveau - ancien) < tol
            break
        end
    end
    
    return nouveau
end
