using Distributed
using MPI
using LinearAlgebra
using SparseArrays
using Random

MPI.Init()
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
        
        temp = zeros(length(ancien))
        MPI.Allreduce(nouveau, temp, MPI.SUM, comm)
        copy!(ancien, temp)
        
        if norm(nouveau - ancien) < tol
            break
        end
    end
    
    return nouveau
end

N = 3
A = Creer_A(N)
F1 = Creer_F(p, q, N)
F2 = Creer_F(α, β, γ, δ, N)
U1 = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
U2 = [u(α, β, γ, δ, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
x0 = rand(1:9, (N*N, 1))

tabY = [i/(N + 1) for i = 1:N, j = 1:N][:]
tabX = [j/(N + 1) for i = 1:N, j = 1:N][:]

sol1 = A\F1

resM1 = Multigrid(A, F1, N, x0, 1, 1)
p1 = plot(tabX, tabY, [U1, sol1,  resM1], label = ["Solution exacte" "Directe" "Multigrid"])

sol2 = A\F2
resM2 = Multigrid(A, F2, N, x0, 1, 1)
p2 = plot(tabX, tabY, [U2, sol2,resM2], label = ["Solution exacte" "Directe" "Multigrid"])

MPI.Finalize()