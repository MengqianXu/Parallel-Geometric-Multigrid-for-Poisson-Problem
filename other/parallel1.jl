using Distributed
using MPI
using LinearAlgebra
using SparseArrays
using Random
using Plots


function Multigrid(A, F, N::Int64, x0, pre, post)
    if N == 1
        return (1/4)*F
    else
        U = copy(x0)
        #Distributed memory parallel computing
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

#At the end of each iteration, use the MPI.Allreduce function to summarize the nouveau vectors of all processes, 
#add their results, and store the result back to nouveau

function parallel_jacobi(A, b, x0, tol = 1e-15, max_iters = 1000)
    ancien = copy(x0)
    nouveau = copy(b)
    n = length(b)
    iter = 0
    temp = zeros(length(ancien))
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
        temp = MPI.Allreduce(nouveau, MPI.SUM, comm)
        nouveau = copy(temp) 
        if norm(nouveau - ancien) < tol
            break
        end
    end
    
    return nouveau

end


    

MPI.Initialized()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)
N = 3
A = Creer_A(N)
F1 = Creer_F(p, q, N)
U1 = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
x0 = rand(1:9, (N*N, 1))
tabY = [i/(N + 1) for i = 1:N, j = 1:N][:]
tabX = [j/(N + 1) for i = 1:N, j = 1:N][:]
sol1 = A\F1
@benchmark begin
    resM1 = Multigrid(A, F1, N, x0, 1, 1)
end
println(resM1)
p1 = plot(tabX, tabY, [U1, sol1,  resM1], label = ["Solution exacte" "Directe" "Multigrid"])

MPI.Finalize()