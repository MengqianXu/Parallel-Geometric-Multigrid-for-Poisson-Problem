using LinearAlgebra
using SparseArrays
using Random
using MPI

include("problem.jl")

function parallel_jacobi_solver(F, U0, N, MaxIter, tol=1e-15)
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    
    a, b = decomposition(size)
    x = ceil(Int, N / a)
    y = ceil(Int, N / b)
    d = rank % a
    c = floor(Int, rank / a)

    U_local = zeros(y, x)
    F_local = zeros(y, x)

    if rank == 0
        temp_U0 = reshape(U0, (N, N))
        temp_F0 = reshape(F, (N, N))
        temp_U = NaN64*ones(a*x, b*y)
        temp_F = NaN64*ones(a*x, b*y)
        
        for m = 0:(a-1)
            for n = 0:(b-1)
                for i = 1:x
                    colonne = m*x + i
                    for j = 1:y
                        ligne = n*y + j
                        if (ligne <= N) && (colonne <= N)
                            temp_U[(m*x+i) + (n*y+j-1)*(a*x)] = temp_U0[ligne, colonne]
                            temp_F[(m*x+i) + (n*y+j-1)*(a*x)] = temp_F0[ligne, colonne]
                        end
                    end
                end
            end
        end

        temp_U = temp_U[:]
        temp_F = temp_F[:]
    else
        temp_U = zeros(a*x*b*y)
        temp_F = zeros(a*x*b*y)
    end
    
    MPI.Scatter!(temp_U, U_local, 0, comm)
    MPI.Scatter!(temp_F, F_local, 0, comm)
    
    U_local = reshape(U_local, (y, x))
    F_local = reshape(F_local, (y, x))
    U_new = zeros(y, x)
    iteration = 0

    requeteS = MPI.MultiRequest(4)
    requeteR = MPI.MultiRequest(4)

    tabgaucheR = zeros(y)
    tabdroiteR = zeros(y)
    tabhautR = zeros(x)
    tabbasR = zeros(x)

    while iteration < MaxIter
        iteration += 1
        nb = 1

        if c > 0
            tabgaucheS = U_local[:, 1]
            MPI.Isend(tabgaucheS, rank - a, 0, comm, requeteS[nb])
            MPI.Irecv!(tabgaucheR, rank - a, 0, comm, requeteR[nb])
            nb += 1
        end

        if c < (a - 1)
            tabdroiteS = U_local[:, x]
            MPI.Isend(tabdroiteS, rank + a, 0, comm, requeteS[nb])
            MPI.Irecv!(tabdroiteR, rank + a, 0, comm, requeteR[nb])
            nb += 1
        end

        if d > 0
            tabhautS = U_local[1, :]
            MPI.Isend(tabhautS, rank - 1, 0, comm, requeteS[nb])
            MPI.Irecv!(tabhautR, rank - 1, 0, comm, requeteR[nb])
            nb += 1
        end

        if d < (b - 1)
            tabbasS = U_local[y, :]
            MPI.Isend(tabbasS, rank + 1, 0, comm, requeteS[nb])
            MPI.Irecv!(tabbasR, rank + 1, 0, comm, requeteR[nb])
            nb += 1
        end

        MPI.Waitall(requeteS)
        MPI.Waitall(requeteR)
        
        norme_local = 0.0
        for i = 1:x
            colonne = c*x + i 
            for j = 1:y
                ligne = d*y + j
                if (ligne <= N) && (colonne <= N)
                    res = (1/(N + 1)^2) * F_local[j, i]
                    if i > 1
                        res += U_local[j, i - 1]
                    else
                        if c > 0
                            res += tabgaucheR[j]
                        end
                    end
    
                    if (colonne + 1) <= N
                        if i < x
                            res += U_local[j, i + 1]
                        else
                            if c < (a - 1)
                                res += tabdroiteR[j]
                            end
                        end
                    end
    
                    if j > 1
                        res += U_local[j - 1, i]
                    else
                        if d > 0
                            res += tabhautR[i]
                        end
                    end
    
                    if (ligne + 1) <= N
                        if j < y
                            res += U_local[j + 1, i]
                        else
                            if d < (b - 1)
                                res += tabbasR[i]
                            end
                        end
                    end
                    U_new[j, i] = res / 4
                    norme_local += (U_new[j, i] - U_local[j, i])^2
                end
            end
        end
        
        norme_global = MPI.Reduce(norme_local, +, 0, comm)
        norme = sqrt(norme_global)
        
        if norme < tol
            break
        end
        
        U_local .= U_new
    end
    
    U_local = U_local[:]
    temp_U = zeros(a*x*b*y)
    
    MPI.Gather!(U_local, temp_U, 0, comm)
    
    U = zeros(N, N)
    if rank == 0
        for m = 0:(a-1)
            for n = 0:(b-1)
                for i = 1:x
                    colonne = m*x + i
                    for j = 1:y
                        ligne = n*y + j
                        if (ligne <= N) && (colonne <= N)
                            U[ligne, colonne] = temp_U[(m*x+i) + (n*y+j-1)*(a*x)]
                        end
                    end
                end
            end
        end
    end
    
    MPI.Finalize()
    return U
end

function main(N)
    p = 1
    q = 1
    A = Creer_A(N)
    F = Creer_F(p, q, N)
    x0 = ones(N * N)
    MaxIter = 1000

    println("Solving with parallel Jacobi...")
    @time x_parallel_jacobi = parallel_jacobi_solver(F, x0, N, MaxIter)
    x_parallel = reshape(x_parallel_jacobi, (N, N))

    println("Solving with non-parallel Jacobi...")
    @time x_jacobi = Jacobi(A, F, x0, MaxIter)
    x_jacobi_reshaped = reshape(x_jacobi, (N, N))

    println("Difference between solutions: ", maximum(abs.(x_jacobi_reshaped - x_parallel)))
end

N = 7
main(N)
