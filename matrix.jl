#Create two matrices (A&B) and solve linear problems
using LinearAlgebra
using SparseArrays
using IterativeSolvers

function create_B_matrix(n)
    # La fonction spdiagm est utilisée pour créer une matrice diagonale clairsemée
    B = spdiagm(0 => fill(4.0, n),
                -1 => fill(-1.0, n-1),
                1 => fill(-1.0, n-1))
    return B
end

function create_A_matrix(n)
    B = create_B_matrix(n)
    A = kron(B, I(n)) + kron(I(n), B) - 4*I(n^2)
    return A
end

function create_F(n)
    f = zeros(n^2)
    return f
end

#AU=h^2f
function direct_solve(n)
    h = 1/(n+1)
    A = create_A_matrix(n)
    f = create_F(n)
    L = cholesky(A)#A=L×L^T
    y = L \ (f * h^2)
    U = L' \ y
    #LL^T⋅U=h^2⋅f, Ly=h^2⋅f,L^T⋅U=y
    return U
end

function iterative_solver_CG(n)
    A = create_A_matrix(n)
    f = create_F(n)
    h = 1 / (n + 1)
    initial_guess = zeros(length(f))
    U = cg(A, f * h^2, x0=initial_guess)
    return U
end

function calculate_error(n)
    cholesky = direct_solve(n)
    iterative = iterative_solver_CG(n)
    error = norm(cholesky - iterative)
    return error
end



n = 2
A = create_A_matrix(n)
println(A)
U_direct = direct_solve(n)
println("Solution(directe) U pour n = $n:\n", U_direct)
println("Solution(iterative) U pour n = $n:\n", U_iterative)
U_iterative = iterative_solver_CG(n)
error = calculate_error(n)
println("Erreur entre Cholesky et Iterative Solver: ", error)
