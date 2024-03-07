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

function calculate_u(x, y, p, q)
    u = -sin(p*pi*x) * sin(q*pi*y)
    return u
end

function calculate_f(x, y, p, q)
    f = (p*pi + q*pi) * sin(p*pi*x) * sin(q*pi*y)
    return f
end

function create_f_matrix(x, y, p, q)
    n = length(x)
    m = length(y)
    f = [calculate_f(x[i], y[j], p, q) for j in 1:m, i in 1:n]
    return f
end

#AU=h^2f
function direct_solve(n,f_mat)
    h = 1/(n+1)
    A = create_A_matrix(n)
    L = cholesky(A)#A=L×L^T
    y = L \ (f_mat * h^2)
    U = L' \ y
    #LL^T⋅U=h^2⋅f, Ly=h^2⋅f,L^T⋅U=y
    return U
end

function iterative_solver_CG(n,f_mat)
    A = create_A_matrix(n)
    h = 1 / (n + 1)
    initial_guess = zeros(length(f_mat))
    U = cg(A, f_mat * h^2, x0=initial_guess)
    return U
end

function calculate_error(n,f_mat)
    cholesky = direct_solve(n,f)
    iterative = iterative_solver_CG(n,f_mat)
    error = norm(cholesky - iterative)
    return error
end



n = 2
m = 2
p = 3
q = 3

# Calculate u and f
x = range(0, stop=1, length=n)
y = range(0, stop=1, length=m)
u = calculate_u.(x', y, p, q)  
f = calculate_f.(x', y, p, q)
f_mat= create_f_matrix(x, y, p, q)
A = create_A_matrix(n)
println(A)
U_direct = direct_solve(n,f_mat)
println("Solution(directe) U pour n = $n:\n", U_direct)
println("Solution(iterative) U pour n = $n:\n", U_iterative)
U_iterative = iterative_solver_CG(n,ff_mat)
error = calculate_error(n,f_mat)
println("Erreur entre Cholesky et Iterative Solver: ", error)



