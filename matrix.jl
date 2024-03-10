#Create two matrices (A&B) and solve linear problems
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using Plots

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


function calculate_u(x, y, p, q)
    u = -sin(p*pi*x) * sin(q*pi*y)
    return u
end

function calculate_f(x, y, u)
    n = length(x)
    h = x[2] - x[1]
    f = zeros(n^2)
    idx = 1
    for j in 1:n
        for i in 1:n
            if i == 1 || i == n || j == 1 || j == n
                f[idx] = 0  # Gérer les conditions aux bords u=g où g=0
            else
                f[idx] = -(u[i+1, j] - 2u[i, j] + u[i-1, j]) / h^2 - (u[i, j+1] - 2u[i, j] + u[i, j-1]) / h^2
            end
            idx += 1
        end
    end
    return f
end


#AU=h^2f
function direct_solve(n,f_mat)
    h = 1/(n+1)
    A = create_A_matrix(n)
    L = cholesky(A)#A=L×L^T
    #LL^T⋅U=h^2⋅f, Ly=h^2⋅f,L^T⋅U=y
    #println("the sizeof f_mat",size(f_mat))
    U = L\(f_mat*h^2)
    return U
end

function iterative_solver_CG1(n,f_mat)
    A = create_A_matrix(n)
    h = 1 / (n + 1)
    initial_guess = zeros(length(f_mat))
    U = cg(A, f_mat * h^2, x0=initial_guess)
    return U
end

function iterative_solver_CG(n, f_mat)
    A = create_A_matrix(n)
    h = 1 / (n + 1)
    U = cg(A, f_mat * h^2)
    return U
end

function calculate_error(U_direct,U_iterative)
    error = error = norm(U_direct - U_iterative)
    return error
end


n = 3
m = 3
p = 3
q = 3

# simple test
x = range(0, stop=1, length=n)
y = range(0, stop=1, length=n)
u = calculate_u.(x', y, p, q)
f_mat = calculate_f(x, y, u)
#println(f_mat)
A = create_A_matrix(n)
#println(A)
U_direct = direct_solve(n,f_mat)
#println("Solution(directe) U pour n = $n:\n", U_direct)
U_iterative = iterative_solver_CG(n,f_mat)
#println("Solution(iterative) U pour n = $n:\n", U_iterative)
error=calculate_error(U_direct,U_iterative)
#println("Erreur entre Cholesky et Iterative Solver: ", error)

#plot
values = 2:10
function get_error(n, p, q)
    x = range(0, stop=1, length=n)
    y = range(0, stop=1, length=n)
    u = calculate_u.(x', y, p, q)
    f_mat = calculate_f(x, y, u)
    U_direct = direct_solve(n, f_mat)
    U_iterative = iterative_solver_CG(n, f_mat)
    error = norm(U_direct - U_iterative)
    println("Erreur entre Cholesky et Iterative Solver: ", error)
    return error
end

# Calculate errors for each value of n
errors = [get_error(val, p, q) for val in values]
plot(n_values, errors, xlabel="n", ylabel="Error", title="Erreur entre Cholesky et Iterative Solver", legend=false, marker=:circle)


