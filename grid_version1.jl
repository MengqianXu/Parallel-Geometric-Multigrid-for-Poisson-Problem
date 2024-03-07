using LinearAlgebra
using SparseArrays
using IterativeSolvers

#The main idea
#Preprocessing: Perform a certain number of iterations (Jacobi, Gauss-Seidel) on the initial mesh to improve the quality of the initial guess.
#Smoothing: Perform a certain number of iterations at each grid level to reduce the residual.
#Restriction: Restrict the solution to a coarse grid to reduce the size of the problem.
#Solve: Solve smaller problems on a coarse grid.
#Interpolation: Interpolate the solution back onto a fine grid to update a more accurate solution.
#Iteration: Repeat the above steps until convergence.


function create_B_matrix(n::Int)
    B = spdiagm(0 => fill(4.0, n),
                -1 => fill(-1.0, n-1),
                1 => fill(-1.0, n-1))
    return B
end

function create_A_matrix(n::Int)
    B = create_B_matrix(n)
    A = kron(B, I(n)) + kron(I(n), B) - 4*I(n^2)
    return A
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

#step 1
function jacobi(A, b, x0)
    D = Diagonal(A)
    E = -(LowerTriangular(A) - D)
    F = -(UpperTriangular(A) - D)
    M = D
    N = E + F
    Minv = inv(M)
    ancien = x0
    res = Minv * (N*ancien + b)
    norme = norm(res - ancien)
    while(norme >= 1e-15)
        ancien = res
        res = Minv * (N*ancien + b)
        norme = norm(res - ancien)
    end

    return res
end

#step2
function smooth(u, f)
    n = size(u, 1)
    u_new = similar(u)
    for i in 2:n-1
        for j in 2:n-1
            u_new[i, j] = (f[i, j] + u[i-1, j] + u[i+1, j] + u[i, j-1] + u[i, j+1]) / 4
        end
    end
    return u_new
end

#step3
function restrict(u)
    n = Int(sqrt(length(u)))
    u_coarse = similar(zeros(n÷2, n÷2))
    for i in 1:n÷2
        for j in 1:n÷2
            u_coarse[i, j] = u[2i-1, 2j-1] + u[2i, 2j-1] + u[2i-1, 2j] + u[2i, 2j]
        end
    end
    return u_coarse
end

# step 4
function interpolate()
    
end

function v_cycle(n::Int, f::Vector; nu=2)
    if n == 2
        return direct_solve(n, f)
    else
        u = zeros(n, n)
        for i in 1:nu
            u = smooth(u, f)
        end
        residual = f - create_A_matrix(n) * vec(u)
        f_coarse = restrict(residual)
        e_coarse = v_cycle(n÷2, f_coarse; nu=nu)
        e_fine = interpolate(e_coarse)
        u += e_fine
        for i in 1:nu
            u = smooth(u, f)
        end
        return u
    end
end

n = 2
m = 2
p = 2
q = 3


x = range(0, stop=1, length=n)
y = range(0, stop=1, length=m)
f_matrix = create_f_matrix(x', y, p, q)


A = create_A_matrix(n)
println("Matrix A:")
println(A)
