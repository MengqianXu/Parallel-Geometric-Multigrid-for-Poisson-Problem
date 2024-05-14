using LinearAlgebra
using Random
using SparseArrays

include("problem.jl")
include("solvers.jl")

N = 3
# p = rand(-5:5) 
# q = rand(-5:5)
# α = rand(1:10) 
# β = rand(1:10) 
# γ = rand(1:10) 
# δ = rand(1:10)

# A = Creer_A(N)
# F1 = Creer_F(p, q, N)
# F2 = Creer_F(α, β, γ, δ, N)
# U1 = Creer_U(p, q, N)
# U2 = Creer_U(α, β, γ, δ, N)
# x0 = ones(N*N)
# test = [1 for j = 1:N, i = 1:N]
# println(test)
# prolong, N2 = Prolongation(test, N)
# println(reshape(prolong, (N2, N2)))
# restrict, N3 = Restriction(prolong, N2)
# println(reshape(restrict, (N3, N3)))
# resD = A\F1
# resJ = Jacobi(A, F1, x0)
# resGS = GaussSeidel(A, F1, x0)
# resMJ = MultigridJ(F1, N, 1, 1)
# resMGS = MultigridGS(F1, N, 1, 1)

# normeD = norm(U1 - resD)
# normeJ = norm(U1 - resJ)
# normeGS = norm(U1 - resGS)
# normeMJ = norm(U1 - resMJ)
# normeMGS = norm(U1 - resMGS)

# test1 = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
# test2 = Creer_U(p, q, N)

# println(test1)
# println(test2)
test = zeros((N, N))
for i = 1:(N*N)
    test[i] = i
end
println(test)