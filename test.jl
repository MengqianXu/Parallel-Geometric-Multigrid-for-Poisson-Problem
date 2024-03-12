using LinearAlgebra
using Random

N = 3
A = Creer_A(N)
F1 = Creer_F(p, q, N)
F2 = Creer_F(α, β, γ, δ, N)
x0 = rand(1:9, (N*N, 1))
w = 0.5
sol1 = A\F1
sol2 = A\F2

resJ1 = Jacobi(A, F1, x0)
resJ2 = Jacobi(A, F2, x0)
resGS1 = GaussSeidel(A, F1, x0)
resGS2 = GaussSeidel(A, F2, x0)
resS1 = SOR(A, F1, w, x0)
resS2 = SOR(A, F2, w, x0)
resM1 = Multigrid(A, F1, N, x0, 1, 1)
resM2 = Multigrid(A, F2, N, x0, 1, 1)

e1 = norm(sol1 - resJ1)
e2 = norm(sol2 - resJ2)
e1 = norm(sol1 - resGS1)
e2 = norm(sol2 - resGS2)
e1 = norm(sol1 - resS1)
e2 = norm(sol2 - resS2)
e1 = norm(sol1 - resM1)
e2 = norm(sol2 - resM2)