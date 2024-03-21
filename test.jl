using LinearAlgebra
using Random
using Plots


# [(i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]

N = 3
A = Creer_A(N)
F1 = Creer_F(p, q, N)
F2 = Creer_F(α, β, γ, δ, N)
U1 = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
U2 = [u(α, β, γ, δ, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
x0 = rand(1:9, (N*N, 1))
w = 0.5

tabY = [i/(N + 1) for i = 1:N, j = 1:N][:]
tabX = [j/(N + 1) for i = 1:N, j = 1:N][:]

sol1 = A\F1
resJ1 = Jacobi(A, F1, x0)
resGS1 = GaussSeidel(A, F1, x0)
resS1 = SOR(A, F1, w, x0)
resM1 = Multigrid(A, F1, N, x0, 1, 1)

p1 = plot(tabX, tabY, [U1, sol1, resJ1, resGS1, resS1, resM1], label = ["Solution exacte" "Directe" "Jacobi" "Gauss-Seidel" "SOR" "Multigrid"])

sol2 = A\F2
resJ2 = Jacobi(A, F2, x0)
resGS2 = GaussSeidel(A, F2, x0)
resS2 = SOR(A, F2, w, x0)
resM2 = Multigrid(A, F2, N, x0, 1, 1)

p2 = plot(tabX, tabY, [U2, sol2, resJ2, resGS2, resS2, resM2], label = ["Solution exacte" "Directe" "Jacobi" "Gauss-Seidel" "SOR" "Multigrid"])

relativeJ1 = norm(sol1 - resJ1)
relativeGS1 = norm(sol1 - resGS1)
relativeS1 = norm(sol1 - resS1)
relativeM1 = norm(sol1 - resM1)

relativeJ2 = norm(sol2 - resJ2)
relativeGS2 = norm(sol2 - resGS2)
relativeS2 = norm(sol2 - resS2)
relativeM2 = norm(sol2 - resM2)

absolueJ1 = norm(U1 - resJ1)
absolueGS1 = norm(U1 - resGS1)
absolueS1 = norm(U1 - resS1)
absolueM1 = norm(U1 - resM1)

absolueJ2 = norm(U2 - resJ2)
absolueGS2 = norm(U2 - resGS2)
absolueS2 = norm(U2 - resS2)
absolueM2 = norm(U2 - resM2)

