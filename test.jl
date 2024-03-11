using LinearAlgebra
using Random

x0 = rand(1:9, (N*N, 1))
w = 0.5
sol1 = A\F1
sol2 = A\F2

resJ1 = Jacobi2(A, F1, x0)
resJ2 = Jacobi2(A, F2, x0)
resGS1 = GaussSeidel2(A, F1, x0)
resGS2 = GaussSeidel2(A, F2, x0)
resS1 = SOR(A, F1, w, x0)
resS2 = SOR(A, F2, w, x0)

e1 = norm(sol1 - resJ1)
e2 = norm(sol2 - resJ2)
e1 = norm(sol1 - resGS1)
e2 = norm(sol2 - resGS2)
e1 = norm(sol1 - resS1)
e2 = norm(sol2 - resS2)
