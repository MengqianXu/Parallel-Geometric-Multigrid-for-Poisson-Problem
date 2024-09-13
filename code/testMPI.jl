using MPI
using LinearAlgebra
using Random
using SparseArrays

include("problem.jl")
include("solvers.jl")
include("solversMPI.jl")


MPI.Init()
comm = MPI.COMM_WORLD
size = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)
@show N = 63
p = 2 
q = 2
@show p, q

A = Creer_A(N)
F = Creer_F(p, q, N)
U = Creer_U(p, q, N)
x0 = [j + N*(i - 1) for j = 1:N, i = 1:N][:]

resD = A\F
resJ = Jacobi(A, F, x0)
resGS = GaussSeidel(A, F, x0)
resJMPI = Jacobi_MPI(F, x0, N)

@show normeD = norm(U - resD)
@show normeJ = norm(U - resJ)
@show normeGS = norm(U - resGS)
@show normeJMPI = norm(U - resJMPI)
# MPI.Finalize()
