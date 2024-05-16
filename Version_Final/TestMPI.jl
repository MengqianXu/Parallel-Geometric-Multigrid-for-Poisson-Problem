using MPI
using LinearAlgebra
using Random
using SparseArrays

include("Problem.jl")
include("Solvers.jl")
include("Jacobi_MPI.jl")


MPI.Init()
comm = MPI.COMM_WORLD
size = MPI.Comm_size(comm)

N = 3
p = rand(-5:5) 
q = rand(-5:5)
@show p, q

A = Creer_A(N)
F = Creer_F(p, q, N)
U = Creer_U(p, q, N)
x0 = ones(N*N)

resD = A\F
resJ = Jacobi(A, F, x0)
resGS = GaussSeidel(A, F, x0)
resJMPI = Jacobi_MPI(F, x0, N)

normeD = norm(U - resD)
normeJ = norm(U - resJ)
normeGS = norm(U - resGS)
normeJMPI = norm(U - resJMPI)