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
rank = MPI.Comm_rank(comm)

N = 3
p = 2 
q = 2
#@show p, q

A = Creer_A(N)
F = Creer_F(p, q, N)
U = Creer_U(p, q, N)
x0 = ones(N*N)

resD = A\F
resJ = Jacobi(A, F, x0)
resGS = GaussSeidel(A, F, x0)
resJMPI = Jacobi_MPI(F, x0, N)
#print("u_mpi:",resJMPI)
#print("u_mpi:",resJ)
if rank == 0
    normeD = norm(U - resD)
    normeJ = norm(U - resJ)
normeGS = norm(U - resGS)
normeJMPI = norm(U - resJMPI)
normeJ_M = norm(resJ - resJMPI)
    println("norm jacobi $normeJ \nnorm Mpi $normeJMPI \nnorm Mpi_jacobi $normeJ_M")
end