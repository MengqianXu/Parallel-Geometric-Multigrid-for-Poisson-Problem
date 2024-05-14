using Distributed
using MPI
using LinearAlgebra
using BenchmarkTools

function test()
    println("hhhh")
end

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)
@benchmark begin
    test()
end

MPI.Finalize()