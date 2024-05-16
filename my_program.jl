using MPI

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

println("Hello from rank $rank out of $size processes")

MPI.Finalize()