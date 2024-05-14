using LinearAlgebra
using SparseArrays
using Random
using Plots


include("problem.jl")
include("solvers.jl")

α = 6
β = 1
γ = 5
δ = 10
p = 2
q = 2

for n = 1:7
	@show n
	N = 2^n - 1
	A = Creer_A(N)
	F1 = Creer_F(p, q, N)
	F2 = Creer_F(α, β, γ, δ, N)
	U1 = Creer_U(p, q, N)
	U2 = Creer_U(α, β, γ, δ, N)
	x0 = ones(N*N, 1)
	w = 0.5

	tabY = [i/(N + 1) for i = 1:N, j = 1:N][:]
	tabX = [j/(N + 1) for i = 1:N, j = 1:N][:]

	sol1 = A\F1
	resJ1 = Jacobi(A, F1, x0)
	resGS1 = GaussSeidel(A, F1, x0)
	resS1 = SOR(A, F1, w, x0)
	resMJ1 = MultigridJ(F1, N, 1, 1)
	resMGS1 = MultigridGS(F1, N, 1, 1)

	sol2 = A\F2
	resJ2 = Jacobi(A, F2, x0)
	resGS2 = GaussSeidel(A, F2, x0)
	resS2 = SOR(A, F2, w, x0)
	resMJ2 = MultigridJ(F2, N, 1, 1)
	resMGS2 = MultigridGS(F2, N, 1, 1)
	
	nom_fichier = string("data/Version3/Resultats/Resultats_avec_N_=_", string(N), ".txt")
	open(nom_fichier, "w") do file
		for i = 1:(N*N)
			println(file, tabX[i], "\t", tabY[i], "\t", U1[i], "\t", sol1[i], "\t", resJ1[i], "\t", resGS1[i], "\t", resS1[i], "\t", resMJ1[i], "\t", resMGS1[i], "\t", U2[i], "\t", sol2[i], "\t", resJ2[i], "\t", resGS2[i], "\t", resS2[i], "\t", resMJ2[i], "\t", resMGS2[i])
		end
	end
end



# p11 = surface(tabX, tabY, U1, title = "Solution exacte")
# p12 = surface(tabX, tabY, sol1, title = "Directe")
# p13 = surface(tabX, tabY, resJ1,  title = "Jacobi")
# p14 = surface(tabX, tabY, resGS1, title = "Gauss-Seidel")
# p15 = surface(tabX, tabY, resS1, title = "SOR")
# p16 = surface(tabX, tabY, resMJ1, title = "Multigrid avec Jacobi")
# p17 = surface(tabX, tabY, resMGS1, title = "Multigrid avec Gauss-Seidel")

# p2 = surface(tabX, tabY, U2, title = "Solution exacte")
# p2 = surface(tabX, tabY, sol2, title = "Directe")
# p2 = surface(tabX, tabY, resJ2, title = "Jacobi")
# p2 = surface(tabX, tabY, resGS2, title = "Gauss-Seidel")
# p2 = surface(tabX, tabY, resS2, title = "SOR")
# p2 = surface(tabX, tabY, resMJ2, title = "Multigrid_Jacobi")
# p2 = surface(tabX, tabY, resMGS2, title = "Multigrid_GaussSeidel")

# relativeJ1 = norm(sol1 - resJ1)
# relativeGS1 = norm(sol1 - resGS1)
# relativeS1 = norm(sol1 - resS1)
# relativeMJ1 = norm(sol1 - resMJ1)
# relativeMGS1 = norm(sol1 - resMGS1)

# relativeJ2 = norm(sol2 - resJ2)
# relativeGS2 = norm(sol2 - resGS2)
# relativeS2 = norm(sol2 - resS2)
# relativeMJ2 = norm(sol2 - resMJ2)
# relativeMGS2 = norm(sol2 - resMGS2)

# absolueJ1 = norm(U1 - resJ1)
# absolueGS1 = norm(U1 - resGS1)
# absolueS1 = norm(U1 - resS1)
# absolueMJ1 = norm(U1 - resMJ1)
# absolueMGS1 = norm(U1 - resMGS1)

# absolueJ2 = norm(U2 - resJ2)
# absolueGS2 = norm(U2 - resGS2)
# absolueS2 = norm(U2 - resS2)
# absolueMJ2 = norm(U2 - resMJ2)
# absolueMGS2 = norm(U2 - resMGS2)

# L, C, V = findnz(A)
# println(L)
# println(C)
# println(V)