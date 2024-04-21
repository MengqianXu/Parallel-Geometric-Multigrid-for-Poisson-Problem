using LinearAlgebra
using Random
using Plots

function Erreur(A, F, N, p, q, w = 0.5)
	x0 = ones(N*N)
	U = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
	
	println("- Direct")
	resD = A\F
	println("- Jacobi")
	resJ = Jacobi(A, F, x0)
	println("- Gauss Seidel")
	resGS = GaussSeidel(A, F, x0)
	println("- SOR")
	resS = SOR(A, F, w, x0)
	
	absolueD = norm(U - resD)
	absolueJ = norm(U - resJ)
	absolueGS = norm(U - resGS)
	absolueS = norm(U - resS)

	relativeD = norm(F - A*resJ)
	relativeJ = norm(F - A*resJ)
	relativeGS = norm(F - A*resGS)
	relativeS = norm(F - A*resS)

	return absolueD, absolueJ, absolueGS, absolueS, relativeD, relativeJ, relativeGS, relativeS
end

function Erreur(A, F, N, α, β, γ, δ, w = 0.5)
	x0 = ones(N*N)
	U = [u(α, β, γ, δ, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
	
	println("- Direct")
	resD = A\F
	println("- Jacobi")
	resJ = Jacobi(A, F, x0)
	println("- Gauss Seidel")
	resGS = GaussSeidel(A, F, x0)
	println("- SOR")
	resS = SOR(A, F, w, x0)
	
	absolueD = norm(U - resD)
	absolueJ = norm(U - resJ)
	absolueGS = norm(U - resGS)
	absolueS = norm(U - resS)

	relativeD = norm(F - A*resJ)
	relativeJ = norm(F - A*resJ)
	relativeGS = norm(F - A*resGS)
	relativeS = norm(F - A*resS)

	return absolueD, absolueJ, absolueGS, absolueS, relativeD, relativeJ, relativeGS, relativeS
end

# # Convergence des méthodes itératives
# taille = 6
α = 6
β = 1
γ = 5
δ = 10
p = 2
q = 2

# open("data/Version3/Convergence_des_methodes_iteratives.txt", "w") do file
# 	for i in 2:taille
# 		@show i
# 		N = 2^i - 1
# 		H = 1/(N + 1)
# 		H2 = (1/(N + 1))^2
		
# 		A = Creer_A(N)
# 		F1 = Creer_F(p, q, N)
# 		F2 = Creer_F(α, β, γ, δ, N)
# 		CondA = cond(Matrix(A))
		
# 		println("Pour F1")
# 		absolueD1, absolueJ1, absolueGS1, absolueS1, relativeD1, relativeJ1, relativeGS1, relativeS1 = Erreur(A, F1, N, p, q)
# 		println("Pour F2")
# 		absolueD2, absolueJ2, absolueGS2, absolueS2, relativeD2, relativeJ2, relativeGS2, relativeS2 = Erreur(A, F2, N, α, β, γ, δ)
		
# 		println(file, N, "\t", H, "\t", H2, "\t", CondA, "\t", absolueD1, "\t", absolueJ1, "\t", absolueGS1, "\t", absolueS1, "\t", relativeD1, "\t", relativeJ1, "\t", relativeGS1, "\t", relativeS1, "\t", absolueD2, "\t", absolueJ2, "\t", absolueGS2, "\t", absolueS2, "\t", relativeD2, "\t", relativeJ2, "\t", relativeGS2, "\t", relativeS2)
# 	end
# end

# Convergence des méthodes multigrid
function ErreurM(A, F, N, p, q, pre, post)
	x0 = ones(N*N)
	U = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
	
	println("- Jacobi")
	resJ = MultigridJ(A, F, N, x0, pre, post)
	println("- Gauss Seidel")
	resGS = MultigridGS(A, F, N, x0, pre, post)
	
	absolueJ = norm(U - resJ)
	absolueGS = norm(U - resGS)
	
	relativeJ = norm(F - A*resJ)
	relativeGS = norm(F - A*resGS)
	
	return absolueJ, absolueGS, relativeJ, relativeGS
end

function ErreurM(A, F, N, α, β, γ, δ, pre, post)
	x0 = ones(N*N)
	U = [u(α, β, γ, δ, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
	
	println("- Jacobi")
	resJ = MultigridJ(A, F, N, x0, pre, post)
	println("- Gauss Seidel")
	resGS = MultigridGS(A, F, N, x0, pre, post)
	
	absolueJ = norm(U - resJ)
	absolueGS = norm(U - resGS)
	
	relativeJ = norm(F - A*resJ)
	relativeGS = norm(F - A*resGS)
	
	return absolueJ, absolueGS, relativeJ, relativeGS
end

tabPre = [pre for post = 1:10, pre = 1:10][:]
tabPost = [post for post = 1:10, pre = 1:10][:]
taille = length(tabPre)
nom_fichier = "data/Version3/Conclusions.txt"

for m = 2:6
	@show m
	N = 2^m - 1
	A = Creer_A(N)
	F1 = Creer_F(p, q, N)
	F2 = Creer_F(α, β, γ, δ, N)
	x0 = ones(N*N)

	open(nom_fichier, "a") do file
		println(file, "Pour N = ", N, " : ")
	end

	minPreAJ1 = Inf
	minPostAJ1 = Inf
	minEAJ1 = Inf
	
	minPreAGS1 = Inf
	minPostAGS1 = Inf
	minEAGS1 = Inf
	
	minPreRJ1 = Inf
	minPostRJ1 = Inf
	minERJ1 = Inf
	
	minPreRGS1 = Inf
	minPostRGS1 = Inf
	minERGS1 = Inf

	minPreAJ2 = Inf
	minPostAJ2 = Inf
	minEAJ2 = Inf
	
	minPreAGS2 = Inf
	minPostAGS2 = Inf
	minEAGS2 = Inf
	
	minPreRJ2 = Inf
	minPostRJ2 = Inf
	minERJ2 = Inf
	
	minPreRGS2 = Inf
	minPostRGS2 = Inf
	minERGS2 = Inf
	
	# nom_fichier = string("data/Version3/Convergence_des_methodes_multigrid_avec_N_=_", string(N), ".txt")
	# open(nom_fichier, "w") do file
		for i = 1:taille
			pre = tabPre[i]
			post = tabPost[i]
			println("i = ", i, ", pre = ", pre, ", post = ", post)
			
			println("Pour F1")
			absolueJ1, absolueGS1, relativeJ1, relativeGS1 = ErreurM(A, F1, N, p, q, pre, post)
			println("Pour F2")
			absolueJ2, absolueGS2, relativeJ2, relativeGS2 = ErreurM(A, F2, N, α, β, γ, δ,  pre, post)
			
			# println(file, pre, "\t", post, "\t", absolueJ1, "\t", absolueGS1, "\t", relativeJ1, "\t", relativeGS1, "\t", absolueJ2, "\t", absolueGS2, "\t", relativeJ2, "\t", relativeGS2)

			if absolueJ1 < minEAJ1
				minPreAJ1 = pre
				minPostAJ1 = post
				minEAJ1 = absolueJ1
			end

			if absolueJ2 < minEAJ2
				minPreAJ2 = pre
				minPostAJ2 = post
				minEAJ2 = absolueJ2
			end

			if absolueGS1 < minEAGS1
				minPreAGS1 = pre
				minPostAGS1 = post
				minEAGS1 = absolueGS1
			end

			if absolueGS2 < minEAGS2
				minPreAGS2 = pre
				minPostAGS2 = post
				minEAGS2 = absolueGS2
			end

			if relativeJ1 < minERJ1
				minPreRJ1 = pre
				minPostRJ1 = post
				minERJ1 = relativeJ1
			end

			if relativeGS1 < minERGS1
				minPreRGS1 = pre
				minPostRGS1 = post
				minERGS1 = relativeGS1
			end

			if relativeJ2 < minERJ2
				minPreRJ2 = pre
				minPostRJ2 = post
				minERJ2 = relativeJ2
			end

			if relativeGS2 < minERGS2
				minPreRGS2 = pre
				minPostRGS2 = post
				minERGS2 = relativeGS2
			end
		end
	# end

	open(nom_fichier, "a") do file
		println(file, "\tPour F1 : ")
		println(file, "\t\tAvec Jacobi : ")
		println(file, "\t\tErreur absolue : ")
		println(file, "\t\t\t- pre = ", minPreAJ1)
		println(file, "\t\t\t- post = ", minPostAJ1)
		println(file, "\t\t\t- erreur = ", minEAJ1, "\n")
		println(file, "\t\tErreur relative : ")
		println(file, "\t\t\t- pre = ", minPreRJ1)
		println(file, "\t\t\t- post = ", minPostRJ1)
		println(file, "\t\t\t- erreur = ", minERJ1, "\n")
		
		println(file, "\t\tAvec Gauss-Seidel : ")
		println(file, "\t\tErreur absolue : ")
		println(file, "\t\t\t- pre = ", minPreAGS1)
		println(file, "\t\t\t- post = ", minPostAGS1)
		println(file, "\t\t\t- erreur = ", minEAGS1, "\n")
		println(file, "\t\tErreur relative : ")
		println(file, "\t\t\t- pre = ", minPreRGS1)
		println(file, "\t\t\t- post = ", minPostRGS1)
		println(file, "\t\t\t- erreur = ", minERGS1, "\n")
		
		println(file, "\tPour F2 : ")
		println(file, "\t\tAvec Jacobi : ")
		println(file, "\t\tErreur absolue : ")
		println(file, "\t\t\t- pre = ", minPreAJ2)
		println(file, "\t\t\t- post = ", minPostAJ2)
		println(file, "\t\t\t- erreur = ", minEAJ2, "\n")
		println(file, "\t\tErreur relative : ")
		println(file, "\t\t\t- pre = ", minPreRJ2)
		println(file, "\t\t\t- post = ", minPostRJ2)
		println(file, "\t\t\t- erreur = ", minERJ2, "\n")
		
		println(file, "\t\tAvec Gauss-Seidel : ")
		println(file, "\t\tErreur absolue : ")
		println(file, "\t\t\t- pre = ", minPreAGS2)
		println(file, "\t\t\t- post = ", minPostAGS2)
		println(file, "\t\t\t- erreur = ", minEAGS2, "\n")
		println(file, "\t\tErreur relative : ")
		println(file, "\t\t\t- pre = ", minPreRGS2)
		println(file, "\t\t\t- post = ", minPostRGS2)
		println(file, "\t\t\t- erreur = ", minERGS2, "\n")
	end
end

println("\n")
println("α = ", α, "β = ", β, "γ = ", γ, "δ = ", δ)
println("p = ", p, "q = ", q)