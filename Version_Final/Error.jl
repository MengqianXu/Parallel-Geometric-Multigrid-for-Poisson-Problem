using LinearAlgebra
using Random
using Plots

include("Problem.jl")
include("Solvers.jl")

function Erreur(A, F, N, p, q, w = 0.5)
	x0 = ones(N*N)
	U = Creer_U(p, q, N)
	
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
	U = Creer_U(α, β, γ, δ, N)
	
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

# Convergence des méthodes itératives
taille = 7
α = 6
β = 1
γ = 5
δ = 10
p = 2
q = 2

tabN = zeros(taille)
tabN2 = zeros(taille)
tabC = zeros(taille)
tabH = zeros(taille)
tabH2 = zeros(taille)

tabEAD1 = zeros(taille)
tabEAJ1 = zeros(taille)
tabEAGS1 = zeros(taille)
tabEAS1 = zeros(taille)

tabERD1 = zeros(taille)
tabERJ1 = zeros(taille)
tabERGS1 = zeros(taille)
tabERS1 = zeros(taille)

tabEAD2 = zeros(taille)
tabEAJ2 = zeros(taille)
tabEAGS2 = zeros(taille)
tabEAS2 = zeros(taille)

tabERD2 = zeros(taille)
tabERJ2 = zeros(taille)
tabERGS2 = zeros(taille)
tabERS2 = zeros(taille)


for i in 1:taille
	@show i
	N = 2^i - 1
	H = 1/(N + 1)
	N2 = N^2
	H2 = H^2
	
	A = Creer_A(N)
	F1 = Creer_F(p, q, N)
	F2 = Creer_F(α, β, γ, δ, N)
	CondA = cond(Matrix(A))
	
	tabN[i] = N
	tabH[i] = H
	tabN2[i] = N2
	tabH2[i] = H2
	tabC[i] = CondA

	println("Pour F1")
	absolueD1, absolueJ1, absolueGS1, absolueS1, relativeD1, relativeJ1, relativeGS1, relativeS1 = Erreur(A, F1, N, p, q)
	tabEAD1[i] = absolueD1
	tabEAJ1[i] = absolueJ1
	tabEAGS1[i] = absolueGS1
	tabEAS1[i] = absolueS1
	tabERD1[i] = relativeD1
	tabERJ1[i] = relativeJ1
	tabERGS1[i] = relativeGS1
	tabERS1[i] = relativeS1

	println("Pour F2")
	absolueD2, absolueJ2, absolueGS2, absolueS2, relativeD2, relativeJ2, relativeGS2, relativeS2 = Erreur(A, F2, N, α, β, γ, δ)
	tabEAD2[i] = absolueD2
	tabEAJ2[i] = absolueJ2
	tabEAGS2[i] = absolueGS2
	tabEAS2[i] = absolueS2
	tabERD2[i] = relativeD2
	tabERJ2[i] = relativeJ2
	tabERGS2[i] = relativeGS2
	tabERS2[i] = relativeS2
end

nom_fichier = "data/Version3/Comparaisons/Convergence_des_méthodes_itératives.txt"
open(nom_fichier, "w") do file
	for i = 1:taille
		println(file, tabN[i], "\t", tabH[i], "\t", tabN2[i], "\t", tabH2[i], "\t", tabC[i], "\t", tabEAD1[i], "\t", tabEAJ1[i], "\t", tabEAGS1[i], "\t", tabEAS1[i], "\t", tabERD1[i], "\t", tabERJ1[i], "\t", tabERGS1[i], "\t", tabERS1[i], "\t", tabEAD2[i], "\t", tabEAJ2[i], "\t", tabEAGS2[i], "\t", tabEAS2[i], "\t", tabERD2[i], "\t", tabERJ2[i], "\t", tabERGS2[i], "\t", tabERS2[i])
	end
end
 
# p1 = plot(tabN, [tabN, tabN2, tabC], label = ["N N^2 Cond(A)"], title = "Conditionnement de la matrice A", xlabel = "N", ylabel = "Conditionnement", xscale=:log2, yscale=:log2)
# p2 = plot(tabN, [tabH, tabH2, tabEAD1, tabEAJ1, tabEAGS1, tabEAS1], label = ["H H^2 Direct Jacobi Gauss-Seidel SOR"], title = "Erreur absolue pour F1", xlabel = "N", ylabel = "Erreur", xscale=:log2, yscale=:log2)
# p3 = plot(tabN, [tabH, tabH2, tabERD1, tabERJ1, tabERGS1, tabERS1], label = ["H H^2 Direct Jacobi Gauss-Seidel SOR"], title = "Erreur relative pour F1", xlabel = "N", ylabel = "Erreur", xscale=:log2, yscale=:log2)
# p4 = plot(tabN, [tabH, tabH2, tabEAD2, tabEAJ2, tabEAGS2, tabEAS2], label = ["H H^2 Direct Jacobi Gauss-Seidel SOR"], title = "Erreur absolue pour F2", xlabel = "N", ylabel = "Erreur", xscale=:log2, yscale=:log2)
# p5 = plot(tabN, [tabH, tabH2, tabERD2, tabERJ2, tabERGS2, tabERS2], label = ["H H^2 Direct Jacobi Gauss-Seidel SOR"], title = "Erreur relative pour F2", xlabel = "N", ylabel = "Erreur", xscale=:log2, yscale=:log2)
# savefig(p1, "image/Julia/Conditionnement de A (1 à 7 - loglog).png")
# savefig(p2, "image/Julia/Erreur absolue pour F1 (1 à 7 - loglog).png")
# savefig(p3, "image/Julia/Erreur relative pour F1 (1 à 7 - loglog).png")
# savefig(p4, "image/Julia/Erreur absolue pour F2 (1 à 7 - loglog).png")
# savefig(p5, "image/Julia/Erreur relative pour F2 (1 à 7 - loglog).png")


# Convergence des méthodes multigrid
function ErreurM(A, F, N, p, q, pre, post)
	U = Creer_U(p, q, N)
	
	println("- Jacobi")
	resJ = MultigridJ(F, N, pre, post)
	println("- Gauss Seidel")
	resGS = MultigridGS(F, N, pre, post)
	
	absolueJ = norm(U - resJ)
	absolueGS = norm(U - resGS)
	
	relativeJ = norm(F - A*resJ)
	relativeGS = norm(F - A*resGS)
	
	return absolueJ, absolueGS, relativeJ, relativeGS
end

function ErreurM(A, F, N, α, β, γ, δ, pre, post)
	U = Creer_U(α, β, γ, δ, N)
	
	println("- Jacobi")
	resJ = MultigridJ(F, N, pre, post)
	println("- Gauss Seidel")
	resGS = MultigridGS(F, N, pre, post)
	
	absolueJ = norm(U - resJ)
	absolueGS = norm(U - resGS)
	
	relativeJ = norm(F - A*resJ)
	relativeGS = norm(F - A*resGS)
	
	return absolueJ, absolueGS, relativeJ, relativeGS
end

tabPre = [pre for post = 1:10, pre = 1:10][:]
tabPost = [post for post = 1:10, pre = 1:10][:]
taille = length(tabPre)
nom_fichier = "data/Version3/Conclusions erreur.txt"

for m = 1:7
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
	
	nom_fichier = string("data/Version3/Comparaisons/Convergence_des_méthodes_multigrid_avec_N_=_", string(N), ".txt")
	open(nom_fichier, "w") do file
		for i = 1:taille
			pre = tabPre[i]
			post = tabPost[i]
			println("i = ", i, ", pre = ", pre, ", post = ", post)
			
			println("Pour F1")
			absolueJ1, absolueGS1, relativeJ1, relativeGS1 = ErreurM(A, F1, N, p, q, pre, post)
			println("Pour F2")
			absolueJ2, absolueGS2, relativeJ2, relativeGS2 = ErreurM(A, F2, N, α, β, γ, δ,  pre, post)
			
			println(file, pre, "\t", post, "\t", absolueJ1, "\t", absolueGS1, "\t", relativeJ1, "\t", relativeGS1, "\t", absolueJ2, "\t", absolueGS2, "\t", relativeJ2, "\t", relativeGS2)

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
	end

	nom_fichier = "data/Version3/Conclusions erreur.txt"
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