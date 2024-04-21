using LinearAlgebra
using Random
using Plots

α = 6
β = 1
γ = 5
δ = 10
p = 2
q = 2

# # Temps de résolution des méthodes itératives
# taille = 6
# N = 3
# A = Creer_A(N)
# F1 = Creer_F(p, q, N)
# F2 = Creer_F(α, β, γ, δ, N)
# x0 = ones(N*N)
# w = 0.5
# @elapsed A\F1
# @elapsed Jacobi(A, F1, x0)
# @elapsed GaussSeidel(A, F1, x0)
# @elapsed SOR(A, F1, w, x0)

# open("data/Version3/Temps_des_methodes_iteratives.txt", "w") do file
# 	for i in 2:taille
# 		@show i
# 		N = 2^i - 1
# 		N2 = N^2
# 		N3 = N^3

# 		A = Creer_A(N)
# 		F1 = Creer_F(p, q, N)
# 		F2 = Creer_F(α, β, γ, δ, N)
# 		x0 = ones(N*N)
# 		w = 0.5

# 		println("Pour F1 : ")
# 		TD1 = @elapsed A\F1
# 		print("- direct = ")
# 		println(TD1)
		
# 		TJ1 = @elapsed Jacobi(A, F1, x0)
# 		print("- Jacobi = ")
# 		println(TJ1)
		
# 		TGS1 = @elapsed GaussSeidel(A, F1, x0)
# 		print("- Gauss-Seidel = ")
# 		println(TGS1)
		
# 		TS1 = @elapsed SOR(A, F1, w, x0)
# 		print("- SOR = ")
# 		println(TS1, "\n")
		
# 		println("Pour F2 : ")
# 		TD2 = @elapsed A\F2
# 		print("- direct = ")
# 		println(TD2)
		
# 		TJ2 = @elapsed Jacobi(A, F2, x0)
# 		print("- Jacobi = ")
# 		println(TJ2)
		
# 		TGS2 = @elapsed GaussSeidel(A, F2, x0)
# 		print("- Gauss-Seidel = ")
# 		println(TGS2)
		
# 		TS2 = @elapsed SOR(A, F2, w, x0)
# 		print("- SOR = ")
# 		println(TS2, "\n", "\n")

# 		println(file, N, "\t", N2, "\t", N3, "\t", TD1, "\t", TJ1, "\t", TGS1, "\t", TS1, "\t", TD2, "\t", TJ2, "\t", TGS2, "\t", TS2)
# 	end
# end

# Temps de résolution des méthodes multigrid
tabPre = [pre for post = 1:10, pre = 1:10][:]
tabPost = [post for post = 1:10, pre = 1:10][:]
taille = length(tabPre)
nom_fichier = "data/Version3/Conclusions.txt"
open(nom_fichier, "w") do file
	println(file, "Les paramètres sont : ")
	println(file, "\t- Pour F1 : p = ", p, ", q = ", q)
	println(file, "\t- Pour F2 : α = ", α, ", β = ", β, ", γ = ", γ, ", δ = ", δ, "\n")
end

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
	
	minPreJ1 = Inf
	minPostJ1 = Inf
	minTMJ1 = Inf
	
	minPreGS1 = Inf
	minPostGS1 = Inf
	minTMGS1 = Inf
	
	minPreJ2 = Inf
	minPostJ2 = Inf
	minTMJ2 = Inf
	
	minPreGS2 = Inf
	minPostGS2 = Inf
	minTMGS2 = Inf
	
	# nom_fichier = string("data/Version3/Temps_des_methodes_multigrid_avec_N_=_", string(N), ".txt")
	# open(nom_fichier, "w") do file
		for i = 1:taille
			@show i
			pre = tabPre[i]
			post = tabPost[i]
			TMJ1 = @elapsed MultigridJ(A, F1, N, x0, pre, post)
			TMGS1 = @elapsed MultigridGS(A, F1, N, x0, pre, post)
			TMJ2 = @elapsed MultigridJ(A, F2, N, x0, pre, post)
			TMGS2 = @elapsed MultigridGS(A, F2, N, x0, pre, post)

			if TMJ1 < minTMJ1
				minPreJ1 = pre
				minPostJ1 = post
				minTMJ1 = TMJ1
			end

			if TMGS1 < minTMGS1
				minPreGS1 = pre
				minPostGS1 = post
				minTMGS1 = TMGS1
			end

			if TMJ2 < minTMJ2
				minPreJ2 = pre
				minPostJ2 = post
				minTMJ2 = TMJ2
			end

			if TMGS2 < minTMGS2
				minPreGS2 = pre
				minPostGS2 = post
				minTMGS2 = TMGS2
			end

			# println(file, pre, "\t", post, "\t", TMJ1, "\t", TMGS1, "\t", TMJ2, "\t", TMGS2)

		end
	# end

	open(nom_fichier, "a") do file
		println(file, "\tPour F1 : ")
		println(file, "\t\tAvec Jacobi : ")
		println(file, "\t\t- pre = ", minPreJ1)
		println(file, "\t\t- post = ", minPostJ1)
		println(file, "\t\t- temps = ", minTMJ1, "\n")
		
		println(file, "\t\tAvec Gauss-Seidel : ")
		println(file, "\t\t- pre = ", minPreGS1)
		println(file, "\t\t- post = ", minPostGS1)
		println(file, "\t\t- temps = ", minTMGS1, "\n")
		
		println(file, "\tPour F2 : ")
		println(file, "\t\tAvec Jacobi : ")
		println(file, "\t\t- pre = ", minPreJ2)
		println(file, "\t\t- post = ", minPostJ2)
		println(file, "\t\t- temps = ", minTMJ2, "\n")
		
		println(file, "\t\tAvec Gauss-Seidel : ")
		println(file, "\t\t- pre = ", minPreGS2)
		println(file, "\t\t- post = ", minPostGS2)
		println(file, "\t\t- temps = ", minTMGS2, "\n", "\n")
	end
end

println("\n")
println("α = ", α, ", β = ", β, ", γ = ", γ, ", δ = ", δ)
println("p = ", p, ", q = ", q)