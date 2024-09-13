using LinearAlgebra
using SparseArrays
using Random
using Plots

include("problem.jl")
include("solvers.jl")

function decomposition_2_facteurs(n)
	i = floor(Int64, sqrt(n))
	while n%i != 0
		i -= 1
	end
	j = floor(Int64, n/i)
	if i > j
		return i, j
	else
		return j, i
	end
end

α = 6
β = 1
γ = 5
δ = 10
p = 2
q = 2

@elapsed decomposition_2_facteurs(p*q*10)

taille = 7
tabA = zeros(taille)
tabN = zeros(taille)
tabN2 = zeros(taille)
tabN4 = zeros(taille)
tabN16 = zeros(taille)
tabH = zeros(taille)
tabH2 = zeros(taille)

tabED1 = zeros(taille)
tabEJ1 = zeros(taille)
tabEGS1 = zeros(taille)
tabES1 = zeros(taille)
tabEMJ1 = zeros(taille)
tabEMGS1 = zeros(taille)

tabED2 = zeros(taille)
tabEJ2 = zeros(taille)
tabEGS2 = zeros(taille)
tabES2 = zeros(taille)
tabEMJ2 = zeros(taille)
tabEMGS2 = zeros(taille)

tabRJ1 = zeros(taille)
tabRGS1 = zeros(taille)
tabRS1 = zeros(taille)
tabRMJ1 = zeros(taille)
tabRMGS1 = zeros(taille)

tabRJ2 = zeros(taille)
tabRGS2 = zeros(taille)
tabRS2 = zeros(taille)
tabRMJ2 = zeros(taille)
tabRMGS2 = zeros(taille)

tabIJ1 = zeros(taille)
tabIGS1 = zeros(taille)
tabIS1 = zeros(taille)

tabIJ2 = zeros(taille)
tabIGS2 = zeros(taille)
tabIS2 = zeros(taille)

tabTD1 = zeros(taille)
tabTJ1 = zeros(taille)
tabTGS1 = zeros(taille)
tabTS1 = zeros(taille)
tabTMJ1 = zeros(taille)
tabTMGS1 = zeros(taille)

tabTD2 = zeros(taille)
tabTJ2 = zeros(taille)
tabTGS2 = zeros(taille)
tabTS2 = zeros(taille)
tabTMJ2 = zeros(taille)
tabTMGS2 = zeros(taille)

println("Début de la boucle")
for n = 1:taille
	@show n
	N = 2^n - 1
	N2 = N^2
	N4 = N2^2
	N16 = N4^2
	h = 1/N
	h2 = 1/N2
	
	A = Creer_A(N) # Pose un problème à partir de n = 8
	F1 = Creer_F(p, q, N)
	F2 = Creer_F(α, β, γ, δ, N)
	U1 = Creer_U(p, q, N)
	U2 = Creer_U(α, β, γ, δ, N)
	x0 = ones(N*N, 1)
	w = 0.5

	tabY = [i/(N + 1) for i = 1:N, j = 1:N][:]
	tabX = [j/(N + 1) for i = 1:N, j = 1:N][:]

	println("Résolution du problème F1")
    println("\t- Direct")
	sol1 = A\F1
	TD1 = @elapsed A\F1
	
	println("\t- Jacobi")
	resJ1, IJ1 = Jacobi(A, F1, x0)
	TJ1 = @elapsed Jacobi(A, F1, x0)
	
	println("\t- Gauss Seidel")
	resGS1, IGS1 = GaussSeidel(A, F1, x0)
	TGS1 = @elapsed GaussSeidel(A, F1, x0)
	
	println("\t- SOR")
	resS1, IS1 = SOR(A, F1, w, x0)
	TS1 = @elapsed SOR(A, F1, w, x0)
	
	println("\t- Multigrid avec Jacobi")
	resMJ1 = MultigridJ(F1, N, 1, 1)
	TMJ1 = @elapsed MultigridJ(F1, N, 1, 1)
	
	println("\t- Multigrid avec Gauss Seidel")
	resMGS1 = MultigridGS(F1, N, 1, 1)
	TMGS1 = @elapsed MultigridGS(F1, N, 1, 1)

	println("Résolution du problème F2")
	println("\t- Direct")
	sol2 = A\F2
	TD2 = @elapsed A\F2
	
	println("\t- Jacobi")
	resJ2, IJ2 = Jacobi(A, F2, x0)
	TJ2 = @elapsed Jacobi(A, F2, x0)
	
	println("\t- Gauss Seidel")
	resGS2, IGS2 = GaussSeidel(A, F2, x0)
	TGS2 = @elapsed GaussSeidel(A, F2, x0)
	
	println("\t- SOR")
	resS2, IS2 = SOR(A, F2, w, x0)
	TS2 = @elapsed SOR(A, F2, w, x0)
	
	println("\t- Multigrid avec Jacobi")
	resMJ2 = MultigridJ(F2, N, 1, 1)
	TMJ2 = @elapsed MultigridJ(F2, N, 1, 1)
	
	println("\t- Multigrid avec Gauss Seidel")
	resMGS2 = MultigridGS(F2, N, 1, 1)
	TMGS2 = @elapsed MultigridGS(F2, N, 1, 1)


	println("Calcul des erreurs et des résidus")
	ED1 = U1 - sol1
	EJ1 = U1 - resJ1
	EGS1 = U1 - resGS1
	ES1 = U1 - resS1
	EMJ1 = U1 - resMJ1
	EMGS1 = U1 - resMGS1

	ED2 = U2 - sol2
	EJ2 = U2 - resJ2
	EGS2 = U2 - resGS2
	ES2 = U2 - resS2
	EMJ2 = U2 - resMJ2
	EMGS2 = U2 - resMGS2

	RJ1 = F1 - A*resJ1
	RGS1 = F1 - A*resGS1
	RS1 = F1 - A*resS1
	RMJ1 = F1 - A*resMJ1
	RMGS1 = F1 - A*resMGS1

	RJ2 = F2 - A*resJ2
	RGS2 = F2 - A*resGS2
	RS2 = F2 - A*resS2
	RMJ2 = F2 - A*resMJ2
	RMGS2 = F2 - A*resMGS2

    # if n > 1
	# 	println("Plot et sauvegarde des solutions")
	# 	p11 = surface(tabX, tabY, U1, title = string("Exact solution for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p11, string("image/Julia/Exact solution for F1 with N = ", string(N), ".png"))
	# 	p12 = surface(tabX, tabY, sol1, title = string("Direct solver for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p12, string("image/Julia/Direct solver for F1 with N = ", string(N), ".png"))
	# 	p13 = surface(tabX, tabY, resJ1,  title = string("Jacobi for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p13, string("image/Julia/Jacobi for F1 with N = ", string(N), ".png"))
	# 	p14 = surface(tabX, tabY, resGS1, title = string("Gauss Seidel for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p14, string("image/Julia/Gauss Seidel for F1 with N = ", string(N), ".png"))
	# 	p15 = surface(tabX, tabY, resS1, title = string("SOR for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p15, string("image/Julia/SOR for F1 with N = ", string(N), ".png"))
	# 	p16 = surface(tabX, tabY, resMJ1, title = string("Multigrid Jacobi for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p16, string("image/Julia/Multigrid Jacobi for F1 with N = ", string(N), ".png"))
	# 	p17 = surface(tabX, tabY, resMGS1, title = string("Multigrid Gauss Seidel for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p17, string("image/Julia/Multigrid Gauss Seidel for F1 with N = ", string(N), ".png"))

	# 	p21 = surface(tabX, tabY, U2, title = string("Exact solution for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p21, string("image/Julia/Exact solution for F2 with N = ", string(N), ".png"))
	# 	p22 = surface(tabX, tabY, sol2, title = string("Direct solver for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p22, string("image/Julia/Direct solver for F2 with N = ", string(N), ".png"))
	# 	p23 = surface(tabX, tabY, resJ2,  title = string("Jacobi for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p23, string("image/Julia/Jacobi for F2 with N = ", string(N), ".png"))
	# 	p24 = surface(tabX, tabY, resGS2, title = string("Gauss Seidel for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p24, string("image/Julia/Gauss Seidel for F2 with N = ", string(N), ".png"))
	# 	p25 = surface(tabX, tabY, resS2, title = string("SOR for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p25, string("image/Julia/SOR for F2 with N = ", string(N), ".png"))
	# 	p26 = surface(tabX, tabY, resMJ2, title = string("Multigrid Jacobi for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p26, string("image/Julia/Multigrid Jacobi for F2 with N = ", string(N), ".png"))
	# 	p27 = surface(tabX, tabY, resMGS2, title = string("Multigrid Gauss Seidel for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "f(x, y)")
	# 	savefig(p27, string("image/Julia/Multigrid Gauss Seidel for F2 with N = ", string(N), ".png"))


	# 	println("Plot et sauvegarde des erreurs")
	# 	p31 = surface(tabX, tabY, ED1, title = string("Error of direct solver for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p31, string("image/Julia/Error of direct solver for F1 with N = ", string(N), ".png"))
	# 	p32 = surface(tabX, tabY, EJ1, title = string("Error of Jacobi for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p32, string("image/Julia/Error of Jacobi for F1 with N = ", string(N), ".png"))
	# 	p33 = surface(tabX, tabY, EGS1,  title = string("Error of Gauss Seidel for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p33, string("image/Julia/Error of Gauss Seidel for F1 with N = ", string(N), ".png"))
	# 	p34 = surface(tabX, tabY, ES1, title = string("Error of SOR for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p34, string("image/Julia/Error of SOR for F1 with N = ", string(N), ".png"))
	# 	p35 = surface(tabX, tabY, EMJ1, title = string("Error of Multigrid Jacobi for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p35, string("image/Julia/Error of Multigrid Jacobi for F1 with N = ", string(N), ".png"))
	# 	p36 = surface(tabX, tabY, EMGS1, title = string("Error of Multigrid Gauss Seidel for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p36, string("image/Julia/Error of Multigrid Gauss Seidel for F1 with N = ", string(N), ".png"))
		
	# 	p41 = surface(tabX, tabY, ED2, title = string("Error of direct solver for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p41, string("image/Julia/Error of direct solver for F2 with N = ", string(N), ".png"))
	# 	p42 = surface(tabX, tabY, EJ2, title = string("Error of Jacobi for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p42, string("image/Julia/Error of Jacobi for F2 with N = ", string(N), ".png"))
	# 	p43 = surface(tabX, tabY, EGS2,  title = string("Error of Gauss Seidel for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p43, string("image/Julia/Error of Gauss Seidel for F2 with N = ", string(N), ".png"))
	# 	p44 = surface(tabX, tabY, ES2, title = string("Error of SOR for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p44, string("image/Julia/Error of SOR for F2 with N = ", string(N), ".png"))
	# 	p45 = surface(tabX, tabY, EMJ2, title = string("Error of Multigrid Jacobi for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p45, string("image/Julia/Error of Multigrid Jacobi for F2 with N = ", string(N), ".png"))
	# 	p46 = surface(tabX, tabY, EMGS2, title = string("Error of Multigrid Gauss Seidel for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "error")
	# 	savefig(p46, string("image/Julia/Error of Multigrid Gauss Seidel for F2 with N = ", string(N), ".png"))


	# 	println("Plot et sauvegarde des residus")
	# 	p51 = surface(tabX, tabY, RJ1, title = string("Residual of Jacobi for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "residual")
	# 	savefig(p51, string("image/Julia/Residual of Jacobi for F1 with N = ", string(N), ".png"))
	# 	p52 = surface(tabX, tabY, RGS1,  title = string("Residual of Gauss Seidel for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "residual")
	# 	savefig(p52, string("image/Julia/Residual of Gauss Seidel for F1 with N = ", string(N), ".png"))
	# 	p53 = surface(tabX, tabY, RS1, title = string("Residual of SOR for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "residual")
	# 	savefig(p53, string("image/Julia/Residual of SOR for F1 with N = ", string(N), ".png"))
	# 	p54 = surface(tabX, tabY, RMJ1, title = string("Residual of Multigrid Jacobi for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "residual")
	# 	savefig(p54, string("image/Julia/Residual of Multigrid Jacobi for F1 with N = ", string(N), ".png"))
	# 	p55 = surface(tabX, tabY, RMGS1, title = string("Residual of Multigrid Gauss Seidel for F1 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "residual")
	# 	savefig(p55, string("image/Julia/Residual of Multigrid Gauss Seidel for F1 with N = ", string(N), ".png"))
		
	# 	p61 = surface(tabX, tabY, RJ2, title = string("Residual of Jacobi for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "residual")
	# 	savefig(p61, string("image/Julia/Residual of Jacobi for F2 with N = ", string(N), ".png"))
	# 	p62 = surface(tabX, tabY, RGS2,  title = string("Residual of Gauss Seidel for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "residual")
	# 	savefig(p62, string("image/Julia/Residual of Gauss Seidel for F2 with N = ", string(N), ".png"))
	# 	p63 = surface(tabX, tabY, RS2, title = string("Residual of SOR for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "residual")
	# 	savefig(p63, string("image/Julia/Residual of SOR for F2 with N = ", string(N), ".png"))
	# 	p64 = surface(tabX, tabY, RMJ2, title = string("Residual of Multigrid Jacobi for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "residual")
	# 	savefig(p64, string("image/Julia/Residual of Multigrid Jacobi for F2 with N = ", string(N), ".png"))
	# 	p65 = surface(tabX, tabY, RMGS2, title = string("Residual of Multigrid Gauss Seidel for F2 with N = ", string(N)), xlabel = "x", ylabel = "y", zlabel = "residual")
	# 	savefig(p65, string("image/Julia/Residual of Multigrid Gauss Seidel for F2 with N = ", string(N), ".png"))
	# end

	println("Sauvegarde des données dans les tableaux")
	i = n
	# tabA[i] = cond(Matrix(A))
	tabN[i] = N
	tabN2[i] = N2
	tabN4[i] = N4
	tabN16[i] = N16
	tabH[i] = h
	tabH2[i] = h2
	
	println("\t- Erreur pour F1")
	tabED1[i] = norm(ED1)/norm(sol1)
	tabEJ1[i] = norm(EJ1)/norm(resJ1)
	tabEGS1[i] = norm(EGS1)/norm(resGS1)
	tabES1[i] = norm(ES1)/norm(resS1)
	tabEMJ1[i] = norm(EMJ1)/norm(resMJ1)
	tabEMGS1[i] = norm(EMGS1)/norm(resMGS1)

	println("\t- Erreur pour F2")
	tabED2[i] = norm(ED2)/norm(sol2)
	tabEJ2[i] = norm(EJ2)/norm(resJ2)
	tabEGS2[i] = norm(EGS2)/norm(resGS2)
	tabES2[i] = norm(ES2)/norm(resS2)
	tabEMJ2[i] = norm(EMJ2)/norm(resMJ2)
	tabEMGS2[i] = norm(EMGS2)/norm(resMGS2)

	println("\t- Residu pour F1")
	tabRJ1[i] = norm(RJ1)/norm(A*resJ1)
	tabRGS1[i] = norm(RGS1)/norm(A*resGS1)
	tabRS1[i] = norm(RS1)/norm(A*resS1)
	tabRMJ1[i] = norm(RMJ1)/norm(A*resMJ1)
	tabRMGS1[i] = norm(RMGS1)/norm(A*resMGS1)

	println("\t- Residu pour F2")
	tabRJ2[i] = norm(RJ2)/norm(A*resJ2)
	tabRGS2[i] = norm(RGS2)/norm(A*resGS2)
	tabRS2[i] = norm(RS2)/norm(A*resS2)
	tabRMJ2[i] = norm(RMJ2)/norm(A*resMJ2)
	tabRMGS2[i] = norm(RMGS2)/norm(A*resMGS2)

	println("\t- Iterations")
	tabIJ1[i] = IJ1
	tabIGS1[i] = IGS1
	tabIS1[i] = IS1

	tabIJ2[i] = IJ2
	tabIGS2[i] = IGS2
	tabIS2[i] = IS2

	println("\t- Temps pour F1")
	tabTD1[n] = TD1
	tabTJ1[n] = TJ1
	tabTGS1[n] = TGS1
	tabTS1[n] = TS1
	tabTMJ1[n] = TMJ1
	tabTMGS1[n] = TMGS1
	
	println("\t- Temps pour F2")
	tabTD2[n] = TD2
	tabTJ2[n] = TJ2
	tabTGS2[n] = TGS2
	tabTS2[n] = TS2
	tabTMJ2[n] = TMJ2
	tabTMGS2[n] = TMGS2
end
println("Fin de la boucle")
println("Ecriture des données dans un fichier texte")
nom_fichier = "data/Version3/Plots2.txt"
open(nom_fichier, "w") do file
	for i = 2:taille
		println(file, tabA[i], "\t", tabN[i], "\t", tabN2[i], "\t", tabN4[i], "\t", tabN16[i], "\t", tabH[i], "\t", tabH2[i], "\t", tabED1[i], "\t", tabEJ1[i], "\t", tabEGS1[i], "\t", tabES1[i], "\t", tabEMJ1[i], "\t", tabEMGS1[i], "\t", tabED2[i], "\t", tabEJ2[i], "\t", tabEGS2[i], "\t", tabES2[i], "\t", tabEMJ2[i], "\t", tabEMGS2[i], "\t", tabRJ1[i], "\t", tabRGS1[i], "\t", tabRS1[i], "\t", tabRMJ1[i], "\t", tabRMGS1[i], "\t", tabRJ2[i], "\t", tabRGS2[i], "\t", tabRS2[i], "\t", tabRMJ2[i], "\t", tabRMGS2[i], "\t", tabIJ1[i], "\t", tabIGS1[i], "\t", tabIS1[i], "\t", tabIJ2[i], "\t", tabIGS2[i], "\t", tabIS2[i], "\t", tabTD1[i], "\t", tabTJ1[i], "\t", tabTGS1[i], "\t", tabTS1[i], "\t", tabTMJ1[i], "\t", tabTMGS1[i], "\t", tabTD2[i], "\t", tabTJ2[i], "\t", tabTGS2[i], "\t", tabTS2[i], "\t", tabTMJ2[i], "\t", tabTMGS2[i])
	end
end
println("Fin de l'écriture des données dans un fichier texte")

# println("\t- Derniers plots")
# p1 = plot(tabN, [tabN, tabN2, tabA], label = ["N" "N^2" "Cond(A)"], title = "Matrix A conditioning", xlabel = "N = discretization size", ylabel = "Conditioning number", xscale=:log2, yscale=:log2)
# p2 = plot(tabN2, [tabN2, tabN4, tabA], label = ["N" "N^2" "Cond(A)"], title = "Matrix A conditioning", xlabel = "N = domaine Ω size", ylabel = "Conditioning number", xscale=:log2, yscale=:log2)
# p3 = plot(tabN16, [tabN4, tabN16, tabA], label = ["N" "N^2" "Cond(A)"], title = "Matrix A conditioning", xlabel = "N = matrix size", ylabel = "Conditioning number", xscale=:log2, yscale=:log2)

# p4 = plot(tabN[2:taille], [tabH[2:taille], tabH2[2:taille], tabED1[2:taille], tabEJ1[2:taille], tabEGS1[2:taille], tabES1[2:taille], tabEMJ1[2:taille], tabEMGS1[2:taille]], label = ["H" "H^2" "Direct" "Jacobi" "Gauss-Seidel" "SOR" "Multigrid_Jacobi" "Multigrid_Gauss-Seidel"], title = "Error for F1", xlabel = "N", ylabel = "Error", xscale=:log2, yscale=:log2)
# p5 = plot(tabN[2:taille], [tabH[2:taille], tabH2[2:taille], tabRJ1[2:taille], tabRGS1[2:taille], tabRS1[2:taille], tabRMJ1[2:taille], tabRMGS1[2:taille]], label = ["H" "H^2" "Jacobi" "Gauss-Seidel" "SOR" "Multigrid_Jacobi" "Multigrid_Gauss-Seidel"], title = "Residual for F1", xlabel = "N", ylabel = "Residual", xscale=:log2, yscale=:log2)
# p6 = plot(tabN[2:taille], [tabH[2:taille], tabH2[2:taille], tabED2[2:taille], tabEJ2[2:taille], tabEGS2[2:taille], tabES2[2:taille], tabEMJ2[2:taille], tabEMGS2[2:taille]], label = ["H" "H^2" "Direct" "Jacobi" "Gauss-Seidel" "SOR" "Multigrid_Jacobi" "Multigrid_Gauss-Seidel"], title = "Error for F2", xlabel = "N", ylabel = "Error", xscale=:log2, yscale=:log2)
# p7 = plot(tabN[2:taille], [tabH[2:taille], tabH2[2:taille], tabRJ2[2:taille], tabRGS2[2:taille], tabRS2[2:taille], tabRMJ2[2:taille], tabRMGS2[2:taille]], label = ["H" "H^2" "Jacobi" "Gauss-Seidel" "SOR" "Multigrid_Jacobi" "Multigrid_Gauss-Seidel"], title = "Residual for F2", xlabel = "N", ylabel = "Residual", xscale=:log2, yscale=:log2)

# p8 = plot(tabN[2:taille], [tabN[2:taille], tabN2[2:taille], tabIJ1[2:taille], tabIGS1[2:taille], tabIS1[2:taille]], label = ["N" "N^2" "Jacobi" "Gauss-Seidel" "SOR"], title = "Number of iterations depending on N for F1", xlabel = "N", ylabel = "Number of iterations")
# p9 = plot(tabN[2:taille], [tabN[2:taille], tabN2[2:taille], tabIJ2[2:taille], tabIGS2[2:taille], tabIS2[2:taille]], label = ["N" "N^2" "Jacobi" "Gauss-Seidel" "SOR"], title = "Number of iterations depending on N for F2", xlabel = "N", ylabel = "Number of iterations")

# p10 = plot(tabN[2:taille], [tabN[2:taille], tabN2[2:taille], tabTD1[2:taille], tabTJ1[2:taille], tabTGS1[2:taille], tabTS1[2:taille], tabTMJ1[2:taille], tabTMGS1[2:taille]], label = ["N" "N^2" "Jacobi" "Gauss-Seidel" "SOR" "Multigrid_Jacobi" "Multigrid_Gauss-Seidel"], title = "Execution time depending on N for F1", xlabel = "N", ylabel = "Execution time (s)")
# p11 = plot(tabN[2:taille], [tabN[2:taille], tabN2[2:taille], tabTD2[2:taille], tabTJ2[2:taille], tabTGS2[2:taille], tabTS2[2:taille], tabTMJ2[2:taille], tabTMGS2[2:taille]], label = ["N" "N^2" "Jacobi" "Gauss-Seidel" "SOR" "Multigrid_Jacobi" "Multigrid_Gauss-Seidel"], title = "Execution time depending on N for F2", xlabel = "N", ylabel = "Execution time (s)")

# println("\t- Dernieres sauvegardes")
# savefig(p1, "image/Julia/Conditionnement de A (N).png")
# savefig(p2, "image/Julia/Conditionnement de A (N2).png")
# savefig(p3, "image/Julia/Conditionnement de A (N4).png")
# savefig(p4, "Erreur pour F1.png")
# savefig(p5, "Residu pour F1.png")
# savefig(p6, "Erreur pour F2.png")
# savefig(p7, "Residu pour F2.png")
# savefig(p8, "Iterations pour F1.png")
# savefig(p9, "Iterations pour F2.png")
# savefig(p10, "Temps pour F1.png")
# savefig(p11, "Temps pour F2.png")
# [2:taille]