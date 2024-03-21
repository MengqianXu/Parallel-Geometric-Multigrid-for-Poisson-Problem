using LinearAlgebra
using Random
using Plots

function Erreur(A, F, N, p, q, w = 0.5)
	x0 = rand(1:9, (N*N, 1))
	sol_direct = A\F
	vrai_sol = [u(p, q, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
	
	resJ = Jacobi(A, F, x0)
	resGS = GaussSeidel(A, F, x0)
	resS = SOR(A, F, w, x0)
	
	relativeJ = norm(sol_direct - resJ)
	relativeGS = norm(sol_direct - resGS)
	relativeS = norm(sol_direct - resS)

	absolueJ = norm(vrai_sol - resJ)
	absolueGS = norm(vrai_sol - resGS)
	absolueS = norm(vrai_sol - resS)

	return absolueJ, absolueGS, absolueS, relativeJ, relativeGS, relativeS
end

function Erreur(A, F, N, α, β, γ, δ, w = 0.5)
	x0 = rand(1:9, (N*N, 1))
	sol_direct = A\F
	vrai_sol = [u(α, β, γ, δ, i/(N + 1), j /(N + 1)) for j = 1:N, i = 1:N][:]
	
	resJ = Jacobi(A, F, x0)
	resGS = GaussSeidel(A, F, x0)
	resS = SOR(A, F, w, x0)
	
	relativeJ = norm(sol_direct - resJ)
	relativeGS = norm(sol_direct - resGS)
	relativeS = norm(sol_direct - resS)

	absolueJ = norm(vrai_sol - resJ)
	absolueGS = norm(vrai_sol - resGS)
	absolueS = norm(vrai_sol - resS)

	return absolueJ, absolueGS, absolueS, relativeJ, relativeGS, relativeS
end

taille = 6
tabEAJ1 = zeros(taille)
tabEAGS1 = zeros(taille)
tabEAS1 = zeros(taille)
tabERJ1 = zeros(taille)
tabERGS1 = zeros(taille)
tabERS1 = zeros(taille)
tabEAJ2 = zeros(taille)
tabEAGS2 = zeros(taille)
tabEAS2 = zeros(taille)
tabERJ2 = zeros(taille)
tabERGS2 = zeros(taille)
tabERS2 = zeros(taille)
tabN = zeros(taille)

for i in 0:5
	@show i
	# index = Int64(i/10)
	index = i + 1
	N = 2^i
	tabN[index] = N
    
	A = Creer_A(N)
    F1 = Creer_F(p, q, N)
    F2 = Creer_F(α, β, γ, δ, N) 
    
	absolueJ1, absolueGS1, absolueS1, relativeJ1, relativeGS1, relativeS1 = Erreur(A, F1, N, p, q)
    absolueJ2, absolueGS2, absolueS2, relativeJ2, relativeGS2, relativeS2 = Erreur(A, F2, N, α, β, γ, δ)
    
	tabEAJ1[index] = absolueJ1
	tabEAGS1[index] = absolueGS1
	tabEAS1[index] = absolueS1
	
	tabERJ1[index] = relativeJ1
	tabERGS1[index] = relativeGS1
	tabERS1[index] = relativeS1
	
	tabEAJ2[index] = absolueJ2
	tabEAGS2[index] = absolueGS2
	tabEAS2[index] = absolueS2
	
	tabERJ2[index] = relativeJ2
	tabERGS2[index] = relativeGS2
	tabERS2[index] = relativeS2
	
end

p1 = plot(tabN, [tabEAJ1, tabEAGS1, tabEAS1], label = ["Jacobi" "Gauss-Seidel" "SOR"])
p2 = plot(tabN, [tabERJ1, tabERGS1, tabERS1], label = ["Jacobi" "Gauss-Seidel" "SOR"])
p3 = plot(tabN, [tabEAJ2, tabEAGS2, tabEAS2], label = ["Jacobi" "Gauss-Seidel" "SOR"])
p4 = plot(tabN, [tabERJ2, tabERGS2, tabERS2], label = ["Jacobi" "Gauss-Seidel" "SOR"])

savefig(p1, "Erreur absolue F1.png")
savefig(p2, "Erreur relative F1.png")
savefig(p3, "Erreur absolue F2.png")
savefig(p4, "Erreur relative F2.png")


i = 3
A = Creer_A(i)
F1 = Creer_F(p, q, i)
F2 = Creer_F(α, β, γ, δ, i) 
absolueJ1, absolueGS1, absolueS1, relativeJ1, relativeGS1, relativeS1 = Erreur(A, F1, i, p, q)
absolueJ2, absolueGS2, absolueS2, relativeJ2, relativeGS2, relativeS2 = Erreur(A, F2, i, α, β, γ, δ)

i = 10
A = Creer_A(i)
F1 = Creer_F(p, q, i)
F2 = Creer_F(α, β, γ, δ, i) 
absolueJ1, absolueGS1, absolueS1, relativeJ1, relativeGS1, relativeS1 = Erreur(A, F1, i, p, q)
absolueJ2, absolueGS2, absolueS2, relativeJ2, relativeGS2, relativeS2 = Erreur(A, F2, i, α, β, γ, δ)