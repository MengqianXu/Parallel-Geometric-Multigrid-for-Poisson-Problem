using LinearAlgebra
using Random
using Plots

function Convergence(A, F, N, p, q, w = 0.5)
	x0 = rand(1:9, (N*N, 1))
	sol = A\F 
	
	resJ = Jacobi(A, F, x0)
	resGS = GaussSeidel(A, F, x0)
	resS = SOR(A, F, w, x0)
	
	normeJ = norm(sol - resJ)
	normeGS = norm(sol - resGS)
	normeS = norm(sol - resS)

	eJ = 0
	eGS = 0
	eS = 0

	for i = 1:N
		x = i/(N + 1)
		for j = 1:N
			y = j/(N + 1)
			eJ += abs(u(p, q, x, y) - resJ[(i - 1)*N + j])^2
			println(eJ)
			eGS += abs(u(p, q, x, y) - resGS[(i - 1)*N + j])^2
			println(eGS)
			eS += abs(u(p, q, x, y) - resS[(i - 1)*N + j])^2
			println(eS)
		end
	end

	return eJ, eGS, eS, normeJ, normeGS, normeS
end

tabJ1 = zeros(9)
tabGS1 = zeros(9)
tabS1 = zeros(9)

tabNJ1 = zeros(9)
tabNGS1 = zeros(9)
tabNS1 = zeros(9)

tabJ2 = zeros(9)
tabGS2 = zeros(9)
tabS2 = zeros(9)

tabNJ2 = zeros(9)
tabNGS2 = zeros(9)
tabNS2 = zeros(9)

tabN = zeros(9)
for i in 2:10
	tabN[i - 1] = i
    A = Creer_A(i)
    F1 = Creer_F(p, q, i)
    F2 = Creer_F(α, β, γ, δ, i) 
    eJ1, eGS1, eS1, normeJ1, normeGS1, normeS1 = Convergence(A, F1, i, p, q)
    eJ2, eGS2, eS2, normeJ2, normeGS2, normeS2 = Convergence(A, F2, i, p, q)
    
	tabJ1[i - 1] = eJ1
    tabJ2[i - 1] = eJ2
    tabGS1[i - 1] = eGS1
    tabGS2[i - 1] = eGS2
    tabS1[i - 1] = eS1
    tabS2[i - 1] = eS2

	tabNJ1[i - 1] = normeJ1
	tabNJ2[i - 1] = normeJ2
	tabNGS1[i - 1] = normeGS1
	tabNGS2[i - 1] = normeGS2
	tabNS1[i - 1] = normeS1
	tabNS2[i - 1] = normeS2
	
end

plot(tabN, [tabJ1, tabGS1, tabS1], label = ["Jacobi" "Gauss-Seidel" "SOR"])
plot(tabN, [tabNJ1, tabNGS1, tabNS1], label = ["Jacobi" "Gauss-Seidel" "SOR"])

i = 10
A = Creer_A(i)
F1 = Creer_F(p, q, i)
F2 = Creer_F(α, β, γ, δ, i) 
eJ1, eGS1, eS1, normeJ1, normeGS1, normeS1 = Convergence(A, F1, i, p, q)
eJ2, eGS2, eS2, normeJ2, normeGS2, normeS2 = Convergence(A, F2, i, p, q)