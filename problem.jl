using LinearAlgebra
using Random

function u(p, q, x, y)
	function g(x, y)
		return sin(p*pi*x)*sin(q*pi*y)
	end
	return g(x, y)
end

function u(α, β, γ, δ, x, y)
	function g(x, y)
		return (x^α)*((1 - x)^β)*(y^γ)*((1 - y)^δ)
	end
	return g(x, y)
end

function Creer_A(N)
	return kron(I(N), Tridiagonal(-ones(N, N)) + 5I(N)) + kron(Tridiagonal(-ones(N, N)) + I(N), I(N))
end

function Creer_F(p, q, N)
	F = zeros(N*N)
	for i = 1:N
		x = i/(N + 1)
		for j = 1:N
			y = j/(N + 1)
			F[N*(i - 1) + j] = u(p, q, x, y)
		end
	end
	return ((1/(N + 1))^2)F
end

function Creer_F(α, β, γ, δ, N)
	F = zeros(N*N)
	for i = 1:N
		x = i/(N + 1)
		for j = 1:N
			y = j/(N + 1)
			F[N*(i - 1) + j] = u(α, β, γ, δ, x, y)
		end
	end
	return ((1/(N + 1))^2)F
end

N = 3
p = rand(-5:5) 
q = rand(-5:5)
α = rand(1:10) 
β = rand(1:10) 
γ = rand(1:10) 
δ = rand(1:10)

A = Creer_A(N)
F1 = Creer_F(p, q, N)
F2 = Creer_F(α, β, γ, δ, N)
