using LinearAlgebra
using SparseArrays
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

function Creer_U(p, q, N)
	U = zeros(N*N)
	for i = 1:N
		x = i/(N + 1)
		for j = 1:N
			y = j/(N + 1)
			U[N*(i - 1) + j] = u(p, q, x, y)
		end
	end
	return U
end

function Creer_U(α, β, γ, δ, N)
	U = zeros(N*N)
	for i = 1:N
		x = i/(N + 1)
		for j = 1:N
			y = j/(N + 1)
			U[N*(i - 1) + j] = u(α, β, γ, δ, x, y)
		end
	end
	return U
end

function Creer_A(N)
	h = 1/(N + 1)
	return ((1/h)^2)sparse(kron(I(N), Tridiagonal(-ones(N, N)) + 5I(N)) + kron(Tridiagonal(-ones(N, N)) + I(N), I(N)))
end

function Creer_F(p, q, N)
	F = zeros(N*N)
	for i = 1:N
		x = i/(N + 1)
		for j = 1:N
			y = j/(N + 1)
			F[N*(i - 1) + j] = u(p, q, x, y)*((pi^2)*((p^2) + (q^2)))
		end
	end
	return F
	# return ((1/(N + 1))^2)F
end

function Creer_F(α, β, γ, δ, N)
	F = zeros(N*N)
	for i = 1:N
		x = i/(N + 1)
		for j = 1:N
			y = j/(N + 1)
			F[N*(i - 1) + j] = 2*α*β*(x^(α - 1))*(y^γ)*((1 - x)^(β - 1))*((1 - y)^δ) - α*(x^(α - 2))*(y^γ)*(α - 1)*((1 - x)^β)*((1 - y)^δ) - β*(x^α)*(y^γ)*(β - 1)*((1 - x)^(β - 2))*((1 - y)^δ) - γ*(x^α)*(y^(γ - 2))*(γ - 1)*((1 - x)^β)*((1 - y)^δ) + 2*γ*δ*(x^α)*(y^(γ - 1))*((1 - x)^β)*((1 - y)^(δ - 1)) - δ*(x^α)*(y^γ)*(δ - 1)*((1 - x)^β)*((1 - y)^(δ - 2))
		end
	end
	return F
	# return ((1/(N + 1))^2)F
end