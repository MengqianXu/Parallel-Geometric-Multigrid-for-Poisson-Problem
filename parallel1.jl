using Distributed
using LinearAlgebra

function Multigrid(A, F, N::Int64, x0, pre, post)
    if N == 1
        return (1/4)*F
    else
        U = copy(x0)

        @sync begin
            @distributed for i = 1:pre
                U = Jacobi(A, F, U)
            end
        end

        newN = floor(Int64, N/2)
        newA = Creer_A(newN)
        newF = zeros(newN*newN)
        newU = zeros(newN*newN)

        for i = 1:newN
            for j = 1:newN
                newF[j + N*(i - 1)] = F[2j + N*(2i - 1)]
                newU[j + N*(i - 1)] = U[2j + N*(2i - 1)]
            end
        end

        newU = Multigrid(newA, newF, newN, newU, pre, post)

        for i = 1:newN
            for j = 1:newN
                U[2j + N*(2i - 1)] = newU[j + N*(i - 1)]
            end
        end

        @sync begin
            @distributed for i = 1:pre
                U = Jacobi(A, F, U)
            end
        end

        return U
    end
end


%Allocate each Jacobi iteration task to multiple processors for parallel execution
function Jacobi(A, b, x0, tol = 10^(-15))
    ancien = copy(x0)
    nouveau = copy(b)
    n = length(b)
    
    for i = 1:n
        for j = 1:n
            if j != i 
                nouveau[i] -= A[i, j]*ancien[j]
            end
        end
        nouveau[i] /= A[i, i]
    end
    
    norme = norm(nouveau - ancien)
    while norme >= tol
        ancien = copy(nouveau)
        nouveau = copy(b)
        for i = 1:n
            for j = 1:n
                if j != i 
                    nouveau[i] -= A[i, j]*ancien[j]
                end
            end
            nouveau[i] /= A[i, i]
        end
        
        norme = norm(nouveau - ancien)
    end

    return nouveau
end
