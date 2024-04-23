using LinearAlgebra


function find_factor(n)
    for i = 2:n
        if n % i == 0
            return i
        end
    end
    return n  
end

function isprime(n)
    if n <= 1
        return false
    end
    if n == 2
        return true
    end
    if n % 2 == 0
        return false
    end
    for i = 3:2:isqrt(n)
        if n % i == 0
            return false
        end
    end
    return true
end

function calculate_optimal_process(N)
    is_prime = isprime(N)
    if !is_prime
        factor = find_factor(N)
        npx = factor
        npy = N ÷ factor
    else
        npx = N
        npy = 1  
    end
    return npx,npy
end

N = 31
println("the result :",calculate_optimal_process(N))