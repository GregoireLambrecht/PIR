function Largeur(T,σ)
    M = maximum(T)
    n = argmax(T)
    k1 = 1
    while (abs(T[k1] - M/2) > σ)
        k1 += 1
    end
    k2 = n+1
    while (abs(T[k2] - M/2) > σ)
        k2 += 1
    end
    k2 - k1
end
