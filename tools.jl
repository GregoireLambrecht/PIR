function ValInd(T,v,σ)
    k1 = 1
    k2 = 1
    while (abs(T[k1] - v) > σ)
        k1 += 1
    end
    k2 = k1+1
    while (abs(T[k2] - v) > σ)
        k2 += 1
    end
    (k1,k2)
end


function Largeur(T,σ)
    M = maximum(T)
    (k1,k2) = ValInd(T,M/2,σ)
    Γ = k2 - k1
end
