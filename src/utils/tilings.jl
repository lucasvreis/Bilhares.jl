function kr(nlados,div)
    α = big(π)/div
    γ = 2big(π)/nlados
    hr = acosh((cos(α)^2 + cos(γ)) / sin(α)^2)
    
    1 - 1/(exp(hr))
end


