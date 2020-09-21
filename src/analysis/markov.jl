
"""
    prob_estacionária(g::SparseMatrixCSC)

Probabilidade estacionária associada à matriz estocástica esparsa `g`,
calculada como um autovetor associado ao autovalor `1` de `g'`, cujas
entradas somam `1`. Caso a matriz não seja irredutível, a probabilidade
pode não ser única. 
"""
prob_estacionária(g::SparseMatrixCSC) = normalize(eigs(g';nev=1,which=:LR)[2][:,1],1) .|> abs


"""
    mat_estocástica(maps::AbstractMatrix)

Retorna uma cópia da matriz `A` com as linhas normalizadas de forma
que as somas de cada linha vale `1`, ou `0` se a linha inteira for nula.
"""
function mat_estocástica(A::AbstractMatrix) 
    P = copy(A)
    for m in eachrow(P)
        if sum(m) != 0
            m ./= sum(m)
        else
            m .*= 0
        end
    end
    P
end


"""
       map_markov(bd::Billiard{T}, Nbounces::Int,
                  np::Int, Nξ, Nφ=Nξ;
                  intervals=arcintervals(bd)) where {T}

Estima um processo de Markov gerado pela partição do bilhar `bd` em
`Nξ` comprimentos de arco iguais e `Nφ` intervalos de ângulos iguais,
resultando em `Nξ*Nφ` estados. As probabilidades de transição são estimadas
simulando `np` partículas cada uma por `Nbounces` reflexões no bordo.
"""
function map_markov(P::Type{<:AbstractParticle}, bd::Billiard{T}, Nbounces::Int,
                           np::Int, Nξ, Nφ=Nξ;
                           intervals=arcintervals(bd)) where {T}
    δξ = totallength(bd) / Nξ # Tamanho dos compartimentos-arco
    δφ = 2 / Nφ               # Tamanho dos compartimentos-ângulo

    # As partículas aleatórias
    ps = randominside(P,bd,np)

    # Vetor para armazenar os dados de cada partícula em paralelo.
    pP = [spzeros(Int, Nξ * Nφ, Nξ * Nφ) for i in 1:np]

    # Transforma (arco,ângulo) em índice.
    vectoind(ξ,sφ) = floor(Int, (sφ + 1) / δφ) * Nξ +
                     floor(Int, ξ / δξ) + 1

    prog = Progress(np) # Barra de progresso

    # Itera pelas partículas em paralelo
    Threads.@threads for ip in 1:np
        # Partícula atual
        p = ps[ip]
        # Índice da reflexão anterior
        lastind = nothing
        for N in 1:Nbounces
        i, = bounce!(p,bd) # Reflete a partícula na próxima borda

        ξ, sφ = to_bcoords(p, bd[i]) # Coordenadas
        ξ += intervals[i] # Transforma em coordenadas globais

        ind = vectoind(ξ,sφ) # Índice
        
        if N > 1
            # Registra que a partícula `ip` saiu da região `lastind`
            # e foi para a região `ind` mais uma vez.
            pP[ip][lastind,ind] += 1
        end
        lastind = ind
        end
        next!(prog) # atualiza a barra de progresso
  end

  mat_estocástica(sum(pP)./1) # transforma em matriz estocástica
end

function média_temporal(bd::Billiard{T}, Nbounces::Int,
                       np::Int, Nξ, Nφ=Nξ;
                       intervals=arcintervals(bd)) where {T}
    δξ = totallength(bd) / Nξ # Tamanho dos compartimentos-arco
    δφ = 2 / Nφ               # Tamanho dos compartimentos-ângulo

    # As partículas aleatórias
    ps = [randominside(bd) for _ in 1:np] 

    # Vetor para armazenar os dados de cada partícula em paralelo.
    pP = [spzeros(Int, Nξ * Nφ) for _ in 1:np]
    firsti = zeros(Int, np)

    # Transforma (arco,ângulo) em índice.
    vectoind(ξ,sφ) = floor(Int, (sφ + 1) / δφ) * Nξ +
                     floor(Int, ξ / δξ) + 1

    prog = Progress(np) # Barra de progresso
    
    # Itera pelas partículas em paralelo
    Threads.@threads for ip in 1:np
        # Partícula atual
        p = ps[ip]
        for N in 1:Nbounces
            i, = bounce!(p,bd) # Reflete a partícula na próxima borda

            ξ, sφ = to_bcoords(p, bd[i]) # Coordenadas
            ξ += intervals[i] # Transforma em coordenadas globais

            ind = vectoind(ξ,sφ) # Índice
            
            pP[ip][ind] += 1

            if N == 1
                firsti[ip] = ind
            end
        end
        next!(prog) # atualiza a barra de progresso
    end
    
    tret = zeros(Nξ * Nφ)
    mret = zeros(Int, Nξ * Nφ)
    for (p,fi) in zip(pP,firsti)
        tret[fi] += Nbounces/p[fi]
        mret[fi] += 1
    end

    vcat(pP'...)./Nbounces, [mret[i] ≠ 0 ? tret[i]/mret[i] : 0 for i in 1:(Nξ * Nφ)]
end
