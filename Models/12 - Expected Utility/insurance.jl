using JuMP
using Ipopt
using Complementarity

GU = Dict(:pi => 0.01, :L => 0.5, :gamma => 0.02, :sigma => 0.5)

GU[:rho] = 1 - 1/GU[:sigma]

function insurance(GU)

    pi = GU[:pi]
    L = GU[:L]
    γ = GU[:gamma]
    σ = GU[:sigma]
    ρ = GU[:rho]

    m = Model(Ipopt.Optimizer)

    @variables(m,begin
        C_G, (start=1,)
        C_B, (start=1,)
        K, (start=1,)
    end)

    @NLobjective(m,Max,(1-pi)* C_G^ρ/ρ+ pi*C_B^ρ/ρ)

    @NLconstraints(m,begin
        budget_G, C_G == 1-γ*K
        budget_B, C_B == 1-L + (1-γ)*K
    end)


    return m
end



function insurance_ev(GU)
    pi = GU[:pi]
    L = GU[:L]
    γ = GU[:gamma]
    σ = GU[:sigma]
    ρ = GU[:rho]

    m = Model(Ipopt.Optimizer)

    @variables(m,begin
        C_G, (start=1,)
        C_B, (start=1,)
        K, (start=1,)
    end)

    @NLobjective(m,Max,100* (( (1-pi)*C_G^ρ + pi * C_B^ρ)^(1/ρ) - 1))

    @NLconstraints(m,begin
        budget_G, C_G == 1-γ*K
        budget_B, C_B == 1-L + (1-γ)*K
    end)


    return m

end

function equilibrium(GU)
    pi = GU[:pi]
    L = GU[:L]
    γ = GU[:gamma]
    σ = GU[:sigma]
    ρ = GU[:rho]

    m = MCPModel()

    @variables(m,begin
        EU, (start=1,)
        EV, (start=1,)
        C_G, (start=1,)
        C_B, (start=1,)
        K, (start=1,)
    end)

    @mapping(m,eudef,EU - ((1-pi)* C_G^(ρ)/ρ+ pi*C_B^(ρ)/ρ))

    @mapping(m,evdef, EV - 100* (( (1-pi)*C_G^ρ + pi * C_B^ρ)^(1/ρ) - 1))

    @mapping(m,budget_G, C_G - (1-γ*K))

    @mapping(m,budget_B, C_B - (1-L + (1-γ)*K))

    @mapping(m,coverage, γ*((1-pi)*C_G^(ρ-1) + pi*C_B^(ρ-1)) - pi*C_B^(ρ-1))
    
    @complementarity(m,eudef,EU)
    @complementarity(m,evdef,EV)
    @complementarity(m,budget_G,C_G)
    @complementarity(m,budget_B,C_B)
    @complementarity(m,coverage,K)

    return m
end

#I = insurance(GU)
#set_silent(I)
#optimize!(I)

I_ev = insurance_ev(GU)
set_silent(I_ev)
optimize!(I_ev)


c = equilibrium(GU)
set_silent(c)

solveMCP(c);

for var in [:C_G,:C_B,:K]
    @assert isapprox(value(I_ev[var]),result_value(c[var]),atol=1e-6) "Variable, $var, doesn't match in insurance.jl"
end

@assert round(result_value(c[:EU]),digits=4) == -1.0083 "Variable, EU, doesn't match in insurance.jl"
@assert round(result_value(c[:EV]),digits=4) == -0.8274 "Variable, EU, doesn't match in insurance.jl"