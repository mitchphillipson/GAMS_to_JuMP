using JuMP
using Ipopt

GU = Dict(:pi => 0.01, :L=> 0.5, :sigma => 2, :kother=>.5)
GU[:rho] = 1 - 1/GU[:sigma]

function nash(GU)
    
    Pi = GU[:pi]
    L = GU[:L]
    σ = GU[:sigma]
    kother = GU[:kother]
    ρ = GU[:rho]    

    m = Model(Ipopt.Optimizer)
    @variables(m,begin
        C_G>=1e-5, (start = 1,)
        C_B>=1e-5, (start = 1,)
        Gamma>=0, (start = Pi,)
        K, (start = 1,)
        EU >= (1-Pi)*1^ρ/ρ + Pi* *(1-L)^ρ/ρ
    end)

    @NLobjective(m,Max, (Gamma-Pi)*(K-kother))

    @NLconstraints(m,begin
        eudef, EU == (1-Pi)*C_G^ρ/ρ + Pi * C_B^ρ/ρ
        budget_G, C_G == 1-Gamma*K
        budget_B, C_B == 1-L + (1-Gamma)*K
        coverage, Gamma*((1-Pi)*C_G^(ρ-1) + Pi*C_B^(ρ-1)) >= Pi*C_B^(ρ-1)
    end)


    return m
end

## Solve loop

GU[:kother] = 0

log = Dict()
for n∈1:5
    dev = 1
    tmp = 0
    for iter∈1:25
        if isapprox(dev,0,atol=1e-4)
            break
        end
        Nash = nash(GU)
        set_silent(Nash)
        optimize!(Nash)

        log[n,iter] = Dict(
            :dev => dev,
            :K => value(Nash[:K]),
            :Gamma => value(Nash[:Gamma]),
            :Profit => objective_value(Nash),
            :C_G => value(Nash[:C_G]),
            :C_B => value(Nash[:C_B])
        )

        dev = abs(GU[:kother] - value(Nash[:K])*((n-1)/n))
        GU[:kother] = value(Nash[:K])*((n-1)/n)

        tmp = value(Nash[:K])*n/(n+1)
    end
    GU[:kother] = tmp
end


# Load GAMS solutions

using CSV

gams_log = Dict()
for row in CSV.File("iterlog.csv")
    gams_log[row[:n],row[:iter]] = Dict(
        :dev => row[:dev],
        :K => row[:K],
        :Gamma => row[:GAMMA],
        :Profit => row[:PROFIT],
        :C_B => row[:C_B],
        :C_G => row[:C_G]
    )
end


#Compare solutions. These don't match when n=3

for (n,iter)∈ sort(collect(intersect(keys(log),keys(gams_log))))
    for key in keys(log[n,iter])
        if !isapprox(log[n,iter][key],gams_log[n,iter][key],atol=1e-3)
            println("$n, $iter, $key, $(log[n,iter][key]-gams_log[n,iter][key])")
        end
    end
end