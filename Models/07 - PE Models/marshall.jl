using JuMP
using Ipopt
using Complementarity
using NamedArrays


# Data

T = [:a,:b,:c,:d]

tech = NamedArray([2 2;5 2; 7 4; 10 1e7],(T,[:cost,:cap]))

c = tech[:,:cost]
k = tech[:,:cap]


### NLP Formulation
nlp_model = Model(Ipopt.Optimizer)

@variables(nlp_model,begin
    P>=0
    PS>=0
    CS>=0
    S>=0 #I believe S is demand, or D
    0<=Q[t=T]<=k[t]
end)

@objective(nlp_model,Max, CS+PS)

@constraints(nlp_model,begin
    price, P == 10-S*10/6
    supply, S == sum(Q[t] for t∈T)
    psurplus, PS == sum((P-c[t])*Q[t] for t∈T)
    csurplus, CS == (10-P) * S/2
end)

optimize!(nlp_model)


### MCP Formulation
mcp_model = MCPModel()

@variables(mcp_model,begin
    Q[t=T]>=0
    RK[t=T]>=0
    D
    P
end)

@mapping(mcp_model, profit[t=T],
    RK[t]+c[t] - P
)

@mapping(mcp_model, capacity[t=T],
    k[t] - Q[t]
)

@mapping(mcp_model, demand, 
    D - (6*(1-P/10))
)

@mapping(mcp_model, market,
    sum(Q[t] for t∈T) - D
)

@complementarity(mcp_model,profit,Q)
@complementarity(mcp_model,capacity,RK)
@complementarity(mcp_model,demand,D)
@complementarity(mcp_model,market,P)

solveMCP(mcp_model)


@assert all(isapprox.(value.(nlp_model[:Q]),result_value.(mcp_model[:Q]);atol = 1e-7)) "Q variable solutions do not match"

@assert value(nlp_model[:P]) ≈ result_value(mcp_model[:P]) "P does not match"

@assert value(nlp_model[:S]) ≈ result_value(mcp_model[:D]) "Demand doesn't match"

@assert result_value.(mcp_model[:RK]) .≈ [3,0,0,0] "RK doesn't match"

println("GAMS gives RK as all zeros")