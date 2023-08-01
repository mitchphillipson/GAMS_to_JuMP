using JuMP
using Ipopt
using Complementarity
using NamedArrays

S = Symbol.("s",1:9) #Load Segements
J = [:Nuclear,:Coal,:Gas,:Diesel] #Generating units
I = [:residential,:commercial,:industrial] #Demand categories

load_domain = [:qref,:hours,:residential,:commercial,:industrial] #Domain of Load

supply_domain = [:mc,:cap] #Domain of supply

load = NamedArray( [1.00 310 0.7 0.2 0.1;
                    0.79 750 0.5 0.2 0.3;
                    0.59 1020 0.4 0.3 0.3;
                    0.50 2460 0.4 0.3 0.3;
                    0.44 1910 0.3 0.3 0.4;
                    0.40 580 0.2 0.3 0.5;
                    0.35 580 0.1 0.3 0.6;
                    0.28 600 0.1 0.3 0.6;
                    0.23 540 0.1 0.3 0.6], 
                (S,load_domain)
)

supply = NamedArray([   4 0.30;
                        6 0.30;
                        7 0.20;
                        9 0.30],
                    (J,supply_domain)
)


epsilon = NamedArray([.1,.2,.5],I)

h = load[S,:hours]
mc = supply[J,:mc]
cap = supply[J,:cap]
qref = load[S,:qref]


reference_prices = Model(Ipopt.Optimizer)

@variable(reference_prices,Y[J,S]>=0)

for s∈S
    set_upper_bound.(Y[J,s],cap.array)
end

@objective(reference_prices,Min, sum(mc[j]*Y[j,s] for s∈S,j∈J))

@constraint(reference_prices,demand[s=S],
    sum(Y[j,s] for j∈J) == qref[s]
)

set_silent(reference_prices)

optimize!(reference_prices)


pref = dual.(reference_prices[:demand])

dref = NamedArray(zeros(length(I),length(S)),(I,S))
for i∈I,s∈S
    dref[i,s] = qref[s] * load[s,i]
end


function calibrated_model_mcp()

    calibrated_mcp = MCPModel()

    @variables(calibrated_mcp,begin
        Y[J,S]>=0
        D[i=I,s=S]>=0, (start = dref[i,s],)
        PI[j=J,s=S]>=0, (start = -dual.(UpperBoundRef.(reference_prices[:Y][j,s])),)
        P[s=S]>=0, (start = pref[s],)
    end)

    @mapping(calibrated_mcp,aggdemand[i=I,s=S],
        D[i,s] - dref[i,s] * (1 - epsilon[i] * (P[s]/pref[s]-1))
    )

    @mapping(calibrated_mcp,supplydemand[s=S],
        sum(Y[j,s] for j∈J) - sum(D[i,s] for i∈I)
    )

    @mapping(calibrated_mcp,profit[j=J,s=S],
        mc[j] + PI[j,s] - P[s] 
    )

    @mapping(calibrated_mcp,capacity[j=J,s=S],
        cap[j] - Y[j,s]
    )   


    @complementarity(calibrated_mcp,aggdemand,D)
    @complementarity(calibrated_mcp,supplydemand,P)
    @complementarity(calibrated_mcp,profit,Y)
    @complementarity(calibrated_mcp,capacity,PI)

    return calibrated_mcp

end
   # solveMCP(calibrated_mcp,cumulative_iteration_limit=1)



function equilibrium_allocation_nlp()
    m = Model(Ipopt.Optimizer)

    @variables(m,begin
        Y[J,S]>=0
        D[i=I,s=S]>=0, (start = dref[i,s],)
    #    PI[j=J,s=S]>=0, (start = -dual.(UpperBoundRef.(reference_prices[:Y][j,s])),)
    #    P[s=S]>=0, (start = pref[s],)
        K[J]>=0
    end)


    @objective(m,Max, 
        sum(pref[s]*D[i,s] * (1 + 1/epsilon[i] * (1-D[i,s]/(2*dref[i,s]))) for i∈I,s∈S)
        - sum(mc[j]*Y[j,s] for j∈J,s∈S)
    )


    @constraint(m,supplydemand[s=S],
        sum(Y[j,s] for j∈J) == sum(D[i,s] for i∈I)
    )

    #@constraint(m,profit[j=J,s=S],
    #    mc[j] + PI[j,s] >= P[s] 
    #)

    @constraint(m,capacity[j=J,s=S],
        cap[j] >= Y[j,s]
    )

    fix.(K[J],supply[J,:cap].array,force=true)

    return m
end


mc[:Coal] = 6

bench_nlp = equilibrium_allocation_nlp()
optimize!(bench_nlp)

bench_mcp = calibrated_model_mcp()
solveMCP(bench_mcp)


mc[:Coal] = 12

shock_nlp = equilibrium_allocation_nlp()
optimize!(shock_nlp)

shock_mcp = calibrated_model_mcp()
solveMCP(shock_mcp)


bench_PI = NamedArray(  [5 3 2 2 2 2 2 0 0; 
                        3 1 0 0 0 0 0 0 0;
                        2 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0],
                    (J,S)
)

shock_PI = NamedArray(  [8 5 5 3 3 3 3 0 0; 
                         0 0 0 0 0 0 0 0 0;
                         5 2 2 0 0 0 0 0 0;
                         3 0 0 0 0 0 0 0 0],
                    (J,S)
)

bench_P = NamedArray([9, 7, 6, 6, 6, 6, 6, 4, 4],S)
shock_P = NamedArray([12, 9, 9, 7, 7, 7, 7, 4, 4],S)


matched_variables = [:Y,:D]
for var in matched_variables
    @assert all(isapprox.(value.(bench_nlp[var]),result_value.(bench_mcp[var]),atol=1e-6)) "Variable, $var, mismatch in qpdispatch"
    @assert all(isapprox.(value.(shock_nlp[var]),result_value.(shock_mcp[var]),atol=1e-6)) "Variable, $var, mismatch in qpdispatch"
end

@assert all(isapprox.(result_value.(bench_mcp[:PI]), bench_PI.array,atol=1e-6)) "Variable, PI, mismatch in qpdispatch"
@assert all(isapprox.(result_value.(shock_mcp[:PI]), shock_PI.array,atol=1e-6)) "Variable, PI, mismatch in qpdispatch"


@assert all(isapprox.(result_value.(bench_mcp[:P]), bench_P.array,atol=1e-6)) "Variable, P, mismatch in qpdispatch"
@assert all(isapprox.(result_value.(shock_mcp[:P]), shock_P.array,atol=1e-6)) "Variable, P, mismatch in qpdispatch"