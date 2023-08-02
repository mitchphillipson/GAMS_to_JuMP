using Complementarity
using NamedArrays
using CSV

using Suppressor

function test_solution(model,path_to_sols)
    vars = [(:K,1),
            (:M,1),
            (:X,2),
            (:PE,2),
            (:PR,1),
            (:PX,1),
            (:PEU,1),
            (:PC,2),
            (:EU,1),
            (:C,2)]

    #path = "gams_solutions/exogenous_maintenance/"

    for (var,n) in vars
        df = CSV.File("$path_to_sols$var.csv",stringtype=String)

        for row in df
            domain = Symbol.([row[i] for i∈1:n])
            value = row[n+1]
            @assert isapprox(result_value(model[var][domain...]),value,atol=1e-5) "Variable, $var, doesn't match in solution $path_to_sols"
        end
    end
    
end

##########
## Data ##
##########

R = [:ercot,:ferc] #Regions
S = [:summer,:winter,:vortex] #States of nature


σ = .25 #Degree of relative risk aversion
ESUB = .25 #Elasticity of substitution between electricity and other goods

pi = NamedArray([300,64,1],S) #State probabilities (read as days per year)
C0 = 31 #Reference aggregate expenditure (electricity plus other goods) /31/,

λ = NamedArray(zeros(length(R),length(S)),(R,S))
λ0 = NamedArray([1 1 1.5;1 1.2 1.2],(R,S))

w = NamedArray(zeros(length(R),length(S)),(R,S))
w0 = NamedArray([0 0 .25;0 0.2 0.25],(R,S))

pi = pi./365

m0 = NamedArray([.4,.8],R)

θe = 1/C0

α = 0
γ = 0


###########
## Model ##
###########

function ERCOT()

    model = MCPModel()

    @variables(model,begin
        K[R]>=0,   (start=1,)   #Capacity
        0<=M[R]<=.99,(start=.5,)                 #Maintenance

        PE[R,S]>=0,(start=1,)   #Wholesale price
        PR[R]>=0,  (start=1,)   #Electricity generation resource

        X[R,S], (start=1,)   #Sales to other region
        PX[S]>=0,  (start=1,)   #Traded electricity price

        PEU[R]>=0, (start=1,)   #Shadow price of expected utility
        PC[R,S]>=0,(start=1,)   #State-contingent price
        EU[R]>=0,  (start=1,)   #Expected utility
        C[R,S]>=0, (start=1,)   #Consumption;
    end)


    @NLexpressions(model,begin
        PC_def[r=R,s=S], PC[r,s] - (θe*(PE[r,s]*λ[r,s])^(1-ESUB) + 1 - θe)^(1/(1-ESUB))
        PEU_def[r=R], PEU[r] - (sum(pi[s]*PC[r,s]^(1-σ) for s∈S)^(1/(1-σ)))
        
        EU_def[r=R], EU[r] - 1/PEU[r]
        
        C_def[r=R,s=S], C[r,s] - EU[r]*(PEU[r]/PC[r,s])^σ
        
        profit_K[r=R], sqrt(PR[r]) - sum( pi[s]*PE[r,s]*(1-w[r,s]*(1-M[r])) for s∈S)
        profit_M[r=R], α/(1-M[r])^γ - sum( pi[s]*PE[r,s]*w[r,s] for s∈S)
        profit_X[r=R,s=S], PE[r,s] - PX[s]
        
        market_PE[r=R,s=S], K[r]*(1-w[r,s]*(1-M[r])) - ( C[r,s] * λ[r,s] * (PC[r,s]/(PE[r,s]*λ[r,s]))^ESUB + X[r,s])
        market_PR[r=R], .5 - .5*K[r]*1/sqrt(PR[r])
        market_PX[s=S], sum(X[r,s] for r∈R) - 0
    end)


    @complementarity(model,profit_K,K)
    @complementarity(model,profit_M,M)
    @complementarity(model,profit_X,X)
    @complementarity(model,market_PE,PE)
    @complementarity(model,market_PR,PR)
    @complementarity(model,market_PX,PX)
    @complementarity(model,PEU_def,PEU)
    @complementarity(model,PC_def,PC)
    @complementarity(model,EU_def,EU)
    @complementarity(model,C_def,C)

    return model


end


##############
## Model 01 ##
##############

λ.=1
w.=0


model_01 = ERCOT()

fix.(model_01[:X],0,force=true)
fix.(model_01[:PX],1,force=true)
fix.(model_01[:M],0,force=true)


#println(start_value.(ignore_weather[:PC]))

out = @capture_out solveMCP(model_01,cumulative_iteration_limit=0)

test_solution(model_01,"gams_solutions/01/")



##############
## Model 02 ##
##############

w .= w0
λ .= λ0

model_02 = ERCOT()

fix.(model_02[:X],0,force=true)
fix.(model_02[:PX],1,force=true)
fix.(model_02[:M],m0[R].array,force=true)



out = @capture_out solveMCP(model_02)


test_solution(model_02,"gams_solutions/02/")


#	Calibrate parameters of the maitenance cost function which are consistent
#	with the assumed maintenance levels in ercot and ferc:


γ = ( log(sum( pi[s] * result_value(model_02[:PE][:ercot,s])*w[:ercot,s] for s∈S)) - 
      log(sum( pi[s] * result_value(model_02[:PE][:ferc,s]) *w[:ferc,s]  for s∈S))
    ) / (log( 1-result_value(model_02[:M][:ferc])) - log( 1-result_value(model_02[:M][:ercot])))


α = sum( pi[s]*result_value(model_02[:PE][:ercot,s])*w[:ercot,s] for s∈S) * (1-result_value(model_02[:M][:ercot]))^γ



##############
## Model 03 ##
##############
w .= w0
λ .= λ0

model_03 = ERCOT()

fix.(model_03[:X],0,force=true)
fix.(model_03[:PX],1,force=true)

set_start_value.(all_variables(model_03),result_value.(all_variables(model_02)))

out = @capture_out solveMCP(model_03,cumulative_iteration_limit=0)

test_solution(model_03,"gams_solutions/03/")


k0 = result_value.(model_03[:K])



##############
## Model 04 ##
##############

model_04 = ERCOT()

out = @capture_out solveMCP(model_04)

test_solution(model_04,"gams_solutions/04/")



##############
## Model 05 ##
##############
model_05 = ERCOT()

fix.(model_05[:X],0,force=true)
fix.(model_05[:PX],1,force=true)
fix.(model_05[:K],k0,force=true)

out = @capture_out solveMCP(model_05)


test_solution(model_05,"gams_solutions/05/")


##############
## Model 06 ##
##############
model_06 = ERCOT()

fix.(model_06[:K],k0,force=true)

out = @capture_out solveMCP(model_06)


#print(generate_report(model_06))

test_solution(model_06,"gams_solutions/06/")


##############
## Model 07 ##
##############
pi[S] = [300,55,10]./365

model_07 = ERCOT()

fix.(model_07[:X],0,force=true)
fix.(model_07[:PX],1,force=true)

out = @capture_out solveMCP(model_07)

test_solution(model_07,"gams_solutions/07/")



##############
## Model 08 ##
##############
model_08 = ERCOT()

out = @capture_out solveMCP(model_08)

test_solution(model_08,"gams_solutions/08/")