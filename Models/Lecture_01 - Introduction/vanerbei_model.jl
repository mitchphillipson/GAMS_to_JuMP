using JuMP
using Clp
using GamsStructure

GU = GamsUniverse()

@GamsSets(GU,"julia_data",begin
    :d, "Dates", false
end);

@GamsParameters(GU,"julia_data",begin
    :day, (:d,), "Day (integer) associated with dates"
    :ave, (:d,), "Average temeratures"
end)


m = Model(Clp.Optimizer)

day = GU[:day]
ave = GU[:ave]

@variables(m,begin
    DEV[GU[:d]] >= 0 #Deviation
    T[GU[:d]]
    X0
    X1
    X2
    X3
    X4
    X5
end)

@objective(m,Min, sum(DEV[d] for d∈GU[:d]))

@constraint(m,tdef[d∈GU[:d]],
    T[d] == X0 + X1*day[[d]] + X2*cos(2*pi*day[[d]]/365.25)
                             + X3*sin(2*pi*day[[d]]/365.25)
                             + X4*cos(2*pi*day[[d]]/(10.7*365.25))
                             + X5*sin(2*pi*day[[d]]/(10.7*365.25))
)

@constraint(m,devlb[d∈GU[:d]],
    DEV[d] >= ave[[d]] - T[d]
)

@constraint(m,devub[d∈GU[:d]],
    DEV[d] >= T[d] - ave[[d]]
)

optimize!(m)


for var in all_variables(m)[end-5:end]
    println("$var => $(value(var))")
end

println("Warming = $(value(m[:X1])*365.25*100)")



