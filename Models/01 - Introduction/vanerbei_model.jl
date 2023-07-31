using JuMP
using NamedArrays
using Ipopt


function run_vanerbei_model(data_path = "julia_data")
    data = Dict()
    data[:d] = [Symbol(r[1]) for r∈CSV.File("$data_path/d.csv")]

    day = [r[2] for r∈ CSV.File("julia_data/day.csv")]
    data[:day] = NamedArray(day,data[:d])

    val = [r[2] for r∈ CSV.File("julia_data/ave.csv")]
    data[:ave] = NamedArray(val,data[:d])

    model = vanerbei_model(data)

    set_silent(model)

    optimize!(model)

    return model

end


function vanerbei_model(data)


    D = data[:d]
    day = data[:day]
    ave = data[:ave]

    m2 = Model(Ipopt.Optimizer)


    @variables(m2,begin
        DEV[D] >= 0 #Deviation
        T[D]
        X0
        X1
        X2
        X3
        X4
        X5
    end)

    @objective(m2,Min, sum(DEV[d] for d∈D))

    @constraint(m2,tdef[d∈D],
        T[d] == X0 + X1*day[d] + X2*cos(2*pi*day[d]/365.25)
                                + X3*sin(2*pi*day[d]/365.25)
                                + X4*cos(2*pi*day[d]/(10.7*365.25))
                                + X5*sin(2*pi*day[d]/(10.7*365.25))
    )

    @constraint(m2,devlb[d∈D],
        DEV[d] >= ave[d] - T[d]
    )

    @constraint(m2,devub[d∈D],
        DEV[d] >= T[d] - ave[d]
    );

    return m2

end