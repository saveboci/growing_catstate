println("BEGIN")
flush(stdout)

using JLD2
using DataFrames
using ITensors

theta = parse(Float64, ARGS[1])

#=
##### TRANSISTOR
Jx = 0.7
Jy = 0.0
Jz = 0.0
h = 0.0
dt = 0.01;

Ls = [160]; bDs = [240, 300]; elle = length(Ls); subsys = Ls[elle]; which_model = "quantum_transistor/standard/L$(Ls[end])_Jx$(Jx)_Jy$(Jy)_Jz$(Jz)_h$(h)"
=#



#=
##### SYM EAST MODEL
J = 0.25
g = 0.22
dt = 0.01;

Ls = [80]; bDs = [20, 40, 80, 160, 240]; elle = length(Ls); subsys = Ls[elle]; which_model= "Symmetric_East_Models/dual_ising/L$(Ls[end])_J$(J)"; 
#Ls = [80]; bDs = [80, 120, 160]; elle = length(Ls); subsys = Ls[elle]; which_model = "Symmetric_East_Models/dual_ising_withX/L$(Ls[end])_J$(J)_g$(g)"; 
#Ls = [140]; bDs = [160, 240]; elle = length(Ls); subsys = Ls[elle]; which_model = "Symmetric_East_Models/dual_ising_withY/L$(Ls[end])_J$(J)_g$(g)"; 
=#


#=
##### WEIRD MODEL

α = 0.7
Δ = 0.5
dt = 0.01;

Ls = [140]; bDs = [100, 160, 240]; elle = length(Ls); subsys = Ls[elle]; which_model = "XtimesXXpYYtimesX/withXX/L$(Ls[end])_resc1_alpha$(α)_Delta$(Δ)"; 
#Ls = [140]; bDs = [100, 160, 240]; elle = length(Ls); subsys = Ls[elle]; which_model = "XtimesXXpYYtimesX/withX/L$(Ls[end])_resc1_alpha$(α)_Delta$(Δ)"; 
=#



#=
##### WEIRD MODEL WITH TRANSISTOR
Jx = 0.3;
Jy = -0.09
#Jz = 0.09
h  = 0.09
Δ  = 0.0
Jw =  -0.15
dt = 0.01

Ls = [100]; bDs = [350, 500]; elle = length(Ls); subsys = Ls[elle]; which_model = "XtimesXXpYYtimesX/withTRANS/L$(Ls[end])_Jx$(Jx)_Jy$(Jy)_h$(h)_Delta$(Δ)_Jw$(Jw)"; 
#Ls = [80]; bDs = [200, 300]; elle = length(Ls); subsys = Ls[elle]; which_model = "XtimesXXpYYtimesX/withTRANS/L$(Ls[end])_Jx$(Jx)_Jy$(Jy)_h$(h)_Delta$(Δ)_Jw$(Jw)"; 
=#



#=
##### TRANSISTOR FULL GLORY
J = 1
γ = 0.5
w = 0.7
Δ = 0.0
Dz = 0.6
hz = 0.0
dt = 0.01

Ls = [100]; bDs = [200, 300]; elle = length(Ls); subsys = Ls[elle]; which_model = "quantum_transistor/fullglory/L$(Ls[end])_J$(J)_gamma$(γ)_w$(w)_Delta$(Δ)_Dz$(Dz)_hz$(hz)"; 
=#


#=
##### LONG RANGE XXZ
n = 2
Jx = 0.25
Jz = 0.1
g = 0.38
dt = 0.01

Ls = [120]; bDs = [240, 300]; elle = length(Ls); subsys = Ls[elle]; which_model = "longrange/XXZ/n$(n)_L$(Ls[end])_Jx$(Jx)_Jz$(Jz)_g$(g)"; 
=#






# for laptop:
#wheretosave = "/Users/saveriobocini/numerics_PhD/full_counted_statistics_CAT/$(which_model)";
#wheretoload = "/Users/saveriobocini/numerics_PhD/$(which_model)";

# for desktop
wheretosave = "/home/sbocini/Numerics_LabDesktop/full_counted_statistics_CAT/$(which_model)";
wheretoload = "/home/sbocini/Numerics_LabDesktop/$(which_model)";


println(which_model)

times = [
    load("$(wheretoload)_BondDim$(bDs[1])_dt$(dt).jld2"
    , "times") for L in Ls
];

states = [
    [
        [
            load("$(wheretoload)_BondDim$(bDs[j])_dt$(dt).jld2"
                , "t$(round(times[k][i], digits=3))")
            for i=1:length(times[1])
        ] for j=1:length(bDs)
    ] for k=1:length(Ls)
];




function eiM(θ, psi, ell) 
    siti = siteinds(psi)
    lunghezza = length(psi)
    outpu = ITensor[]
    for j in siti[- div(ell, 2) + 1 + div(lunghezza, 2) : div(ell, 2) + div(lunghezza, 2)]
        ope = op("Sz", j)
        push!(outpu, exp(2im * pi * θ * ope / (ell + 1)))
    end
    outpu
end

function G(state, ell) 
    [inner(state, apply(eiM(k, state, ell), state)) for k=-div(ell, 2):div(ell, 2)]
end

function full_counting_ofZ(state, ell)
    Gk = G(state, ell)
    outpu = ComplexF64[]
    rangeh = div(ell, 2)
    for m = -rangeh:rangeh
        expophase = [exp(- 2im * pi * k * m / (ell + 1)) for k=-rangeh:rangeh]
        push!(outpu, sum(expophase .* Gk) / (ell + 1))
    end
    outpu
end






println("Full counting of the MAGNETIZATION in the ROTATED CASE started.")
flush(stdout)

full_countedZs = [
    [
        real(full_counting_ofZ(sin(theta/2) * states[elle][bd][i] + cos(theta/2) * productMPS(siteinds(states[elle][bd][i]), n -> "Up"), subsys)) 
        for i=1:length(times[elle])
    ] for bd = 1:length(bDs)
];
        

println("Full counting ended, just need to save it now.")



for bonddim = 1:length(bDs)
    df = DataFrame([full_countedZs[bonddim][i][j] for j=1:length(full_countedZs[bonddim][1]), i=1:length(full_countedZs[bonddim])], :auto)
    rename!(df, Symbol.(times[elle]))
    save("$(wheretosave)_dt$(dt)_bD$(bDs[bonddim]).csv", df)
end

println("THE END")
