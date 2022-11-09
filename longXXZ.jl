using ITensors
using JLD2

println("BEGIN")

nmax = parse(Int, ARGS[1])

#gn = [1/n^2 for n=1:nmax]
gn = [parse(Float64, ARGS[6]), parse(Float64, ARGS[7])]

Jx = 0.25
Jz = round(0.25 * parse(Float64, ARGS[4]), digits=2)
g = round(0.25 * parse(Float64, ARGS[5]), digits=2)

Jxs = Jx .* gn
Jzs = Jz .* gn

bondD = parse(Int, ARGS[3])
dt = 0.01
L = parse(Int, ARGS[2]) 
sites = siteinds("S=1/2", L);

println("Long range XXZ with L=$(L), bD=$(bondD), n=$(nmax), Jx=$(Jx), Jz=$(Jz), g=$(g), gn = $(gn), dt=$(dt)")

# here I define time evolution-----------------------------
# it works only for nmax = 1

function TEBD1(psi, cutoff, bond_dim_max, timestep, ttotal)
        
    tau=timestep
    
    siti = siteinds(psi)
    LL = length(psi)
        
    Nsteps = round(Int, ttotal/timestep) 
    gates = ITensor[] 
    
    for j=2:LL-nmax-1
        ss = siti[j:j+nmax]
        X1 = 4Jxs[1] * op("Sx",ss[1]) * op("Sx",ss[2])
        Y1 = 4Jxs[1] * op("Sy",ss[1]) * op("Sy",ss[2])
        Z1 = 4Jzs[1] * op("Sz",ss[1]) * op("Sz",ss[2])
        Ej = 4g * op("Sx",ss[1]) * op("Sz",ss[2])
        Fj = -4g * op("Sz",ss[1]) * op("Sx",ss[2])
        hj = X1 + Y1 + Z1 + Ej + Fj
        Gj = exp(-1.0im * tau/2 * hj)
        push!(gates,Gj)
    end
    
    #LEFT EDGE:
    ss = siti[1:1+nmax]
    X1 = 4Jxs[1] * op("Sx",ss[1]) * op("Sx",ss[2])
    Y1 = 4Jxs[1] * op("Sy",ss[1]) * op("Sy",ss[2])
    Z1 = 4Jzs[1] * op("Sz",ss[1]) * op("Sz",ss[2])
    Fj = -4g * op("Sz",ss[1]) * op("Sx",ss[2])
    hj = X1 + Y1 + Z1 + Fj
    Gj = exp(-1.0im * tau/2 * hj)
    push!(gates,Gj)
    
    #RIGHT EDGE:
    ss = siti[L-nmax:L]
    X1 = 4Jxs[1] * op("Sx",ss[1]) * op("Sx",ss[2])
    Y1 = 4Jxs[1] * op("Sy",ss[1]) * op("Sy",ss[2])
    Z1 = 4Jzs[1] * op("Sz",ss[1]) * op("Sz",ss[2])
    Ej = 4g * op("Sx",ss[1]) * op("Sz",ss[2])
    hj = X1 + Y1 + Z1 + Ej
    Gj = exp(-1.0im * tau/2 * hj)
    push!(gates,Gj)
    
    append!(gates, reverse(gates));
    
    for step in 1:Nsteps
        psi = apply(gates, psi; cutoff=cutoff, maxdim = bond_dim_max) 
    end
    
    psi
end

# here I define time evolution
# it works only for nmax = 2

function TEBD2(psi, cutoff, bond_dim_max, timestep, ttotal)
        
    tau=timestep
    
    siti = siteinds(psi)
    LL = length(psi)
        
    Nsteps = round(Int, ttotal/timestep) 
    gates = ITensor[] 
    
    for j=2:LL-nmax-1
        ss = siti[j:j+nmax]
        X1 = 4Jxs[1] * op("Sx",ss[1]) * op("Sx",ss[2]) * op("Id",ss[3])
        Y1 = 4Jxs[1] * op("Sy",ss[1]) * op("Sy",ss[2]) * op("Id",ss[3])
        Z1 = 4Jzs[1] * op("Sz",ss[1]) * op("Sz",ss[2]) * op("Id",ss[3])
        X2 = 4Jxs[2] * op("Sx",ss[1]) * op("Id",ss[2]) * op("Sx",ss[3])
        Y2 = 4Jxs[2] * op("Sy",ss[1]) * op("Id",ss[2]) * op("Sy",ss[3])
        Z2 = 4Jzs[2] * op("Sz",ss[1]) * op("Id",ss[2]) * op("Sz",ss[3])
        Ej = 4g * op("Sx",ss[1]) * op("Sz",ss[2]) * op("Id",ss[3])
        Fj = -4g * op("Sz",ss[1]) * op("Sx",ss[2]) * op("Id",ss[3])
        hj = X1 + X2 + Y1 + Y2 + Z1 + Z2 + Ej + Fj
        Gj = exp(-1.0im * tau/2 * hj)
        push!(gates,Gj)
    end
    
    #LEFT EDGE:
    ss = siti[1:1+nmax]
    X1 = 4Jxs[1] * op("Sx",ss[1]) * op("Sx",ss[2]) * op("Id",ss[3])
    Y1 = 4Jxs[1] * op("Sy",ss[1]) * op("Sy",ss[2]) * op("Id",ss[3])
    Z1 = 4Jzs[1] * op("Sz",ss[1]) * op("Sz",ss[2]) * op("Id",ss[3])
    X2 = 4Jxs[2] * op("Sx",ss[1]) * op("Id",ss[2]) * op("Sx",ss[3])
    Y2 = 4Jxs[2] * op("Sy",ss[1]) * op("Id",ss[2]) * op("Sy",ss[3])
    Z2 = 4Jzs[2] * op("Sz",ss[1]) * op("Id",ss[2]) * op("Sz",ss[3])
    Fj = -4g * op("Sz",ss[1]) * op("Sx",ss[2]) * op("Id",ss[3])
    hj = X1 + X2 + Y1 + Y2 + Z1 + Z2 + Fj
    Gj = exp(-1.0im * tau/2 * hj)
    push!(gates,Gj)
    
    #RIGHT EDGE:
    ss = siti[L-nmax:L]
    X1 = 4Jxs[1] * op("Sx",ss[1]) * op("Sx",ss[2]) * op("Id",ss[3])
    Y1 = 4Jxs[1] * op("Sy",ss[1]) * op("Sy",ss[2]) * op("Id",ss[3])
    Z1 = 4Jzs[1] * op("Sz",ss[1]) * op("Sz",ss[2]) * op("Id",ss[3])
    X2 = 4Jxs[2] * op("Sx",ss[1]) * op("Id",ss[2]) * op("Sx",ss[3])
    Y2 = 4Jxs[2] * op("Sy",ss[1]) * op("Id",ss[2]) * op("Sy",ss[3])
    Z2 = 4Jzs[2] * op("Sz",ss[1]) * op("Id",ss[2]) * op("Sz",ss[3])
    Ej = 4g * op("Sx",ss[1]) * op("Sz",ss[2]) * op("Id",ss[3])
    hj = X1 + X2 + Y1 + Y2 + Z1 + Z2 + Ej
    Gj = exp(-1.0im * tau/2 * hj)
    push!(gates,Gj)
    
    append!(gates, reverse(gates));
    
    for step in 1:Nsteps
        psi = apply(gates, psi; cutoff=cutoff, maxdim = bond_dim_max) 
    end
    
    psi
end

function TEBD(n, psi, cutoff, bond_dim_max, timestep, ttotal)
    if n==1
        TEBD1(psi, cutoff, bond_dim_max, timestep, ttotal)
    elseif n==2
        TEBD2(psi, cutoff, bond_dim_max, timestep, ttotal)
    else
        println("Houston, we have a problem!")
    end
end

#------------- time evolution defined.


state_all_up = productMPS(sites, n -> "Up");
vmax = 3;

timeSTART = 0.2;
timeSTEPS = 11;

timeMAX = (L/2 - 6) / vmax;
timeSTEP = (timeMAX - timeSTART) / timeSTEPS;

print("initializing step zero\n")
flush(stdout)

psi0 = apply(2 * op("Sx", sites[round(Int,L/2)]), state_all_up);
evolved = TEBD(nmax, psi0, 1E-12, bondD, dt, timeSTART);
provisory = evolved
 
print(0,"/",timeSTEPS," time steps done\n")
flush(stdout)

times = [timeSTART + (j-1) * timeSTEP for j in 1:timeSTEPS+1];
lege = reshape([round(timeSTART + (i-1)*timeSTEP, digits=3) for i in 1:timeSTEPS+1],(1,timeSTEPS+1));

for i in 1:timeSTEPS
    global provisory = TEBD(nmax, provisory, 1E-12, bondD, dt, timeSTEP);
    global evolved = vcat(evolved, provisory)

    save("/home/sbocini/Numerics_LabDesktop/longrange/XXZ/n$(nmax)_L$(L)_Jx$(Jx)_Jz$(Jz)_g$(g)_g1$(gn[1])_g2$(gn[2])_BondDim$(bondD)_dt$(dt).jld2"
    #save("/Users/saveriobocini/numerics_PhD/longrange/XXZ/n$(nmax)_L$(L)_Jx$(Jx)_Jz$(Jz)_g$(g)_g1$(gn[1])_g2$(gn[2])_BondDim$(bondD)_dt$(dt).jld2"
        , "times", times
        , "t$(lege[1])", evolved[1]
        , "t$(lege[2])", evolved[2]
        , "t$(lege[3])", (i > 1 ? evolved[3] : "not here yet")
        , "t$(lege[4])", (i > 2 ? evolved[4] : "not here yet")
        , "t$(lege[5])", (i > 3 ? evolved[5] : "not here yet")
        , "t$(lege[6])", (i > 4 ? evolved[6] : "not here yet")
        , "t$(lege[7])", (i > 5 ? evolved[7] : "not here yet")
        , "t$(lege[8])", (i > 6 ? evolved[8] : "not here yet")
        , "t$(lege[9])", (i > 7 ? evolved[9] : "not here yet")
        , "t$(lege[10])", (i > 8 ? evolved[10] : "not here yet")
    , "t$(lege[11])", (i > 9 ? evolved[11] : "not here yet")
    , "t$(lege[12])", (i > 10 ? evolved[12] : "not here yet")
 )	

    print(i,"/",timeSTEPS," time steps done\n")
    flush(stdout)
end

println("THE END")
