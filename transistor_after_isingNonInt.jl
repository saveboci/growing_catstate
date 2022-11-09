# requires as imput:
# 1) chain length
# 2) bond dim max
# 3) Jy
# 4) h 
# 5) how long we apply ising for

using ITensors
# using Plots
using JLD2

#gr()

println("BEGIN")

rescale = .7;
Jx = round(1 * rescale, digits = 2);
Jy = round(parse(Float64, ARGS[3]) * rescale, digits = 2);
Jz = round(0 * rescale, digits = 2);
h = round(parse(Float64, ARGS[4]) * rescale, digits = 2);

hx_ising = 0.7
hz_ising = 0.6
J_ising = 1 

L = parse(Int64, ARGS[1]); 
bondD = parse(Int64, ARGS[2]);
dt = 0.01

timeSTART = parse(Float64, ARGS[5]);
timeMAX = 20;

println("We take L = $(L), bond dim = $(bondD), dt = $(dt). We first evolve the state allup with an TFIC hamiltonian with overall constant $(J_ising) and external fields hz=$(hz_ising) and hx=$(hx_ising) for a time $(timeSTART). Then we evolve with the quantum transistor with Jx = $(Jx), Jy = $(Jy), Jz = $(Jz), h = $(h) for $(timeMAX) timeunits.")

sites = siteinds("S=1/2", L);




############## HERE WE DEFINE EVOLUTION ##############


# Compute the real time evolution with an Ising Hamiltonian with external field along Z.
# The output is the the MPS describing the state at the end of the evolution

function TEBD_via_Ising(psi, cutoff, bond_dim_max, timestep, ttotal, mpo = false)
        
    tau=timestep
    
    siti = siteinds(psi)
    LL = length(psi)
    
    Nsteps = round(Int, ttotal/timestep) 
    gates = ITensor[] 
    for j=1:LL-2
        s1 = siti[j]
        s2 = siti[j+1]
        Aj = - 4 * J_ising * op("Sx",s1) * op("Sx",s2)
        Bj = - 2 * J_ising * hz_ising * op("Sz",s1) * op("Id",s2) 
	Cj = - 2 * J_ising * hx_ising * op("Sx",s1) * op("Id",s2) 
        hj = Aj + Bj + Cj
        Gj = exp(-1.0im * tau/2 * hj)
        push!(gates,Gj)
    end
    s1 = siti[LL-1]
    s2 = siti[LL]
    A_end = - 4 * J_ising * op("Sx",s1) * op("Sx",s2)
    B_end = - 2 * J_ising * hz_ising * op("Sz",s1) * op("Id",s2) - 2 * J_ising * hz_ising * op("Id",s1) * op("Sz",s2) 
    C_end = - 2 * J_ising * hx_ising * op("Sx",s1) * op("Id",s2) - 2 * J_ising * hx_ising * op("Id",s1) * op("Sx",s2) 
    h_end = A_end + B_end + C_end
    G_end = exp(-1.0im * tau/2 * h_end)
    push!(gates, G_end)
    append!(gates,reverse(gates));
    
    for step in 1:Nsteps
        psi = apply(gates, psi; cutoff=cutoff, maxdim = bond_dim_max, apply_dag = mpo) 
        # When you use `apply(U, A; apply_dag=true)`, it performs the operation `U * A * U^dag`, 
        # where `U^dag` both applies `dag` to the gates as well as swaps the prime levels of the gates 
        # to perform the transpose
    end
    
    psi
end


function TEBD(psi, cutoff, bond_dim_max, timestep, ttotal)
        
    tau=timestep
    
    siti = siteinds(psi)
    LL = length(psi)
        
    Nsteps = round(Int, ttotal/timestep) 
    gates = ITensor[] 
    for j=2:LL-2
        s0 = siti[j-1]
        s1 = siti[j]
        s2 = siti[j+1]
        Aj =  2Jx * op("Sx",s0) * op("Id",s1) * op("Sx",s2) + 2Jy * op("Sy",s0) * op("Id",s1) * op("Sy",s2) + 2Jz * op("Sz",s0) * op("Id",s1) * op("Sz",s2)
        Bj = -4Jx * op("Sx",s0) * op("Sz",s1) * op("Sx",s2) - 4Jy * op("Sy",s0) * op("Sz",s1) * op("Sy",s2) - 4Jz * op("Sz",s0) * op("Sz",s1) * op("Sz",s2)
        Cj =  -2h * op("Sz",s0) * op("Id",s1) * op("Id",s2)
        hj = Aj + Bj + Cj
        Gj = exp(-1.0im * tau/2 * hj)
        push!(gates,Gj)
    end
    
    s0 = siti[LL-2]
    s1 = siti[LL-1]
    s2 = siti[LL]
    Aj =  2Jx * op("Sx",s0) * op("Id",s1) * op("Sx",s2) + 2Jy * op("Sy",s0) * op("Id",s1) * op("Sy",s2) + 2Jz * op("Sz",s0) * op("Id",s1) * op("Sz",s2)
    Bj = -4Jx * op("Sx",s0) * op("Sz",s1) * op("Sx",s2) - 4Jy * op("Sy",s0) * op("Sz",s1) * op("Sy",s2) - 4Jz * op("Sz",s0) * op("Sz",s1) * op("Sz",s2)
    Cj =  -2h * op("Sz",s0) * op("Id",s1) * op("Id",s2) - 2h * op("Id",s0) * op("Sz",s1) * op("Id",s2) - 2h * op("Id",s0) * op("Id",s1) * op("Sz",s2)
    hj = Aj + Bj + Cj
    Gj = exp(-1.0im * tau/2 * hj)
    push!(gates,Gj)
    
    append!(gates,reverse(gates));
    
    for step in 1:Nsteps
        psi = apply(gates, psi; cutoff=cutoff, maxdim = bond_dim_max) 
        # When you use `apply(U, A; apply_dag=true)`, it performs the operation `U * A * U^dag`, 
        # where `U^dag` both applies `dag` to the gates as well as swaps the prime levels of the gates 
        # to perform the transpose
    end
    
    psi
end


############## HERE WE PERFORM TIME EVOLUTION ##############


state_all_up = productMPS(sites, n -> "Up");
vmax = 1.5;

timeSTEPS = 9;

timeSTEP = (timeMAX - timeSTART) / timeSTEPS;

print("initializing step zero\n")
flush(stdout)


evolved = TEBD_via_Ising(state_all_up, 1E-12, bondD, dt, timeSTART);
evolved = apply(2 * op("Sx",sites[round(Int,L/2)]), evolved);
provisory = evolved
 
print(0,"/",timeSTEPS," time steps done\n")
flush(stdout)

times = [timeSTART + (j-1) * timeSTEP for j in 1:timeSTEPS+1];
lege = reshape([round(timeSTART + (i-1)*timeSTEP, digits=3) for i in 1:timeSTEPS+1],(1,timeSTEPS+1));

for i in 1:timeSTEPS
    global provisory = TEBD(provisory, 1E-12, bondD, dt, timeSTEP);
    global evolved = vcat(evolved, provisory)

    save("/home/sbocini/Numerics_LabDesktop/quantum_transistor/after_isingNonInt/L$(L)_Jx$(Jx)_Jy$(Jy)_Jz$(Jz)_h$(h)_Deltat$(timeSTART)_hZising$(hz_ising)_hXising$(hx_ising)_tTOT$(timeMAX)_BondDim$(bondD)_dt$(dt).jld2"
    #save("/Users/saveriobocini/numerics_PhD/quantum_transistor/after_isingNonInt/L$(L)_Jx$(Jx)_Jy$(Jy)_Jz$(Jz)_h$(h)_Deltat$(timeSTART)_hZising$(hz_ising)_hXising$(hx_ising)_tTOT$(timeMAX)_BondDim$(bondD)_dt$(dt).jld2"
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
	)

    print(i,"/",timeSTEPS," time steps done\n")
    flush(stdout)
end




############## HERE WE MAKE A FEW PLOTS ##############

#= ixes = [[(i-L/2)/times[j] for i in 1:L] for j in 1:length(evolved)];
ixes_sqrt = [[(i-L/2)/sqrt(times[j]) for i in 1:L] for j in 1:length(evolved)];


howmanyplots = 5;
start = timeSTEPS - howmanyplots + 2;

ypsilons = [expect(evolved[i], "Sz") for i=start:length(evolved)];

lege = reshape([round(timeSTART + (i-1)*timeSTEP, digits=3) for i in start:length(evolved)],(1,length(evolved)-start+1))

plot(ixes[start:length(evolved)], ypsilons
#plot(ypsilons
    , xlim = (-2,2)
    , label=lege
    , xlabel = "l"
    , ylabel = "<S(l)>"
    , legend = :bottomright) =#








println("THE END")
