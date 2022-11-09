using ITensors
# using Plots
#using LinearAlgebra
using JLD2

#gr()

println("BEGIN")
println("Quantum transistor, starting with all up and a spin flip.")

J  = 1
γ  = parse(Float64, ARGS[3])
w  = parse(Float64, ARGS[4])
Δ  = parse(Float64, ARGS[5])
Dz = parse(Float64, ARGS[6])
hz = parse(Float64, ARGS[7])

L = parse(Int64, ARGS[1]); 
bondD = parse(Int64, ARGS[2]);
dt = 0.01

println("L = $(L), J = $(J), gamma = $(γ), w = $(w), Delta = $(Δ), Dz = $(Dz), hz = $(hz), bond dim = $(bondD), dt = $(dt)")

sites = siteinds("S=1/2", L);


############## HERE WE PERFORMTIME EVOLUTION ##############

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

        # terms without sigmaZ in the projector
        term1a = 0.25J * (1+γ) * op("Sx",s0) * op("Id",s1) * op("Sx",s2) + 0.25J * (1-γ) * op("Sy",s0) * op("Id",s1) * op("Sy",s2) 
        term2a = 0.5J * Δ * op("Sz",s0) * op("Id",s1) * op("Sz",s2)
        term3a = 0.5J * w * op("Sx",s0) * op("Id",s1) * op("Sy",s2) + 0.5J * w * op("Sy",s0) * op("Id",s1) * op("Sx",s2) 
        term4a = 0.5Dz * op("Sx",s0) * op("Id",s1) * op("Sy",s2) - 0.5Dz * op("Sy",s0) * op("Id",s1) * op("Sx",s2) 

        # terms with sigmaZ in the projector
        term1b = -0.5J * (1+γ) * op("Sx",s0) * op("Sz",s1) * op("Sx",s2) - 0.5J * (1-γ) * op("Sy",s0) * op("Sz",s1) * op("Sy",s2) 
        term2b = -J * Δ * op("Sz",s0) * op("Sz",s1) * op("Sz",s2)
        term3b = -J * w * op("Sx",s0) * op("Sz",s1) * op("Sy",s2) - J * w * op("Sy",s0) * op("Sz",s1) * op("Sx",s2) 
        term4b = -Dz * op("Sx",s0) * op("Sz",s1) * op("Sy",s2) + Dz * op("Sy",s0) * op("Sz",s1) * op("Sx",s2)     
        
        # external field
        term5 = -hz * op("Id",s0) * op("Sz",s1) * op("Id",s2)

        #everything together
        hj = term1a + term2a + term3a + term4a + term1b + term2b + term3b + term4b + term5
        Gj = exp(-1.0im * tau/2 * hj)
        push!(gates,Gj)
    end
    

    append!(gates,reverse(gates));
    
    for step in 1:Nsteps
        psi = apply(gates, psi; cutoff=cutoff, maxdim = bond_dim_max) 
        # When you use `apply(U, A; apply_dag=true)`, it performs the operation `U * A * U^dag`, 
        # where `U^dag` both applies `dag` to the gates as well as swaps the prime levels of the gates 
        # to perform the transpose
    end
    
    psi
end

state_all_up = productMPS(sites, n -> "Up");
vmax = 1.5;

timeSTART = 0.2;
timeSTEPS = 12;

timeMAX = (L/2 - 6) / vmax;
timeSTEP = (timeMAX - timeSTART) / timeSTEPS;

print("initializing step zero\n")
flush(stdout)

psi0 = apply(2 * op("Sx",sites[round(Int,L/2)]), state_all_up);
evolved = TEBD(psi0, 1E-12, bondD, dt, timeSTART);
provisory = evolved
 
print(0,"/",timeSTEPS," time steps done\n")
flush(stdout)

times = [timeSTART + (j-1) * timeSTEP for j in 1:timeSTEPS+1];
lege = reshape([round(timeSTART + (i-1)*timeSTEP, digits=3) for i in 1:timeSTEPS+1],(1,timeSTEPS+1));

for i in 1:timeSTEPS
    global provisory = TEBD(provisory, 1E-12, bondD, dt, timeSTEP);
    global evolved = vcat(evolved, provisory)

    savingarray = []
    for k = 1:timeSTEPS+1
        savingarray = vcat(savingarray, "t$(lege[k])")
        savingarray = vcat(savingarray, i+1 >= k ? evolved[k] : "not here yet")
    end

    # first line for Desktop, second line for laptop
    save("/home/sbocini/Numerics_LabDesktop/quantum_transistor/fullglory/L$(L)_J$(J)_gamma$(γ)_w$(w)_Delta$(Δ)_Dz$(Dz)_hz$(hz)_BondDim$(bondD)_dt$(dt).jld2"
    #save("/Users/saveriobocini/numerics_PhD/quantum_transistor/fullglory/L$(L)_J$(J)_gamma$(γ)_w$(w)_Delta$(Δ)_Dz$(Dz)_hz$(hz)_BondDim$(bondD)_dt$(dt).jld2"
        , "times", times
        , savingarray...
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
