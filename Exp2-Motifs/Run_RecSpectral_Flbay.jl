using MAT
using SparseArrays
using LinearAlgebra
include("../src/recursive_spectral.jl")
mat = matread("../data/Florida_Bay_Dataset.mat")
A125 = mat["A125"]                    # This is the processed dataset--a directed who-eats-who graph
n = size(A125,1)
C = round.(Int64,mat["C"])            # These are four diffrent ways to classify the species
label = C[:,3]                        # These is the labels of Li et al. NeurIPs 2017, ignoring "predator" and "no predator" tags

## Need to form different graphs based on different reductions
include("../src/Motif_Adj_Mat.jl")
AFan = ExtractFanMat(A125)      # So we can run Lambda-Louvain on the bi-fan hypergraph
ASim = turntounweighted(A125)   # So we can run Louvain on an unweighted version

# Scale the graphs so that their total volume is the same, for easier comparison
mSim = sum(ASim)
mFan = sum(AFan)
AFan = mSim/mFan*AFan
m = mSim
n = size(ASim,1)

numpoints = 20
mins = 10 .^range(log10(5),log10(75),length=numpoints)
times = zeros(20,2)

SpecFan = zeros(n,numpoints)
SpecSim = zeros(n,numpoints)
SpecFari = zeros(numpoints)
SpecSari = zeros(numpoints)

for i = 1:length(mins)
    minsize = mins[i]

    # Run recursive spectral on the original graph
    start = time()
    println("Minsize = $minsize")
    sets = recursive_spectral(ASim,minsize)
    c = set2clusvec(sets,n)
    times[i,1] = time()-start
    SpecSim[:,i] .= c
    SpecSari = ari(round.(Int64,c),label)

    # Run it on the motif adjacency matrix
    start = time()
    sets = recursive_spectral(AFan,minsize)
    c = set2clusvec(sets,n)
    times[i,2] = time()-start
    SpecFan[:,i] .= c
    SpecFari[i] = ari(round.(Int64,c),label)
end

matwrite("Output/RecursiveSpectral_Flbay.mat", Dict("times"=>times,
"SpecFan"=>SpecFan, "SpecSim"=>SpecSim))
