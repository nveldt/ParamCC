using MAT
using SparseArrays
using LinearAlgebra
include("../src/recursive_spectral.jl")
mat = matread("../data/emailEUcore.mat")
A = mat["A"]            # Largest connected component of email EU
ASim = mat["Asim"]      # Undirected version of same graph: ASim = spones(A+A')
m = sum(A)
label = vec(round.(Int64,mat["labels"]))

include("../src/Motif_Adj_Mat.jl")
Ata = ExtractTriMat(ASim)
Acc, p = largest_component(Ata)
core_inds = findall(x->x==true,p)
Ata = Ata[p,p]      # Just extracting triangles
ASim = ASim[p,p]    # Symmetric version of the graph
label = label[p]
n = size(ASim,1)

numpoints = 20
mins = 10 .^range(log10(5),log10(500),length=numpoints)
times = zeros(20,2)

SpecTri = zeros(n,numpoints)
SpecSim = zeros(n,numpoints)

for i = 1:length(mins)
    minsize = mins[i]

    # Run recursive spectral on the original graph
    start = time()
    @show minsize
    sets = recursive_spectral(ASim,minsize)
    c = set2clusvec(sets,n)
    times[i,1] = time()-start
    SpecSim[:,i] .= c

    # Run it on the motif adjacency matrix
    start = time()
    sets = recursive_spectral(Ata,minsize)
    c = set2clusvec(sets,n)
    times[i,2] = time()-start
    SpecTri[:,i] .= c

end

matwrite("Output/RecursiveSpectral_Email.mat", Dict("times"=>times,
"SpecTri"=>SpecTri, "SpecSim"=>SpecSim))
