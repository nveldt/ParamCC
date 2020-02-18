using MAT
using LinearAlgebra
using SparseArrays
using MatrixNetworks
mat = matread("emailEUcore.mat")
A = mat["A"]            # Largest connected component of email EU
ASim = mat["Asim"]      # Undirected version of same graph: ASim = spones(A+A')
m = sum(A)
label = vec(round.(Int64,mat["labels"]))
include("../src/Motif_Adj_Mat.jl")

## Performs weighted triangle clique expansion on an already undirected graph
Ata = ExtractTriMat(ASim)

## Here, the motif is a complete triangle, with 6 directed edges. More restrictive
Atb = ExtractAllTriMat(A)

## The clique expanded graph isn't connected, so extract the largest connected component
Acc, p = largest_component(Ata)
core_inds = findall(x->x==true,p)
matwrite("EmailTriCore.mat", Dict("core_inds"=>core_inds))
Ata = Ata[p,p]
Atb = Atb[p,p]
ASim = ASim[p,p]
label = label[p]

n = size(ASim,1)

# Get total weight of graphs
mta = sum(Ata)
mtb = sum(Atb)
mSim = sum(A)

# Normalize so that total waith is the same for all graphs
Atb = mta/mtb*Atb
ASim = mta/mSim*ASim

## Get weighted degree vectors
dta = vec(sum(Ata,dims = 2))
dtb = vec(sum(Atb,dims = 2))
dSim = vec(sum(ASim,dims = 2))

include("../src/Graph_Louvain.jl")
include("../src/ClusterValidation.jl")

##
numpoints = 20

# For each motif adjacency matrix, store the
# ari, nmi, purity, probility valitation, number of clusters
SimLouv = zeros(numpoints,5)
taLouv = zeros(numpoints,5)
tbLouv = zeros(numpoints,5)
Lams = range(0.05/m,3/m,length=numpoints)

# Store clusterings
SimC  = zeros(n,numpoints)
taC = zeros(n,numpoints)
tbC = zeros(n,numpoints)
##
maxits = 1000
numtimes = 5            # Number of times to run Lambda-Lovain for each
for i = 1:length(Lams)
    lam = Lams[i]

    cSim, BestObj = Many_Louvain(ASim,dSim,lam,numtimes,false,maxits)
    SimC[:,i] = cSim
    v1, v2, v3, v4 = AllValidations(cSim,label)
    SimLouv[i,:] = [v1, v2, v3, v4, maximum(cSim)]

    cta, BestObj = Many_Louvain(Ata,dta,lam,numtimes,false,maxits)
    taC[:,i] = cta
    v1, v2, v3, v4 = AllValidations(cta,label)
    taLouv[i,:] = [v1, v2, v3, v4, maximum(cta)]

    println("One more step done")
end

matwrite("Output/Simple_Email_clusterings.mat", Dict("tbC" => tbC,"taC" => taC,"SimC" => SimC,
     "taLouv" => taLouv, "tbLouv" => tbLouv, "SimLouv" => SimLouv, "Lams" => collect(Lams)))
