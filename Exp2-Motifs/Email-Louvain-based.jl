using MAT
using LinearAlgebra
using SparseArrays
using MatrixNetworks
include("../src/Graph_Louvain.jl")
mat = matread("../data/emailEUcore.mat")
A = mat["A"]            # Largest connected component of email EU
ASim = mat["Asim"]      # Undirected version of same graph: ASim = spones(A+A')
m = sum(A)
label = vec(round.(Int64,mat["labels"]))

## Performs weighted triangle clique expansion on an already undirected graph
include("../src/Motif_Adj_Mat.jl")
Ata = ExtractTriMat(ASim)

## The clique expanded graph isn't connected, so extract the largest connected component
Acc, p = largest_component(Ata)
core_inds = findall(x->x==true,p)
matwrite("../data/EmailTriCore.mat", Dict("core_inds"=>core_inds))
Ata = Ata[p,p]
ASim = ASim[p,p]    # Symmetric version of the graph
label = label[p]

matwrite("../data/EmailTriangle.mat",Dict("Ata"=>Ata))
n = size(ASim,1)

# Get total weight of graphs
mta = sum(Ata)
mSim = sum(A)

# Normalize so that total weight is the same for both graphs
ASim = mta/mSim*ASim

## Get weighted degree vectors
dta = vec(sum(Ata,dims = 2))
dSim = vec(sum(ASim,dims = 2))

numpoints = 20
Lams = range(0.05/m,3/m,length=numpoints)

# Store clusterings
SimC  = zeros(n,numpoints)
taC = zeros(n,numpoints)    # triangle motifs

##
maxits = 1000           # Louvain maximum iterations
numtimes = 1            # Number of times to run Lambda-Lovain for each
totaltimes = zeros(length(Lams),2)
for i = 1:length(Lams)
    lam = Lams[i]

    println("Lambda = $lam")
    start = time()
    cSim, BestObj = Many_Louvain(ASim,dSim,lam,numtimes,maxits)
    totaltimes[i,1] = time()-start
    SimC[:,i] = cSim

    start = time()
    cta, BestObj = Many_Louvain(Ata,dta,lam,numtimes,maxits)
    totaltimes[i,2] = time()-start
    taC[:,i] = cta
end

matwrite("Output/Louvain_Email_clusterings.mat", Dict("taC" => taC,
    "SimC" => SimC,"Lams" => collect(Lams),"totaltimes"=>totaltimes))
