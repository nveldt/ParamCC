##
using Plots
using Statistics
using SparseArrays
include("../src/Graph_Louvain.jl")

using MAT
num = 9
mat = matread("../data/Amazon_$num.mat")
H = mat["H"]
m,n = size(H)
NodeLabels = mat["labels"]
order = round.(Int64,vec(sum(H,dims=2)))

## Three different ways of converting a graph to hypergraph
d = vec(sum(H,dims=1))
A = [spzeros(n,n) sparse(H'); H spzeros(m,m)]
w = [d; zeros(m)]
A3 = SimpleCliqueExp(H)
M = sum(d)
A2 = WeightedCliqueExp(H,order)
A2 = M/sum(A2)*A2
A3 = M/sum(A3)*A3

w2 = vec(sum(A2,dims=2))
w3 = vec(sum(A3,dims=2))

matwrite("Matrices/Amazon_($num).mat", Dict("A"=>A, "A2"=>A2, "A3"=>A3))
##
numtimes = 5
numpoints = 10
lin_ari = zeros(numpoints)
wce_ari = zeros(numpoints)
sce_ari = zeros(numpoints)

lin_nc = zeros(numpoints)
wce_nc = zeros(numpoints)
sce_nc = zeros(numpoints)

lin_t = zeros(numpoints)
wce_t = zeros(numpoints)
sce_t= zeros(numpoints)

LinC = zeros(n,numpoints)
WC = zeros(n,numpoints)
SC = zeros(n,numpoints)

Mstart = 4
Mend = 4
Lams = range(1/(Mstart*M),Mend/M,length=numpoints)
numits = 1000
for i = 1:length(Lams)


    lam = round(Lams[i],digits = 5)

    # Run the method with a linear cut
    s = time()
    c, BestObj = Many_Louvain(A,w,lam,numtimes,numits)
    lin_t[i] = time()-s
    cLin = c[1:n]
    LinC[:,i] = cLin
    lin_ari[i] = ari(cLin,NodeLabels)
    lin_nc[i] = length(unique(cLin))

    s = time()
    cW, BestObj = Many_Louvain(A2,w2,lam,numtimes,numits)
    wce_t[i] = time()-s
    WC[:,i] = cW
    wce_ari[i] = ari(cW,NodeLabels)
    wce_nc[i] = length(unique(cW))

    s = time()
    cS, BestObj = Many_Louvain(A2,w3,lam,numtimes,numits)
    sce_t[i] = time()-s
    SC[:,i] = cS
    sce_ari[i] = ari(cS,NodeLabels)
    sce_nc[i] = length(unique(cS))

    println("$(lin_ari[i]) \t $(wce_ari[i]) \t $(sce_ari[i]) ")

end

##
matwrite("Output/Amazon_$(num)_($numtimes)_$(Mstart)_$(Mend).mat", Dict(
"lin_ari"=>lin_ari, "lin_t"=>lin_t, "lin_nc"=>lin_nc, "LinC"=>LinC,
"wce_ari"=>wce_ari, "wce_t"=>wce_t, "wce_nc"=>wce_nc, "WC"=>WC,
"sce_ari"=>sce_ari, "sce_t"=>sce_t, "sce_nc"=>sce_nc,   "SC"=>SC))
