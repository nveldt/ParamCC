using MAT
using LinearAlgebra
using SparseArrays
using MatrixNetworks
using Clustering

mat = matread("data/Florida_Bay_Dataset.mat")
A125 = mat["A125"]                    # This is the processed dataset--a directed who-eats-who graph
n = size(A125,1)
C = round.(Int64,mat["C"])            # These are four diffrent ways to classify the species
label = C[:,3]                        # These is the labels of Li et al. NeurIPs 2017, ignoring "predator" and "no predator" tags
Li_c = vec(round.(Int64,mat["Li_c"])) # This is the clustering obtained by Li et al.

##
include("../src/Graph_Louvain.jl")
for i = 1:4
    li_ari = ari(Li_c,C[:,i])
    println("$i $li_ari")
end

## Need to form different graphs based on different reductions
include("../src/Motif_Adj_Mat.jl")
AFan = ExtractFanMat(A125)      # So we can run Lambda-Louvain on the bi-fan hypergraph
ASim = turntounweighted(A125)   # So we can run Louvain and Graclus on an unweighted version

# Scale the graphs so that their total volume is the same, for easier comparison
mSim = sum(ASim)
mFan = sum(AFan)
AFan = mSim/mFan*AFan
m = mSim

# Get degree vectors for each of them
dSim = vec(sum(ASim,dims = 2))
dFan = vec(sum(AFan,dims = 2))

include("../src/Graph_Louvain.jl")

##
numpoints = 25

# For each motif adjacency matrix, store the ari scores and number of clusters
SimLouv = zeros(numpoints,2)
FanLouv= zeros(numpoints,2)

Lams = range(0.05/m,5/m,length=numpoints)
SimC  = zeros(n,numpoints)
FanC = zeros(n,numpoints)

numtimes = 25           # number of times to run Louvain.

for i = 1:length(Lams)
    lam = Lams[i]

    # Run the standard Louvain algorithm (i.e., the motif is an edge)
    cSim, BestObj = Many_Louvain(ASim,dSim,lam,numtimes)
    SimC[:,i] = cSim
    sim_ari = ari(cSim,label)
    SimLouv[i,:] = [sim_ari length(unique(cSim))]

    # Run the Lambda-Louvain algorithm for the bi-fan motif
    cFan, BestObj = Many_Louvain(AFan,dFan,lam,numtimes)
    FanC[:,i] = cFan
    fan_ari = ari(cFan,label)
    FanLouv[i,:] = [fan_ari length(unique(cFan))]
    println("One more through.")
end

matwrite("Output/Florida_Bay_125_M6_Bifan.mat", Dict("FanC" => FanC,
    "FanLouv" => FanLouv,"Lams" => collect(Lams)))

## Extract the graclus clusterings

num = size(GracC,2)
Grac = zeros(num,2)
G = matread("Output/Graclus_2to50.mat")
GracC = round.(Int64,G["GracC"])
for i = 1:num
    cG = vec(GracC[:,i])
    Grac[i,1] = ari(cG,label)
    Grac[i,2] = length(unique(cG))
end


## Plot some results
using Plots
p1 = sortperm(SimLouv[:,2])
p2 = sortperm(FanLouv[:,2])
##

lw = 3
msw = 0
s1 = 300
s2 = 250
ms = 5
lw = 2

rightmax = 51

val = 1
bb = 17
p = plot(FanLouv[p2[1:22],2], FanLouv[p2[1:22],val],color=:green,linewidth= lw, label = "")
scatter!(FanLouv[p2[1:22],2],FanLouv[p2[1:22],val],color=:green,markerstrokewidth = msw, markersize = ms, label = "HyperLam")
plot!(SimLouv[p1[1:bb],2], SimLouv[p1[1:bb],val],color=:black,linewidth= lw, xlim = [0, rightmax],grid = false, label = "",xaxis = "Number of Clusters")
scatter!(SimLouv[p1[1:bb],2],SimLouv[p1[1:bb],val],color=:black,markerstrokewidth = msw, markersize = ms, label = "LamLouvain",markershape= :circle)
plot!(Grac[:,2], Grac[:,1],color=:blue,linewidth= lw, label = "",ylabel = "ARI Scores",size = (s1,s2))
scatter!(Grac[:,2], Grac[:,1],color=:blue,markerstrokewidth = msw, markersize = ms, label = "Graclus")
savefig(p,"plots/Florida_Bay_125_ARI_ByClusNum.pdf")

##
val = 1
p = plot(Lams*m*2, SimLouv[:,val],color=:black,linewidth= lw, grid = false, label = "",legend = false,xaxis = "vol\\(G \\) * \\lambda")
scatter!(Lams*m*2,SimLouv[:,val],color=:black,markerstrokewidth = msw, markersize = ms, label = "Edge LamCC",markershape= :circle)
plot!(Lams*m*2, FanLouv[:,val],color=:green ,linewidth= lw,label = "")
scatter!(Lams*m*2,FanLouv[:,val],color=:green,markerstrokewidth = msw, markersize = ms, label = "Fan Motif LamCC")
savefig(p,"plots/Florida_Bay_125_ARI_ByLam.pdf")
