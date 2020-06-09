using MAT
using LinearAlgebra
using SparseArrays
using MatrixNetworks
using Clustering

## Load Data
mat = matread("../data/Florida_Bay_Dataset.mat")
A125 = mat["A125"]                    # This is the processed dataset--a directed who-eats-who graph
n = size(A125,1)
C = round.(Int64,mat["C"])            # These are four diffrent ways to classify the species
label = C[:,3]                        # These is the labels of Li et al. NeurIPs 2017, ignoring "predator" and "no predator" tags
Li_c = vec(round.(Int64,mat["Li_c"])) # This is the clustering obtained by Li et al.

## Form different graphs based on different reductions
include("../src/Motif_Adj_Mat.jl")
AFan = ExtractFanMat(A125)      # Motif-adjacency matrix
ASim = turntounweighted(A125)   # Symmetrized version of directed graph

# Save, to run other algorithms requiring MATLAB interface
matwrite("../data/Fl_bay_Afan.mat",Dict("AFan"=>AFan))

# Scale the graphs so that their total volume is the same
mSim = sum(ASim)
mFan = sum(AFan)
AFan = mSim/mFan*AFan
m = mSim

# Get weighted degree vectors for each graph
dSim = vec(sum(ASim,dims = 2))
dFan = vec(sum(AFan,dims = 2))


## Run Louvain-based algorithms
include("../src/Graph_Louvain.jl")
numpoints = 20  # number of lambda values to try

# For each graph, store the ari scores and number of clusters
SimLouv = zeros(numpoints,2)
FanLouv= zeros(numpoints,2)

Lams = range(0.05/m,5/m,length=numpoints)
SimC  = zeros(n,numpoints)
FanC = zeros(n,numpoints)

numtimes = 1          # number of times to run Louvain.
totaltimes = zeros(length(Lams),2)  # store runtimes

for i = 1:length(Lams)
    lam = Lams[i]
    println(" Lambda = $lam")
    # Run the standard Lambda-Louvain algorithm (i.e., the motif is an edge)
    start = time()
    cSim, BestObj = Many_Louvain(ASim,dSim,lam,numtimes)
    totaltimes[i,1] = time()-start
    SimC[:,i] = cSim
    sim_ari = ari(cSim,label)
    SimLouv[i,:] = [sim_ari length(unique(cSim))]

    # Run the Lambda-Louvain algorithm for the bi-fan motif
    start = time()
    cFan, BestObj = Many_Louvain(AFan,dFan,lam,numtimes)
    totaltimes[i,2] = time()-start
    FanC[:,i] = cFan
    fan_ari = ari(cFan,label)
    FanLouv[i,:] = [fan_ari length(unique(cFan))]
end

p1 = sortperm(SimLouv[:,2])
p2 = sortperm(FanLouv[:,2])
SimLouv = SimLouv[p1,:]
FanLouv = FanLouv[p2,:]
matwrite("Output/Florida_Bay_125_output.mat", Dict("FanC" => FanC,
    "FanLouv" => FanLouv,"Lams" => collect(Lams),"totaltimes"=>totaltimes))

## Extract the graclus clusterings
G = matread("Output/Graclus_Flbay.mat")
GracC = round.(Int64,G["GracFan"])
num = size(GracC,2)
Grac = zeros(num,2)
GracRun = G["times"][:,2]
for i = 1:num
    cG = vec(GracC[:,i])
    Grac[i,1] = ari(cG,label)
    Grac[i,2] = length(unique(cG))
end

## Extract the Metis clusterings
G = matread("Output/Metis_Flbay.mat")
MetC = round.(Int64,G["MetFan"])
num = size(MetC,2)
Met = zeros(num,2)
MetRun = G["times"][:,2]
for i = 1:num
    cM = vec(MetC[:,i])
    Met[i,1] = ari(cM,label)
    Met[i,2] = length(unique(cM))
end

## Extract the hMetis clusterings
G = matread("Output/hmetis_2to50_10.mat")
hMetC = round.(Int64,G["hmetC"])
num = size(hMetC,2)
hMet = zeros(num,2)
hMetRun = G["times"]
for i = 1:num
    cM = vec(hMetC[:,i])
    hMet[i,1] = ari(cM,label)
    hMet[i,2] = length(unique(cM))
end
# order by cluster size--this isn't always automatic
p = sortperm(hMet[:,2])
hMet = hMet[p,:]
## Extract the Recursive Spectral Clusterings
G = matread("Output/RecursiveSpectral_Flbay.mat")
SpecC = round.(Int64,G["SpecFan"])
num = size(MetC,2)
Spec = zeros(num,2)
SpecRun = G["times"][:,2]
for i = 1:num
    cM = vec(SpecC[:,i])
    Spec[i,1] = ari(cM,label)
    Spec[i,2] = length(unique(cM))
end

p = sortperm(Spec[:,2])
Spec = Spec[p,:]
## Plot results
using Plots

lw = 3
msw = 0
s1 = 300
s2 = 250
ms = 5
lw = 2

rightmax = 51

p = plot(FanLouv[:,2], FanLouv[:,1],color=:green,markerstrokecolor=:green,
markerstrokewidth = msw, markersize = ms,markershape = :circle,linewidth= lw, label = "HypLam")

plot!(SimLouv[:,2], SimLouv[:,1],color=:black,markerstrokecolor=:black,
markerstrokewidth = msw, markersize = ms,markershape = :circle,linewidth= lw, xlim = [0, rightmax],grid = false, label = "LamLouv",xaxis = "Number of Clusters")

plot!(Grac[:,2], Grac[:,1],color=:blue,markerstrokecolor=:blue,
markerstrokewidth = msw, markersize = ms,markershape = :circle,linewidth= lw,ylabel = "ARI Scores",size = (s1,s2), label = "Graclus")

plot!(Met[:,2], Met[:,1],color=:orange,markerstrokecolor=:orange,
markerstrokewidth = msw, markersize = ms,markershape = :circle,linewidth= lw,ylabel = "ARI Scores",size = (s1,s2), label = "Metis")

plot!(hMet[:,2], hMet[:,1],color=:purple,markerstrokewidth = msw, markerstrokecolor=:purple,
markersize = ms,markershape = :circle,linewidth= lw,ylabel = "ARI Scores",size = (s1,s2), label = "hMetis")

plot!(Spec[:,2], Spec[:,1],color=:brown,markerstrokewidth = msw, markersize = ms,markershape = :circle,linewidth= lw,ylabel = "ARI Scores",size = (s1,s2),label = "Spectral",
ylim = [0,.5],foreground_color_legend = nothing,background_color_legend=nothing,markerstrokecolor=:brown)

savefig(p,"plots/Florida_Bay_ari.pdf")

## Runtime plots
p = plot(FanLouv[:,2], totaltimes[:,2],color=:green,linewidth= lw, label = "")
scatter!(FanLouv[:,2],totaltimes[:,2],color=:green,markerstrokecolor=:green,markerstrokewidth = msw, markersize = ms, label = "HyperLam")

plot!(SimLouv[:,2], totaltimes[:,1],color=:black,linewidth= lw, xlim = [0, rightmax],grid = false, label = "",xaxis = "Number of Clusters")
scatter!(SimLouv[:,2],totaltimes[:,1],color=:black,markerstrokecolor=:black,markerstrokewidth = msw, markersize = ms, label = "LamLouvain",markershape= :circle)

plot!(Grac[:,2], GracRun,color=:blue,linewidth= lw, label = "",ylabel = "Runtime (seconds)",size = (s1,s2))
scatter!(Grac[:,2], GracRun,color=:blue,markerstrokecolor=:blue,markerstrokewidth = msw, markersize = ms, label = "Graclus")

plot!(Met[:,2], MetRun,color=:orange,linewidth= lw, label = "",size = (s1,s2),legend = false)
scatter!(Met[:,2], MetRun,color=:orange,markerstrokecolor=:orange,markerstrokewidth = msw, markersize = ms, label = "Metis")

plot!(hMet[:,2], hMetRun,color=:purple,linewidth= lw, label = "",size = (s1,s2),yaxis = :log10)
scatter!(hMet[:,2], hMetRun,color=:purple,markerstrokewidth = msw,
markersize = ms, label = "hMetis",markerstrokecolor=:purple)

plot!(Spec[:,2], SpecRun,color=:brown,linewidth= lw, label = "",size = (s1,s2),yaxis = :log10)
scatter!(Spec[:,2], SpecRun,color=:brown,markerstrokewidth = msw, markersize = ms, label = "Spectral",
foreground_color_legend = nothing,background_color_legend=nothing,markerstrokecolor=:brown)

savefig(p,"plots/Florida_Bay_runtime.pdf")


## ARI/Runtime Tradeoff Plots
ms = 5
p = scatter(totaltimes[:,2],FanLouv[:,1],color=:green,linewidth= lw,markerstrokecolor=:green,
markerstrokewidth = msw, markersize = ms, label = "HyperLam", legend = false)

scatter!(totaltimes[:,1],SimLouv[:,1],color=:black,linewidth= lw,
grid = false,xaxis = "Runtime", yaxis = "ARI Scores",markerstrokecolor=:black,
markerstrokewidth = msw, markersize = ms, label = "LamLouvain",markershape= :circle)

scatter!(GracRun,Grac[:,1],color=:blue,linewidth= lw,markerstrokecolor=:blue,
size = (s1,s2),markerstrokewidth = msw, markersize = ms, label = "Graclus")

scatter!(MetRun,Met[:,1],color=:orange,linewidth= lw,size = (s1,s2),markerstrokecolor=:orange,
legend = false,markerstrokewidth = msw, markersize = ms, label = "Metis")

scatter!(hMetRun,hMet[:,1],color=:purple,linewidth= lw,markerstrokecolor=:purple,
markerstrokewidth = msw, markersize = ms, label = "hMetis")

scatter!(SpecRun, Spec[:,1], color=:brown,linewidth= lw,size = (s1,s2),xaxis = :log10,
markerstrokewidth = msw, markersize = ms, label = "Spectral",markerstrokecolor=:brown,
foreground_color_legend = nothing,background_color_legend=nothing)


## Increase size of largest ARI dot
ms = 10
best = findmax(FanLouv[:,1])[2]
scatter!([totaltimes[best,2]],[FanLouv[best,1]],color=:green,linewidth= lw,markerstrokecolor=:green,
markerstrokewidth = msw, markersize = ms, label = "HyperLam", legend = false)

best = findmax(SimLouv[:,1])[2]
scatter!([totaltimes[best,1]],[SimLouv[best,1]],color=:black,linewidth= lw,
grid = false,xaxis = "Runtime", yaxis = "ARI Scores",markerstrokecolor=:black,
markerstrokewidth = msw, markersize = ms, label = "LamLouvain",markershape= :circle)

best = findmax(Grac[:,1])[2]
scatter!([GracRun[best]],[Grac[best,1]],color=:blue,linewidth= lw,markerstrokecolor=:blue,
size = (s1,s2),markerstrokewidth = msw, markersize = ms, label = "Graclus")

best = findmax(Met[:,1])[2]
scatter!([MetRun[best]],[Met[best,1]],color=:orange,linewidth= lw,size = (s1,s2),markerstrokecolor=:orange,
legend = false,markerstrokewidth = msw, markersize = ms, label = "Metis")

best = findmax(hMet[:,1])[2]
scatter!([hMetRun[best]],[hMet[best,1]],color=:purple,linewidth= lw,markerstrokecolor=:purple,
markerstrokewidth = msw, markersize = ms, label = "hMetis")

best = findmax(Spec[:,1])[2]
scatter!([SpecRun[best]], [Spec[best,1]], color=:brown,linewidth= lw,size = (s1,s2),xaxis = :log10,
markerstrokewidth = msw, markersize = ms, label = "Spectral",markerstrokecolor=:brown,
foreground_color_legend = nothing,background_color_legend=nothing)

savefig(p,"plots/Florida_Bay_tradeoff.pdf")
