using MAT
using Plots
using Statistics
include("../src/Graph_Louvain.jl")
num = 9
numtimes = 5
Mstart = 4
Mend = 4
mat = matread("Output/Amazon_$(num)_($numtimes)_$(Mstart)_$(Mend).mat")

lin_ari = mat["lin_ari"]
wce_ari = mat["wce_ari"]
sce_ari = mat["sce_ari"]
SC = mat["SC"]
WC = mat["WC"]
LinC = mat["LinC"]
numtimes = length(lin_t)
d = vec(sum(H,dims=1))
vol = sum(d)

##
M = matread("Amazon_$num.mat")
names = M["names"]
labels = M["labels"]
uLabels = [1 3 2 12 18 25 17 24 15]
uLs = uLabels[1:num]

H = M["H"]
order = vec(sum(H,dims = 2))
m,n = size(H)

## Confirm the hypergraph is one connecte component
A = H'*H
using MatrixNetworks
Acc, p = largest_component(A)

## Compute what percentage of nodes are completely contained in a cluster
mu = 0.75
High = 0
All = 0
for e = 1:m
    global mu, High, All
    edge = findall(x->x>0, H[e,:])
    clusters = labels[edge]
    mm = mode(clusters)
    prcnt = length(findall(x->x==mm,clusters))/length(clusters)
    if prcnt == 1
        All += 1
    elseif prcnt > mu
        High += 1
    end
end

@show All/m, High/m

##
using Plots

lw = 3
msw = 0
s1 = 300
s2 = 250
ms = 5
lw = 2

lin = mat["lin_ari"]
wce = mat["wce_ari"]
Lams = range(1/(Mstart*vol),Mend/vol,length=numpoints)

## Begin by plotting ARI scores for detecting different category-clusters
plot(vol*Lams, lin,color=:purple,linewidth= lw,grid = false,legend = :bottom,
xaxis = "\\lambda * vol\\(G\\)", yaxis= "ARI Scores", label = "HyperLam (Star) ",
markerstrokewidth = msw, markersize = ms, markershape = :circle)

plot!(vol*Lams, wce ,color=:green,linewidth= lw, label = "HyperLam (Clique) ",
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))
savefig("Plots/Amazon_aris.pdf")


## Now track individual clusters
Fl = zeros(num,numtimes)
Rl = zeros(num,numtimes)
for i = 1:numtimes
    c = round.(Int64,LinC[:,i])
    F1s, PRs, REs =  ClusterTracker(c,labels, uLs)
    Fl[:,i] = F1s
    Rl[:,i] = REs
end
F = zeros(num,numtimes)
R = zeros(num,numtimes)
for i = 1:numtimes
    c = round.(Int64,WC[:,i])
    F1s, PRs, REs =  ClusterTracker(c,labels, uLs)
    F[:,i] = F1s
    R[:,i] = REs
end

##

plot()
#"Software",
lbs = ["Pantry", "I&S"]
clus = [8 9]
colors = [:red, :black, :blue]
for i = 2:-1:1
    l = clus[i]
    plot!(Lams*vol, Fl[l,:], grid = false, linewidth = lw, markershape = :circle,
    label = lbs[i],legend = :topleft,markerstrokewidth = 0,  size = (s1,s2),
    markersize = ms, color = colors[i])
end
for i = 2:-1:1
    l = clus[i]
    plot!(Lams*vol, F[l,:], grid = false, linewidth = lw, markershape = :circle,
    label = "", legend = :bottomleft,linestyle = :dash, markerstrokewidth = 0,  size = (s1,s2),
    markersize = ms, color = colors[i],xaxis = "\\lambda * vol\\(G\\)", yaxis= "F1 Scores")
end
savefig("Plots/TwoClusters_F1.pdf")


##

plot()

lbs = ["Soft.", "Pantry", "I&S"]
clus = [6 8 9]
colors = [:red, :black, :blue]
for i = 1:3
    l = clus[i]
    plot!(Lams*vol, Fl[l,:], grid = false, linewidth = lw, markershape = :circle,
    label = lbs[i],legend = :topleft,markerstrokewidth = 0,  size = (s1,s2),
    markersize = ms, color = colors[i])
end
for i = 1:3
    l = clus[i]
    plot!(Lams*vol, F[l,:], grid = false, linewidth = lw, markershape = :circle,
    label = "", legend = :topleft,linestyle = :dash, markerstrokewidth = 0,  size = (s1,s2),
    markersize = ms, color = colors[i],xaxis = "\\lambda * vol\\(G\\)", yaxis= "F1 Scores")
end
savefig("Plots/ThreeClusters_F1.pdf")
