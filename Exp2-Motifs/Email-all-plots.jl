## Plot results from email experiment
using MAT
mat = matread("../data/emailEUcore.mat")
label = vec(round.(Int64,mat["labels"]))

mat3  = matread("../data/EmailTriCore.mat")
coreinds = mat3["core_inds"]
label = label[coreinds]

mat = matread("Output/Louvain_Email_clusterings.mat")
Lams = mat["Lams"]

## Get runtimes and clusterings
# Get clusters from running triangle motif clustering...
triC = mat["taC"]
tt = mat["totaltimes"]
TriRun = tt[:,2]

# And standard Louvain clustering with resolution parameter lambda
SimC = mat["SimC"]
SimRun = tt[:,1]

# Results from running other algorithms
grac = matread("Output/Graclus_clusterings_2to340.mat")
GracC = grac["GracCtri"]
GracRun = grac["times"][:,2]

met = matread("Output/Metis_2to340.mat")
MetC = met["MetCtri"]
MetRun = met["times"][:,2]

hmet = matread("Output/hmetis_2to340.mat")
hMetC = hmet["hmetC"]
hMetRun = hmet["times"][:,1]

spec = matread("Output/RecursiveSpectral_Email.mat")
specC = spec["SpecTri"]
SpecRun = spec["times"][:,2]

# Check ARI scores
numpoints = length(Lams)
SimARI = zeros(numpoints)
TriARI = zeros(numpoints)
SimSize = zeros(numpoints)
TriSize = zeros(numpoints)
GracARI= zeros(numpoints)
GracSize = zeros(numpoints)
MetARI= zeros(numpoints)
MetSize = zeros(numpoints)
hMetARI= zeros(numpoints)
hMetSize = zeros(numpoints)
SpecARI= zeros(numpoints)
SpecSize = zeros(numpoints)

include("../src/Graph_Louvain.jl")
for i = 1:length(Lams)
    lam = Lams[i]
    cSim = round.(Int64,SimC[:,i])
    SimARI[i] = ari(cSim,label)
    SimSize[i] = maximum(cSim)

    cTri = round.(Int64,triC[:,i])
    TriARI[i] = ari(cTri,label)
    TriSize[i] = maximum(cTri)

    cGrac = round.(Int64,GracC[:,i])
    GracARI[i] = ari(cGrac,label)
    GracSize[i] = maximum(cGrac)

    cMet = round.(Int64,MetC[:,i])
    MetARI[i] = ari(cMet,label)
    MetSize[i] = maximum(cMet)

    chMet = round.(Int64,hMetC[:,i])
    hMetARI[i] = ari(chMet,label)
    hMetSize[i] = maximum(chMet)

    cSpec = round.(Int64,specC[:,i])
    SpecARI[i] = ari(cSpec,label)
    SpecSize[i] = maximum(cSpec)

end

# After checking, rec-spec is the only alg that
# produces cluster sizes in different orders.
# We adjust for clearest comparison in plots
p = sortperm(SpecSize)
SpecSize = SpecSize[p]
SpecARI = SpecARI[p]
SpecRun = SpecRun[p]

## Make Plots
using Plots

lw = 3
msw = 0
s1 = 300
s2 = 250
ms = 5
lw = 2

## Runtime vs. ARI
p=scatter(TriRun, TriARI,color=:green,linewidth= lw,label = "HypLam",markerstrokecolor = :green,
markerstrokewidth = msw, markersize = ms, markershape = :circle, legend = false)

scatter!(SimRun, SimARI,color=:black,linewidth= lw,grid = false, xscale = :log10,
xaxis = "Runtime", yaxis= "ARI Scores", label = "LamLouv",markerstrokecolor = :black,
markerstrokewidth = msw, markersize = ms, markershape = :circle)

scatter!(GracRun, GracARI,color=:blue,linewidth= lw, label = "Graclus",markerstrokecolor = :blue,
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))

scatter!(MetRun, MetARI,color=:orange,linewidth= lw, label = "Metis",markerstrokecolor = :orange,
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))

scatter!(hMetRun, hMetARI,color=:purple,linewidth= lw, label = "hMetis",markerstrokecolor = :purple,
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))

scatter!(SpecRun, SpecARI,color=:brown,linewidth= lw, label = "Spectral",markerstrokecolor = :brown,
markerstrokewidth = msw, markersize = ms,markershape = :circle, size = (s1,s2),
foreground_color_legend = nothing,background_color_legend=nothing)

##
ms = 10
scatter!([TriRun[findmax(TriARI)[2]]], [TriARI[findmax(TriARI)[2]]],color=:green,linewidth= lw,label = "HypLam",
markerstrokewidth = msw, markersize = ms, markershape = :circle,markerstrokecolor = :green)

scatter!([SimRun[findmax(SimARI)[2]]], [SimARI[findmax(SimARI)[2]]],color=:black,linewidth= lw,grid = false, xscale = :log10,
xaxis = "Runtime", yaxis= "ARI Scores", label = "LamLouv",markerstrokecolor = :black,
markerstrokewidth = msw, markersize = ms, markershape = :circle)

scatter!([GracRun[findmax(GracARI)[2]]], [GracARI[findmax(GracARI)[2]]],color=:blue,linewidth= lw, label = "Graclus",
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2),markerstrokecolor = :blue)

scatter!([MetRun[findmax(MetARI)[2]]], [MetARI[findmax(MetARI)[2]]],color=:orange,linewidth= lw, label = "Metis",
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2),markerstrokecolor = :orange)

scatter!([hMetRun[findmax(hMetARI)[2]]], [hMetARI[findmax(hMetARI)[2]]],color=:purple,linewidth= lw, label = "hMetis",
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2),markerstrokecolor=:purple)

scatter!([SpecRun[findmax(SpecARI)[2]]], [SpecARI[findmax(SpecARI)[2]]],color=:brown,linewidth= lw, label = "Spectral",
markerstrokewidth = msw, markersize = ms,markershape = :circle, size = (s1,s2),markerstrokecolor=:brown)

savefig(p,"Plots/Email-tradeoff.pdf")

## ARI vs. ClusterNumbers
ms = 5
rightmax = minimum( [maximum(SimSize), maximum(TriSize)])
p=plot(TriSize, TriARI,color=:green,markerstrokecolor = :green,linewidth= lw,label = "HypLam",
markerstrokewidth = msw, markersize = ms, markershape = :circle, legend = :topright)

plot!(SimSize, SimARI,color=:black,markerstrokecolor = :black,linewidth= lw,
xlim = [0, rightmax],grid = false,
xaxis = "Number of Clusters", yaxis= "ARI Scores", label = "LamLouv",
markerstrokewidth = msw, markersize = ms, markershape = :circle)

plot!(GracSize, GracARI,color=:blue,markerstrokecolor = :blue,linewidth= lw, label = "Graclus",
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))

plot!(MetSize, MetARI,color=:orange,markerstrokecolor = :orange,
linewidth= lw, label = "Metis",
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))

plot!(hMetSize, hMetARI,color=:purple,markerstrokecolor = :purple,
linewidth= lw, label = "hMetis",
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))

plot!(SpecSize, SpecARI,color=:brown,markerstrokecolor = :brown,
linewidth= lw, label = "Spectral",
markerstrokewidth = msw, markersize = ms,markershape = :circle, size = (s1,s2),
foreground_color_legend = nothing,background_color_legend=nothing)

savefig(p,"Plots/Email-ari.pdf")

## Runtime vs cluster-number

p=plot(TriSize, TriRun,color=:green,markerstrokecolor = :green,
linewidth= lw,label = "HyperLam",yscale = :log10,
markerstrokewidth = msw, markersize = ms, markershape = :circle, legend = :topright)

plot!(SimSize, SimRun,color=:black,linewidth= lw,
xlim = [0, rightmax],grid = false,markerstrokecolor = :black,
xaxis = "Number of Clusters", yaxis= "Runtime (seconds)", label = "LamLouv",
markerstrokewidth = msw, markersize = ms, markershape = :circle)

plot!(GracSize, GracRun,color=:blue,markerstrokecolor = :blue,
linewidth= lw, label = "Graclus",
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))

plot!(MetSize, MetRun,color=:orange,markerstrokecolor = :orange,
linewidth= lw, label = "Metis",
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))

plot!(hMetSize, hMetRun,color=:purple,markerstrokecolor = :purple,
linewidth= lw, label = "hMetis", legend = false,
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))

plot!(SpecSize, SpecRun,color=:brown,markerstrokecolor = :brown,
linewidth= lw, label = "Spectral",
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))

savefig(p,"Plots/Email-runtime.pdf")
