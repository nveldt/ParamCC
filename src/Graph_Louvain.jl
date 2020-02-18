using Random
using SparseArrays

# From the adjacency matrix, build an adjacency list for the graph
function ConstructAdj(C::SparseMatrixCSC,n::Int64)
    rp = C.rowval
    ci = C.colptr
    Neighbs = Vector{Vector{Int64}}()
    d = zeros(Int64,n)
    for i = 1:n
        # chop up the rp vector and put it in Neighbs
        push!(Neighbs,rp[ci[i]:ci[i+1]-1])
        d[i] = ci[i+1]-ci[i]
    end

    # d is the number of neighbors. This is the unweighted degree,
    # but note importantly that if the original graph is weighted this is
    # not the same as the degree vector d we will sometimes use
    return Neighbs, d
end

function LamCCobj(A::SparseMatrixCSC,c,w,lam)

    w_volA = sum(w)
    # d = sum(A,dims = 1)
    # volA = sum(d)
    obj = (lam*(w_volA)^2-lam*sum(w.^2))/2

    for i = 1:maximum(c)
        S = findall(x->x==i,c)
        AS = A[:,S]
        vol = sum(AS.nzval);
        SAS = AS[S,:]
        edges = sum(SAS.nzval);
        cut = vol-edges
        w_volS = sum(w[S])
        obj += 1/2*(cut - lam*w_volS*(w_volA-w_volS))
    end
    return obj
end


# Take a clustering vector that may have gaps in the numbers, and re-name clusters
function renumber(c::Vector{Int64},Clusters::Vector{Vector{Int}})

    n = length(c)
    map = unique(c)
    cnew = zeros(Int64,n)

    Clusters = Clusters[map]

    # Rename the clusters
    for i = 1:n
        newClus = findfirst(x->x == c[i],map)
        cnew[i] = newClus
        push!(Clusters[newClus],i)
    end

    return cnew, Clusters

end


function renumber(c::Vector{Int64})

    n = length(c)
    map = sort(unique(c))
    cnew = zeros(Int64,n)

    # Rename the clusters
    for i = 1:n
        newClus = findfirst(x->x == c[i],map)
        cnew[i] = newClus
    end

    return cnew

end


function Many_Louvain(A::SparseMatrixCSC{Float64,Int64},w::Vector{Float64},lam::Float64,numtimes::Int64,maxits::Int64=10000)

    n = size(A,1)
    BestObj = Inf
    cBest = collect(1:n)

    for k = 1:numtimes
        Cs = LambdaLouvain(A,w,lam,true,maxits)
        c = Cs[:,end]
        obj = LamCCobj(A,c,w,lam)
        if obj < BestObj
            BestObj = obj
            cBest = c
        end
    end

    return cBest, BestObj
end


function LambdaLouvain(A::SparseMatrixCSC{Float64,Int64},w::Vector{Float64},lam::Float64,randflag::Bool=false,maxits::Int64=10000)

    @assert(issymmetric(A))
    n = size(A,1)
    # Run a first step
    c, improved = LambdaLouvain_Step(A,w,lam,randflag,maxits)

    # Keep track of every improvement you find
    if improved
        Cs = c
        c_old = copy(c)
    else
        Cs = c
    end

    while improved
        # Collapse the clustering and run again
        Anew, wnew = collapse_clustering(A,w,c_old)
        cSuper, improved = LambdaLouvain_Step(Anew,wnew,lam,randflag,maxits)
        N = length(wnew)
        if improved
            # Extract what that new clustering is
            c_new = zeros(Int64,n)
            # For each supernode, place all the nodes that make it up into a cluster
            # of the supernodes label
            for i = 1:N

                # What cluster is supernode i in?
                SuperI_Cluster = cSuper[i]

                # What individual nodes are in supernode i?
                SuperI_nodes = findall(x->x==i,c_old)
                c_new[SuperI_nodes] .= SuperI_Cluster
            end
            Cs = [Cs c_new]
            c_old = copy(c_new)
        end
    end

    return Cs
end

function LambdaLouvain_Step(A::SparseMatrixCSC{Float64,Int64},w::Vector{Float64},lam::Float64,randflag::Bool=false,maxits::Int64=Inf)
    @assert(issymmetric(A))
    n = size(A,1)
    # println("Merging $n Communities")

    if randflag
        p = randperm(n)
        A = A[p,p]
        undop = sortperm(p)
        w = w[p]
    end

    c = collect(1:n)
    @assert(size(w,1) == n)
    improving = true
    nextclut = n+1
    its = 0
    Neighbs, degvec = ConstructAdj(A,n)

    # The find function is inefficient for getting nodes in a cluster.
    # We will store a vector of n clusters (by the end, many will be empty clusters)
    Clusters = Vector{Vector{Int64}}()
    for v = 1:n
        push!(Clusters, Vector{Int64}())
        push!(Clusters[v],v)
    end

    # As long as we can keep improving, continue
    while improving && its < maxits
        # println("Iter $its")
        its += 1
        # println("Iteration $its")
        improving = false
        nextclus = 1

        for i = 1:n

            # println("\n\t NODE $i")

            # Cluster index
            Ci_ind = c[i]

            # Get the indices of nodes in i's cluster
            # Ci = findall(x->x == Ci_ind,c)
            Ci = Clusters[Ci_ind]

            # Get the indices of i's neighbors
            # Ni = findall(x->x>0, A[i,:])
            Ni = Neighbs[i]

            # Find the indices of nodes that neighbor i and are in c[i]
            # ci_neighbs = intersect(Ni,Ci)

            # Weight of negative mistakes currently at i
            neg_inner = w[i]*(sum(w[Ci]) - w[i])

            # Weight of positive mistakes if we remove i from Ci
            # pos_inner = sum(A[i,ci_neighbs])
            pos_inner = sum(A[i,Ci])

            # Increase in mistakes if we remove i from Ci
            total_inner = pos_inner - lam*neg_inner

            BestImprove = 0
            BestC = Ci_ind

            # Get the neighboring clusters of i
            NC = unique(c[Ni])
            #NC = Ni

            # @show NC
            # See if it's better to move to a nearby cluster, Cj
            for j = 1:length(NC)
                Cj_ind = NC[j]

                # Check how much it would improve to move i to to cluster j
                if Cj_ind == Ci_ind

                    change = 0

                else

                    # Cj = findall(x->x == Cj_ind,c)
                    Cj = Clusters[Cj_ind]

                    # Find the neighbors of i in Cj
                    # cj_neighbs = intersect(Ni,Cj)

                    # Moving i from Ci to Cj adds negative mistakes
                    neg_outer = w[i]*(sum(w[Cj]))

                    # Moving i from Ci to Cj decreases positive mistakes
                    # pos_outer = sum(A[i,cj_neighbs])
                    pos_outer = sum(A[i,Cj])

                    total_outer = lam*neg_outer - pos_outer

                    # This is the overall change in objective if we move i to Cj
                    change = total_outer + total_inner

                end

                # Check if this would be a more beneficial change than other options
                if change < BestImprove
                    BestImprove = change
                    BestC = Cj_ind
                    improving = true
                end
            end

            # Move i to the best cluster to move it to
            if BestImprove < 0
                ci_old = c[i]
                # if its > 100
                #     println("moved node $i from $ci_old to $BestC")
                #     obj = LamCCobj_slow(A,c,w,lam)
                #     println("obj = $obj")
                # end
                c[i] = BestC

                # Remove i from its old cluster
                Clusters[ci_old] = setdiff(Clusters[ci_old],i)

                # And add it to its new cluster
                push!(Clusters[BestC],i)

                improving = true
            end

        end
    end
    if its == 1
        improved = false
    else
        improved = true
        c, Clusters = renumber(c,Clusters)
    end

    if randflag
        c = c[undop]
    end

    return c, improved

end

# Collapse a clustering into a new network of supernodes and weighted edges
function collapse_clustering(A::SparseMatrixCSC{Float64,Int64},w::Vector{Float64},c::Vector{Int64})

    n = size(A,1)
    # Fill cluster arrays (faster than working with the "find" function)
    Clusters = Vector{Vector{Int64}}()
    for v = 1:maximum(c)
        push!(Clusters, Vector{Int64}())
    end

    for i = 1:n
        push!(Clusters[c[i]],i)
    end

    # Number of supernodes to form
    N = round(Int64,maximum(c))

    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    wnew = zeros(N)
    Anew = zeros(N,N)
    start = time()
    # Construct a new sparse matrix with new node weights
    for i = 1:N
        #Ci = findall(x->x == i,c)
        Ci = Clusters[i]
        wnew[i] = sum(w[Ci])
        ACi = A[Ci,:]
        for j = i+1:N
            Cj = Clusters[j]
            Eij = sum(ACi[:,Cj])
            Anew[i,j] = Eij
        end
    end
    getedges = time()-start

    start = time()
    # Anew = sparse(I,J,V,N,N)
    Anew = Anew+Anew'
    Anew = sparse(Anew)

    return Anew, wnew
end

function SimpleCliqueExp(H::SparseMatrixCSC{Float64,Int64})

    A = H'*H
    for i = 1:size(A,1)
        A[i,i] = 0.0
    end
    dropzeros!(A)
    return A
end

function WeightedCliqueExp(H::SparseMatrixCSC{Float64,Int64}, order::Vector{Int64})

    m,n = size(H)
    I = Vector{Int64}()
    J = Vector{Int64}()
    vals = Vector{Float64}()
    Hyperedges = incidence2elist(H)
    for e = 1:m
        Edge = Hyperedges[e]
        Ord = order[e]
        for ii = 1:length(Edge)
            for jj = ii+1:length(Edge)
                i = Edge[ii]
                j = Edge[jj]
                push!(I,i); push!(J,j); push!(vals,1/(Ord-1))
            end
        end
    end

    A = sparse(I,J,vals,n,n)
    A = sparse(A+A')
    return A
end


## Sort sizes: Return a set of cluster sizes, arranged in order
function ClusterSizes(c)
    c = renumber(c)

    sizes = Vector{Int64}()
    for i = 1:maximum(c)
        inds = findall(x->x==i,c)
        push!(sizes,length(inds))
    end

    return sort(sizes,rev=true)
end


# return precision, recall, and f1 scores for how well
# the clusters in c "track" the ground truth in labels
function ClusterTracker(c,labels, uLab)

    # uLab = unique(labels)
    nl = length(uLab)

    Clusters = Vector{Vector{Int64}}()
    for j = 1:length(unique(c))
        push!(Clusters,Vector{Int64}())
    end
    for i = 1:length(c)
        push!(Clusters[c[i]],i)
    end

    ## For each ground truth cluster, find the best
    F1s = zeros(nl)
    PRs = zeros(nl)
    REs = zeros(nl)
    for i = 1:nl
        label = uLab[i]

        clus = findall(x->x==label,labels)
        f1 = 0; pr = 0; re = 0
        # @show length(clus)
        # For each cluster, track the best pr, re, and f1
        for j = 1:length(Clusters)
            myclus = Clusters[j]
            p,r,f = PRF(clus,myclus)
            # println("Cluster $i, mycluster $j : $p \t $r \t $f")
            if p > pr
                pr = p
            end
            if r > re
                re = r
            end
            if f > f1
                f1 = f
            end
        end
        F1s[i] = f1
        PRs[i] = pr
        REs[i] = re
    end
    return F1s, PRs, REs
end


# This computes the precision, recall, and F1 score for a set Returned
# compared against a Target set
function PRF(Target,Returned)

    if length(Returned) == 0
        pr = 0; re = 0; F1 = 0
    else
        TruePos = intersect(Returned,Target)
        pr = length(TruePos)/length(Returned)
        re = length(TruePos)/length(Target)
        F1 = 2*(pr*re)/(pr+re)

        if length(TruePos) == 0
            F1 = 0
        end
    end

    return pr, re, F1

end

# ARI scores
using Statistics
using Clustering

function ari(x,y)
    evaluations = randindex(x, y)
    ari = evaluations[1]
    return ari
end
