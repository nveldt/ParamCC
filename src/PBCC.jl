# Functions for Parametric Bipartite Correlation Clustering
using JuMP
using Gurobi
using SparseArrays
using LinearAlgebra

gurobi_env = Gurobi.Env()


"""
Returns the linear programming relaxation of the Parametric Bipartite
Correlation Clustering objective on a bipartite graph H with parameters
mu and beta.
"""
function PBCC_LP(H::SparseMatrixCSC,mu::Float64,beta::Float64,outputflag::Bool=true,optimalflag::Bool=false)

    m = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=outputflag, gurobi_env))

    l,r = size(H)           # number of nodes on each side
    n = l+r                 # total number of nodes in graph
    L = collect(1:l)        # Left indices
    R = collect(l+1:(n))    # Right indices, R[1] is the global index of 1st node on right side

    Hbeta = Matrix(H)-beta*ones(l,r)
    A = [-mu*ones(l,l) Hbeta; Matrix(Hbeta') -mu*ones(r,r)]


    Al = H*H'       # Gives neighbors of left side: Al[i,j] = 1 if nodes i and j, on same side, share a neighbor on the other side
    Ar = H'*H       # Gives neighbors or right side

    if optimalflag
        @variable(m, x[1:n,1:n],Bin)
    else
        @variable(m, x[1:n,1:n])
        @constraint(m,x .<= ones(n,n))
        @constraint(m,x .>= zeros(n,n))
    end

    @objective(m, Min, sum((A[i,j])*x[i,j] for i=1:n-1 for j = i+1:n))

    # println("Symmetry constraints")
    for i = 1:n
        for j = 1:n
            @constraint(m, x[i,j] == x[j,i])
        end
    end
    # println("First solve")
    # One run of method
    JuMP.optimize!(m)

    while true
         D = JuMP.value.(x)

         # Store violating tuples in a vector
         violations = Vector{Tuple{Int,Int,Int}}()
         find_violations!(D,violations)

         # Iterate through and add constraints
         numvi = size(violations,1)
         # println("\t\t\t\t\tViolations = $numvi")
         for v in violations
             @constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
         end

         if numvi == 0
             break
         end
         JuMP.optimize!(m)
     end

     ## Output the answer
     X = JuMP.value.(x)
     LPval= JuMP.objective_value(m)

     LPbound = 0.0
     for i = 1:n
         for j = i+1:n
            aij = A[i,j]
            if aij < 0
                LPbound += aij*(X[j,i]-1)
            else
                LPbound += aij*X[j,i]
            end
         end
     end

    return X, LPbound
end


## Round vector X. There are far more efficient ways, but this won't
# be the bottleneck in the code anyways.
function RoundX(X::Matrix{Float64}, delta::Float64)
    Clusters = Vector{Vector{Int64}}()
    n = size(X,1)
    clustered = 0
    Unclustered = collect(1:n)
    label = 1
    Labels = zeros(Int64,n)
    while clustered < n
        pivot = rand(Unclustered)
        neighbs = findall(x->x<delta, X[:,pivot])
        cluster = intersect(neighbs,Unclustered)
        Unclustered = setdiff(Unclustered,cluster)
        push!(Clusters,cluster)
        clustered += length(cluster)
        Labels[cluster] .= label
        label += 1
    end

    C = spzeros(n,length(Clusters))
    for j = 1:length(Clusters)
        C[Clusters[j],j] .= 1
    end
    return C, Labels
end


## Round vector X. There are far more efficient ways, but this won't
# be the bottleneck in the code anyways.
function ManyRoundX(H::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64}, delta::Float64,mu::Float64,beta::Float64,k::Int64)

    best = Inf
    bestC = spzeros(1,1)
    bestL = Vector{Int64}()

    for i = 1:k
        C, Labels = RoundX(X,delta)
        obj = PBCC_Objective(H,C,mu,beta)
        if obj < best
            best = obj
            bestC = C
            bestL = Labels
        end
    end

    return bestC, bestL, best
end

## Round vector X with a bunch of delta values and take the best
function RangeRoundX(H::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64},mu::Float64,beta::Float64,k::Int64)

    best = Inf
    bestC = spzeros(1,1)
    bestL = Vector{Int64}()
    bestD = 0

    for delta = 0.05:0.05:0.95
        C, Labels, obj = ManyRoundX(H,X,delta,mu,beta,k)
        if obj < best
            best = obj
            bestC = C
            bestL = Labels
            bestD = delta
        end
    end

    return bestC, bestL, best, bestD
end


# Given a clustering of the matrix A = [0 H, H' 0], output a clustering
# of just H, i.e. a set of two vectors, x and y, so that x gives the cluster
# labels for nodes on one side of bipartite matrix H, and y gives the cluster
# labels for the other side
function  Biclustering(H::SparseMatrixCSC{Float64,Int64},L::Vector{Int64})

    l,r = size(H)
    x = L[1:l]
    y = L[l+1:l+r]

    return x,y
end

# Objective for a clustering, a super inefficient way to do things
function PBCC_Objective(H::SparseMatrixCSC{Float64,Int64},
    C::SparseMatrixCSC{Float64,Int64},mu::Float64,beta::Float64)

    n = size(C,1)

    l,r = size(H)
    Hbeta = Matrix(H)-beta*ones(l,r)
    A = [-mu*ones(l,l) Hbeta; Matrix(Hbeta') -mu*ones(r,r)]

    D = ones(n,n)  # dense, inefficient
    for i = 1:n
        cluster_ind = findall(x->x==1,C[i,:])   # cluster index
        cluster = findall(x->x==1,C[:,cluster_ind])
        D[i,cluster] .= 0.0
    end

    obj = 0.0
    for i = 1:n
        for j = i+1:n
           aij = A[i,j]
           if aij < 0
               obj += aij*(D[j,i]-1)
           else
               obj += aij*D[j,i]
           end
        end
    end
    return obj
end

"""
Returns the linear programming relaxation of the Bicluster Deletion Problem.
"""
function BiclusterDeletion_LP(H::SparseMatrixCSC,outputflag::Bool=true,optimalflag::Bool=false)

    m = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=outputflag, gurobi_env))

    l,r = size(H)           # number of nodes on each side
    n = l+r                 # total number of nodes in graph
    L = collect(1:l)        # Left indices
    R = collect(l+1:(n))    # Right indices, R[1] is the global index of 1st node on right side

    A = [spzeros(l,l) H; sparse(H') spzeros(r,r)]

    # Al = Project(H,1)       # Gives neighbors of left side: Al[i,j] = 1 if nodes i and j, on same side, share a neighbor on the other side
    # Ar = Project(H,2)       # Gives neighbors or right side

    Al = H*H'       # Gives neighbors of left side: Al[i,j] = 1 if nodes i and j, on same side, share a neighbor on the other side
    Ar = H'*H       # Gives neighbors or right side

    if optimalflag
        @variable(m, x[1:n,1:n],Bin)
    else
        @variable(m, x[1:n,1:n])
        @constraint(m,x .<= ones(n,n))
        @constraint(m,x .>= zeros(n,n))
    end

    @objective(m, Min, sum((A[i,j])*x[i,j] for i=1:n-1 for j = i+1:n))

    # Make the matrix symmetric
    println("First constraints")
    for i = 1:n
        for j = 1:n
            @constraint(m, x[i,j] == x[j,i])
        end
    end

    println("Second constraints")
    # Two nodes on opposite sides can't be together if they don't share an edge
    for i = 1:l
        for j = 1:r
            if H[i,j] == 0
                @constraint(m,x[i,R[j]] == 1)
            end
        end
    end

    # println("Third constraints")
    # Two nodes on the same side can't be together if they don't share a neighbor
    for i = 1:l-1
        for j = i+1:l
            if Al[i,j] == 0
                # No shared neighbor means can't be clustered together
                @constraint(m, x[i,j] == 1)
            end
        end
    end

    # println("Fourth constraints")
    # Do the same for the right side, careful to adjust indices using R
    for i = 1:r-1
        for j = i+1:r
            if Ar[i,j] == 0
                # No shared neighbor means can't be clustered together
                @constraint(m, x[R[i],R[j]] == 1)
            end
        end
    end

    # println("First solve")
    # One run of method
    JuMP.optimize!(m)

    while true
         D = JuMP.value.(x)

         # Store violating tuples in a vector
         violations = Vector{Tuple{Int,Int,Int}}()
         find_violations!(D,violations)

         # Iterate through and add constraints
         numvi = size(violations,1)
         # println("\t\t\t\t\tViolations = $numvi")
         for v in violations
             @constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
         end

         if numvi == 0
             break
         end
         JuMP.optimize!(m)
     end

     ## Output the answer
     X = JuMP.value.(x)
     LPval= JuMP.objective_value(m)

    return X, LPval
end


# find_violations
# Given a candidate distance matrix D, iterate through all 3-tuples of nodes
# and return the tuples where triangle inequality constraints have been violated.
#
# Output is stored in vector 'violations'
#
# Heuristic speed ups:
#       - Only grab Dij, Dik, Djk once per tuple
#       - if a - b - c > 0, then a>b and a>c. By checking these first
#           we rule out a lot of situations where constraints do not need to be
#           be added, so in all test cases this gave a speed up
#
# Note that we want D to be lower triangular here, though if it is symmetric
# this will work fine as well. We just need to make sure the distance information
# in D is not just stored in the upper triangular portion
function find_violations!(D::Matrix{Float64}, violations::Vector{Tuple{Int,Int,Int}},epsi::Float64=1e-8)
  n = size(D,1)

  @inbounds for i = 1:n-2
       for j = i+1:n-1
          a = D[j,i]
           for k = j+1:n
              b = D[k,i]
              c = D[k,j]
        if a - b > epsi && a - c > epsi && a-b-c > epsi
            push!(violations, (i,j,k))
            # @constraint(m, x[i,j] - x[i,k] - x[j,k] <= 0)
        end
        if b - a > epsi && b - c > epsi && b-a-c > epsi
            push!(violations, (i,k,j))
            # @constraint(m, x[i,k] - x[i,j] - x[j,k] <= 0)
        end

        if c - a > epsi && c-b>epsi && c-a-b > epsi
            push!(violations, (j,k,i))
            # @constraint(m, x[j,k] - x[i,k] - x[i,j] <= 0)
        end
      end
    end
  end
end
