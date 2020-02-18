# Form clique-reduction graph
function ExtractTriMat(A::SparseMatrixCSC{Float64,Int64})
    @assert(issymmetric(A))
    n = size(A,1)
    Am = zeros(n,n)
    for i = 1:n-2
        for j = i+1:n-1
            Aji = A[j,i]
            if Aji == 0
                continue
            end

            for k = j+1:n

                if A[k,j]+A[k,i] == 2
                    Am[j,i]+=1/2
                    Am[k,i]+=1/2
                    Am[k,j]+=1/2
                end
            end
        end
    end
    Am = sparse(Am + Am')
    return Am
end

# Form clique-reduction graph
function ExtractAllTriMat(A::SparseMatrixCSC{Float64,Int64})
    n = size(A,1)
    Am = zeros(n,n)
    for i = 1:n-2
        for j = i+1:n-1
            Aji = A[j,i]
            Aij = A[i,j]
            if Aji == 0 || Aij == 0
                continue
            end

            for k = j+1:n

                if A[k,j]+A[k,i]+A[j,k]+A[i,k] == 4
                    Am[j,i]+=1/2
                    Am[k,i]+=1/2
                    Am[k,j]+=1/2
                end
            end
        end
    end
    Am = sparse(Am + Am')
    return Am
end

function ExtractFanMat(A::SparseMatrixCSC{Float64,Int64})
    n = size(A,1)
    Am = zeros(n,n)
    for i = 1:n
        for j = i+1:n
            for k = j+1:n
                for l = k+1:n
                    # Yuck, quad loop.
                    ps= [i j k l; i k j l; i l j k; j k i l; j l i k; k l i j]
                    for w = 1:size(ps,1)
                        a = ps[w,1]; b = ps[w,2]; c = ps[w,3]; d = ps[w,4]
                        if A[a,c] > 0 && A[b,c] > 0 && A[a,d] > 0 && A[b,d] > 0 && A[c,a] ==0 && A[d,a] ==0 && A[c,b] ==0 && A[d,b] ==0
                            Am[a,b] += 1/2
                            Am[b,a] += 1/2
                            Am[c,d] += 1/2
                            Am[d,c] += 1/2
                        end
                    end
                end
            end
        end
    end
    Am = sparse(Am+Am')
    return Am
end


function turntounweighted(Aorig::SparseMatrixCSC{Float64,Int64})
    A = copy(Aorig)
    n = size(A,1)
    for i = 1:n
        for j= 1:n
            if A[i,j] > 0
                A[i,j] = 1
                A[j,i] = 1
            end
        end
    end
    return A
end
