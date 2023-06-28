module EulerAngles
using LinearAlgebra
export Angles

mutable struct Angles{T, N, VT<:AbstractVector{T}}
    θs::VT
    r::T
    n::NTuple{N, Int}
end

function Angles(v::AbstractVector)
    n = norm(v)
    return Angles(angles(v), n, (length(v),))
end
# function angles(v::AbstractVector)
#     v = normalize(v)
#     if length(v) < 2
#         throw(ArgumentError("Length of the vector needs to be larger than 2"))
#     end
#     θs = similar(v, length(v) - 1)
#     angles!(θs, v)
#     return θs
# end
function angles(v::AbstractVector)
    v = normalize(v)
    if length(v) < 2
        throw(ArgumentError("Length of the vector needs to be larger than 2"))
    end
    n = length(v)
    angles = Float64[]
    for i = 1:n-1
        v[1] = clamp(v[1], -1, 1)

        θ = asin(v[1])
        if length(v) == 2 && v[2] < 0.0
            if θ > 0
                θ = π - θ
            elseif θ < 0
                θ = - π - θ
            end
        end

        push!(angles, θ)
        v = v[2:end]./cos(θ)
    end
    return Angles(angles, 1.0, (n,n))
end
function angles!(θs, v)
    n = length(v)
   
    for i = 1:n - length(θs)-1
        push!(θs, 0)
    end
    @inbounds begin    
        θs[n-1] = atan(v[end-1], v[end])
        for i = n-2:-1:1
            t = sin(θs[i+1])
            if t == 0
                θs[i] = 0
            else
                θs[i] = atan(v[i], v[i+1]/t)
            end
        end
    end
end

function angles!(e::Angles, v::AbstractVector)
    r = norm(v)
    angles!(e.θs, v./r)
    e.r = r
    e.n = (length(v), )
end

function angles!(e::Angles, m::AbstractMatrix)
    r = tr(m)
    angles!(e.θs, m)
    e.r = r
    e.n = size(m)
end

function Angles(m::AbstractMatrix)
    n = size(m, 1)
    e = Angles(zeros(div(n * (n - 1),2)), 1.0, (n,n))
    matcache = similar(m)
    curid = 1
    t = copy(m)
    mulcache = similar(m)
    ret = Float64[]
    for i in 2:n
        # Lets take vectors starting from the last column on ortho_matrix
        an = m[:, end]
        # print('Vector: \n%s' % an)
        as = angles(an)
        # print('Angles: \n%s' % angles)
        a = rotmat(as)
        # print('Matrix a: \n%s' % a)
        m = transpose(m)* a
        # print('New matrix b: \n%s' % b)
        m = m[1:end-1, 1:end-1]
        # print('New matrix b (fixed): \n%s' % b)
        append!(ret, as.θs[1:end-1])
    end
    
    return Angles(ret, 1.0, (n,n))
end

function rotmat(a::Angles{T, 2}) where T
    push!(a.θs, π/2)
    n = a.n[1]
    out =diagm(0 => [T(1) for i = 1:n])
    for i in 1:n-1
        out[i,i] = cos(a.θs[i])
        out[i,n] = tan(a.θs[i])
        for j = 1:i+1
            out[i, n] *= cos(a.θs[j])
        end
    end

    for i in 1:n
        for k in 1:n-1
            if i > k
                out[i,k] = -tan(a.θs[i]) * tan(a.θs[k])
                for l in k+1:i
                    out[i, k] *= cos(a.θs[l])
                end
            end
        end
    end

    out[n, n] = tan(a.θs[n])
    for j in 1:n
        out[n, n] *= cos(a.θs[j])
    end
    return out
    
end

    
# function Angles(m::AbstractMatrix)
#     n = size(m, 1)
#     e = Angles(zeros(div(n * (n - 1),2)), 1.0, (n,n))
#     matcache = similar(m)
#     curid = 1
#     t = copy(m)
#     mulcache = similar(m)
#     @inbounds for i = n:-1:2
#         finalid = curid+i-2
#         te = Angles(view(e.θs, curid:finalid), 1.0, (n,))
#         @show e.θs
#         angles!(te.θs, view(t, 1:i, i))
#         mul!(mulcache, transpose(t), fill!(matcache, te))
#         curid = finalid+1
#         t, mulcache = mulcache, t
#     end
#     return e
# end

function Base.Vector(e::Angles{T}) where {T}
    out = ones(T, e.n)
    @inbounds for i = 1:e.n - 1
        out[i] *= sin(e.θs[i])
        for j = i+1:n
            out[j] *= cos(e.θs[i])
        end
    end
    return out * e.r
end
    
@inline function Base.fill!(out::Matrix, e::Angles{T,1}) where {T}
    @inbounds begin
        n = length(e.θs)+1
        ang = e.θs
        for i = 1:n-1
            out[i, i] = cos(ang[i])
        end
        # Padding for if matrix is bigger than required
        out[n, n] = 0.0
        for i = n+1:e.n[1]
            for j = 1:e.n[1]
                out[j, i] = 0.0
            end
        end
        for j = 1:e.n[1]
            for i = n+1:e.n[1]
                out[i, j] = i == j ? 1.0 : 0.0
            end
        end
        # Region 3
        for k=1:n-1
            s = - sin(ang[k])
            for i = k+1:n
                t = i == n ? s : sin(ang[i]) * s
                for l in k+1:i-1
                    t *= out[l, l]
                end
                out[i, k] = t
            end
        end
        
        # Region 2
        out[1, n] = 1.0
        out[2, n] = out[1, 1]
        for i = 3:n
            out[i, n] = out[i-1, n] * out[i-1, i - 1]
        end
        for i in 1:n-1
            out[i, n] *= sin(ang[i])
        end
        # Region 4
        for i in 1:n-1
            for j in i+1:n-1
                out[i, j] = 0.0
            end
        end

        return out
    end
end

function Base.Matrix(e::Angles)
    out = similar(e.θs, e.n[1], e.n[1])
    return fill!(out, e)
end

function Base.Matrix(e::Angles{T, 2}) where {T}
    n = e.n[1]
    norm = e.r/n
    out = Matrix(diagm(0 => ones(T,n)))
    matcache = similar(out)
    mulcache = copy(out)
    ang = e.θs
    curid = length(ang)
    @inbounds for i = 2:n
        finalid = curid+i-2
        te = Angles(view(e.θs, curid:finalid), 1.0, (n,))
        mul!(mulcache, fill!(matcache, te), transpose(out))
        out, mulcache = mulcache, out
        curid -= i
    end
    return out
end


end # module EulerAngles
