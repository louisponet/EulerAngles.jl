module EulerAngles
using LinearAlgebra
export Angles

# mutable struct Angles{T, N, VT<:AbstractVector{T}}
#     θs::VT
#     r::T
#     n::NTuple{N, Int}
# end
#
# function Angles(v::AbstractVector)
#     n = norm(v)
#     return Angles(angles(v), n, (length(v),))
# end

struct Angles{T, VVT<:Vector{Vector{T}}}
    θs::VVT
    flat::Vector{T}
end

function toVVA(θs::Vector{Vector{T}}) where T
    n = length(θs)
    flat = Vector{T}(undef, Int(n*(n+1)//2))
    head = 1
    for (i, len) = zip(1:n, n:-1:1)
        flat[head:head+len-1] = θs[i]
        head += len
    end
    Angles(θs, flat)
end

function Angles(flat::Vector{T}) where T
    L = size(flat, 1)
    n = Int((sqrt(8*L+1) - 1) * 0.5)
    θs::Vector{Vector{T}} = []
    head = 1
    for i =n:-1:1
        push!(θs, flat[head:head+i-1])
        head += i
    end
    Angles(θs, flat)
end

# function angles_t(t::Matrix{T}) where T
#     ν = size(t, 1)
#     agls = zeros(T, ν)
#     # last angle is always π/2
#     agls[ν] = π/2
#     agls[ν-1] = atan(t[ν-1, ν], t[ν, ν])
#     if t[ν, ν] / sin(agls[ν-1]) ≈ 1
#         @info "k = ν-1:"  agls
#         return agls
#     end
#     for k = ν-2:-1:1
#         if t[k, ν] ≈ 0.
#             agls[k] = 0 
#             continue
#         elseif abs(t[k, ν]) ≈ 1. 
#             agls[k] = sign(t[k, ν]) * π/2 
#             return agls
#         end
#         agls[k] = atan(t[k, ν] / t[k+1, ν] * sin(agls[k+1]))
#     end
#     return agls
# end

function angles_t(t::Matrix{T}) where T
    ν = size(t, 1)
    tν = t[:, end]
    agls = zeros(T, ν)
    # last angle is always π/2
    agls[ν] = π/2
    agls[1] = asin(tν[1])
    Σcos::T = 1.
    for k = 2:ν-1
        Σcos *= cos(agls[k-1])
        agls[k] = asin(tν[k] / Σcos)
        if agls[k] ≈ π/2 || cos(agls[k]) ≈ 0.
            return agls
        end
    end
    agls[ν-1] = atan(t[ν-1, ν], t[ν, ν])
    return agls
end

function matrix_a(agls::Vector{T}) where T
    n = size(agls, 1)
    ma = diagm(ones(T, n))
    for i = 1:n-1
        ma[i, i] = cos(agls[i])
        ma[i, n] = tan(agls[i])
        for j = 1:i
            ma[i, n] *= cos(agls[j])
        end
    end

    for i = 1:n
        for k = 1:n-1
            if i > k
                ma[i, k] = -tan(agls[i]) * tan(agls[k])
                for l = k:i
                    ma[i, k] *= cos(agls[l])
                end
            end
        end
    end
    
    ma[n, n] = tan(agls[n])
    for j = 1:n
        ma[n, n] *= cos(agls[j])
    end
    return ma
end
    
function Angles(t::Matrix{T}) where T
    @assert [norm(t[:,i]) for i in 1:size(t, 2)] ≈ ones(size(t, 1))
    n = size(t, 1)
    res::Vector{Vector{T}} = []
    mat_t = deepcopy(t)
    for _ = 2:n
        agls = angles_t(mat_t)
        a = matrix_a(agls) 
        tmp = similar(mat_t)
        mul!(tmp, a', mat_t) 
        mat_t = tmp[1:end-1, 1:end-1]
        # remove the last angle π/2
        push!(res, agls[1:end-1])
    end
    return Angles(res)
end
    
function Base.Matrix(angles::Angles{T}) where T
    agls = angles.θs
    n = size(agls, 1) + 1
    b_mat = diagm(ones(T, n))
    for i = 1:n-1
        size_a = size(agls[i], 1) + 1
        tmp_agl = deepcopy(agls[i])
        push!(tmp_agl, π/2)
        padded = Matrix{T}(I, n, n)
        padded[1:size_a,1:size_a] = matrix_a(tmp_agl)
        b_mat = b_mat * padded
    end
    return b_mat
end

end # module EulerAngles

