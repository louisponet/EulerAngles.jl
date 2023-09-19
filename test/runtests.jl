using EulerAngles
using LinearAlgebra
using Test

@testset "Back and forth from random orthonomal matrix" begin
    for n_dim = 1:10
        for _ in 1:100
            # construct real symmetric matrix,
            # so the eigenvalues are real
            m = rand(n_dim,n_dim)
            m = (m + m') /2
            vals, vecs = eigen(m)
            angles = Angles(vecs)
            vecs_re = Matrix(angles)
            m_re = vecs_re * diagm(vals) * vecs_re' 
            # @test isapprox(vecs, vecs_re; rtol=1e-8)
            @test isapprox(m_re, m; rtol=1e-10)
        end
    end
end

@testset "Back and forth from random angles" begin
    for n_dim = 1:10
        for _ in 1:100
            angles = Angles(n_dim)
            mat = Matrix(angles)
            angles_re = Angles(mat)
            @test isapprox(angles_re.θs, angles.θs;rtol=1e-4)
        end
    end
end

@testset "Convert Identiy matrix" begin
    for n_dim = 1:20
        n_angles::Int = n_dim*(n_dim-1)/2
        # Matrix to Angles
        m = diagm(ones(n_dim))
        @test Angles(m).θs == zeros(n_angles)

        # Angles to Matrix
        an = Angles(zeros(n_angles))
        @test Matrix(an) == m
    end
end

@testset "Unit eigen vector" begin
    for n_dim = 1:20
        # construct random symmetry matrix
        m = rand(n_dim,n_dim)
        m = (m + m') /2
        for i = 1:n_dim
            m[:, i] = zeros(n_dim)
            m[i, :] = zeros(n_dim)
            m[i, i] = 0.5
            vals, vecs = eigen(m)
            agls = Angles(vecs)
            re_vecs = Matrix(agls)
            # @test re_vecs ≈ vecs
            re_m = Matrix(Eigen(vals, re_vecs))
            @test re_m ≈ m
        end
    end
end
