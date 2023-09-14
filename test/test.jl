using EulerAngles
using LinearAlgebra
using Test

@testset "back and forth from random orthonomal matrix" begin
    for N = 1:20
        m = rand(N,N)
        m = (m + m') /2
        vals, vecs = eigen(m)
        angles = Angles(vecs)
        vecs_re = Matrix(angles)
        m_re = vecs_re * diagm(vals) * vecs_re' 
        @test isapprox(vecs, vecs_re; rtol=1e-8)
        @test isapprox(m_re, m; rtol=1e-10)
    end
end

@testset "back and forth from random angles" begin
    for N = 1:20
        angles = Angles([rand(i) for i = N:-1:1])
        mat = Matrix(angles)
        angles_re = Angles(mat)
        @test isapprox(angles_re.θs, angles.θs;rtol=1e-10)
    end
end

@testset "test Identiy matrix" begin
    for N = 1:20
        m = diagm(ones(5))
        @test Angles(m).flat == zeros(10)

        an = Angles(zeros(10))
        @test Matrix(an) == m
    end
end

@testset "test corner cases" begin
    for N = 1:20
        m = rand(N,N)
        m = (m + m') /2
        for i = 1:N
            m[:, i] = zeros(N)
            m[i, :] = zeros(N)
            m[i, i] = 0.8
            vals, vecs = eigen(m)
            agls = Angles(vecs)
            re_vecs = Matrix(agls)
            # @test re_vecs ≈ vecs
            re_m = Matrix(Eigen(vals, re_vecs))
            @test re_m ≈ m
        end
    end
end
