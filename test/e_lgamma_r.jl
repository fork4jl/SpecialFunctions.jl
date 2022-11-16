# From openlibm/test/libm-test-ulps.h, openlibm/test/libm-test.c

# lgamma_test block
# for (T, lgamma) in ((Float64, _lgamma_r), (Float32, _lgammaf_r))
for T in (Float64, Float32)
    @testset "lgamma_test, $T" begin
        @test logabsgamma(T(Inf))[1] === T(Inf)
        @test logabsgamma(T(0))[1] === T(Inf)
        @test logabsgamma(T(NaN))[1] === T(NaN)

        @test logabsgamma(T(-3))[1] === T(Inf)
        @test logabsgamma(T(-Inf))[1] === T(Inf)

        # logabsgamma(1) == 0, lgamma (1) sets signgam to 1
        y, signgam = logabsgamma(T(1))
        @test y === T(0.0)
        @test signgam == 1

        # logabsgamma(3) == log(2), lgamma (3) sets signgam to 1
        y, signgam = logabsgamma(T(3))
        @test y === log(T(2.0))
        @test signgam == 1

        # logabsgamma(0.5) == log(sqrt(pi)), logabsgamma(0.5) sets signgam to 1
        y, signgam = logabsgamma(T(0.5))
        @test y === T(0.5log(π))
        @test signgam == 1

        # logabsgamma(-0.5) == log(2sqrt(pi)), logabsgamma(-0.5) sets signgam to -1
        y, signgam = logabsgamma(T(-0.5))
        @test y === T(0.5log(4π))
        @test signgam == -1

        # In the two "broken" tests, an exact match not possible, even
        # in Float64, thus, we check for as close a tolerance as
        # possible.

        # logabsgamma(0.7) == 0.26086724653166651439, logabsgamma(0.7) sets signgam to 1
        y, signgam = logabsgamma(T(0.7))
        # @test_broken y === 0.26086724653166651439
        if T === Float64
            @test y ≈ 0.26086724653166651439 atol=6e-17
        else
            @test y ≈ 0.26086724653166651439 atol=3e-8
        end
        @test signgam == 1

        # logabsgamma(1.2) == -0.853740900033158497197e-1, logabsgamma(1.2) sets signgam to 1
        y, signgam = logabsgamma(T(1.2))
        # @test_broken y === -0.853740900033158497197e-1
        if T === Float64
            @test y ≈ -0.853740900033158497197e-1 atol=2e-17
        else
            @test y ≈ -0.853740900033158497197e-1 atol=2e-8
        end
        @test signgam == 1
    end
end

# Comparison against the "gold standard" (in context of SpecialFunctions.jl)
using OpenLibm_jll

function openlibm_logabsgamma(x::Float64)
    signp = Ref{Int32}()
    y = ccall((:lgamma_r,libopenlibm),  Float64, (Float64, Ptr{Int32}), x, signp)
    return y, Int(signp[])
end
function openlibm_logabsgamma(x::Float32)
    signp = Ref{Int32}()
    y = ccall((:lgammaf_r,libopenlibm),  Float32, (Float32, Ptr{Int32}), x, signp)
    return y, Int(signp[])
end

meetstol(x, atol) = isapprox(openlibm_logabsgamma(x)[1], logabsgamma(x)[1], atol=atol)

@testset "logabsgamma validation against OpenLibm, Float64" begin
    @test all(x -> meetstol(x, 1e-13), -50:1e-4:50)
    @test all(x -> meetstol(x, 1e-12), 50:1e-4:500)
    @test all(x -> meetstol(x, 1e-11), 500:1e-3:1000)
    @test all(x -> meetstol(x, 1e-10), 1000:1e-1:10000)
end

@testset "logabsgamma validation against OpenLibm, Float32" begin
    @test all(x -> meetstol(x, 1f-5), -0.0f0:1f-4:25.0f0)
    @test all(x -> meetstol(x, 1f-4), -25.0f0:1f-4:0.0f0)
    @test all(x -> meetstol(x, 1f-4), 25.0f0:1f-4:100.0f0)
    @test all(x -> meetstol(x, 1f-3), 100.0f0:1f-4:1000.0f0)
    @test all(x -> meetstol(x, 1f-1), 1000.0f0:1f-1:10000.0f0)
end
