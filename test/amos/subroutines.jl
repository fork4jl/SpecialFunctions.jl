
"Special Float32"
const SPECIAL_FLOAT32 = Float64[
    0.0,        # +0.0
    0x1p-149,   # Min Subnormal Float
    0x1p-127,   # mid Subnormal Float
    0x1.fffffcp-127,  # Max Subnormal Float
    0x1p-126,   # Min Normalized Float
    0x1p+1,     # 2
    0x1.fffffep+127,  # Max Normalized Float
    Inf,        # Inf
    reinterpret(Float32, 0x7f800001),  # NaN (+min)
    NaN,        # NaN (+min)
    reinterpret(Float32, 0x7fffffff),  # NaN (+max)

    # eps( 10^0 ~ 10^16 )
    [ eps(10.0^i) for i in 0:16 ]...,
] # SPECIAL_FLOAT32
# TOOD: add Float64 const


"""
input:
    zarr[ (+) ]

output:
    zarr[
        (+), (-),
        inv.(+), inv.(-),
    ]
"""
function gen_neg_inv(zarr::Vector{Float64})
    vcat(
        zarr, -1.0*zarr,
        inv.(zarr), -1.0*inv.(zarr),
    )
end

"""
input:
    zarr[ (+,+) ]

output:
    zarr[
        (-,+),  (+,+),
        (-,-),  (+,-),
    ]
"""
function gen_phase4(zarr::Vector{ComplexF64})
    vcat(
        1.0im * zarr,
        zarr
        -1.0 * zarr,
        -1.0im * zarr,
    )
end
function gen_phase4(zvv::Vector{Vector{ComplexF64}})
    gen_phase4.(zvv)
end

"""test if ComplexF64 contains NaN or Inf in real or imag parts.
"""
function contains_inf_nan(z::Vector{ComplexF64})
    z |> real .|> isinf |> any ||
        z |> imag .|> isinf |> any ||
        z |> real .|> isnan |> any ||
        z |> imag .|> isnan |> any
end


# ==== Test Set ====#

@testset "AMOS.uchk" begin
    special_inputs = [
        #= ( yr,yi, ascle, tol ) =#
        # min > ascle
        ( 2.0,3.0, 1e-6, 1.0),
        # min <= ascle && max < (min/tol)
        ( 1e-8,1e-9, 1e-6, 1e-3),
        # min <= ascle && max >= (min/tol)
        ( 1e-8,1e-13, 1e-6, 1e-3),
    ]
    for (yr,yi,ascle,tol) in special_inputs
        for y in [
            complex(-yr, yi), complex( yr, yi),
            complex(-yr,-yi), complex( yr,-yi),
        ]
            @test AMOS.uchk(y,ascle,tol) == AMOS._uchk(y,ascle,tol)
        end
    end

    test_y = [
        SPECIAL_FLOAT32...,
        [ 10.0^i for i in 0:16 ]...,
    ]
    test_ascle = [
        [ 10.0^i for i in 0:-1:-16 ]...,
    ]
    test_tol = [
        0.0,  # test divide by zero
        [ 10.0^i for i in 6:-1:-6 ]...,
    ]

    for tol in test_tol,
        ascle in test_ascle,
        yr in test_y
        for y in [
            complex(yr, 2.),
            complex(2., yr),
            complex(yr, yr),
        ]
            @test AMOS.uchk(y,ascle,tol) == AMOS._uchk(y,ascle,tol)
        end
    end
end


"Warp `gammalog`, return for `z <= 0.0`"
function gammalog(z::Float64)
    if z > 0.0
        return loggamma(z)
    end

    return NaN
end

@testset "AMOS.gammaln" begin
    special_inputs = Float64[
        #= if z > 0.0 =#
        #  z <= 0.0
        -1.0, 0.0,
    
        #  z > 0.0
        #= if z <= 101.0 =#
        #  z in Int(0, 100]
        1:100...,
        #  z in Float64(100, 101]
        100 + rand(),

        #= if z < zmin,  zmin=7.0 =#
        #  z < zmin
        3.14,
        #  z >= zmin
        7.0, 7.1,

        #= if zp >= wdtol,  wdtol=2.2e-16 =#
        1e-16, 1e-10,
        
        #= if 0.0 == zinc =#
        #  when z > zmin(=7.0), zinc=0
        8.1,
        π * 100,
        7 + rand()*10.0^rand(0:300),  # (7, 1e300]
        #  when z in (0.0, 7.0], zinc != 0 
        rand() * 7,
    ]
    for z in special_inputs
        @test AMOS.gammaln(z) === AMOS._gammaln(z)
    end

    test_y = [
        SPECIAL_FLOAT32...,
        # [0, 1)
        # rand(Float64, 1000_00)...,
        # 1e+300 ~ 1e-16
        [ 10.0^i for i in 300:-1:-16 ]...,
    ]

    broken = []
    for y in (test_y)
        ref_jl = gammalog(y)
        ref = AMOS._gammaln(y)  # call AMOS, baseline
        res = AMOS.gammaln(y)   # fortran translation

        if isnan(ref)
            @test ref === res
            if isinf(ref_jl)
                # broken on Inf
                @test_broken ref_jl == res
            else
                @test ref_jl === res
            end
        else
            @test ref ≈ res
            @test ref == res
            if ref != res
                @info "ref != res" y
                push!(broken, y)
            end
            
            @test ref_jl ≈ res
        end
    end
    
    println("broken = $broken")
    
    for y in broken
        if AMOS.gammaln(y) === loggamma(y)
            println("amos_impl_bad_case1_ignore +=  $y,")
        end
    end

    """
    The output of these tests is inconsistent on all three implementations.

    The inconsistency between the julia impl (`AMOS.gammaln`) and the fortran 
        reference impl (`AMOS._gammaln`) has been analyzed and determined to be a flaw 
        in the libm library used by the fortran reference impl.
    Specifically, this is due to inaccuracies in the output of the log function in libm.
        Therefore these tests will not be fixed by this rewrite.

    The inconsistency between the julia implementation and the loggamma output will be 
        investigated further after the rewrite is complete.
    """
    test_broken_wontfix = [
        # err: last step: log(zp)
        0.02479404681512587,
        # err: tlg = log(zdmy);
        0.11466853192189908,
        # err: tlg = log(zdmy);
        0.20595890907476233,
        # err: tlg = log(zdmy);
        0.30508527867637447,
        # err: tlg = log(zdmy);
        0.9041842141513168
    ]
    for y in test_broken_wontfix
        # Won't fix
        @test AMOS._gammaln(y) !== AMOS.gammaln(y)
        @test AMOS._gammaln(y) ≈ AMOS._gammaln(y)
        # TODO: fix after rewrite
        @test_broken loggamma(y) === AMOS._gammaln(y)
        @test loggamma(y) ≈ AMOS._gammaln(y)
    end

    """
    The AMOS implementation is incorrect, 
        currently AMOS.gammaln is implemented correctly.
    Keep these tests for now and remove them when the rewrite is complete.
    """
    amos_impl_bad_case1_ignore = [
        # Reference results are generated using Matlab R2023b
        #   fprintf('(%.16e, %.16e),\n', f, gammaln(f))
        (4.8005616413904217e-03, 5.3362703098007724e+00),
        (6.5484884510764729e-03, 5.0247762942629715e+00),
        (9.5307120639248621e-03, 4.6478089333560009e+00),
        (1.0437445211128327e-02, 4.5564199285322493e+00),
        (1.5161032651903383e-02, 4.1804632661918237e+00),
        (1.7681503535653342e-02, 4.0252850865635805e+00),
        (3.9183257716550957e-02, 3.2181278030939726e+00),
        (4.5328177330525743e-02, 3.0693159532467242e+00),
        (5.6114036963397118e-02, 2.8505009301958908e+00),
        (6.2753362300222215e-02, 2.7354647445616491e+00),
        (8.7850721600321258e-02, 2.3874984673865742e+00),
        (9.3047948621315602e-02, 2.3277486080640775e+00),
        (1.0785155408548408e-01, 2.1738439502350180e+00),
        (1.1878681593647111e-01, 2.0728423113377392e+00),
        (1.9463793225911241e-01, 1.5528082232364799e+00),
        (2.4022730231778189e-01, 1.3301780924913191e+00),
        (2.5274438262508425e-01, 1.2764850947623523e+00),
        (2.9137968135096370e-01, 1.1264540568414718e+00),
        (3.6482220971878199e-01, 8.9153332152762865e-01),
        (4.9546148733747319e-01, 5.8132744574887640e-01),
        (6.1087936625438222e-01, 3.8168597849277025e-01),
        (6.4073094472553438e-01, 3.3839219410920057e-01),
    ]
    for (y, mat_ref) in amos_impl_bad_case1_ignore
        @test mat_ref === AMOS.gammaln(y)
        @test mat_ref === loggamma(y)
        # Won't fix
        @test mat_ref !== AMOS._gammaln(y)
        @test mat_ref ≈ AMOS._gammaln(y)
    end

    # This may be a problem with the original fortran implementation
    max_subnormal = 0x1.fffffcp-127
    @test AMOS.gammaln(max_subnormal) == AMOS._gammaln(max_subnormal)
    # Confirm in Matlab: `gammaln(1.1754942106924411e-38) == 87.3365448697624`
    #   and wolframcloud: `LogGamma[1.1754942106924411e-38] == 87.3365448697624`
    @test loggamma(max_subnormal) == 87.3365448697624
    @test_broken AMOS.gammaln(max_subnormal) == loggamma(max_subnormal)
end

@testset "AMOS.kscl!" begin
    function tkscl!(
        y::Vector{ComplexF64},
        zr::ComplexF64,
        fnu::Float64=2.0,
        n::Int=2,
        rz::ComplexF64=0.5+0.5*im,
        ascle::Float64=1e-30,
        tol::Float64=1e-15,
        elim::Float64=20.,
    )
        # @info "input" y zr fnu n rz ascle tol elim
        y_ref = copy(y)
        y_res = copy(y)
        ref = AMOS._kscl!(y_ref,zr,fnu,n,rz,ascle,tol,elim)
        res = AMOS.kscl!(y_res,zr,fnu,n,rz,ascle,tol,elim)

        @test ref == res
        if ref != res
            @info "ref != res" ref res
            @info "params" zr fnu n rz ascle elim
            @show y y_ref y_res
        end

        if contains_inf_nan(y_ref)
            # use === to test NaN/Inf
            @test all(y_ref .=== y_res)
        else
            @test isapprox(y_ref, y_res)
            # @info "isapprox"  ref zr y y_ref y_res
            # @show y_ref
            # @show y_res
            if !isapprox(y_ref, y_res)
                @info "!isapprox"  ref zr y y_ref y_res
                @info "params" zr fnu n rz ascle elim
                @show y y_ref y_res
            end
        end
    end

    test_y = [
        # ComplexF64[],  # TODO: n==0
        # n=1
        ComplexF64[0.0],
        ComplexF64[1.0],
        ComplexF64[pi],
        [ rand(ComplexF64, 1) for _ in 1:10 ]...,
        # n=2
        ComplexF64[0.0, 0.0],
        ComplexF64[1., 1.],
        ComplexF64[1., 0.],
        ComplexF64[0., 1.],
        ComplexF64[1., 2.],
        [ rand(ComplexF64, 2) for _ in 1:10 ]...,
        # n=5
        ComplexF64[ complex(i,i) for i in 1:5 ],
    ]
    test_zr = [
        SPECIAL_FLOAT32...,
        # [0, 1)
        rand(Float64, 10)...,
        # 1e-0 ~ 1e-16
        [ 10.0^i for i in 0:-1:-16 ]...,
        1.0,
        pi, ℯ,
    ]

    for y in gen_phase4(test_y),
        zr in gen_neg_inv(test_zr)
        # TODO: @test_broken
        isnan(zr) && continue
        isinf(zr) && continue
        0x1.fffffep+127==zr && continue

        tkscl!(y, complex(zr), 2.0, length(y))
    end
end

"test if one of input is NaN"
contains_nan(ref, impl) = isnan(ref) || isnan(impl)
contains_nan(ref::Ref{T}, impl::Ref{T}) where {T} =
    contains_nan(ref[], impl[])
contains_nan(ref::Vector{T}, impl::Vector{T}) where {T} =
    any(contains_nan.(ref, impl))

naninf(f) = isnan(f) || isinf(f)

@testset "AMOS.s1s2!" begin
    function test_s1s2(
        zr::ComplexF64,
        s1::ComplexF64,
        s2::ComplexF64,
        ascle::Float64=1e-30,
        alim::Float64=1e-6,
        iuf::Int=0
    )
        s1_ref, s1_impl = Ref(s1), Ref(s1)
        s2_ref, s2_impl = Ref(s2), Ref(s2)
        
        nz_ref, iuf_ref = AMOS._s1s2!(zr, s1_ref, s2_ref, ascle, alim, iuf)
        nz_impl, iuf_impl = AMOS.s1s2!(zr, s1_impl, s2_impl, ascle, alim, iuf)
        
        @test nz_ref == nz_impl
        @test iuf_ref == iuf_impl
        if contains_nan(s1_ref, s1_impl)
            @test s1_ref[] === s1_impl[]
        else
            @test s1_ref[] ≈ s1_impl[]
        end
        if contains_nan(s2_ref, s2_impl)
            @test s2_ref[] === s2_impl[]
        else
            @test s2_ref[] ≈ s2_impl[]
        end
    end

    test_zr = [
        SPECIAL_FLOAT32...,
    ]
    s1s2 = [
        SPECIAL_FLOAT32...,
    ] |> complex
    lim_scale = [
        # 1e-0 ~ 1e-8
        [ 10.0^i for i in 0:-1:-8 ]...,
    ]
    iuf_test = [ -1, 0, 1, 2 ]
    
    for zr in test_zr,
        s1 in s1s2,
        s2 in reverse(s1s2),
        ascle in lim_scale,
        alim in lim_scale,
        iuf in iuf_test
        # TODO
        naninf(zr) && continue
        naninf(s1) && continue
        naninf(s2) && continue

        test_s1s2(complex(zr), s1, s2, ascle, alim, iuf)
    end
end

@testset "AMOS.asyi!" begin
    function test_asyi(
        y::Vector{ComplexF64},
        z::ComplexF64, 
        fnu::Float64, 
        kode::Int, 
        n::Int, 
        rl::Float64, 
        tol::Float64,
        elim::Float64, 
        alim::Float64
    )
        y_ref, y_impl = copy(y), copy(y)

        nz_ref = AMOS._asyi!(y_ref, z, fnu, kode, n, rl, tol, elim, alim)
        nz_impl = AMOS.asyi!(y_impl, z, fnu, kode, n, rl, tol, elim, alim)
        
        @test nz_ref == nz_impl
        if contains_nan(y_ref, y_impl)
            cmp = all(y_ref .=== y_impl)
            @test cmp
            if !cmp
                @info "contains_nan" y_ref y_impl
                @show y z fnu kode n rl tol elim alim
            end
        else
            @test y_ref ≈ y_impl
        end
    end

    test_y = [
        ComplexF64[],  # TODO: n==0
        # n=1
        ComplexF64[0.0],
        ComplexF64[1.0],
        ComplexF64[pi],
        [ rand(ComplexF64, 1) for _ in 1:10 ]...,
        # n=2
        ComplexF64[0.0, 0.0],
        ComplexF64[1., 1.],
        ComplexF64[1., 0.],
        ComplexF64[0., 1.],
        ComplexF64[1., 2.],
        [ rand(ComplexF64, 2) for _ in 1:10 ]...,
        # n=5
        ComplexF64[ complex(i,i) for i in 1:5 ],
    ]
    test_z = [
        SPECIAL_FLOAT32...,
    ] |> complex

    for y in test_y,
        z in test_z,
        fnu in [ 0., 1. ],
        kode in [ 1, 2 ]
        naninf(z) && continue
        z==0.0 && continue  # div by zero

        test_asyi(y, z, fnu, kode, length(y), 3.1, 1e-6, eps(), 2*eps())
    end
end


@testset "AMOS.mlri!" begin
    function test_mlri(
        y::Vector{ComplexF64},
        z::ComplexF64, 
        fnu::Float64,
        kode::Int, 
        n::Int, 
        tol::Float64
    )
        y_ref, y_impl = copy(y), copy(y)

        nz_ref = AMOS._mlri!(y_ref, z, fnu, kode, n, tol)
        nz_impl = AMOS.mlri!(y_impl, z, fnu, kode, n, tol)
        
        @test nz_ref == nz_impl
        if contains_nan(y_ref, y_impl)
            cmp = all(y_ref .=== y_impl)
            @test cmp
            if !cmp
                @info "contains_nan" y_ref y_impl
                @show y z fnu kode n tol
            end
        else
            @test y_ref ≈ y_impl
        end
    end

    test_y = [
        # ComplexF64[],  # TODO: n==0
        # n=1
        ComplexF64[0.0],
        ComplexF64[1.0],
        ComplexF64[pi],
        [ rand(ComplexF64, 1) for _ in 1:10 ]...,
        # n=2
        ComplexF64[0.0, 0.0],
        ComplexF64[1., 1.],
        ComplexF64[1., 0.],
        ComplexF64[0., 1.],
        ComplexF64[1., 2.],
        [ rand(ComplexF64, 2) for _ in 1:10 ]...,
        # n=5
        ComplexF64[ complex(i,i) for i in 1:5 ],
    ]
    test_z = [
        SPECIAL_FLOAT32...,
    ] |> complex

    for y in test_y,
        z in test_z,
        fnu in [ 0., 1. ],
        kode in [ 1, 2 ]
        naninf(z) && continue
        z==0.0 && continue  # div by zero
        abs(z)>typemax(Int) && continue  # InexactError: Int64

        test_mlri(y, z, fnu, kode, length(y), eps())
    end
end

@testset "AMOS.seri!" begin
    function test_seri(
        y::Vector{ComplexF64},
        z::ComplexF64, 
        fnu::Float64,
        kode::Int, 
        n::Int, 
        tol::Float64,
        elim::Float64, 
        alim::Float64
    )
        y_ref, y_impl = copy(y), copy(y)

        nz_ref = AMOS._seri!(y_ref, z, fnu, kode, n, tol, elim, alim)
        nz_impl = AMOS.seri!(y_impl, z, fnu, kode, n, tol, elim, alim)
        
        @test nz_ref == nz_impl
        if contains_nan(y_ref, y_impl)
            cmp = all(y_ref .=== y_impl)
            @test cmp
            if !cmp
                @info "contains_nan" y_ref y_impl
                @show y z fnu kode n tol elim alim
            end
        else
            @test y_ref ≈ y_impl
        end
    end

    test_y = [
        # n=1
        ComplexF64[0.0],
        ComplexF64[1.0],
        ComplexF64[pi],
        [ rand(ComplexF64, 1) for _ in 1:10 ]...,
        # n=2
        ComplexF64[0.0, 0.0],
        ComplexF64[1., 1.],
        ComplexF64[1., 0.],
        ComplexF64[0., 1.],
        ComplexF64[1., 2.],
        [ rand(ComplexF64, 2) for _ in 1:10 ]...,
        # n=5
        ComplexF64[ complex(i,i) for i in 1:5 ],
    ]
    test_z = [
        SPECIAL_FLOAT32...,
    ] |> complex

    for y in test_y,
        z in test_z,
        fnu in [ 0., 1. ],
        kode in [ 1, 2 ]
        naninf(z) && continue
        z==0.0 && continue  # div by zero
        abs(z)>typemax(Int) && continue  # InexactError: Int64

        test_seri(y, z, fnu, kode, length(y), 1e-6, eps(), 2*eps())
    end
end

# TODO: improve code coverage
@testset "AMOS.bknu!" begin
    function test_bknu(
        y::Vector{ComplexF64},
        z::ComplexF64, 
        fnu::Float64,
        kode::Int, 
        n::Int, 
        tol::Float64,
        elim::Float64, 
        alim::Float64
    )
        y_ref, y_impl = copy(y), copy(y)

        nz_ref = AMOS._bknu!(y_ref, z, fnu, kode, n, tol, elim, alim)
        nz_impl = AMOS.bknu!(y_impl, z, fnu, kode, n, tol, elim, alim)
        
        @test nz_ref == nz_impl
        if contains_nan(y_ref, y_impl)
            cmp = all(y_ref .=== y_impl)
            @test cmp
            if !cmp
                @info "contains_nan" y_ref y_impl
                @show y z fnu kode n tol elim alim
            end
        else
            @test y_ref ≈ y_impl
        end
    end

    test_y = [
        # n=1
        ComplexF64[0.0],
        ComplexF64[1.0],
        ComplexF64[pi],
        [ rand(ComplexF64, 1) for _ in 1:10 ]...,
        # n=2
        ComplexF64[0.0, 0.0],
        ComplexF64[1., 1.],
        ComplexF64[1., 0.],
        ComplexF64[0., 1.],
        ComplexF64[1., 2.],
        [ rand(ComplexF64, 2) for _ in 1:10 ]...,
        # n=5
        ComplexF64[ complex(i,i) for i in 1:5 ],
    ]
    test_z = [
        SPECIAL_FLOAT32...,
    ] |> complex

    for y in test_y,
        z in test_z,
        fnu in [ 0., 1. ],
        kode in [ 1, 2 ]
        naninf(z) && continue
        z==0.0 && continue  # div by zero
        abs(z)>typemax(Int) && continue  # InexactError: Int64

        test_bknu(y, z, fnu, kode, length(y), 1e-6, eps(), 2*eps())
    end
end


@testset "AMOS.acai!" begin
    function test_acai(
        y::Vector{ComplexF64},
        z::ComplexF64, 
        fnu::Float64, 
        kode::Int,
        mr::Int, 
        n::Int, 
        rl::Float64, 
        tol::Float64,
        elim::Float64, 
        alim::Float64
    )
        y_ref, y_impl = copy(y), copy(y)

        nz_ref = AMOS._acai!(y_ref, z, fnu, kode, mr, n, rl, tol, elim, alim)
        nz_impl = AMOS.acai!(y_impl, z, fnu, kode, mr, n, rl, tol, elim, alim)
        
        if contains_nan(y_ref, y_impl)
            cmp = all(y_ref .=== y_impl)
            if !cmp
                # @info "contains_nan" y_ref y_impl
                # @show y z fnu kode mr n rl tol elim alim
                @test_broken cmp
            else
                @test nz_ref == nz_impl
                @test cmp
            end
        else
            @test nz_ref == nz_impl
            @test y_ref ≈ y_impl
        end
    end

    test_y = [
        # n=1
        ComplexF64[0.0],
        ComplexF64[1.0],
        ComplexF64[pi],
        [ rand(ComplexF64, 1) for _ in 1:10 ]...,
        # n=2
        ComplexF64[0.0, 0.0],
        ComplexF64[1., 1.],
        ComplexF64[1., 0.],
        ComplexF64[0., 1.],
        ComplexF64[1., 2.],
        [ rand(ComplexF64, 2) for _ in 1:10 ]...,
        # n=5
        ComplexF64[ complex(i,i) for i in 1:5 ],
    ]
    test_z = [
        SPECIAL_FLOAT32...,
    ] |> complex

    for y in test_y,
        z in test_z,
        fnu in [ 0., 1. ],
        kode in [ 1, 2 ]
        naninf(z) && continue
        z==0.0 && continue  # div by zero

        test_acai(y, z, fnu, kode, 1, length(y), 3.1, 1e-6, eps(), 2*eps())
    end
end


@testset "AMOS.airy" begin
    function test_airy(
        z::ComplexF64,
        id::Int,
        kode::Int,
    )
        y_ref = SpecialFunctions._airy(z, Int32(id), Int32(kode))
        y_impl = AMOS.airy(z, id, kode)

        if contains_nan(y_ref, y_impl)
            cmp = all(y_ref .=== y_impl)
            @test cmp
            if !cmp
                @info "contains_nan" y_ref y_impl
                @show z id kode
            end
        else
            @test y_ref ≈ y_impl
            if !isapprox(y_ref, y_impl)
                @info "!isapprox" y_ref y_impl
                @show z id kode
            end
        end
    end

    test_z = Float64[
        SPECIAL_FLOAT32...,
        [ rand() for _ in 1:10 ]...,
        pi,
        [ 10.0^i for i in 1:6 ]...,
    ] |> complex    
    test_z = [
        test_z...,
        [ rand(ComplexF64) for _ in 1:10 ]...,
        
        0.1 + 0.2im,
        3.14 + 2.7im,
        2.7 + 3.14im
    ]

    # TODO: gen_phase4
    for z in (test_z),
        id in [ 0, 1 ],
        kode in [ 1, 2 ]
        isnan(z) && continue  # Int64(NaN)
        abs(z)>1.3e6 && continue  # too large arg

        test_airy(z, id, kode)
    end
end
