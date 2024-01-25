
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
