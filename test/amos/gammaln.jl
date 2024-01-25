
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
        # 1e+300 ~ 1e-16
        [ 10.0^i for i in 300:-1:-16 ]...,
    ]

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
            if !isapprox(ref, res)
                @warn y
            end

            @test ref == res
            @test ref_jl ≈ res
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
    test_wontfix_in_rewrite = [
        # err: last step: log(zp)
        (1.7140765668835524e-32, 7.3143848485170963e+01),
        (2.4794046815125870e-02, 3.6833397905797782e+00),

        # err: tlg = log(zdmy);
        (1.1466853192189908e-01, 2.1097745674575785e+00),
        (2.0595890907476233e-01, 1.4930045975739470e+00),
        (3.0508527867637447e-01, 1.0781433866410493e+00),
        (9.0418421415131678e-01, 6.3234254316328384e-02),
        (1.0030054521742426e+00, -1.7273757975118590e-03),

        # err: tlg = log(z);
        (1.1275345484474153e+02, 4.1858665390175162e+02),
    ]
    for (y, mat_ref) in test_wontfix_in_rewrite
        # Won't fix
        @test AMOS._gammaln(y) !== AMOS.gammaln(y)
        @test AMOS._gammaln(y) ≈ AMOS.gammaln(y)
        # TODO: fix after rewrite
        @test_broken mat_ref === AMOS.gammaln(y)
        @test mat_ref ≈ AMOS.gammaln(y)
    end

    """
    The current impl matches the output of the fortran reference impl.

    However, it does not matches the output of Matlab.
    """
    test_not_match_mat_todo = [
        (1.4012984643248171e-45, 1.0327892990343184e+02),
        # Max Subnormal: 0x1.fffffcp-127
        # Confirm with Matlab: `gammaln(1.1754942106924411e-38) == 87.3365448697624`
        #   and wolframcloud: `LogGamma[1.1754942106924411e-38] == 87.3365448697624`
        (1.1754942106924411e-38, 8.7336544869762406e+01),
        (1.0000438400824574e-16, 3.6841317648783225e+01),
        (1.0000327253233265e-08, 1.8420648013392153e+01),
        (1.0000020646814845e-04, 9.2102805938354670e+00),
        (1.0000000000000002e-02, 4.5994798780420210e+00),
        (1.0000019272127880e+00, -1.1124143560825666e-06),
        (1.0000181598531222e+02, 3.5914255921250549e+02),
        (1.0000244355122401e+04, 8.2101968081058963e+04),
    ]
    for (y, mat_ref) in test_not_match_mat_todo
        # Match Fortran Ref impl
        @test AMOS._gammaln(y) === AMOS.gammaln(y)
        # NOT match Matlab
        # TODO: fix after rewrite
        @test_broken mat_ref === AMOS.gammaln(y)
        @test mat_ref ≈ AMOS.gammaln(y)
    end

    """
    The AMOS impl is not precise; current `AMOS.gammaln` impl is precise.
    Keep these tests for now and remove them when the rewrite is complete.
    """
    test_remove_after_rewrite = [
        # Reference results are generated using Matlab R2023b
        #   fprintf('(%.16e, %.16e),\n', f, gammaln(f))
        #
        # ( input_y, reference_output ),
        (1.0938769749613498e-16, 3.6751633244559891e+01),
        (1.3528203464811494e-08, 1.8118489177381669e+01),
        (1.0443000379825174e-04, 9.1669332614365313e+00),
        (1.0278570979004344e-02, 4.5718475469812878e+00),
        (2.5334589554522413e+00, 3.0848277234626065e-01),
        (1.0796531772254237e+02, 3.9608659648662461e+02),
        (1.2208401549613756e+04, 1.0266740169949362e+05),
        (1.2179734219771166e+08, 2.1458096221985698e+09),
        (1.1586095324839718e+16, 4.1696714276010867e+17),
        (1.0107989710858223e+32, 7.3478478670872164e+33),
        (1.3254736271879463e+64, 1.9437701704564903e+66),
        (1.9336249168162124e+128, 5.6924039680253297e+130),
        (4.5221368056396525e+256, 2.6679285261833809e+259),
        (3.4425072585673266e+300, 2.3788128892517655e+303),

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
    for (y, mat_ref) in test_remove_after_rewrite
        @test mat_ref === AMOS.gammaln(y)
        # Won't fix: AMOS impl is not precise
        @test mat_ref !== AMOS._gammaln(y)
        @test mat_ref ≈ AMOS._gammaln(y)
    end

end
