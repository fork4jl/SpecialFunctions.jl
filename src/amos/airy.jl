# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""
    airy(z::ComplexF64, id::Int, kode::Int)

COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z.

ON KODE=1, ZAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
-PI/3.LT.ARG(Z).LT.PI/3 AND THE EXPONENTIAL GROWTH IN
PI/3.LT.ABS(ARG(Z)).LT.PI WHERE ZTA=(2/3)*Z*CSQRT(Z).

WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
MATHEMATICAL FUNCTIONS (REF. 1).


# Arguments
- `z` Complex Number
- `id`: ORDER OF DERIVATIVE
    id=0, or id=1
- `kode`: INDICATE THE SCALING OPTION
    kode=1, 
        AI=AI(Z)                ON ID=0 OR
        AI=DAI(Z)/DZ            ON ID=1
    kode=2, 
        AI=CEXP(ZTA)*AI(Z)      ON ID=0 OR
        AI=CEXP(ZTA)*DAI(Z)/DZ  ON ID=1 WHERE
        ZTA=(2/3)*Z*CSQRT(Z)

# Return
COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND KODE


# Examples



# Extended help

## AMOS Prologue
* BEGIN PROLOGUE  ZAIRY
* DATE WRITTEN   830501   (YYMMDD)
* REVISION DATE  890801   (YYMMDD)
* CATEGORY NO.  B5K
* KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
* AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES

## Long Description
AI AND DAI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE K BESSEL
FUNCTIONS BY

AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
                C=1.0/(PI*SQRT(3.0))
                ZTA=(2/3)*Z**(3/2)

WITH THE POWER SERIES FOR CABS(Z).LE.1.0.

IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
FLAG IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
MACHINES.

THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
OR -PI/2+P.

## References
- [`openspecfun/amos/zairy.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zairy.f)
- [`scipy/scipy/special/_amos.c:amos_airy`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L650-L994)
- HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
    AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
    COMMERCE, 1955.
- COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
- A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
    1018, MAY, 1985
- A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
    MATH. SOFTWARE, 1986

## Implementation
ROUTINES CALLED  ZACAI,ZBKNU,AZEXP,AZSQRT,I1MACH,D1MACH
"""
function airy(z::ComplexF64, id::Int, kode::Int)
    # aa, ad, ak, alim, atrm, az, az3, bk, ck, dig, dk, d1, d2, elim, fid, fnu, rl, r1m5, sfac, tol,  bb, alaz = 0.0
    TTH = 2.0 / 3.0
    GAMMA_C1 = 0.35502805388781723926  # 1/(Gamma(2/3) * 3**(2/3))
    GAMMA_C2 = 0.25881940379280679841  # 1/(Gamma(1/3) * 3**(1/3))
    COEF = 0.18377629847393068317  # 1 / (sqrt(3) * PI)

    cy = [ 0.0 + 0.0im ]
    zta = 0.0 + 0.0im
    csq = 0.0 + 0.0im
    
    # TODO: return or throw
    ierr = 0
    nz = 0
    ai = 0.0 + 0.0im
    if (id < 0) || (id > 1)
        ierr = 1
    end
    if (kode < 1) || (kode > 2)
        ierr = 1
    end
    if ierr != 0
        # TODO: throw(AmosException(ierr)
        return 0.0 + 0.0im
    end

    az = abs(z)
    tol = max(D1_MACH[4], 1e-18)
    fid = float(id)

    if az <= 1.0
        #
        # POWER SERIES FOR ABS(Z) <= 1.
        #
        s1 = 1.0 + 0.0im
        s2 = 1.0 + 0.0im
        if az < tol
            #= 170 =#
            aa = 1.0e3 * D1_MACH[1]
            s1 = 0.0 + 0.0im
            if id != 1
                if az > aa
                    s1 = GAMMA_C2 * z
                end
                #= 180 =#
                ai = GAMMA_C1 - s1
                # @info "180: ai" nz ierr
                # @show csq s1 ai
                return ai
            end
        
            #= 190 =#
            ai = complex(-GAMMA_C2)
            aa = sqrt(aa)
            if az > aa
                s1 = z * z * 0.5
            end
            #= 200 =#
            ai += s1 * GAMMA_C1
            # @info "200: ai" nz ierr
            # @show csq s1 ai
            return ai
        end
    
        aa = az * az
        if aa >= (tol / az)
            trm1 = 1.0 + 0.0im
            trm2 = 1.0 + 0.0im
            atrm = 1.0
            z3 = z * z * z
            az3 = az * aa
            ak = 2.0 + fid
            bk = 3.0 - fid - fid
            ck = 4.0 - fid
            dk = 3.0 + fid + fid
            d1 = ak * dk
            d2 = bk * ck
            ad = min(d1, d2)
            ak = 24.0 + 9.0 * fid
            bk = 30.0 - 9.0 * fid
            for k = 1:25
                trm1 *= z3 / d1
                s1 += trm1
                trm2 *= z3 / d2
                s2 += trm2
                atrm *= az3 / ad
                d1 += ak
                d2 += bk
                ad = (d1 > d2 ? d2 : d1)
                if atrm < tol * ad
                    break
                end
                ak += 18.0
                bk += 18.0
            end
        end #= 30 =#
        #= 40 =#
        if id != 1
            ai = s1*GAMMA_C1 - z*s2*GAMMA_C2
            if kode == 1
                return ai
            end
        
            zta = z * sqrt(z) * TTH
            ai *= exp(zta)
            # @info "40: ai" nz ierr
            # @show zta
            return ai
        end
        #= 50 =#
        ai = -s2 * GAMMA_C2
        if az > tol
            ai += z * z * s1 * GAMMA_C1 / (1.0 + fid)
        end
        if kode == 1
            return ai
        end
    
        zta = z * sqrt(z) * TTH
        # @info "50: ai * exp(zta)" nz ierr
        # @show csq s1 ai
        return ai * exp(zta)
    end

    #= 70 =#
    #
    # CASE FOR CABS(Z).GT.1.0
    #
    fnu = (1.0 + fid) / 3.0
    #
    # SET PARAMETERS RELATED TO MACHINE CONSTANTS.
    # TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    # ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
    # EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
    # EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
    # UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
    # RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
    # DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
    #
    k1 = trunc(Int, I1_MACH[15])
    k2 = trunc(Int, I1_MACH[16])
    r1m5 = D1_MACH[5]
    k = min(abs(k1), abs(k2))
    elim = 2.303 * (k * r1m5 - 3.0)
    k1 = trunc(Int, I1_MACH[14]) - 1
    aa = r1m5 * k1
    dig = min(aa, 18.0)
    aa *= 2.303
    alim = elim + max(-aa, -41.45)
    rl = 1.2 * dig + 3.0
    alaz = log(az)

    #
    # TEST FOR PROPER RANGE
    #
    aa = 0.5 / tol
    bb = I1_MACH[9] * 0.5
    aa = min(aa, bb)
    aa = aa^(TTH)
    if az > aa
        #= 260 =#
        ierr = 4
        nz = 0
        # TODO: throw(AmosException(ierr)
        # @info "260: 0.0im" nz ierr
        # @show az aa
        return 0.0 + 0.0im
    end

    aa = sqrt(aa)
    if az > aa
        ierr = 3
    end
    csq = sqrt(z)
    zta = z * csq * TTH

    #
    # RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
    #
    iflag = 0
    sfac = 1.0
    ak = imag(zta)
    if real(z) < 0.0
        bk = real(zta)
        ck = -abs(bk)
        zta = complex(ck, ak)
    end
    #= 80 =#
    if (imag(z) == 0.0) && (real(z) <= 0.0)
        zta = complex(0.0, ak)
    end
    #= 90 =#
    aa = real(zta)
    if (aa < 0.0) || (real(z) <= 0.0)
        if kode != 2
            #
            # OVERFLOW TEST
            #
            if aa <= -alim
                aa = -aa + 0.25 * alaz
                iflag = 1
                sfac = tol
                if aa > elim
                    #= 270 =#
                    nz = 0
                    ierr = 2
                    # TODO: throw(AmosException(ierr)
                    # @info "90: ai" nz ierr iflag
                    # @show aa sfac
                    return ai
                end
            end
        end

        #= 100 =#
        #
        # CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
        #
        mr = 1
        if imag(z) < 0.0
            mr = -1
        end
        nn = acai!(cy, zta, fnu, kode, mr, 1, rl, tol, elim, alim)
        @assert length(cy) == 1
        if nn < 0
            #= 280 =#
            if nn == -1
                nz = 1
                # @info "280: 0.0" nz ierr iflag
                return 0.0
            else
                nz = 0
                ierr = 5
                # TODO: throw(AmosException(ierr)
                # @info "208-else: 0.0" nz ierr iflag
                return 0.0
            end
        end
    
        nz += nn
    else
        #= 110 =#
        if kode != 2
            #
            # OVERFLOW TEST
            #
            if aa >= alim
                aa = -aa - 0.25 * alaz
                iflag = 2
                sfac = 1.0 / tol
                if aa < -elim
                    nz = 1
                    # @info "110: 0.0" nz ierr iflag
                    return 0.0
                end
            end
        end
        #= 120 =#
        nz = bknu!(cy, zta, fnu, kode, 1, tol, elim, alim)
        @assert length(cy) == 1
    end

    #= 130 =#
    s1 = cy[1] * COEF

    # 0: normal;  3: underflow
    @assert ierr==0 || ierr==3
    if iflag == 0
        if id != 1
            # @info "130: csq * s1" nz ierr iflag
            # @show csq s1
            return csq * s1
        end
        #= 140 =#
        # @info "140: csq * s1" nz ierr iflag
        # @show csq s1
        return (-z * s1)
    end

    #= 150 =#
    s1 *= sfac
    if id != 1
        s1 *= csq
        # @info "150: s1 / sfac" nz ierr iflag
        # @show csq s1
        return (s1 / sfac)
    end
    
    #= 160 =#
    s1 *= -z
    # @info "160: s1 / sfac" nz ierr iflag
    # @show csq s1
    return (s1 / sfac)
end
