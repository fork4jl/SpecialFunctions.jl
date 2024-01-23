# SPDX-License-Identifier: BSD-3-Clause OR MIT

"""
    airy(z::Complex{Float64}, id::Int, kode::Int)

Compute Airy functions `Ai(z)` or its derivative `dAi(z)` for complex `z`.

On `kode=1`, airy computes the complex Airy function `Ai(z)` or
its derivative `dAi(z)/dz` on `id=0` or `id=1` respectively. 

On `kode=2`, a scaling option `exp(ζ) * Ai(z)` or 
`exp(ζ) * dAi(z)/dz` is provided to remove the exponential decay in
`-π/3 < arg(z) < π/3` and the exponential growth in
`π/3 < abs(arg(z)) < π` where `ζ=(2/3) * z * sqrt(z)`

While the Airy functions `Ai(z)` and `dAi(z)/dz` are analytic in
the whole `Z` plane, the corresponding scaled functions defined
for `kode=2` have a cut along the negative real axis.
Definitions and notation are found in the NBS Handbook of
Mathematical Functions (Abramowitz and Stegun, 1964).

# Arguments
## Input
- `z::Complex{Float64}`: Complex number to compute.
- `id::Int`: Order of Derivative.
    - `id=0`: Compute Airy functions, `Ai(z)`
    - `id=1`: Compute first derivative of Airy function, `dAi(z)/dz`
    - Other values: Return error
- `kode::Int`: Indicate the scaling option.
    - `kode=1`: No scaling.
    - `kode=2`: Apply scaling factor `exp(ζ)` to output,
        where `ζ = (2/3) * z * sqrt(z)`
    - Other values: Return error

## Output
Complex answer depending on the choices for id and kode:
- `kode=1`
    - `id=0`: return `Ai(z)`
    - `id=1`: return `dAi(z)/dz`
- `kode=2`, where `ζ = (2/3) * z * sqrt(z)`
    - `id=0`: return `exp(ζ) * Ai(z)`
    - `id=1`: return `exp(ζ) * dAi(z)/dz`


# Examples
```jldoctest
julia> AMOS.airy(0.0 + 0.0im, 0, 1)
0.3550280538878172

julia> AMOS.airy(0.0 + 0.0im, 1, 1)
-0.2588194037928068

julia> AMOS.airy(0.0 + 0.0im, 0, 2)
0.3550280538878172

julia> AMOS.airy(0.0 + 0.0im, 1, 2)
-0.2588194037928068
```


# Extended help

## AMOS Prologue
- Date Written: 1983/05/01 (YYYY/MM/DD)
- Revision Date: 1989/08/01 (YYYY/MM/DD)
- Category No.: B5k
- Keywords: Airy Function, Bessel Functions Of Order One Third
- Author: Amos, Donald E., Sandia National Laboratories

## Long Description
`Ai(z)` and `dAi(z)` are computed for `abs(z) > 1.0` from the K bessel
functions by

     Ai(z) =  C * sqrt(z) * K(1/3,zta),
    dAi(z) = -C * z * K(2/3,zta)
              C = 1.0 / (pi*sqrt(3.0))
            zta = (2/3) * z^(3/2)

with the power series for `abs(z) <= 1.0`.

In most complex variable computation, one must evaluate elementary
functions.  When the magnitude of `z` is large, losses
of significance by argument reduction occur.  Consequently, if
the magnitude of `zeta = (2/3) * z^1.5` exceeds `u1 = sqrt(0.5/ur)`,
then losses exceeding half precision are likely and an error
flag `ierr=3` is triggered where `ur = dmax1(d1mach(4),1.0d-18)` is
double precision unit roundoff limited to 18 digits precision.
Also, if the magnitude of `zeta` is larger than `u2 = 0.5/ur`, then
all significance is lost and `ierr=4`.  In order to use the int
function, zeta must be further restricted not to exceed the
largest integer, `u3=i1mach(9)`.  Thus, the magnitude of zeta
must be restricted by `min(u2,u3)`.  On 32 bit machines, `u1`, `u2`,
and `u3` are approximately `2.0e+3`, `4.2e+6`, `2.1e+9` in single
precision arithmetic and `1.3e+8`, `1.8e+16`, `2.1e+9` in double
precision arithmetic respectively.  This makes `u2` and `u3` limiting
in their respective arithmetics.  This means that the magnitude
of `z` cannot exceed `3.1e+4` in single and `2.1e+6` in
double precision arithmetic.  This also means that one can
expect to retain, in the worst cases on 32 bit machines,
no digits in single precision and only 7 digits in double
precision arithmetic.  Similar considerations hold for other
machines.

The approximate relative error in the magnitude of a complex
bessel function can be expressed by `p*10^s` where `p = max(unit
roundoff,1.0e-18)` is the nominal precision and `10^s` represents
the increase in error due to argument reduction in the
elementary functions.  Here, `s = max(1,abs(log10(abs(z)))`,
`abs(log10(fnu)))` approximately (i.e. `s=max(1,abs(exponent of
abs(z),abs(exponent of fnu))` ).  However, the phase angle may
have only absolute accuracy.  This is most likely to occur when
one component (in absolute value) is larger than the other by
several orders of magnitude.  If one component is `10^k` larger
than the other, then one can expect only `max(abs(log10(p))-k,0)`
significant digits; or, stated another way, when `k` exceeds
the exponent of `p`, no significant digits remain in the smaller
component.  However, the phase angle retains absolute accuracy
because, in complex arithmetic with precision `p`, the smaller
component will not (as a rule) decrease below `p` times the
magnitude of the larger component.  In these extreme cases,
the principal phase angle is on the order of `+p`, `-p`, `pi/2-p`,
or `-pi/2+p`.

## References
- OpenSpecfun (fortran): [`openspecfun/amos/zairy.f`](https://github.com/JuliaMath/openspecfun/blob/v0.5.6/amos/zairy.f)
- Scipy (C): [`scipy/scipy/special/_amos.c:amos_airy`](https://github.com/scipy/scipy/blob/b882f1b7ebe55e534f29a8d68a54e4ecd30aeb1a/scipy/special/_amos.c#L650-L994)
- Abramowits, M., & Stegun, I. A. (1964). Handbook of Mathematical Functions AMS55. National Bureau of Standards.
- Amos, D. E. (1983). Computation of Bessel functions of complex argument and large order (No. SAND-83-0643). Sandia National Lab., Albuquerque, NM (United States).
- Amos, D. E. (1985). Subroutine package for Bessel functions of a complex argument and nonnegative order (No. SAND-85-1018). Sandia National Labs., Albuquerque, NM (United States).
- Amos, D. E. (1986). Algorithm 644: A portable package for Bessel functions of a complex argument and nonnegative order. ACM Transactions on Mathematical Software (TOMS), 12(3), 265-273.

## Implementation
- Routines called: `acai`, `bknu`
- Constants used: `D1_MACH`, `I1_MACH`
"""
function airy(z::ComplexF64, id::Int, kode::Int)
    TTH = 2.0 / 3.0
    "1/(Gamma(2/3) * 3^(2/3))"
    GAMMA_C1 = 0.35502805388781723926
    "1/(Gamma(1/3) * 3^(1/3))"
    GAMMA_C2 = 0.25881940379280679841
    "1 / (sqrt(3) * PI)"
    COEF = 0.18377629847393068317

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
        """
        Power series for `abs(z) <= 1.0`
        """
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
    """
    Case for `abs(z) > 1.0`
    """
    fnu = (1.0 + fid) / 3.0
    """
    Set parameters related to machine constants.

    - `tol` is the approximate unit roundoff limited to `1.0e-18`.
    - `elim` is the approximate exponential over- and underflow limit.
            exp(-elim) < exp(-alim) = exp(-elim)/tol    and
            exp( elim) > exp( alim) = exp( elim)*tol      
        are intervals near underflow and overflow limits
        where scaled arithmetic is done.
    - `rl` is the lower boundary of the asymptotic expansion for large `z`.
    - `dig` = number of base 10 digits in `tol = 10^(-dig)`.
    """
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

    """
    Test for proper range
    """
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

    """
    `real(zta) <= 0` when `real(z) < 0`, especially when `imag(z)` is small
    """
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
            """
            Overflow test
            """
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
        """
        `bknu` and `acai` return `exp(zta) * K(fnu,zta)` on `kode=2`
        """
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
            """
            Overflow test
            """
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
