module joukowskiDynamics

cd("/Users/blake/Dropbox/CMU/julia/swim")
include("generateStructs.jl")

using LinearAlgebra, PolynomialRoots, DifferentialEquations
using .generateStructs
using Debugger

export joukowskiDyn


function joukowskiDyn(dx,x,p,t)
    vp, fp, sp = p;
    # du = zeros(sp.N,1);

    # define the shape change and time-variations of the shape variables
    rc = fp.rc; # radius of the circle in zeta plane. Considered to be 1
    ua = fp.ua; # amplitude of flapping
    uf = fp.uf; # frequency of flapping
    zcy = ua*sin(2*pi*uf*t);
    zcydot = ua*(2*pi*uf)*cos(2*pi*uf*t);


    # ak is a constant that defines the shape of the foil with constant area.
    ak = 0.85^2/(rc^2 - 0.15^2); # variable 'k' Equation 8.2, Page# 85 (Chapter of Model Reduction) from Hailong's Thesis document
    zcx = (1-ak)/(1+ak)*sqrt(rc^2 - zcy^2); # This is equivalent to Zeta_x from Thesis document; % variable 'Zeta_x' -- Equation 8.2, Page# 85 (Chapter of Model Reduction) from Hailong's Thesis document
    a = -2*ak/(1+ak)*sqrt(rc^2 - zcy^2); # variable 'a' Equation 8.2, Page# 85 (Chapter of Model Reduction) from Hailong's Thesis document

    # zcx is expressed in terms of zcy (Zeta_y)-- Single degree of freedom
    zcxdot = (1-ak)*(-zcy*zcydot)/(1+ak)/(zcx-a); # Zeta_xdot
    adot = -2*ak/(1+ak)*(-zcy*zcydot)/(zcx-a);

    fp.x = x[1]; fp.y = x[2]; fp.theta = x[3];

    # Update location of vortices
    for j = 1:vp.maxNumOfVortices
        vp.mvortex[j,1] = vp.mvortex0[j,1] + x[2+j*2];
        vp.mvortex[j,2] = vp.mvortex0[j,2] + x[3+j*2];
    end
    # indices = 1:vp.maxNumOfVortices;
    # vp.mvortex[:,1] = vp.mvortex0[:,1] + x[2 .+ 2*indices];
    # vp.mvortex[:,2] = vp.mvortex0[:,2] + x[3 .+ 2*indices];

    (velocZeta,vf) = computeVelocities(t,zcx,zcy,rc,a,zcxdot,zcydot,adot,ak,fp,vp,sp);

    dx[1] = fp.U*cos(x[3]) - fp.V*sin(x[3]); # Page # 19 from Dr. Xiong thesis document
    dx[2] = fp.U*sin(x[3]) + fp.V*cos(x[3]); # Page # 19 from Dr. Xiong thesis document
    dx[3] = fp.Omega; # Page # 19 from Dr. Xiong thesis document

    for j = 1:vf
        dx[2+j*2] = real(velocZeta[j]); # It starts with 4 as 1,2 & 3 are used for location of foil-fixed frame
        dx[3+j*2] = imag(velocZeta[j]);
    end
    # vfInd = 1:vf;
    # dx[2 .+ 2*vfInd] = real.(velocZeta[vfInd]);
    # dx[3 .+ 2*vfInd] = imag.(velocZeta[vfInd]);

    # dx[sp.N] = sqrt(fp.U^2 + fp.V^2);
    # print("Working...")
end

function computeVelocities(t,zcx,zcy,rc,a,zcxdot,zcydot,adot,ak,fp,vp,sp)

    # Compute how many vortices are in the flow at time t
    # This loop is probably unnecessary. What is a good way to identify which j is zero? Shouldn't be difficult...
    vp.vortexFlag = 0;
    # vfIndCheck = 1:vp.maxNumOfVortices;

    for j = 1:vp.maxNumOfVortices
        if vp.mvortex[j,3] != 0.0
            vp.vortexFlag = j; # vortexFlag is the number of vortices at time t
        end
    end

    # Compute simplified inertia matrix
    zc = zcx + im*zcy;
    conjzc = zcx - im*zcy;
    zcdot = zcxdot + im*zcydot;

    NewI11 = 2*pi*(rc^2 - a^2); # Elements of Simplified Inertial Matrix, refer page 86/87 from thesis document (approximation of Momenta)
    NewI22 = 2*pi*(rc^2 + a^2);
    NewI13 = real(2*pi*im*(rc^2*zc + a^2*conjzc - a^2*ak*zc));
    NewI23 = imag(2*pi*im*(rc^2*zc + a^2*conjzc - a^2*ak*zc));
    NewI33 = pi*rc^4*(1-ak^4)/2 + pi*a^2*(1+ak^2)*(zc^2+conjzc^2)+2*pi*(rc^2-ak*a^2)*zc*conjzc + 2*pi*a^4;
    NewI33 = abs(NewI33);

    Imatrix = [NewI11 0 NewI13; 0 NewI22 NewI23; NewI13 NewI23 NewI33]; #Simplified Inertia matrix for foil-fluid system Page # 87, Equation # 8.10 (thesis document)

    #= Compute parameters defining the complex potential due to changes in the
    foil's shape with constant area (page 85, 86 in Hailong thesis)=#
    us1 = -rc^2; us2 = a^2;
    vs1 = -im*rc^2; vs2 = -im*a^2;
    omegas1 = -im*rc^2*(zc + a^2/zc);
    omegas2  = -im*a^2*(-1)*(ak*zc + a^2/ak/zc);

    NewA1 = -rc^2*(zcdot + 2*a*adot/zc - a^2*zcdot/zc^2); # A1, page #86
    NewA2 = (rc^2*zcdot*(a^2/zc^2 - ak^2) + conj(zcdot)*(a^2 - ak^2*zc^2) - 2*a*adot*ak*zc);  # A2, page #86
    NewA3 = rc^2*(2*a*adot/zc^2 - 2*a^2*zcdot/zc^3);  # A3, page #86

    # Cast inerital position into body-fixed frame
    # Rmatrix converts vectors in spatially-fixed frame to vectors in foil-fixed frame
    Rmatrix = [cos(fp.theta) sin(fp.theta); -sin(fp.theta) cos(fp.theta)];

    Lb = Rmatrix*[fp.x, fp.y];
    xbody = Lb[1]; ybody = Lb[2];

    deltaTForNewVortex = 0.005;

    # If it is good time to shed a new vortex
    # if (abs(round(t/deltaTForNewVortex) - t/deltaTForNewVortex) < 1e-10) && t != vp.timeOfVortexSheddingEvent
    if (t - vp.timeOfVortexSheddingEvent > deltaTForNewVortex && t != vp.timeOfVortexSheddingEvent && vp.vortexFlag < vp.maxNumOfVortices)
        # shedNewVortex()
        z0 = (a-zc); # This is on Zeta Plane
        zeta0 = (a-zc);
        te = zeta0 + zc + a^2/(zeta0+zc); # location of the trailing edge
        beta = atan(zcy,(zcx - a)); # Double check which atan representation to use
        Zstart = 1.2*zeta0 + zc + a^2/(1.2*zeta0+zc);

        if  vp.vortexFlag == sp.npv # npv = number of initially already present vortices
            oldpv = Zstart;
        else
            oldpvzeta = vp.mvortex[vp.vortexFlag, 1] + im*vp.mvortex[vp.vortexFlag, 2];
            oldpv = oldpvzeta + zc + (a^2)/(oldpvzeta+zc);
        end

        # Hailong's thesis page #43, Modeled from Mason / Streitlein method for vortex shedding
        dto = oldpv - te;
        mason = (dto*dto)/abs(dto*dto)*exp(-im*4*beta);

        if (zcy < 0 && angle(dto) < pi/2+2*beta && angle(dto) > 0)
            newpv = 2*a + (oldpv - 2*a)/(1+mason^(1/3)*exp(-im*2*pi/3) + mason^(2/3)*exp(-im*4*pi/3));
            # println("First case")
        elseif (zcy < 0 && angle(dto) > -( pi/2-2*beta) && angle(dto) < 0)
            newpv = 2*a + (oldpv - 2*a)/(1+mason^(1/3)*exp(im*2*pi/3) + mason^(2/3)*exp(im*4*pi/3));
            # println("Second case")
        elseif (zcy > 0 && angle(dto) < pi/2+2*beta && angle(dto) > 0)
            newpv = 2*a + (oldpv - 2*a)/(1+mason^(1/3)*exp(-im*2*pi/3) + mason^(2/3)*exp(-im*4*pi/3));
            # println("Third case")
        elseif (zcy > 0 && angle(dto) > -( pi/2-2*beta) && angle(dto) < 0)
            newpv = 2*a + (oldpv - 2*a)/(1+mason^(1/3)*exp(im*2*pi/3) + mason^(2/3)*exp(im*4*pi/3));
            # println("Fourth case")
        else
            newpv = 2*a + (oldpv - 2*a)/(1+mason^(1/3) + mason^(2/3));
            # println("Final case")
        end

        #= This is finding the roots of the equation BELOW equation 4.6 (page #43 in Hailong thesis).
        Note that this is the conformal mapping between z-plane and ζ-plane. =#
        A = 1; B = -newpv; C = a^2;
        posRoots = (-B + sqrt(B^2 - 4*A*C))/(2*A) - zc;
        negRoots = (-B - sqrt(B^2 - 4*A*C))/(2*A) - zc;
        pzeta = [posRoots, negRoots];
        testCond = broadcast(abs, pzeta);
        newvortex = pzeta[testCond .> rc][end];
        z1 = newvortex;

        # % Complex potentials due to translation and rotation
        # %         z0= Zeta
        # %         zc= Zeta_c
        # %         zcx= Zeta_x
        # %         zcy= Zeta_y
        # %         us1 = -rc^2
        # %         us2 = a^2
        # %         vs1 = -i*rc^2
        # %         vs2 = -i*a^2
        # % Page # 86

        # Solve for Udot, Vdot, Omegadot, and Gammadot (page 95 in Hailong thesis)
        # to update U, V, and Omega according to newly-shed vortex
        dw1dz = -us1/z0^2 - us2 /(z0+zc)^2; # w1
        dw2dz = -vs1/z0^2 - vs2 /(z0+zc)^2; # w2
        dw3dz = -omegas1/z0^2 - omegas2 /(z0+zc)^2; # w3
        dwsdz = -NewA1/z0^2 - NewA2/(z0+zc)^2 + NewA3*(1/(z0+zc) - 1/z0); # ws(zeta)
        dwdz = fp.U*dw1dz + fp.V*dw2dz + fp.Omega*dw3dz + dwsdz; # Equation 8.22, Page # 91
        dwdz = dwdz + im*fp.gammac/z0;

        # This loop is probably not necessary (currently working on a more vectorized version of this code)
        if vp.vortexFlag != 0
            for pn = 1:vp.vortexFlag
                strengthj = vp.mvortex[pn,3];
                zj = vp.mvortex[pn,1] + im*vp.mvortex[pn,2];
                dwdz = dwdz + strengthj*im*(1/(z0-zj) - 1/(z0 - rc^2/conj(zj))); # Page # 48
            end
        end

        C = zeros(4,4) .+ 0*im;
        C[1,1] = dw1dz;
        C[1,2] = dw2dz;
        C[1,3] = dw3dz;
        C[1,4] = im*(1/(z0-z1) - 1/(z0 - rc^2/conj(z1)));
        C[2:4,1:3] = Imatrix;

        # Location of the k-th vortex?
        zk = z1;
        Zk =  zk + zc + (a^2)/(zk+zc);
        Zkx = real(Zk);
        Zky = imag(Zk);
        (impulse, impulseCouple) = computeImpulseAndImpulseCoupleForVortex(zk, rc, zc, a);

        B = [real(impulse), imag(impulse), impulseCouple];
        Bsh = [-Zky, Zkx, (Zkx^2+Zky^2)/2];
        C[2:4,4] = (B .+ Bsh*2*pi);
        posAndStrength = -[real(dwdz), 0, 0, 0] .- im*[imag(dwdz), 0, 0, 0];
        ss = C\posAndStrength;

        # @bp
        deltaU = real(ss[1]);
        deltaV = real(ss[2]);
        deltaOmega = real(ss[3]);
        newStrength = real(ss[4]);

        # Change in U, V, and Omega in foil fixed frame
        fp.U = fp.U + deltaU;
        fp.V = fp.V + deltaV;
        fp.Omega = fp.Omega + deltaOmega;

        vp.timeOfVortexSheddingEvent = t;

        # Update vortex flag since we have shed another vortex
        vp.vortexFlag = vp.vortexFlag + 1;
        println(vp.vortexFlag)

        # Store information about vortex in vp.mvortex and vp.mvortex0
        vf = vp.vortexFlag;
        vp.mvortex[vf,1] = real(newvortex);
        vp.mvortex[vf,2] = imag(newvortex);
        vp.mvortex[vf,3] = newStrength;

        # Initial locations of vortex, location  of foil fixed frame, and time of generating new vortex
        vp.mvortex0[vf,1] = real(newvortex);
        vp.mvortex0[vf,2] = imag(newvortex);
        vp.mvortex0[vf,3] = fp.x;
        vp.mvortex0[vf,4] = fp.y;
        vp.mvortex0[vf,5] = fp.theta;
        vp.mvortex0[vf,6] = t;
    end

    # Loop version
    velocZeta = zeros(vp.maxNumOfVortices,1) .+ 0*im;
    if vp.vortexFlag != 0
        for j = 1:vp.vortexFlag
            # z0 is the location of jth vortex
            z0 = vp.mvortex[j,1] + im*vp.mvortex[j,2];
            strength = vp.mvortex[j,3];

            proxTol = 0.00001;
            if (abs(z0) - rc < proxTol) # if it is too close to the foil, destroy this vortex
                # println("destroying vortex...")
                 vp.mvortex[j,3] = 0.0;
                 velocZeta[j] = 0.0;
            else
                dw1dz = -us1/z0^2 - us2 /(z0+zc)^2;
                dw2dz = -vs1/z0^2 - vs2 /(z0+zc)^2;
                dw3dz = -omegas1/z0^2 - omegas2 /(z0+zc)^2;
                dwsdz = -NewA1/z0^2 - NewA2/(z0+zc)^2 + NewA3*(1/(z0+zc) - 1/z0);
                dwdz = fp.U*dw1dz + fp.V*dw2dz + fp.Omega*dw3dz + dwsdz;
                dwdz = dwdz + im*fp.gammac/z0 +  strength*im*(-1/(z0 - rc^2/conj(z0)));

                for pn = 1:vp.vortexFlag
                    if pn != j
                       strengthj = vp.mvortex[pn,3];
                       zj = vp.mvortex[pn,1] + im*vp.mvortex[pn,2];
                       dwdz = dwdz + strengthj*im*(1/(z0-zj) - 1/(z0 - rc^2/conj(zj)));
                     end
                end

                dFdz = 1 - (a^2)/(z0+zc)^2;  # Page# 15, Equation 3.15
                dF2dz2 = 2*(a^2)/(z0+zc)^3;
                Veloc0 = dwdz/dFdz - im*strength*dF2dz2/2/dFdz^2; # Page# 61, Equation 5.62

                # Page# 61 Equation 5.64
                Z0 = (z0 + zc + (a^2)/(z0+zc));
                veloc0 = conj(Veloc0) - (fp.U + im*fp.V + im*fp.Omega*Z0);

                dFdzc = dFdz;
                dFda = 2*a/(z0+zc);

                velocZeta[j] =  (veloc0 - dFdzc*zcdot - dFda*adot)/dFdz;
            end
        end
    end

    #= Compute the updated momenta of the foil-fluid system according to conservation
    of momentum (equation 3.64 on page 34 in Hailong thesis)=#
    NewLf1 = -pi*(NewA1 - a^2/rc^2*conj(NewA1) + NewA2 - rc^2*ak^2/a^2*conj(NewA2)); # Lsf1, Page# 88 Eqn 8.13
    NewLf2 = -pi*(NewA3*zc - conj(NewA3)*ak*conj(zc)); # Lsf2, Page# 88 Eqn 8.13
    Lf = NewLf1 + NewLf2; # Lsf = Lsf1 + Lsf2, Page# 88 Eqn 8.13 Impulse of the fluid due to shape change

    Lfx = real(Lf);
    Lfy = imag(Lf);

    NewAf1 = real(2*pi*im*(1+ak^2)*conjzc*NewA1 + 2*pi*im*(ak^3*rc^2*conjzc/a^2 + conjzc - ak*zc)*NewA2 + 2*pi*im/rc^2*(ak^2*zc*(2*rc^2-zc*conjzc)-a^2*conjzc)*conj(NewA1)
    + 4*pi*im*rc^2 *zc*ak^3/a^2*conj(NewA2))/2;
    NewAf2 = real(2*pi*im*zc*conjzc*NewA3 + 2*pi*im*zc*conjzc*ak^2*conj(NewA3))/2;

    Af = NewAf1 + NewAf2;
    # Page 88, Equation 8.12, Linear and Angular momentum components of the deforming foil from shape change
    Newm2 = pi*rc^2*(1-2*ak^2+ak^4*conjzc^2/a^2); # m2
    Newm3 = pi*(2*a*conjzc - 2*a*ak*zc + 2*rc^2*ak^3*conjzc/a); # m3
    NewJ2 = -pi*rc^2*im*(conjzc + 2*ak^3*zc - ak^5*conjzc*(rc^2+zc*conjzc)/a^2); # J2
    NewJ3 = -pi*im*(2*a^3 + a*conjzc^2 + a*ak^2*zc^2 - rc^2*ak^4*(rc^2+2*zc*conjzc)/a); # J3
    Lb = Newm2*zcdot + Newm3*adot; # Lsb
    Ab = real(zcdot*NewJ2 + adot*NewJ3); # Psb

    Lbx = real(Lb);
    Lby = imag(Lb);

    P1 = [(Lbx + Lfx), (Lby + Lfy), (Ab + Af)];
    # [(Linear momentum of body in x + impulse of fluid in x dirtn); (Linear momentum of body in y dirtn + impulse of fluid in y dirtn);(Angular momentum of body + impulse couple of fluid)]

    Pvortex = [0.0, 0.0, 0.0];

    # Loop version
    if vp.vortexFlag != 0
        for pn = 1:vp.vortexFlag
            # @bp
            strength = vp.mvortex[pn,3];
            zk = vp.mvortex[pn,1] + im*vp.mvortex[pn,2];
            Zk = (zk + zc + (a^2)/(zk+zc));
            Zkx = real(Zk);
            Zky = imag(Zk);

            (impulse, impulseCouple) = computeImpulseAndImpulseCoupleForVortex(zk, rc, zc, a);

            B = [real(impulse), imag(impulse), impulseCouple];
            Bsh = [-Zky, Zkx, (Zkx^2+Zky^2)/2];
            Pvortex = Pvortex .+ B*strength .+ Bsh*strength*2*pi;
        end
    end

    # impulse0 and impulseCouple0 are the impulse and impulse couple from gammac

    if fp.gammac != 0.0
       impulse0 = fp.gammac*2*pi*im*zc + fp.gammac*2*pi*im*(xbody + im*ybody);
       impulseCouple0 = fp.gammac*2*pi*(-imag((-im)*(rc^2 + zc*conjzc + (a^2)*(a^2)/(rc^2 - zc*conjzc))/2)) + 1/2*(-2*pi*fp.gammac)*(xbody^2+ybody^2);
       Pvortex = Pvortex .+ [real(impulse0), imag(impulse0), impulseCouple0];
    end

    # @bp
    P = P1 + Pvortex;
    LMb = Rmatrix*[sp.Pinitial[1], sp.Pinitial[2]];
    AMb = sp.Pinitial[3] - (fp.x*sp.Pinitial[2] - fp.y*sp.Pinitial[1]);
    Mb = [LMb[1], LMb[2], AMb];

    # @bp
    Vbody = Imatrix\(Mb - P); # equation  3.64, Page 34
    fp.U = Vbody[1]; fp.V = Vbody[2]; fp.Omega = Vbody[3];
    # @bp
    return velocZeta, vp.vortexFlag
end

function computeImpulseAndImpulseCoupleForVortex(zk, rc, zc, a)    # page 86-90 in Hailong thesis
    # Page 86-90, Fluid momentum due to point vortices
    conjzc = conj.(zc);
    conjzk = conj.(zk);
    Zk = (zk .+ zc) .+ (a^2)./(zk .+zc);
    impulse = (-2*pi*im).*(Zk .- zk  .+ rc^2 ./conj.(zk) );
    w3zk = (-im/2).*(rc^2 + zc.*conjzc + a^4 ./(rc^2 .-zc.*conjzc) .+ 2*rc^2 .*(zc .+ a^2 ./zc)./zk .- 2*a^2 .*(rc^2 .- zc.*conjzc)./zc./(zk.+zc) .- 2*zc*a^4 ./(rc^2 .- zc.*conjzc)./(zk.+zc));
    impulseCouple = 2*pi*imag.(w3zk);
    return impulse, impulseCouple
end

end
