# Make an animation of the Joukowski foil shedding vortices
plotlyjs()

xplot = plot(sol.t,[sol[1,:], sol[2,:], sol[3,:]], xaxis = ("t", (0.0,T)), label=["x" "y" "θ"],linewidth=2.0)
# plot([0,1,2],[0,1,2], seriestype=:scatter, zcolor=[1,2,3], markerstrokewidth=0, markersize=2.0, c=:redsblues, legend=false, colorbar=false)

solArray = convert(Array,sol);
fname = "video/foilVortices";
Tspaniter = 0.0:deltaT:T;
numDigits = ndigits(length(sol.t));
i = 1;
maxStrength = maximum(abs.(vp.mvortex[:,3]));
(finalFoilData, finalVortexData, finalColorData, finalAlphaData) = getFoilAndVortexPlotData(sol,sol.t[end],maxStrength);

# currentR = 0; currentG = 0; currentB = 1;
# tint_factor = 0.0;
# newR = floor(Int8, currentR + (255 - currentR) * tint_factor)
# newG = floor(Int8, currentG + (255 - currentG) * tint_factor)
# newB = floor(Int8, currentB + (255 - currentB) * tint_factor)
# testcolor = RGB(newR/255.0,newG/255.0,newB/255.0);
# testcolors1 = [RGB(0,0,1) for i in 1:100];
# testcolors2 = [RGB(1,0,0) for i in 101:length(finalStrengthData)];
# testcolors = [testcolors1; testcolors2];
# alphas = [0.5 for i in 1:length(finalStrengthData)];

plot(finalFoilData, xaxis = ("x", (-5,15)), yaxis=("y",(-3,3)), color=:black, linewidth=1.5, legend=false,aspect_ratio=0.85)
plot!(real.(finalVortexData), imag.(finalVortexData), seriestype=:scatter, markersize=0.8, markerstrokewidth=0, markeralpha = finalAlphaData, color=finalColorData, legend=false, colorbar=false, grid=false, showaxis=false, xlabel=false,ylabel=false)
savefig("joukowskinew2.pdf")

for t in sol.t
    (foilData, vortexData, colorData, alphaData) = getFoilAndVortexPlotData(sol,fp,t,maxStrength);
    foilPlot = plot(foilData, xaxis = ("x", (-5,15)), yaxis=("y",(-5,5)), color=:black, linewidth=1.5, legend=false, aspect_ratio=0.85);

    if isempty(vortexData) != 1
        plot!(real.(vortexData), imag.(vortexData), seriestype=:scatter, markersize=0.8, markerstrokewidth=0, markeralpha = alphaData, color=colorData, clims=(0.0,1.0), legend=false, colorbar=false, grid=false, showaxis=false, xlabel=false,ylabel=false);
    end
    numZeros = numDigits - ndigits(i);
    firstPart = join(["0" for j in 1:numZeros]);
    imagei = firstPart*string(i);
    fn = fname*imagei;
    savefig(fn)
    i = i + 1;
end

function getFoilAndVortexPlotData(sol,t,maxStrength)
    rc = fp.rc; # radius of the circle in zeta plane. Considered to be 1
    ua = fp.ua; # amplitude of flapping
    uf = fp.uf; # frequency of flapping
    zcy = ua*sin(2*pi*uf*t);

    # ak is a constant that defines the shape of the foil with constant area.
    ak = 0.85^2/(rc^2 - 0.15^2); # variable 'k' Equation 8.2, Page# 85 (Chapter of Model Reduction) from Hailong's Thesis document
    zcx = (1-ak)/(1+ak)*sqrt(rc^2 - zcy^2); # This is equivalent to Zeta_x from Thesis document; % variable 'Zeta_x' -- Equation 8.2, Page# 85 (Chapter of Model Reduction) from Hailong's Thesis document
    a = -2*ak/(1+ak)*sqrt(rc^2 - zcy^2);

    th = range(0,stop=2*pi,length=100);

    # zeta = x .+ y;
    zetac = zcx + im*zcy;
    zeta = rc*exp.(im.*th);

    foilLocationConfMap = sol(t)[1] .+ sol(t)[2]*im .+ exp(im*sol(t)[3]).*(zeta .+ zetac .+ a^2 ./ (zeta .+ zetac));
    # foilLocationConfMap = solArray[1,1,i] .+ solArray[2,1,i]*im .+ exp(im*solArray[3,1,i]).*(zeta .+ zetac .+ a^2 ./ (zeta .+ zetac));

    # Get vortex locations at time t
    # xindices = 4:2:(vp.vortexFlag);
    #
    # # xindices = xindices[broadcast(abs, sol(t)[xindices]) .> 0.0];
    # xindices = xindices[broadcast(abs, solArray[xindices,1,i]) .> 0.0];
    # xindices = xindices[broadcast(abs, vp.mvortex[xindices,3]) .> 0.0];
    # yindices = xindices .+ 1;
    #
    # # whichVorticesToPlotX = filter!(x -> vp.mvortex[]  whichVorticesToPlotX)
    # # xpv = vp.mvortex0[xindices .- 3,1] .+ sol(t)[xindices];
    # xpv = vp.mvortex0[xindices .- 3,1] .+ solArray[xindices,1,i];
    # # ypv = vp.mvortex0[yindices .- 3,2] .+ sol(t)[yindices];
    # ypv = vp.mvortex0[yindices .- 3,2] .+ solArray[yindices,1,i];
    #
    # vortexLocation = xpv + ypv.*im;
    # # vortexLocConfMap = exp(im*sol(t)[3]).*(vortexLocation .+ zetac .+ a^2 ./(vortexLocation .+ zetac)) .+ sol(t)[1] .+ sol(t)[2]*im;
    # vortexLocConfMap = exp(im*solArray[3,1,i]).*(vortexLocation .+ zetac .+ a^2 ./(vortexLocation .+ zetac)) .+ solArray[1,1,i] .+ solArray[2,1,i]*im;
    (vortexLocConfMap, colorData, alphaData) = getVortexPlotData(sol,vp,fp,t,maxStrength);

    # colors = abs.(vp.mvortex[xindices .- 3,3]);
    return foilLocationConfMap, vortexLocConfMap, colorData, alphaData;
end

function getVortexPlotData(sol,vp,fp,t,maxStrength)
    b = 0.0;
    rc = fp.rc; # radius of the circle in zeta plane. Considered to be 1
    ua = fp.ua; # amplitude of flapping
    uf = fp.uf; # frequency of flapping
    zcy = ua*sin(2*pi*uf*t);

    # ak is a constant that defines the shape of the foil with constant area.
    ak = 0.85^2/(rc^2 - 0.15^2); # variable 'k' Equation 8.2, Page# 85 (Chapter of Model Reduction) from Hailong's Thesis document
    zcx = (1-ak)/(1+ak)*sqrt(rc^2 - zcy^2); # This is equivalent to Zeta_x from Thesis document; % variable 'Zeta_x' -- Equation 8.2, Page# 85 (Chapter of Model Reduction) from Hailong's Thesis document
    a = -2*ak/(1+ak)*sqrt(rc^2 - zcy^2);
    zc = zcx + im*zcy;

    vortexData = [0.0+0.0*im];
    colorData = [RGB(1,1,1) for i in 1:vp.vortexFlag]; alphaData = zeros(vp.vortexFlag);
    # vortexData = [0.0+0.0*im for i in 1:vp.vortexFlag];

    for jj = 1:vp.vortexFlag
        # jj = index(N);
        vortexlocationx = sol(t)[2*jj+2];
        vortexlocationy = sol(t)[2*jj+3];
        if vortexlocationx != 0.0

            vlocationx = vortexlocationx + vp.mvortex0[jj,1];
            vlocationy = vortexlocationy + vp.mvortex0[jj,2];
            p_zeta = vlocationx + im*vlocationy;
            p_z = (p_zeta + zc + (a^2+b^2)/(p_zeta + zc) - (a*b)^2/3/(p_zeta+zc)^3 )*exp(im*sol(t)[3]) + sol(t)[1] + im*sol(t)[2];

            if abs(vp.mvortex[jj,3]) > 0.0
                if jj == 1
                    vortexData[jj] = p_z;
                else
                    append!(vortexData, p_z);
                    # vortexData[jj] = p_z;
                    scalingFactor = 5.0;
                    alphaVal = vp.mvortex[jj,3]/maxStrength;
                    alphaData[jj] = scalingFactor*abs(alphaVal);
                    # append!(alphaData, abs(alphaVal));
                    if alphaVal > 0.0
                        # append!(colorData, RGB(1,0,0) );
                        colorData[jj] = RGB(1,0,0);
                    else
                        # append!(colorData, RGB(0,0,1) );
                        colorData[jj] = RGB(0,0,1);
                    end
                end
            end
        end
    end
    return vortexData, colorData, alphaData
end
