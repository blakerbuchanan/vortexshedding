# This file is common for both Joukowski and Von Mises foil... Hence the variable 'b' is initialized to zero.
# Positive strength-- Red Vortex, clockwise
# Negative strength-- Blue Vortex, Counter clockwise
rc = 1;
zcy = 0.0;
alpha = 7*pi/24;
ak = 0.85^2/(rc^2 - 0.15^2);
zcx = (1-ak)/(1+ak)*sqrt(rc^2 - zcy^2);
a = -2*ak/(1+ak)*sqrt(rc^2 - zcy^2);
b = 0;

zc = zcx + im*zcy;

conjzc = conj(zc);

t = 0;

conjb = conj(b);

# This is where the initial velocity of the foils is given as input
U0 = 1;  V0 = 0; Omega0 = 0;
# To decide the distribution of initial ambient vortices
startp = 1.5;

npv = 0;

if npv > 0
    for k = 1:npv

    # zeta(k) = 1.5 + (-1)^(k-1)*0.6*i;
    # z(k) = zeta(k) + zc + (a^2+b^2)/(zeta(k)+zc) - (a*b)^2/3/(zeta(k)+zc)^3;
    zeta[k] = startp + 4*0 + (k-1)*2 + (-1)^(k-1)* 0.75*im;
    z[k] = zeta[k] + zc + (a^2+b^2)/(zeta[k]+zc) - (a*b)^2/3/(zeta[k]+zc)^3;
    gamma[k] = (1)^k*1;

   end
end

us0 = zc;
us1 = -rc^2;
us2 = a^2+b^2;
us3 = 0;
us4 = - (a*b)^2/3;

vs0 = -i*zc;
vs1 = -i*rc^2;
vs2 = -i*(a^2+b^2);
vs3 = 0;
vs4 = i*(a*b)^2/3;


omegas0 = -i*(rc^2 + zc*conjzc + (a^2+b^2)*(a^2+ conjb^2)/(rc^2 - zc*conjzc) - 2*(a^2+b^2)*(a*conjb)^2*zc^2/3/(rc^2 - zc*conjzc)^3 ...
                   + (a*b)^2*(a*conjb)^2/9 *( (rc^4 + 4*rc^2*zc*conjzc + zc^2*conjzc^2)/(rc^2 - zc*conjzc)^5 ) )/2;
omegas1 = -i*(2*rc^2*zc + 2*(a^2+b^2)*rc^2/zc  - 2*(a*b)^2/3*rc^2/zc^3 )/2;
omegas2 = -i*(-2*(a^2+b^2)*(rc^2 - zc*conjzc)/zc + (a^2+b^2)*(a^2+conjb^2)*(-2*zc)/(rc^2 - zc*conjzc) + 2*(a*b)^2/3*(rc^2/zc^3) ...
                   + 2*(a^2+conjb^2)*(a*b)^2/3*rc^2*conjzc/(rc^2 - zc*conjzc)^3 + 2*(a^2+b^2)*(a*conjb)^2/3*zc^3/(rc^2 - zc*conjzc)^3 ...
                   + (a*b)^2*(a*conjb)^2/9*( -6*zc*conjzc*zc/(rc^2 - zc*conjzc)^4 - 3*zc/(rc^2 - zc*conjzc)^3 - 6*zc*rc^2*(rc^2+zc*conjzc)/(rc^2-zc*conjzc)^5 )  )/2;
omegas3 = -i*(-2*(a*b)^2/3*(-rc^2/zc^2)  - 2*(a^2+conjb^2)*(a*b)^2/3*rc^2/(rc^2 - zc*conjzc)^2 ...
                   + (a*b)^2*(a*conjb)^2/9*( 3*zc*conjzc*zc^2/(rc^2-zc*conjzc)^4+3*rc^2*zc^2/(rc^2-zc*conjzc)^4 + 3*zc^2/(rc^2 - zc*conjzc)^3  ) )/2;
omegas4 = -i*(2*(a*b)^2/3*((rc^2 - zc*conjzc)/zc )  + 2*(a^2+conjb^2)*(a*b)^2/3*zc/(rc^2 - zc*conjzc)  ...
                   +  (a*b)^2*(a*conjb)^2/9*( -2*zc^3/(rc^2 - zc*conjzc)^3    ) )/2 ;


[imU, imcU] = GetImpulseAndImpulseCouple0(us0, us1, us2, us3, us4, rc, zc, a, b);

[imV, imcV] = GetImpulseAndImpulseCouple0(vs0, vs1, vs2, vs3, vs4, rc, zc, a, b);

[imOmega, imcOmega] = GetImpulseAndImpulseCouple0(omegas0, omegas1, omegas2, omegas3, omegas4, rc, zc, a, b);


Ifluid = [real(imU)    real(imV)    real(imOmega);
             imag(imU)  imag(imV)   imag(imOmega);
             imcU          imcV          imcOmega];

S = pi*rc^2* ...
    (1 - (a^2+b^2)*(a^2 + conjb^2)/(rc^2 - zc*conjzc)^2  ...
    + (a^2*conjb^2*(a^2+b^2)*zc^2 + a^2*b^2*(a^2+ conjb^2)*conjzc^2)/(rc^2 - zc*conjzc)^4 ...
    - (a*b)^2*(a*conjb)^2*(rc^4 + 6*rc^2*zc*conjzc + 3*zc^2*conjzc^2)/3/(rc^2 - zc*conjzc)^6  );

S = real(S);

m1 = i*pi*rc^2* ...
(zc +  (a^2+conjb^2)*conjzc*((a^2+b^2)^2 - 2*(a*b)^2/3)/(rc^2 - zc*conjzc)^3 ...
+  4*(a^2+b^2)*(a^2+conjb^2)*(-(a*b)^2/3)*conjzc^3/(rc^2 - zc*conjzc)^5  ...
+  3*((a^2+b^2)^2 - 2*(a*b)^2/3)*(-(a*conjb)^2/3)*zc*(rc^2 + zc*conjzc)/(rc^2 - zc*conjzc)^5 ...
+ 12*conjzc*(a*b)^2*(a*conjb)^2/9*(a^2+b^2)*(rc^4 + 3*rc^2*zc*conjzc + zc^2*conjzc^2)/(rc^2 - zc*conjzc)^7 ...
+ 3*conjzc^5*(a*b)^4/9*(a^2+conjb^2)/(rc^2 - zc*conjzc)^7 ...
+ 3*conjzc^3*(a*b)^4/9*(-(a*conjb)^2/3)*(10*rc^4 + 15*rc^2*zc*conjzc + 3*zc^2*conjzc^2)/(rc^2 - zc*conjzc)^9);


Iomega = -1/2/(rc^2 - conjzc*zc)^12* ...
    (pi*rc^2*(3*(-(a*conjb)^2/3)^2*((a^2+b^2)^2*zc^4*(rc^2 - conjzc*zc)^4*(5*rc^2 + 2*conjzc*zc) ...
    + (-(a*b)^2/3)^2*(rc^10 + 30*conjzc*rc^8*zc + 150*conjzc^2*rc^6*zc^2 + 200*conjzc^3*rc^4*zc^3 ...
    + 75*conjzc^4*rc^2*zc^4 + 6*conjzc^5*zc^5) ...
    + 2*(-(a*b)^2/3)*zc^2*(rc^2 - conjzc*zc)^2*(zc^2*(rc^2 - conjzc*zc)^2*(5*rc^2 + 2*conjzc*zc) ...
    + 2*(a^2+b^2)*(5*rc^6 + 20*conjzc*rc^4*zc + 15*conjzc^2*rc^2*zc^2 +  2*conjzc^3*zc^3)))...
    + (rc^2 - conjzc*zc)^4*(-(rc^2 - conjzc*zc)^8*(rc^2 + 2*conjzc*zc) + (a^2+conjb^2)^2*(((a^2+b^2)^2 + 2*(-(a*b)^2/3))*rc^10 ...
    - 2*conjzc*((a^2+b^2)^2 + 2*(-(a*b)^2/3))*rc^8*zc + 2*conjzc^2*rc^6*(6*(a^2+b^2)*(-(a*b)^2/3) - ((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^2) + 8*conjzc^3*rc^4*zc*(-2*(a^2+b^2)*(-(a*b)^2/3) + ((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^2) ...
    + conjzc^4*rc^2*(15*(-(a*b)^2/3)^2 - 4*(a^2+b^2)*(-(a*b)^2/3)*zc^2 - 7*((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^4) + 2*conjzc^5*(3*(-(a*b)^2/3)^2*zc + 4*(a^2+b^2)*(-(a*b)^2/3)*zc^3 + ((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^5))) ...
    + 2*(-(a*conjb)^2/3)*(rc^2 - conjzc*zc)^2*((rc^2 - conjzc*zc)^2*(((a^2+b^2)^2 + 2*(-(a*b)^2/3))*rc^10 -  2*conjzc*((a^2+b^2)^2 + 2*(-(a*b)^2/3))*rc^8*zc ...
    + 2*conjzc^2*rc^6*(6*(a^2+b^2)*(-(a*b)^2/3) - ((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^2) +  8*conjzc^3*rc^4*zc*(-2*(a^2+b^2)*(-(a*b)^2/3) + ((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^2) ...
    + conjzc^4*rc^2*(15*(-(a*b)^2/3)^2 - 4*(a^2+b^2)*(-(a*b)^2/3)*zc^2 - 7*((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^4) + 2*conjzc^5*(3*(-(a*b)^2/3)^2*zc + 4*(a^2+b^2)*(-(a*b)^2/3)*zc^3 + ((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^5)) ...
    + 2*(a^2+conjb^2)*(20*conjzc^3*(-(a*b)^2/3)*rc^4*zc*(3*(-(a*b)^2/3) - 2*(a^2+b^2)*zc^2) + 10*conjzc*rc^8*zc*(2*(a^2+b^2)*(-(a*b)^2/3) - ((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^2) ...
    + rc^10*(2*(a^2+b^2)*(-(a*b)^2/3) + 3*((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^2) + 5*conjzc^4*rc^2*zc^2*(9*(-(a*b)^2/3)^2 + 4*(a^2+b^2)*(-(a*b)^2/3)*zc^2 - ((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^4) ...
    + 5*conjzc^2*rc^6*(3*(-(a*b)^2/3)^2 - 2*(a^2+b^2)*(-(a*b)^2/3)*zc^2 + 2*((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^4) + 2*conjzc^5*(3*(-(a*b)^2/3)^2*zc^3 + 4*(a^2+b^2)*(-(a*b)^2/3)*zc^5 + ((a^2+b^2)^2 + 2*(-(a*b)^2/3))*zc^7)))));


Iomega = real(Iomega);


Ibody = [S    0   real(m1);
             0    S   imag(m1);
            real(m1)  imag(m1)  Iomega];


Imatrix = Ifluid + Ibody;




P1 = Imatrix*[U0 V0 Omega0]';



if npv == 0

    % initially there is no vortex

    Pic = P1;




    fid = fopen('icdata.m','w');

    fprintf(fid, '%12.8f\n', U0);
    fprintf(fid, '%12.8f\n', V0);
    fprintf(fid, '%12.8f\n', Omega0);

    fprintf(fid, '%12.8f\n', npv);

    fprintf(fid, '%12.8f %12.8f %12.8f\n', Pic);

    fclose(fid);



 else

    Pvortex = 0;


         for   k = 1 : npv

                zk = zeta(k);
                Zk = z(k);
                strength = gamma(k);
                Zkx = real(Zk);
                Zky = imag(Zk);


                [im, imc] = GetVortexImpulse0(zk, rc, zc, a, b);

                Impulsevortex = im;
                Imcouplevortex = imc;

                B = [real(Impulsevortex) imag(Impulsevortex) Imcouplevortex]';
                Bsh = [-Zky Zkx (Zkx^2+Zky^2)/2]';

                Pvortex = Pvortex + B*strength + Bsh*strength*2*pi;



          end




         Pic = P1 + Pvortex;

         pvlocation = [real(zeta); imag(zeta);gamma];

         fid = fopen('icdata.m','w');

         fprintf(fid, '%12.8f\n', U0);
         fprintf(fid, '%12.8f\n', V0);
         fprintf(fid, '%12.8f\n', Omega0);

         fprintf(fid, '%12.8f\n', npv);

         fprintf(fid, '%12.8f %12.8f %12.8f\n', pvlocation);

         fprintf(fid, '%12.8f %12.8f %12.8f\n', Pic);



         fclose(fid);



end


tt = [0:0.01:2*pi];

    for k = 1:length(tt)

       circle(k) = rc*exp(j*tt(k));
       foil(k) = (circle(k) + zc + (a^2+b^2)/(circle(k) + zc) - (a*b)^2/3/(circle(k)+zc)^3 );


   end


if npv > 0

    plot(real(foil), imag(foil), real(z), imag(z),'o')
    axis equal
    axis([-3 20 -2 2])
    grid on

end
