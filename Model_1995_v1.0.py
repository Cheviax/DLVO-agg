import math

pi      =   3.141592653
kk      =   1.38e-23

hc      =   1.25
uc_kt   =   6.77
fc      =   15.7
k0      =   13.8

hc      *=  1e-9
fc      *=  1e-12

temp    =   298.15
uc      =   uc_kt*kk*temp
zz      =   5
aa      =   17.5e-9
dl      =   1.3
df      =   1.6
miu     =   1.65
miu     *=  1e-3
phi_p   =   0.0139
phi_m   =   0.6
phi_ref =   1.1

tau_0   =   (6*pi*miu*(aa**3))/(kk*temp)

tau_sic =   (tau_0**-1)*(k0**0.5)*((uc*aa)/(kk*temp*hc))*math.exp((-zz*uc)/(kk*temp))
tau_si  =   tau_sic ** -1

for i in  range(10):
    gm      =   i*0.5-1.5
    gama    =   10**gm
    qq_min = aa
    qq_max = aa * 100

    while True:
        qq = (qq_max + qq_min) / 2
        phi_int = (qq / aa) ** (df - 3)
        phi_a = phi_p / phi_int
        if  phi_a > phi_m*0.99/(phi_ref**3):
            phi_a = phi_m*0.99/(phi_ref**3)
            phi_int = phi_p / phi_a
            qq = aa*(phi_int**(1/(df-3)))

        tau_sc = (tau_0 ** -1) * (k0 ** 0.5) * (uc_kt ** 1.5) \
                 * ((aa / hc) ** 2) * (aa / qq) * math.exp(-uc_kt)
        tau_s = tau_sc ** -1
        tau_rc = (k0 ** 0.5) * gama * (aa / hc) * ((qq / aa) ** -dl)
        tau_r = tau_rc ** -1
        tau_sr = tau_0 * qq / aa
        tau_rs = tau_si * ((aa / qq) ** dl)

        alpha = (tau_s * (tau_r + tau_rs)) / ((tau_s * tau_rs) + (tau_sr * (tau_r + tau_rs)))

        eta1 = miu * ((1 - (phi_a / phi_m)*(phi_ref**3)) ** (-2))
        eta2 = ((alpha * fc * hc * (k0 ** -0.5) * (phi_a ** 2))
                / ((aa ** 3) * gama)) * ((phi_p/phi_a) ** ((5 - 2 * df - dl)/(3 - df)))
        eta = eta1 + eta2
        etar = eta / miu
        gama_c = (2 * fc * aa) / (5 * pi * eta * (qq ** 3))
        gama_if = gama_c / gama

        if gama_if > 1.0001:
            qq_min = qq

        elif gama_if < 0.9999:
            qq_max = qq

        else:
            qa = qq / aa
            line    =   'q/a: '+str('%.3g' % qa)+' '+str('%.3g' % gama)+' '+str('%.3g' % eta1)\
                        +' '+str('%.3g' % eta2)+' '+str('%.3g' % eta)+' '+str('%.3g' % etar)\
                        +' phi_a: '+str('%.3g' % phi_a)
            print(line)
            break

