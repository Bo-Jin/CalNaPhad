import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#necessary input variables 
Vm = float(eval(input("Molar volume (cm3/mol) = ")))
print('\n')

Tm = float(eval(input("Melting temperature (K) = ")))
print('\n')

DeltaSm = float(eval(input("Melting entropy (J/(mol·K))= ")))
print('\n')

Sigma_sl = float(eval(input("Interfacial energy of solid/liquid (J/m2) = ")))
print('\n')

Sigma_lg = float(eval(input("Interfacial energy of liquid/gas (J/m2) = ")))
print('\n')

#fitting curve
def func(Tdata, xi, DeltaSigma):
    return xi*np.log((DeltaSigma*Vm*1000)/(xi*DeltaSm*(Tm-Tdata)))

xi = []

d0 = []

rcritical = []

DeltaSigma = []

Sigma_sg = []

N = int(input("The number of surfaces (1 or 2 or 3): "))
print('\n')

for i in range(0,N,1):
    #read experimental data
    filename = input("Please enter text name of surface melting data for macro-crystals (including file type, columns seperated by space): ")
    Tdata = []
    deqdata = []
    Tdata, deqdata = np.loadtxt(filename, usecols=(0,1), dtype=float, unpack=True)
    plt.plot(Tdata, deqdata, 'o', color='b', label='data')

    Tdatamin = min(Tdata) - 50

    print('\n')
    print("**********Calculate surface melting**********")
    print("**********Please wait**********")

    #fit for the parameters
    popt, pcov = curve_fit(func, Tdata, deqdata)
    #print(popt)
    print("The covariance =")
    print(pcov)

    xi.append(popt[0])
    DeltaSigma.append(popt[1])

    print('\n')
    print("The Interaction length (nm) = %.6f" % xi[i])

    print('\n')
    print("The Difference among interfacial energies (J/m2) = %.6f" % DeltaSigma[i])

    print('\n')
    input("Press <enter> to continue")

    print('\n')
    d0.append(float(eval(input("Minimum thickness (nm) = "))))

    Sigma_sg.append(Sigma_sl + Sigma_lg + DeltaSigma[i])

    print('\n')
    print("**********Calculate premelting low limit and high limit**********")
    print("**********Please wait**********")


    #surface premelting low limit
    Tsmlow = Tm-(DeltaSigma[i]*Vm*10**(-6))/(xi[i]*10**(-9)*DeltaSm)

    #surface premelting high limit
    Tsmhigh = Tm - (DeltaSigma[i]*Vm*10**(-6))/(xi[i]*10**(-9)*DeltaSm)*math.exp(-d0[i]/xi[i])

    print('\n')
    print("The surface premelting low limit = %.6f" % Tsmlow)

    print('\n')
    print("The surface premelting high limit = %.6f" % Tsmhigh)

    print('\n')
    input("Press <enter> to continue")


    print('\n')
    print("**********Plot surface melting diagram**********")
    print("**********Please wait**********")

    Tval = np.arange(Tdatamin, Tm + 200, 1)
    deqval = func(Tval, xi[i], DeltaSigma[i])
    plt.plot(Tval, deqval, color='r', label='fit')
    plt.xlabel('Temperature, K')
    plt.ylabel('Thickness, nm')
    plt.axvline(Tsmlow, ls='--', c="green")
    plt.axvline(Tsmhigh, ls='--', c="green")
    plt.axvline(Tm, ls='--', c="green")
    plt.legend()
    plt.savefig('Surface melting_%i.ps' % (i+1), format='ps',dpi=1000)
    plt.show()

    print('\n')
    input("Press <enter> to continue")



    print('\n')
    print("**********Calculate critical size**********")
    print("**********Please wait**********")


    #critical size
    n = 30 #change n to increase the accuracy of critical size r
    yliquidus0 = []
    ysolidus0 = []
    V_r = []
    r = 0.1
    while r<=100:
        ysolidus = (1-d0[i]/r)**3
        y0 = []
        Gsl0 = []
        Gl0 = []
        for yy in range(1,1001,1):
            y = float(yy/1000)
            Teq = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
            Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
            Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg
            y0.append(y)
            Gsl0.append(Gsl)
            Gl0.append(Gl)
        M = len(y0)
        for k in range(0,M-1,1):
            if (Gsl0[k]-Gl0[k])*(Gsl0[k+1]-Gl0[k+1])<0:
                yliquidus0.append((y0[k]+y0[k+1])/2)
                ysolidus0.append(ysolidus)
                V_r.append(r)
                break
        r = r + 0.1

    M = len(V_r)
    for k in range(0,M-1,1):
        if (yliquidus0[k]-ysolidus0[k])*(yliquidus0[k+1]-ysolidus0[k+1])<0:
            r0 = V_r[k]
            r1 = V_r[k+1]

    ###############
    r = r0
    ysolidus = (1-d0[i]/r)**3
    y0 = 10**(-6)
    y1 = 1
    y = y0
    Teq = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
    Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
    Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg
    Y0 = Gsl - Gl

    y = y1
    Teq = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
    Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
    Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg
    Y1 = Gsl - Gl

    if Y0*Y1<0:
        for t in range(0,n,1):
            y = (y0 + y1)/2
            Teq = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
            Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
            Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg
            Y = Gsl - Gl

            if Y*Y0>0:
                y0 = y
                Y0 = Y
            else:
                y1 = y
                Y1 = Y
        yliquidus = (y0 + y1)/2

    R0 = ysolidus - yliquidus

    ####################
    r = r1
    ysolidus = (1-d0[i]/r)**3
    y0 = 10**(-6)
    y1 = 1
    y = y0
    Teq = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
    Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
    Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg
    Y0 = Gsl - Gl

    y = y1
    Teq = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
    Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
    Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg
    Y1 = Gsl - Gl

    if Y0*Y1<0:
        for t in range(0,n,1):
            y = (y0 + y1)/2
            Teq = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
            Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
            Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg
            Y = Gsl - Gl

            if Y*Y0>0:
                y0 = y
                Y0 = Y
            else:
                y1 = y
                Y1 = Y
        yliquidus = (y0 + y1)/2

    R1 = ysolidus - yliquidus

    if R1*R0<0:
        for h in range(0,n,1):
            r = (r0 + r1)/2
            ysolidus = (1-d0[i]/r)**3
            y0 = 10**(-6)
            y1 = 1
            y = y0
            Teq = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
            Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
            Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg
            Y0 = Gsl - Gl

            y = y1
            Teq = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
            Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
            Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg
            Y1 = Gsl - Gl
            if Y0*Y1<0:
                for t in range(0,n,1):
                    y = (y0 + y1)/2
                    Teq = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                    Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                    Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Teq)+Sigma_lg
                    Y = Gsl - Gl

                    if Y*Y0>0:
                        y0 = y
                        Y0 = Y
                    else:
                        y1 = y
                        Y1 = Y
                yliquidus = (y0 + y1)/2
            R = ysolidus - yliquidus
            if R0*R>0:
                r0 = r
                R0 = R
            else:
                r1 = r
                R1 = R

        rcritical.append(round((r0 + r1)/2, 6))

    print('\n')
    print("The critical size (nm) = %.6f" % rcritical[i])

    #print("{:.4}".format(rcritical))

    print('\n')
    input("Press <enter> to continue")

    print('\n')
    print("**********Plot intersection diagram**********")
    print("**********Please wait**********")

    #intersection diagram
    #liquidus_y_T
    yliquidus0 = []
    Tliquidus0 = []
    r0 = []
    n = 30
    r = 0.5*rcritical[i]
    while r<=100:
        r0.append(r)
        y0 = 10**(-6)
        y1 = 1
        y = y0
        Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
        Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
        Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
        Y0 = Gsl-Gl

        y = y1
        Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
        Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
        Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
        Y1 = Gsl-Gl

        if Y0*Y1<0:
            for t in range(0,n,1):
                y = (y0 + y1)/2
                Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                Y = Gsl - Gl

                if Y*Y0>0:
                    y0 = y
                    Y0 = Y
                else:
                    y1 = y
                    Y1 = Y
            yliquidus = (y0 + y1)/2
        yliquidus0.append(yliquidus)
        r = r + 0.01

    j = 0
    r = 0.5*rcritical[i]
    while r<=100:
        Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(yliquidus0[j]**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(yliquidus0[j]**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(yliquidus0[j]**(1/3)-1))
        Tliquidus0.append(Tliquidus)
        j = j+1
        r = r + 0.01


    #solidus_y_T
    ysolidus0 = []
    Tsolidus0 = []
    r = 0.5*rcritical[i]
    while r<=100:
        ysolidus = (1-d0[i]/r)**3
        ysolidus0.append(ysolidus)
        Tsolidus = Tm-2*1000*Vm/(r*DeltaSm*(ysolidus**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1))
        Tsolidus0.append(Tsolidus)
        r = r + 0.01


    #plot radius - solid phase fraction
    plt.xlabel('Radius, nm')
    plt.ylabel('Solid phase fraction')
    plt.plot(r0,yliquidus0, linewidth=2.0)
    plt.plot(r0,ysolidus0, linewidth=2.0)
    plt.xlim(0, 100)
    plt.ylim(0, 1)
    plt.savefig('Intersection_%i.ps'% (i+1), format='ps',dpi=1000)
    plt.show()

    print('\n')
    input("Press <enter> to continue")

    print('\n')
    print("**********Plot nano-phase diagram**********")
    print("**********Please wait**********")
    print('\n')

    #below critical size
    ycritical = (1-d0[i]/rcritical[i])**3
    Tcritical = Tm-2*1000*Vm/(rcritical[i]*DeltaSm*(ycritical**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1))
    fcor = (Tm-Tcritical)*(rcritical[i]*DeltaSm)/(3*1000*Vm*(Sigma_sg[i]-Sigma_lg))
    Tnano0 = []
    rnano = []
    r = 0.000001
    while r<=rcritical[i]:
        rnano.append(r)
        Tnano = Tm-3*1000*Vm/(r*DeltaSm)*(Sigma_sg[i]-Sigma_lg)*fcor
        Tnano0.append(Tnano)
        r = r + 0.000001


    #above critical size
    l = len(r0)
    r1 = []
    Tliquidus = []
    Tsolidus = []
    for j in range(0,l,1):
        if r0[j]<rcritical[i]:
            j = j + 1
        else:
            r1.append(r0[j])
            Tliquidus.append(Tliquidus0[j])
            Tsolidus.append(Tsolidus0[j])

    #read experimental data
    filename = input("Please enter text name of melting data for nano-crystal (including file type, columns seperated by space) or press <enter> to continue: ")
    print('\n')
    rex = []
    Tex = []
    if filename != "":
        rex, Tex = np.loadtxt(filename, usecols=(0,1), dtype=float, unpack=True)
        plt.plot(rex, Tex, 'o', color='b')


    #plot nano-phase diagram
    plt.xlabel('Radius, nm')
    plt.xlabel('Radius, nm')
    plt.ylabel('Temperature, K')
    plt.plot(rnano,Tnano0, linewidth=2.0)
    plt.plot(r1,Tliquidus, linewidth=2.0)
    plt.plot(r1,Tsolidus, linewidth=2.0)
    plt.xlim(0, 100)
    plt.ylim(0, Tm+100)
    plt.savefig('Nano-phase diagram_%i.ps'% (i+1), format='ps',dpi=1000)
    plt.show()



#minimum method
if N == 2:
    print("**********Plot nano-phase diagram by Gibbs energy minimum method**********")
    print("**********Please wait**********")
    
    rcritical_min_index = rcritical.index(min(rcritical))
    rcritical_max_index = rcritical.index(max(rcritical))

    rcritical_sort = np.sort(rcritical)
    n = 30
    
    rcritical_min = rcritical_sort[0]
    rcritical_max = rcritical_sort[1]

    #from maximum critical size to 100 nm
    r0 = []
    Tliqexact0 = []
    Tsolexact0 = []
    r = rcritical_max
    while r <= 100:
        r0.append(r)
        #liquidus
        Tliquidus0 = []
        Gsl0 = []
        #solidus
        Tsolidus0 = []
        Gs0 = []
        for i in range(0,N,1):
            y0 = 10**(-6)
            y1 = 1
            y = y0
            Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
            Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
            Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
            Y0 = Gsl-Gl

            y = y1
            Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
            Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
            Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
            Y1 = Gsl-Gl

            if Y0*Y1<0:
                for t in range(0,n,1):
                    y = (y0 + y1)/2
                    Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                    Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                    Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                    Y = Gsl - Gl

                    if Y*Y0>0:
                        y0 = y
                        Y0 = Y
                    else:
                        y1 = y
                        Y1 = Y
                yliquidus = (y0 + y1)/2
            Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
            Tliquidus0.append(Tliquidus)
            Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
            Gsl0.append(Gsl)

            ysolidus = (1-d0[i]/r)**3
            Tsolidus = Tm-2*1000*Vm/(r*DeltaSm*(ysolidus**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1))
            Tsolidus0.append(Tsolidus)
            Gs = 10**(-3)*(1-ysolidus)*r/(3*Vm)*DeltaSm*(Tm-Tsolidus)+Sigma_lg+ysolidus**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1)))
            Gs0.append(Gs)
        minindex = Gsl0.index(min(Gsl0))
        Tliqexact = Tliquidus0[minindex]
        Tliqexact0.append(Tliqexact)

        minindex = Gs0.index(min(Gs0))
        Tsolexact = Tsolidus0[minindex]
        Tsolexact0.append(Tsolexact)

        r = r + 0.001
    #plot nano-phase diagram
    plt.xlabel('Radius, nm')
    plt.ylabel('Temperature, K')
    plt.plot(r0,Tliqexact0, linewidth=2.0)
    plt.plot(r0,Tsolexact0, linewidth=2.0)
    plt.xlim(0, 100)
    plt.ylim(0, Tm+100)


    #from minimum critical size to maximum critical size
    r0 = []
    Tliqexact0 = []
    Tsolexact0 = []
    r = rcritical_min
    while r <= rcritical_max:
        r0.append(r)
        #liquidus
        Tliquidus0 = []
        Gsl0 = []
        #solidus
        Tsolidus0 = []
        Gs0 = []
        for i in range(0,N,1):
            if i == rcritical_min_index:
                y0 = 10**(-6)
                y1 = 1
                y = y0
                Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                Y0 = Gsl-Gl

                y = y1
                Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                Y1 = Gsl-Gl

                if Y0*Y1<0:
                    for t in range(0,n,1):
                        y = (y0 + y1)/2
                        Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                        Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                        Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                        Y = Gsl - Gl

                        if Y*Y0>0:
                            y0 = y
                            Y0 = Y
                        else:
                            y1 = y
                            Y1 = Y
                    yliquidus = (y0 + y1)/2
                Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                Tliquidus0.append(Tliquidus)
                Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                Gsl0.append(Gsl)

                ysolidus = (1-d0[i]/r)**3
                Tsolidus = Tm-2*1000*Vm/(r*DeltaSm*(ysolidus**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1))
                Tsolidus0.append(Tsolidus)
                Gs = 10**(-3)*(1-ysolidus)*r/(3*Vm)*DeltaSm*(Tm-Tsolidus)+Sigma_lg+ysolidus**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1)))
                Gs0.append(Gs)
            else:
                ycritical = (1-d0[i]/rcritical[i])**3
                Tcritical = Tm-2*1000*Vm/(rcritical[i]*DeltaSm*(ycritical**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1))
                fcor = (Tm-Tcritical)*(rcritical[i]*DeltaSm)/(3*1000*Vm*(Sigma_sg[i]-Sigma_lg))
                Tnano = Tm-3*1000*Vm/(r*DeltaSm)*(Sigma_sg[i]-Sigma_lg)*fcor

                Tliquidus0.append(Tnano)
                Tsolidus0.append(Tnano)

                Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tnano)+Sigma_lg
                Gsl0.append(Gl)
                Gs0.append(Gl)

        minindex = Gsl0.index(min(Gsl0))
        Tliqexact = Tliquidus0[minindex]
        Tliqexact0.append(Tliqexact)

        minindex = Gs0.index(min(Gs0))
        Tsolexact = Tsolidus0[minindex]
        Tsolexact0.append(Tsolexact)

        r = r + 0.001

    #plot nano-phase diagram
    plt.plot(r0,Tliqexact0, linewidth=2.0)
    plt.plot(r0,Tsolexact0, linewidth=2.0)

    #below minimum critical size
    r0 = []
    Tliqexact0 = []
    Tsolexact0 = []
    r = 0.001
    while r <= rcritical_min:
        r0.append(r)
        #liquidus
        Tliquidus0 = []
        Gsl0 = []
        #solidus
        Tsolidus0 = []
        Gs0 = []
        for i in range(0,N,1):
            ycritical = (1-d0[i]/rcritical[i])**3
            Tcritical = Tm-2*1000*Vm/(rcritical[i]*DeltaSm*(ycritical**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1))
            fcor = (Tm-Tcritical)*(rcritical[i]*DeltaSm)/(3*1000*Vm*(Sigma_sg[i]-Sigma_lg))
            Tnano = Tm-3*1000*Vm/(r*DeltaSm)*(Sigma_sg[i]-Sigma_lg)*fcor

            Tliquidus0.append(Tnano)
            Tsolidus0.append(Tnano)

            Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tnano)+Sigma_lg
            Gsl0.append(Gl)
            Gs0.append(Gl)

        minindex = Gsl0.index(min(Gsl0))
        Tliqexact = Tliquidus0[minindex]
        Tliqexact0.append(Tliqexact)

        minindex = Gs0.index(min(Gs0))
        Tsolexact = Tsolidus0[minindex]
        Tsolexact0.append(Tsolexact)

        r = r + 0.001
    #plot nano-phase diagram
    plt.plot(r0,Tliqexact0, linewidth=2.0)
    plt.plot(r0,Tsolexact0, linewidth=2.0)

    #read experimental data
    print('\n')
    filename = input("Please enter text name of melting data for nano-crystal (including file type, columns seperated by space) or press <enter> to continue: ")
    rex = []
    Tex = []
    if filename != "":
        rex, Tex = np.loadtxt(filename, usecols=(0,1), dtype=float, unpack=True)
        plt.plot(rex, Tex, 'o', color='b')
    plt.savefig('Nano-phase diagram.ps', format='ps',dpi=1000)
    plt.show()

    print('\n')
    input("Press <enter> to close program")
    
elif N == 3:
    print("**********Plot nano-phase diagram by Gibbs energy minimum method**********")
    print("**********Please wait**********")

    rcritical_min_index = rcritical.index(min(rcritical))
    rcritical_max_index = rcritical.index(max(rcritical))

    rcritical_sort = np.sort(rcritical)
    n = 30 #increase accuracy
    
    rcritical_min = rcritical_sort[0]
    rcritical_mid = rcritical_sort[1]
    rcritical_max = rcritical_sort[2]

    #from maximum critical size to 100 nm
    r0 = []
    Tliqexact0 = []
    Tsolexact0 = []
    r = rcritical_max
    while r<=100:
        r0.append(r)
        #liquidus
        Tliquidus0 = []
        Gsl0 = []
        #solidus
        Tsolidus0 = []
        Gs0 = []
        for i in range(0,N,1):
            y0 = 10**(-6)
            y1 = 1
            y = y0
            Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
            Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
            Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
            Y0 = Gsl-Gl

            y = y1
            Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
            Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
            Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
            Y1 = Gsl-Gl

            if Y0*Y1<0:
                for t in range(0,n,1):
                    y = (y0 + y1)/2
                    Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                    Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                    Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                    Y = Gsl - Gl

                    if Y*Y0>0:
                        y0 = y
                        Y0 = Y
                    else:
                        y1 = y
                        Y1 = Y
                yliquidus = (y0 + y1)/2
            Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
            Tliquidus0.append(Tliquidus)
            Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
            Gsl0.append(Gsl)

            ysolidus = (1-d0[i]/r)**3
            Tsolidus = Tm-2*1000*Vm/(r*DeltaSm*(ysolidus**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1))
            Tsolidus0.append(Tsolidus)
            Gs = 10**(-3)*(1-ysolidus)*r/(3*Vm)*DeltaSm*(Tm-Tsolidus)+Sigma_lg+ysolidus**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1)))
            Gs0.append(Gs)
        minindex = Gsl0.index(min(Gsl0))
        Tliqexact = Tliquidus0[minindex]
        Tliqexact0.append(Tliqexact)

        minindex = Gs0.index(min(Gs0))
        Tsolexact = Tsolidus0[minindex]
        Tsolexact0.append(Tsolexact)

        r = r + 0.001
    #plot nano-phase diagram
    plt.xlabel('Radius, nm')
    plt.ylabel('Temperature, K')
    plt.plot(r0,Tliqexact0, linewidth=2.0)
    plt.plot(r0,Tsolexact0, linewidth=2.0)
    plt.xlim(0, 100)
    plt.ylim(0, Tm+100)

    
    #from middle critical size to maximum critical size
    r0 = []
    Tliqexact0 = []
    Tsolexact0 = []
    r = rcritical_mid
    while r<=rcritical_max:
        r0.append(r)
        #liquidus
        Tliquidus0 = []
        Gsl0 = []
        #solidus
        Tsolidus0 = []
        Gs0 = []
        for i in range(0,N,1):
            if i == rcritical_max_index:
                ycritical = (1-d0[i]/rcritical[i])**3
                Tcritical = Tm-2*1000*Vm/(rcritical[i]*DeltaSm*(ycritical**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1))
                fcor = (Tm-Tcritical)*(rcritical[i]*DeltaSm)/(3*1000*Vm*(Sigma_sg[i]-Sigma_lg))
                Tnano = Tm-3*1000*Vm/(r*DeltaSm)*(Sigma_sg[i]-Sigma_lg)*fcor

                Tliquidus0.append(Tnano)
                Tsolidus0.append(Tnano)

                Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tnano)+Sigma_lg
                Gsl0.append(Gl)
                Gs0.append(Gl)
            else:
                y0 = 10**(-6)
                y1 = 1
                y = y0
                Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                Y0 = Gsl-Gl

                y = y1
                Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                Y1 = Gsl-Gl

                if Y0*Y1<0:
                    for t in range(0,n,1):
                        y = (y0 + y1)/2
                        Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                        Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                        Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                        Y = Gsl - Gl

                        if Y*Y0>0:
                            y0 = y
                            Y0 = Y
                        else:
                            y1 = y
                            Y1 = Y
                    yliquidus = (y0 + y1)/2
                Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                Tliquidus0.append(Tliquidus)
                Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                Gsl0.append(Gsl)

                ysolidus = (1-d0[i]/r)**3
                Tsolidus = Tm-2*1000*Vm/(r*DeltaSm*(ysolidus**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1))
                Tsolidus0.append(Tsolidus)
                Gs = 10**(-3)*(1-ysolidus)*r/(3*Vm)*DeltaSm*(Tm-Tsolidus)+Sigma_lg+ysolidus**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1)))
                Gs0.append(Gs)


        minindex = Gsl0.index(min(Gsl0))
        Tliqexact = Tliquidus0[minindex]
        Tliqexact0.append(Tliqexact)

        minindex = Gs0.index(min(Gs0))
        Tsolexact = Tsolidus0[minindex]
        Tsolexact0.append(Tsolexact)

        r = r + 0.001

    #plot nano-phase diagram
    plt.plot(r0,Tliqexact0, linewidth=2.0)
    plt.plot(r0,Tsolexact0, linewidth=2.0)


    #from minimum critical size to middle critical size
    r0 = []
    Tliqexact0 = []
    Tsolexact0 = []
    r = rcritical_min
    while r<=rcritical_mid:
        r0.append(r)
        #liquidus
        Tliquidus0 = []
        Gsl0 = []
        #solidus
        Tsolidus0 = []
        Gs0 = []
        for i in range(0,N,1):
            if i == rcritical_min_index:
                y0 = 10**(-6)
                y1 = 1
                y = y0
                Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                Y0 = Gsl-Gl

                y = y1
                Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                Y1 = Gsl-Gl

                if Y0*Y1<0:
                    for t in range(0,n,1):
                        y = (y0 + y1)/2
                        Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                        Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                        Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg
                        Y = Gsl - Gl

                        if Y*Y0>0:
                            y0 = y
                            Y0 = Y
                        else:
                            y1 = y
                            Y1 = Y
                    yliquidus = (y0 + y1)/2
                Tliquidus = Tm-2*1000*Vm/(r*DeltaSm*(y**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1))
                Tliquidus0.append(Tliquidus)
                Gsl = 10**(-3)*(1-y)*r/(3*Vm)*DeltaSm*(Tm-Tliquidus)+Sigma_lg+y**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(y**(1/3)-1)))
                Gsl0.append(Gsl)

                ysolidus = (1-d0[i]/r)**3
                Tsolidus = Tm-2*1000*Vm/(r*DeltaSm*(ysolidus**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1))
                Tsolidus0.append(Tsolidus)
                Gs = 10**(-3)*(1-ysolidus)*r/(3*Vm)*DeltaSm*(Tm-Tsolidus)+Sigma_lg+ysolidus**(2/3)*(Sigma_sl+DeltaSigma[i]*math.exp(r/xi[i]*(ysolidus**(1/3)-1)))
                Gs0.append(Gs)
            else:
                ycritical = (1-d0[i]/rcritical[i])**3
                Tcritical = Tm-2*1000*Vm/(rcritical[i]*DeltaSm*(ycritical**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1))
                fcor = (Tm-Tcritical)*(rcritical[i]*DeltaSm)/(3*1000*Vm*(Sigma_sg[i]-Sigma_lg))
                Tnano = Tm-3*1000*Vm/(r*DeltaSm)*(Sigma_sg[i]-Sigma_lg)*fcor

                Tliquidus0.append(Tnano)
                Tsolidus0.append(Tnano)

                Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tnano)+Sigma_lg
                Gsl0.append(Gl)
                Gs0.append(Gl)

        minindex = Gsl0.index(min(Gsl0))
        Tliqexact = Tliquidus0[minindex]
        Tliqexact0.append(Tliqexact)

        minindex = Gs0.index(min(Gs0))
        Tsolexact = Tsolidus0[minindex]
        Tsolexact0.append(Tsolexact)

        r = r + 0.001

    #plot nano-phase diagram
    plt.plot(r0,Tliqexact0, linewidth=2.0)
    plt.plot(r0,Tsolexact0, linewidth=2.0)

    #below minimum critical size
    r0 = []
    Tliqexact0 = []
    Tsolexact0 = []
    r = 0.001
    while r<=rcritical_min:
        r0.append(r)
        #liquidus
        Tliquidus0 = []
        Gsl0 = []
        #solidus
        Tsolidus0 = []
        Gs0 = []
        for i in range(0,N,1):
            ycritical = (1-d0[i]/rcritical[i])**3
            Tcritical = Tm-2*1000*Vm/(rcritical[i]*DeltaSm*(ycritical**(1/3)))*(Sigma_sl+DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1)))-1000*Vm/(xi[i]*DeltaSm)*DeltaSigma[i]*math.exp(rcritical[i]/xi[i]*(ycritical**(1/3)-1))
            fcor = (Tm-Tcritical)*(rcritical[i]*DeltaSm)/(3*1000*Vm*(Sigma_sg[i]-Sigma_lg))
            Tnano = Tm-3*1000*Vm/(r*DeltaSm)*(Sigma_sg[i]-Sigma_lg)*fcor

            Tliquidus0.append(Tnano)
            Tsolidus0.append(Tnano)

            Gl = 10**(-3)*r/(3*Vm)*DeltaSm*(Tm-Tnano)+Sigma_lg
            Gsl0.append(Gl)
            Gs0.append(Gl)

        minindex = Gsl0.index(min(Gsl0))
        Tliqexact = Tliquidus0[minindex]
        Tliqexact0.append(Tliqexact)

        minindex = Gs0.index(min(Gs0))
        Tsolexact = Tsolidus0[minindex]
        Tsolexact0.append(Tsolexact)

        r = r + 0.001
    #plot nano-phase diagram
    plt.plot(r0,Tliqexact0, linewidth=2.0)
    plt.plot(r0,Tsolexact0, linewidth=2.0)

    #read experimental data
    print('\n')
    filename = input("Please enter text name of melting data for nano-crystal (including file type, columns seperated by space) or press <enter> to continue: ")
    rex = []
    Tex = []
    if filename != "":
        rex, Tex = np.loadtxt(filename, usecols=(0,1), dtype=float, unpack=True)
        plt.plot(rex, Tex, 'o', color='b')
    plt.savefig('Nano-phase diagram.ps', format='ps',dpi=1000)
    plt.show()

    print('\n')
    input("Press <enter> to close program")

elif N == 1:
    print('\n')
    input("Press <enter> to close program")

else:
    print("Sorry, it is just for one surface or two and three surfaces")
    print('\n')
    input("Press <enter> to close program")
