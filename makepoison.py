from dataclasses import replace
from cmath import pi
import os
import math
import random as rnd
import numpy as np
import itertools as it


def scale_write(nPin, solutions):
    gap = 0.009
    clad1 = 0.03
    clad2 = 0.057
    daya = 3.5
    H_f = 200
    wf_u1 = 92.3
    wf_u2 = 94.3
    dens_usi = 11.468
    dens_un = 13
    
    #all the parameter bassed on the previous result
    if nPin == 10 :
        fuel_opt = 'USi'
        clad_opt = 'FeCrAl'
        PD = 1.34
        PinRad = 0.4033
        enr =0.077
    elif nPin == 11 :
        fuel_opt = 'USi'
        clad_opt = 'FeCrAl'
        PD = 1.6
        PinRad = 0.3707
        enr =0.0679
    elif nPin == 12 :
        fuel_opt = 'UN'
        clad_opt = 'FeCrAl'
        PD = 1.81
        PinRad = 0.3621
        enr = 0.0751
    elif nPin == 13 :
        fuel_opt = 'UN'
        clad_opt = 'SiC'
        PD = 1.56
        PinRad = 0.41
        enr =0.0717
    elif nPin == 14 :
        fuel_opt = 'USi'
        clad_opt = 'SiC'
        PD = 1.41
        PinRad = 0.5197
        enr = 0.0639
    elif nPin == 15 :
        fuel_opt = 'USi'
        clad_opt = 'SiC'
        PD = 1.72
        PinRad = 0.4642
        enr = 0.0653
    elif nPin == 16 :
        fuel_opt = 'USi'
        clad_opt = 'FeCrAl'
        PD = 1.62
        PinRad = 0.4543
        enr =0.0667
    else:
        fuel_opt = 'USi'
        clad_opt = 'FeCrAl'
        PD = 1.46
        PinRad = 0.4233
        enr = 0.0679
    
    if solutions[1] == 0: #choose the poison pin type
        readScale = open("InputHomogen.inp", "r")
    else:
        readScale = open("InputRing.inp", "r")
    #name should corelated with solution
    name_map = np.array(solutions[1:],int)
    name = ''.join(map(str,name_map))
    Nama_file = str(solutions[0])+"_"+name+".inp"
    # done!!!
    NewFile = open(Nama_file, "w")
  
    if fuel_opt == 'USi':
        dens = dens_usi
        wf_u = wf_u1
        comp_fuel = "comp c_usi   : WT Si=7.3 c_u=-100"
        mat_fuel = "mat  FUEL.1  : c_usi dens=11.468 temp=900"
        fuel_mat = "c_usi"
    else:
        dens = dens_un
        wf_u = wf_u2
        comp_fuel = "comp c_n     : WT 7015=50 7014=-100\ncomp c_un    : WT N=5.7 c_u=-100"
        mat_fuel = "mat  FUEL.1  : c_un dens=13 temp=900"
        fuel_mat = "c_un"
    
    power = 60000 / 3650
    gap1 = PinRad + gap

    if clad_opt == 'SiC':
        clad = clad2
        Clad_1 = gap1 + clad
        comp_clad = "comp sic     : WT Si=70.08 C=29.92"
        mat_clad = "mat  CLAD.1  : sic dens=2.58  temp=600"
        pin_F = ("pin 0 : "+str("%.4f" %PinRad)+"    "+str("%.4f" %gap1)+"  "+str("%.4f" %Clad_1))
    else:
        clad = clad1
        Clad_1 = gap1 + clad
        comp_clad = "comp fecral  : WT Fe=75 Cr=20 Al=5"
        mat_clad = "mat  CLAD.1  : fecral dens=7.1  temp=600"
        pin_F = ("pin 0 : "+str("%.4f" %PinRad)+"    "+str("%.4f" %gap1)+"  "+str("%.4f" %Clad_1))
    
    ppitch = (PD) * ((PinRad+gap+clad) * 2)
    enric_u = enr*100
    nxn = (str("%.0f" % nPin)+"x"+str("%.0f" % nPin))
    # Map of pins # need more adjustment bassed on the generated configuration
    nPinArray = nPin
    nPin_sym = math.ceil(nPinArray/2)
    map_solution = np.int_(solutions[2:])
    #array for pin map
    confpin = np.array(map_solution).reshape(nPin_sym,nPin_sym)
    pinMap = ""
    for i in range(nPin_sym):
        conf = ' '.join(map(str,confpin[i]))
        pinMap += " "+conf+"\n"
    #done!!!
    #poison
    if solutions[1] == 0: #choose the poison pin type #homogen
        poison_ratio = solutions[0] * 100
        dens_mix = ((1- solutions[0]) * dens) + ( solutions[0] * 7.41)
        comp_mix = ("comp u_gd     : WT GD2O3="+str("%.4f" %poison_ratio)+" "+str(fuel_mat)+"=-100 ")
        mat_mix = ("mat MIX.1   : u_gd dens="+str(dens_mix)+" temp=900")
        pin_P = ("pin 1 : "+str("%.4f" %PinRad)+"    "+str("%.4f" %gap1)+"  "+str("%.4f" %Clad_1))
        for line in readScale:
            NewFile.write(line.replace('nn', nxn).replace('jp', "%.0f" %nPin).replace('pp1', "%.4f" %ppitch).replace('u235', "%.4f" %enric_u)
                .replace('c_usi1', comp_fuel).replace('matf', mat_fuel).replace('comp_cl', comp_clad).replace('mat_cl', mat_clad).replace('comp_ugd', comp_mix).replace('mat_m', mat_mix)
                .replace('pinF', pin_F).replace('pinP', pin_P).replace('pow_1', "%.4f" %power).replace('pmap', pinMap))
    
    else: #ring
        poison_ratio = solutions[0]
        volume_fuel = math.pi * H_f * PinRad ** 2
        inner_fuel = math.sqrt(((1-poison_ratio)*volume_fuel)/(H_f * math.pi))
        pin_P = ("pin 1 : "+str("%.4f" %inner_fuel)+"  "+str("%.4f" %PinRad)+"    "+str("%.4f" %gap1)+"  "+str("%.4f" %Clad_1))
        for line in readScale:
            NewFile.write(line.replace('nn', nxn).replace('jp', "%.0f" %nPin).replace('pp1', "%.4f" %ppitch).replace('u235', "%.4f" %enric_u)
                .replace('c_usi1', comp_fuel).replace('matf', mat_fuel).replace('comp_cl', comp_clad).replace('mat_cl', mat_clad)
                .replace('pinF', pin_F).replace('pinP', pin_P).replace('pow_1', "%.4f" %power).replace('pmap', pinMap))
    
    #list nama file
    
    runFile = open("Nama_file", "a")
    runFile.writelines((Nama_file)+"\n")
    
    return #power

