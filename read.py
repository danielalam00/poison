from calendar import day_name
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

def scale_read(nPin, solutions):

    time = [0, 1, 2, 4, 8, 30, 120, 240, 360, 480, 960, 1400, 2200, 3000, 3140, 3280, 3420, 3560, 3700]
    time_arr = np.array(time)
    kw = 'Transport: k-eff ='
    #open output file
    name_map = np.array(solutions[1:],int)
    name = ''.join(map(str,name_map))
    Nama_file = str(solutions[0])+"_"+name+".out"
    f = open(Nama_file, "r")
    print(Nama_file)
    kval = []
    if nPin == 10 :
        k_asli = 1.41074
    elif nPin == 11 :
        k_asli = 1.48980
    elif nPin == 12 :
        k_asli = 1.46944
    elif nPin == 13 :
        k_asli = 1.45209
    elif nPin == 14 :
        k_asli = 1.48336
    elif nPin == 15 :
        k_asli = 1.56876
    elif nPin == 16 :
        k_asli = 1.50510
    else :
        k_asli = 1.45475
    try:
        for line in f:
            if kw in line:
                keff = line.split()
                kval.append(float(keff[3]))
        k_arr = np.array(kval)
        k_awal = k_arr[0]
        #del_k = k_awal-1
        ord = k_arr[-5:]
        if k_awal > 1:
            del_k = k_asli - k_awal
        else :
            del_k = k_awal - 1
        print('file dapat dibaca')
    except:
        knilai = [0.2, 0.3, 0.4, 0.5, 0.6]
        knilai_arr= np.array(knilai)
        power =[200, 500]
        del_k = -100
        k_awal = 0
        ord = knilai_arr[-5:]
        print("file kosong ")
    #print(ord)
    abs = time_arr[-5:]
    if len(ord) ==0 :
        knilai = [0.2, 0.3, 0.4, 0.5, 0.6]
        knilai_arr= np.array(knilai)
        power =[200, 500]
        del_k = -100
        k_awal = 0
        ord = knilai_arr[-5:]
        print("file kosong ")
    if k_awal < 1:
        del_k = k_awal - 1
        day = 0
        print("to much poison ")
    else:
        ord_log = np.log(ord)
        curve_fit = np.polyfit(abs, ord_log,1)
        c1 = curve_fit[0]
        c2 = np.exp(curve_fit[1])
        func = c2*np.exp(c1*abs)
        day = -1*np.log(c2)/c1
        print("all good ")
        cek_k = k_arr[:-1]
        if any(i <= 1 for i in cek_k):
            index_subs = np.empty((0, 0))
            for i in k_arr :
                if i <= 1 :
                    ind_temp = np.where(k_arr == i)
                    index_subs = np.append(index_subs, ind_temp) 
            ind_first = index_subs[0]
            ind_last = index_subs[-1]
            selisih = ind_last - ind_first
            if len(index_subs) != (selisih + 1):
                day = -1000
                del_k = k_awal - 1
                print("subcritical in the middle ")
            elif len(index_subs) != 0 and ind_last != 18 :
                day = -1000
                del_k = k_awal - 1
                print("subcritical in the middle ")
    fitt = (day * 0.6 ) +(del_k *4000)
    print('fitness:', fitt,   'day:', day,   'delta k:', del_k)
    return [fitt]
