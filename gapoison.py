#optimization of the burnable poison ratio and configuration for fuel assembly

import os
import numpy as np
import pandas as pd
import random
import math

from read import scale_read
#need to make sure the make file is the corenpondent configuration for a particular FA
from makepoison import scale_write

#manually set parameters
pins = 15      #coresponded to the total generated pins type
pin_maps = (math.ceil(pins/2))**2 #need to be adjust bassed on pins number
total_parameters = pin_maps + 2
bin_opt = [0, 1]   # binner option 
#clad_opt = 2   # 2 is FeCrAl and 1 is SiC

#parameters initialization
target_time = 3650 #days
c = 0
rankedsolutions = np.empty([0, (total_parameters + 1)]) #need to be adjust bassed on pins number
mutan = np.empty([0, total_parameters])  #need to be adjust bassed on pins number
solutions = np.empty([0, total_parameters]) #need to be adjust bassed on pins number
runscl = 'runScale'
pathin = 'inputs/genepool.csv'
exist = 0
mutation_rate =  0.5# need to be adjusted according to the requairenment
#mutation_degree = 0.02 # 2% of the binary gene mutated, as for the floating gane it's the degree of change
convergency = 0
population_number = 80
picked_population = int(population_number/2)
ratio_mutation = 0.1

#check if inputs folder exist
pathfol = './inputs'
in_exist =os.path.exists(pathfol)

if in_exist == False :
    os.system('mkdir inputs')

# there is some concern regarding the pin number which affect the generation of solution, need some opinion on the topics
# function
def fitness(w): # Fitness Function
    fit = w
    return fit
def data_skip(database, solutions): #Skip if data exist in Genepool Function
    global fitnessValue
    global fitnessValueExist
    fitnessValueExist = 0
    for i in range(len(database)):
        dtbase = database.iloc[i]
        newdat =np.array(dtbase[1:])
        check_arr = np.array_equal(newdat, solutions)
        if check_arr == 1:
            print("found the same data")
            fitnessValue = dtbase[0]
            fitnessValueExist = 1
            continue
    if fitnessValueExist == 1:
        return fitnessValue
    else:
        return 0

def dataRead(pathcsv):
    solList = np.empty([0, total_parameters])
    exist = 0
    database = 0
    if os.path.isfile(pathcsv) == 1: # Check whether the file of previous solutions exist
        exist = 1
        database = pd.read_csv(pathcsv, index_col=0) # read data from the gene pool
        solSorted = database.sort_values('0', ascending=False)
        del solSorted['0']
        solList = solSorted[:population_number].to_numpy() # Import the set of solutions 
        if solList.size < population_number*total_parameters: # this section aren't really a necessary part considering the gane pool file always written after 2 generation has been done calculated
            for s in range(population_number-(solutions.size/total_parameters)): # Initialize Initial Solutions with the different pin number, need an if function that produced an array with matching column but still has matching row
                 solList_temp = np.empty([0, total_parameters] ) #temporary array to accomodated different pins number
                 for f in range(total_parameters ): 
                    if f == 0:
                        solList_temp = np.append(solList_temp, round(random.uniform(0.0001, 0.1),4)) #generate poison ratio
                    else:
                        solList_temp = np.append(solList_temp, int(random.choice(bin_opt))) # generate pin type and configuration
                 solList = np.append(solList, [solList_temp], axis = 0)
        
    else:
        for s in range(population_number): # Initialize Initial Solutions
             solList_temp = np.empty([0, total_parameters] )
             for f in range(total_parameters ):
                 if f == 0:
                     solList_temp = np.append(solList_temp, round(random.uniform(0.0001, 0.1),4))#generate poison ratio
                 else:
                     solList_temp = np.append(solList_temp, int(random.choice(bin_opt)))# generate pin type and configuration
             solList = np.append(solList, [solList_temp], axis = 0)
    return [solList, exist, database]

def ComSol (lengthnewsol, tem_solutions):
    prev_bestsolutions = np.empty([0, total_parameters ])
    bestsolutions = tem_solutions[:picked_population]
    for k in bestsolutions: # Add the best solutions in 1 generation to generational run pool
        prev_bestsolutions = np.append(prev_bestsolutions, np.array([k]), axis=0)
        prev_bestsolutions = prev_bestsolutions[prev_bestsolutions[:, 0].argsort()] #Sort
        prev_bestsolutions = prev_bestsolutions[::-1] #Reverse
    
    # mutation and crossing section need more thought on this specially on the mutation degree,
    mutan = np.empty([0, total_parameters]) # need some adjustment
    elements = np.empty([0, total_parameters])
    
    for s in prev_bestsolutions[:picked_population]: # Select the inputs from best solutions
        elements_temp = np.empty([0, total_parameters] )
        for f in range(total_parameters ): 
            elements_temp = np.append(elements_temp, s[f])
        elements = np.append(elements, [elements_temp], axis=0)
    parents = np.unique(elements, axis=0)
    print("parent in comsol")
    print(parents)
    newGen = np.empty([0, total_parameters])
    for _ in range(lengthnewsol): # Crossing
        break_point = int(random.uniform(2, total_parameters))
        newGen_temp = np.empty(0)
        a1 = random.choice(np.unique(parents, axis = 0))
        a2 = random.choice(np.unique(parents, axis = 0))
        b1 = a1[0:(break_point)]
        b2 = a2[(break_point):]
        newGen_temp = np.append(newGen_temp, b1)
        newGen_temp = np.append(newGen_temp, b2)
        newGen = np.append(newGen, [newGen_temp], axis=0)
        print("crossing penambahan")
        print(newGen)
        #crossing mechanism using the first option
       #mutation mechanism based on the type of the variable
    for l in range(int(lengthnewsol)):
        rand = random.randint(0, (lengthnewsol-1))
        o1 = round(newGen[rand, 0] * random.uniform(0.95, 1.05), 4) # Poison Ratio
        while o1 > 0.1 or o1 < 0.001:
            o1 = round(newGen[rand, 0] * random.uniform(0.95, 1.05), 4) # Poison Ratio
        newGen[rand, 0] = o1
        #biner bassed variale mutation
        flip_number = random.randint(1, 3)
        for _ in range(int(flip_number)):
            rand_biner1 = random.randint(2, (total_parameters-1))
            o2 = abs(newGen[rand, rand_biner1] - 1) # configuration
            newGen[rand, rand_biner1] = o2
            print("fliping complete")
            print(newGen)
    for m in range(int(lengthnewsol)):
        rand = random.randint(0, (lengthnewsol-1))
        o3 = abs(newGen[rand, 1] - 1) # type
        newGen[rand, 1] = o3
        print("ratio complete")
        print(newGen)
    tem_sol = newGen
    return tem_sol

# this section need more adjustment 
#Generate Solutions
while convergency == 0:
    tempOutFunc = dataRead(pathin)# calling dataRead fuction with output [solution, exist, database]
    solutions = tempOutFunc[0]
    exist = tempOutFunc[1]
    database = tempOutFunc[2]
    
    prev_bestsolutions = np.empty([0, total_parameters + 1])
    for gen in range(2): # Generation Number
        rankedsolutions = np.empty([0, total_parameters + 1])
        skippedData = np.empty([0, total_parameters + 1])
        temp_solutions = np.empty([0, total_parameters + 1])
        #mutationIdentifier = np.empty((0, 1))
        if exist == 1:
            tempData = np.empty([0, total_parameters + 1])
            for i in range(len(solutions)): # Skip Data if Exist
                temp1 = np.empty([0, 1])
                calculatedFitness = data_skip(database, solutions[i])# calling data_skip function with output fitness value if solution match database
                temp1 = np.append(temp1, calculatedFitness)
                temp1 = np.append(temp1, solutions[i])
                tempData = np.append(tempData, [temp1], axis = 0)

            print("Temporary Data")
            print(tempData)
            
            skippedData = tempData[np.where(tempData[:,0] != 0)] #Concantenate Skipped Data
            temp_solutions = solutions[np.where(tempData[:,0] == 0)] #Concantenate solutions
            solutions = np.unique(temp_solutions, axis=0)
            temcomsol = 0
            while temcomsol == 0:
                add_solutions = np.empty([0, total_parameters])
                if (len(solutions) + len(skippedData)) != population_number :
                    new_solutions = np.empty([0, total_parameters])
                    tempData2 = np.empty([0, total_parameters + 1])
                    solfunc = dataRead(pathin)
                    new_sol = tempOutFunc[0]
                    lengthnewsol = population_number - (len(solutions) + len(skippedData))
                    comsolfunc = ComSol(lengthnewsol, new_sol)
                    new_solutions = np.append(new_solutions, comsolfunc, axis=0)
                    for i in range(len(new_solutions)): # Skip Data if Exist
                        temp2 = np.empty([0, 1])
                        calculatedFitness2 = data_skip(database, new_solutions[i])# calling data_skip function with output fitness value if solution match database
                        temp2 = np.append(temp2, calculatedFitness2)
                        temp2 = np.append(temp2, new_solutions[i])
                        tempData2 = np.append(tempData2, [temp2], axis = 0)
                    skippedData1 = tempData2[np.where(tempData2[:,0] != 0)] #Concantenate Skipped Data
                    temp_solutions1 = new_solutions[np.where(tempData2[:,0] == 0)] #Concantenate solutions
                    add_solutions = np.unique(temp_solutions1, axis=0)
                    solutions = np.append(solutions, add_solutions, axis = 0)
                    solutions = np.unique(solutions, axis=0)
                    temcomsol = 0
                else :
                    temcomsol = 1
            print("Skipped Data after concantenated")
            print(skippedData)                             
                                             
            print("Solution after concantenated")
            print(solutions)
        
        if solutions.size != 0:
            polaris_out = np.empty([0, 1])
            
            if os.path.isfile(runscl) == 1: # Check whether the file of previous solutions exist
                os.system("rm runScale")
                os.system("rm Nama_file")
            else:
                print("memulai perhitungan dari awal")
                
            for i in range(len(solutions)):
                create = scale_write(pins, solutions[i]) # concern, the name of the file need to considered 
                #create runScale

            text_file = open("Nama_file", "r")
            lines = text_file.readlines()
            lines = [x.replace('\n','') for x in lines]

            runFile = open("runScale", "a")
            a=[7, 15, 23, 31, 39, 47, 55, 63, 71, 79 ]
            b=[8, 16, 24, 32, 40, 48, 56, 64, 72, 80 ]
            c=(len(solutions) - 1)
            for i in range(len(solutions)):
                if i in a:
                        runFile.writelines("/usr/local/SCALE-6.2.4/bin/scalerte "+(lines[i])+"\n")
                        runFile.writelines("wait \n")
                elif i == c:
                        runFile.writelines("/usr/local/SCALE-6.2.4/bin/scalerte "+(lines[i])+"\n")
                else:
                        runFile.writelines("/usr/local/SCALE-6.2.4/bin/scalerte "+(lines[i])+" &\n")
            if len(solutions) not in b :
                runFile.writelines("wait \n")
            runFile.close()
            print("Segment B2-4")
            
            os.system("chmod u+x runScale")
            #runScale >> isinya perintah untuk run semua case

            os.system("./runScale")

            #Read the Outputs from SCALE Polaris
            
            for i in range(len(solutions)):
                 readOut = scale_read(pins, solutions[i])
                 polaris_out = np.append(polaris_out, readOut)
                 #print(polaris_out)
                 #print(scale_read)
            
            
            rankedsolutions = np.empty((0, (total_parameters + 1)))
            for index, s in enumerate(solutions, start=0):
                rankedsolutions_temp = np.empty([0, (total_parameters + 1)])
                f_out = fitness(polaris_out[index])
                for f in range(total_parameters + 1): 
                    if f == 0:
                        rankedsolutions_temp = np.append(rankedsolutions_temp, f_out)
                    else:
                        rankedsolutions_temp = np.append(rankedsolutions_temp, s[f-1])
                rankedsolutions = np.append(rankedsolutions, [rankedsolutions_temp], axis = 0)
                c+=1 # Solution Number Tracker

            print("Rankedsolutions before append")
            print(rankedsolutions)
        
        if exist == 1:  
            rankedsolutions = np.append(rankedsolutions, skippedData, axis = 0)
        
        print("Rankedsolutions after append")
        print(rankedsolutions)
        
        rankedsolutions = rankedsolutions[rankedsolutions[:, 0].argsort()] #Sort
        rankedsolutions = rankedsolutions[::-1] #Reverse
        
        #print("=== Gen {"+str(genname)+"} best solutions ===")
        print("=== best solutions for this generation ===")
        print(rankedsolutions[0])
    
        bestsolutions = rankedsolutions[:picked_population]
        print("bestsolutions")
        print(bestsolutions)
        for k in bestsolutions: # Add the best solutions in 1 generation to generational run pool
            prev_bestsolutions = np.append(prev_bestsolutions, np.array([k]), axis=0)
            prev_bestsolutions = prev_bestsolutions[prev_bestsolutions[:, 0].argsort()] #Sort
            prev_bestsolutions = prev_bestsolutions[::-1] #Reverse
        
        print(" top solutions :")
        print(rankedsolutions[:picked_population])
        
        # mutation and crossing section need more thought on this specially on the mutation degree,
        mutan = np.empty([0, total_parameters]) # need some adjustment
        elements = np.empty([0, total_parameters])
        
        for s in prev_bestsolutions[:picked_population]: # Select the inputs from best solutions
            elements_temp = np.empty([0, total_parameters] )
            for f in range(total_parameters +1 ): 
                if f != 0:
                    elements_temp = np.append(elements_temp, s[f])
            elements = np.append(elements, [elements_temp], axis=0)
        parents = np.unique(elements, axis=0)
        print("parent")
        print(parents)
        newGen = np.empty([0, total_parameters])
        for _ in range(population_number): # Crossing
            break_point = int(random.uniform(2, total_parameters))
            newGen_temp = np.empty(0)
            a1 = random.choice(np.unique(parents, axis = 0))
            a2 = random.choice(np.unique(parents, axis = 0))
            b1 = a1[0:(break_point)]
            b2 = a2[(break_point):]
            newGen_temp = np.append(newGen_temp, b1)
            newGen_temp = np.append(newGen_temp, b2)
            newGen = np.append(newGen, [newGen_temp], axis=0)
        print("crossing")
        print(newGen)
            #crossing mechanism using the first option
           #mutation mechanism based on the type of the variable
        for t in range(int(np.rint(population_number * mutation_rate))):
            rand = random.randint(0, (population_number-1))
            o1 = round(newGen[rand, 0] * random.uniform(0.95, 1.05), 4) # Poison Ratio
            while o1 > 0.1 or o1 < 0.001:
                o1 = round(newGen[rand, 0] * random.uniform(0.95, 1.05), 4) # Poison Ratio
            newGen[rand, 0] = o1
            #biner bassed variale mutation
            flip_number = random.randint(1, 3)
            for _ in range(flip_number):
                rand_biner1 = random.randint(2, (total_parameters-1))
                o2 = abs(newGen[rand, rand_biner1] - 1) # configuration
                newGen[rand, rand_biner1] = o2
        print("flipping")
        print(newGen)
        for p in range(int(population_number * ratio_mutation)):
            rand = random.randint(0, (population_number-1))
            o3 = abs(newGen[rand, 1] - 1) # type
            newGen[rand, 1] = o3
        print("ratio")
        print(newGen)
        print("New Generations")
        print(newGen)

        solutions = newGen

    pd_pb = pd.DataFrame(data=prev_bestsolutions)
    if exist == 1:
        pd_pb.to_csv(pathin, mode = 'a', header=False)
    else:
        pd_pb.to_csv(pathin, mode = 'a', header=True)
        
    convergence_data = prev_bestsolutions
    avg_data = np.average(convergence_data[:,0])
    top1 = convergence_data[:,0]
    error = np.absolute((top1-avg_data)/top1)
    error_t = np.array([error]).T
    error_t = error_t.T
    
    pd_c = pd.DataFrame(data=error_t)
    if os.path.isfile('inputs/convergence.csv') == 1:
        pd_c.to_csv('inputs/convergence.csv', mode = 'a', header=False)
    else:
        pd_c.to_csv('inputs/convergence.csv', mode = 'a', header=True)
        
    if np.average(error_t) < 0.00001:
        convergency = 1
