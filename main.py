#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 13:39:44 2024

@author: luismoncayo
"""
import sys
from Gurobi_Solution.Implementation_Gurobi import GUROBI_Sol_Class
from QuantumAnnealing_Solution.Implementation_Quantum import Quantum_Sol_Class
# Instance to be solved '/Users/luismoncayo/Dropbox/Python/Line_Balancing/Instances_Data/very large data set_n=1000/instance_n=1000_125.alb'
#--------------------------------------------------------------------
path_to_instance = "very large data set_n=1000/instance_n=1000_125.alb"
#--------------------------------------------------------------------

parts = path_to_instance.split("instance", 1)  # Split into two parts, at the first occurrence of 'fox'
model_name = "None"
if len(parts) > 1:
    result = parts[1].strip()  # Get the part after the keyword
    after_dot = result.split('.', 1)[0]
    #print("instance"+after_dot)  # Output: "jumps over the lazy dog"
    model_name = "instance"+after_dot
else:
    print("instance word is not found")

with open("Instances_Data/"+path_to_instance, 'r') as f:
    lines = [line.rstrip('\n') for line in f]
    #print(lines)
f.close()

number_tasks = int(lines[lines.index('<number of tasks>')+1])
cycle_time =  int(lines[lines.index('<cycle time>')+1])

index_task_time = lines.index('<task times>')
tasks_num = []
tasks_times = []
for i in range(index_task_time+1, index_task_time+number_tasks+1):
    number = int(lines[i].split(" ")[0])
    time = int(lines[i].split(" ")[1])
    tasks_num.append(number)
    tasks_times.append(time)

index_precedence = lines.index('<precedence relations>')
precedence = []
k = index_precedence+1
while len(lines[k]) != 0:
    from_t = int(lines[k].split(",")[0])
    to_t = int(lines[k].split(",")[1])
    precedence.append((from_t, to_t))
    k+=1

with open("Outputs_Gurobi/"+model_name+".txt", "w") as f:
    # Redirect standard output to the file
    sys.stdout = f
    
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("-----> Results for the ", model_name)
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    
    ##### ##### Set the number of stations         ##### #####
    number_stations = 135
    ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

    # Test the data ---------------------------------------------------------------
    print("------------------------------------------------------------------------")
    test_tasks = list(range(1,number_tasks+1))
    if test_tasks==tasks_num:
        print(f"There are {number_tasks} and are enumerated consecutively.")
    else:
        print("The tasks are not enumerated consecutively.")
    print()
    if (number_tasks == len(tasks_num)) and (len(tasks_times) == len(tasks_num)):
        print("The vector of tasks and the tasks times are of the same length.")
    else:
        print("The vector of tasks and the tasks times are not of the same length.")
        print("or/and the number of task quoted in the read file is different.")
    print()
    for i in precedence:
        if (i[0] not in tasks_num) or (i[1] not in tasks_num) or (i[0] > i[1]):
            print("A precedence value is not in the tasks vector or is enumerated wrongly.")
    print("The precedence relationships are correct.") 
    print()
    # -----------------------------------------------------------------------------
    
    tasks = tasks_num
    times = tasks_times
    precedence = precedence
    
    # ------ GUROBI approach ------------------------------------
    solution = GUROBI_Sol_Class(tasks, times, precedence, model_name)
    SALB2 = solution.SALB_2(number_stations)
    #-------------------------------------------------------------
    sys.stdout = sys.__stdout__
print("Gurobi has finisnhed")

# ----------- Quantum Annealing Approach -----------------------
sol_quantum = Quantum_Sol_Class(tasks, times, precedence, model_name)

print = False
solve = True
sol_quantum.Build_Model(number_stations,print,solve)
# --------------------------------------------------------------



