#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:00:03 2024

@author: luismoncayo
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 13:10:44 2022

@author: luismoncayo
"""
import gurobipy as gp
from gurobipy import GRB
import math
import pandas as pd
pd.set_option('display.max_columns', None)

class GurobiSolution:
    
    cycle_C = 0
    work_W = 0
    
    def __init__(self, tasks, processing_times, precedence_relationships, instance_name):
        self.tasks = tasks
        self.processing_times = processing_times
        self.precedence_relationships = precedence_relationships
        self.instance_name = instance_name
        
    def salb_1(self, cycle_time):
        self.cycle_time = cycle_time
        max_stations = math.ceil(sum(self.processing_times)/self.cycle_time)
        
        SALB1_model = gp.Model(self.instance_name)
        y = SALB1_model.addVars(range(1,max_stations+1), vtype=GRB.BINARY, name="y")
        
        SALB1_model.setObjective(sum(y[s] for s in range(1,max_stations+1)), GRB.MINIMIZE)
        
        x = SALB1_model.addVars(range(1,len(self.tasks)+1), range(1,max_stations+1), vtype=GRB.BINARY, name="x")
        
        for s in range(1,max_stations+1):
            SALB1_model.addConstr( sum(self.processing_times[t-1]*x[t,s] for t in self.tasks)<= self.cycle_time*y[s], "Station[%d]" % s )

        for t in self.tasks:
            SALB1_model.addConstr( sum(x[t,s] for s in range(1,max_stations+1)) == 1, "Task[%d]" % t )
            
        for precedence in self.precedence_relationships:
            i = precedence[0]
            j = precedence[1]
            SALB1_model.addConstr( sum(s*x[i,s] for s in range(1,max_stations+1))<= sum(s*x[j,s] for s in range(1,max_stations+1)), "Precedence(%d,%d)" %(i,j) )
            
        # for s in range(max_stations-1):
        #     SALB1_model.addConstr(y[s+1] <= y[s], "Stations(%d,%d)" %(s+1,s) )

        #SALB1_model.write("Outputs_Gurobi/"+self.instance_name+'.lp')
        
        SALB1_model.optimize()
        
        result_SALB_1 = pd.DataFrame(columns=["Task", "Workstation"])
        for t in range(1,len(self.tasks)+1):
            for s in range(1,max_stations+1):
                if x[t,s].X > 0.99:
                    result_SALB_1 = pd.concat([result_SALB_1, pd.DataFrame({"Task": [t], "Workstation": [s]})], ignore_index=True)
                    #print(x[t,s])
        result_SALB_1["Times"] = self.processing_times
        
        global cycle_C
        global work_W
        cycle_C = self.cycle_time
        work_W = max_stations
        
        return result_SALB_1
        

    def salb_2(self, stations):
        self.stations = stations
        SALB2_model = gp.Model(self.instance_name)

        cycle_time = SALB2_model.addVar(name="C")

        SALB2_model.setObjective(cycle_time, GRB.MINIMIZE)

        x = SALB2_model.addVars(range(1,len(self.tasks)+1), range(1,self.stations+1), vtype=GRB.BINARY, name="x")

        # The workstation's time must be less that or equal to the cycle time
        for s in range(1,self.stations+1):
            SALB2_model.addConstr( sum(self.processing_times[t-1]*x[t,s] for t in self.tasks)<= cycle_time, "Station[%d]" % s )

        # each task must be assigned to a station
        for t in self.tasks:
            SALB2_model.addConstr( sum(x[t,s] for s in range(1,self.stations+1)) == 1, "Task[%d]" % t )
            
        # prece
        for precedence in self.precedence_relationships:
            i = precedence[0]
            j = precedence[1]
            SALB2_model.addConstr( sum(s*x[i,s] for s in range(1,self.stations+1))<= sum(s*x[j,s] for s in range(1,self.stations+1)), "Precedence(%d,%d)" %(i,j) )

        #SALB2_model.write("Outputs_Gurobi/"+self.instance_name+'.lp')

        # Solve
        SALB2_model.optimize()
        
        result_SALB_2 = pd.DataFrame(columns=["Task", "Workstation"])
        for t in range(1,len(self.tasks)+1):
            for s in range(1,self.stations+1):
                if x[t,s].X > 0.99:
                    result_SALB_2 = pd.concat([result_SALB_2, pd.DataFrame({"Task": [t], "Workstation": [s]})], ignore_index=True)
                    #print(x[t,s])
        result_SALB_2["Times"] = self.processing_times
        
        global cycle_C
        global work_W
        cycle_C = SALB2_model.ObjVal
        work_W = self.stations
        
        return result_SALB_2

    def print_solution(self, results_table):
        self.results_table = results_table
        print("------------------------------------------")
        print('Cycle Time is: %g' % cycle_C)
        print('Number of stations: %d' % work_W)
        print("------------------------------------------")
        print("\n---> Task per station")
        for w in range(1,work_W+1):
            filtered_df = results_table[results_table["Workstation"] == w]
            cycle = filtered_df["Times"].sum()
            print()
            print('Station %d' % w)
            print("Tasks:", ", ".join(map(str, filtered_df["Task"])))
            print("Times:", ", ".join(map(str, filtered_df["Times"])))
            print(f"The cycle of station {w} is {cycle}")
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            
    # def test_gurobi_solution(self, solution):
    #     if solution.shape[0] != len(self.tasks):
    #         print("Not the same")
    #     else:
    #         print("The same size")
        
        
        