#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 10:22:15 2024

@author: luismoncayo
"""
from dimod import ConstrainedQuadraticModel,BinaryQuadraticModel,Binary,Real,quicksum
from dwave.system import LeapHybridCQMSampler,DWaveSampler,EmbeddingComposite
#from dwave.system import DWaveSampler, EmbeddingComposite
#from dimod import BinaryQuadraticModel, Binary
import dimod
import pandas as pd
pd.set_option('display.max_columns', None)
import pprint
import sys
import math
import copy

class QuantumSolution:
    
    def __init__(self, tasks, processing_times, precedence_relationships, model_name):
        self.tasks = tasks
        self.times = processing_times
        self.precedence = precedence_relationships
        self.model_name = model_name
        #self.print_model = False
        #self.solve_model = False
        
    def quantum_salb_1(self,cycle_time,print_model:bool,solve_model:bool):
        workstations = math.ceil(sum(self.times)/cycle_time)
        
        y_list = [f"y{i}" for i in range(1,workstations+1)]
        x = []
        for t in self.tasks:
            for w in range(1,workstations+1):
                x.append(f't{t}_w{w}')
        
        bqm = BinaryQuadraticModel('BINARY')
        
        for r in range(len(x)):
            bqm.add_variable(x[r])
        
        # Objective function
        for y in y_list:
            bqm.add_variable(v=y,bias=1)
            
        #print("Variables in BQM:", bqm.variables)
        
        c0 = [(y_list[u],1) for u in range(len(y_list))]
        bqm.add_linear_equality_constraint(terms = c0, 
                                           lagrange_multiplier =1000, 
                                           constant = -workstations)
        if print_model:
            print("Objective function")
            print(c0)

        #workstation constraints, less than the cycle time ----------------------------
        for q in range(workstations):
            work = []
            increment = q
            while len(work) < len(self.tasks):
                work.append(x[increment])
                increment+=workstations
            c1 = [(work[a],self.times[a]) for a in range(len(self.tasks))]
            bqm.add_linear_inequality_constraint(c1+[(y_list[q],-cycle_time)], 
                                                 constant = 0,
                                                 ub=0,
                                                 lb=-20,
                                                 lagrange_multiplier = 1000,
                                                 #penalization_method = "slack",
                                                 label = 'c1_time_')
            
            if print_model:
                if q==0:
                    print("Cycle time constraints")
                print(c1+[(y_list[q],-cycle_time)])
        
        #task assigned to one station -----------------------------------
        for w in range(len(self.tasks)):
            e = workstations*w
            task_workstation = x[e:e+workstations]
            c2 = [(task_workstation[r],1) for r in range(workstations)]
            bqm.add_linear_equality_constraint(terms = c2, 
                                               lagrange_multiplier = 2000, 
                                               constant = -1)
            if print_model:
                if w==0:
                    print("One task to one worstation constriants")
                print(c2)
        
        if print_model:
            print("Precedence constaints")
            
        for p in self.precedence:
            first_task = (p[0]-1)*workstations
            second_task = (p[1]-1)*workstations
            to_task = x[first_task:first_task+workstations]
            from_task = x[second_task:second_task+workstations]
            c3_to = [(to_task[a],a+1) for a in range(workstations)]
            c3_from = [(from_task[a],-(a+1)) for a in range(workstations)]
            c3 = c3_to+c3_from
            bqm.add_linear_inequality_constraint(c3, 
                                                 constant =0,
                                                 ub=0,
                                                 lb=-workstations+1,
                                                 lagrange_multiplier = 10000,
                                                 #penalization_method = "slack",
                                                 label = p)
        
            if print_model:
                print(c3)
        
        if solve_model:
            print()
            print("Solving the model")
            print()
            # Create a D-Wave sampler
            sampler = DWaveSampler()
            # Use EmbeddingComposite to automatically handle embedding
            embedding_sampler = EmbeddingComposite(sampler)
            # Sample your problem
            response = embedding_sampler.sample(bqm, num_reads=3600, label=self.model_name)
            
            outputs = pd.DataFrame(response)
            outputs.to_csv("/Users/luismoncayo/Dropbox/Python/Line_Balancing/Outputs_Quantum/"+self.model_name+"_B1.csv")
        
            print()
            print("The solution D-Wave returned is in the file: "+self.model_name+"_B1.csv")
        
        return outputs, workstations
            
    def quantum_salb_2(self,num_stations,print_model:bool,solve_model:bool ):
        self.num_stations = num_stations
        workstations = list(range(1,self.num_stations+1))
        
        cqm = ConstrainedQuadraticModel()

        # Define the variables
        name = []
        for t in self.tasks:
            for w in workstations:
                name.append(f't{t}_w{w}')
        x = [Binary(label=i) for i in name]
        cycle_time = Real(label="cycle_time", upper_bound=sum(self.times))

        # Set the objective function
        cqm.set_objective(cycle_time)

        # workstation constraints
        for i in range(len(workstations)):
            work = []
            increment = i
            while len(work) < len(self.tasks):
                work.append(x[increment])
                increment+=len(workstations)
            cqm.add_constraint(quicksum(self.times[j]*work[j] for j in range(len(self.tasks))) - cycle_time <= 0, label='w'+str(i+1))
            
        # tasks constrains
        for i in range(len(self.tasks)):
            j = len(workstations)*i
            task_to_ = x[j:j+len(workstations)]
            cqm.add_constraint(quicksum(task_to_[i] for i in range(len(task_to_))) == 1, label='t'+str(i+1))
            
        # precedences constraints
        for p in self.precedence:
            first_task = (p[0]-1)*len(workstations)
            second_task = (p[1]-1)*len(workstations)
            to_task = x[first_task:first_task+len(workstations)] 
            from_task = x[second_task:second_task+len(workstations)]
            cqm.add_constraint(quicksum((i+1)*to_task[i] for i in range(len(workstations)))-quicksum((j+1)*from_task[j] for j in range(len(workstations))) <= 0, label=p)
        
        # with open("cqm_model.cqm","wb") as self.f:
        #     cqm.to_file(self.f)
        # Print the model 
        if isinstance(print_model, bool):
            self.print_model = print_model
            if self.print_model:
                with open("Outputs_Quantum/"+self.model_name+'_2.txt', "w") as file:
                    # Redirect sys.stdout to the file
                    sys.stdout = file
                    #print(cqm) # Here the model is printed
                    print("Variables in the model:")
                    for var in cqm.variables:
                        print(var)
                    print("\nObjective Function:")
                    print(cqm.objective)
                    print("\nConstraints in the model:")
                    for constraint in cqm.constraints:
                        print(f"{constraint}: {cqm.constraints[constraint]}")
                    sys.stdout = sys.__stdout__
        else:
            raise ValueError("Input True or False to print the model.")
        
        # Solve the model
        if isinstance(solve_model, bool):
            self.solve_model = solve_model
            if self.solve_model:
                # Using ExactSolver()
                # - Submit to the CQM sampler
                cqm_sampler = LeapHybridCQMSampler()
                sampleset = cqm_sampler.sample_cqm(cqm, label=self.model_name)
                df = pd.DataFrame(sampleset)
                df.to_csv("Outputs_Quantum/"+self.model_name+'_B2.csv', index=False)
                print()
                print("+++++++++++++++++++++++++++++++++++++++++++++++")
                print()
                print("DWave returns a solution")
                print("The solution is stored in the file:"+self.model_name)
        else:
            raise ValueError("Input True or False to solve the model.")
        
        return df
        
    def test_solution(self, outputs, workstations,tag):
        print()
        print('The solution is tested')
        print()
        with open("Outputs_Quantum/"+self.model_name+'_'+tag+'_sol.txt', "w") as file:
            sys.stdout = file
            outputs = outputs.drop(columns=[col for col in outputs.columns if 'slack' in col or 'y' in col or 'Unnamed' in col])
            feasible_sol = []
            for j in range(len(outputs)):
                selected_response = outputs.loc[j]
                task_work_var = selected_response[selected_response == 1].index.tolist()
                if len(task_work_var) == len(self.tasks):
                    sol_data = pd.DataFrame([item.split('_') for item in task_work_var], columns=['Task', 'Workstation'])
                    sol_data['Task'] = sol_data['Task'].str[1:]
                    sol_data['Workstation'] = sol_data['Workstation'].str[1:]
                    sol_data = sol_data.astype(int)
                    any_tasks = sol_data['Task'].isin(self.tasks).any() and all(name in sol_data['Task'].values for name in self.tasks)
                    any_workstation = sol_data['Workstation'].isin(list(range(1, workstations+1))).any() and all(name in sol_data['Workstation'].values for name in list(range(1, workstations+1)))
                    if any_tasks and any_workstation:
                        flag = True
                        for p in self.precedence:
                            initial_task = p[0]
                            arriving_task = p[1]
                            initial_workstation = sol_data.loc[sol_data['Task'] == initial_task, 'Workstation'].iloc[0]
                            final_workstation = sol_data.loc[sol_data['Task'] == arriving_task, 'Workstation'].iloc[0]
                            if initial_workstation > final_workstation :
                                flag = False
                                break
                        if flag:
                            feasible_sol.append(j)
                            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
                            #print(sol_data)
                            df = pd.DataFrame({'Task': self.tasks,'Time': self.times})
                            merged_df = pd.merge(df, sol_data, on='Task')
                            print("Solution: "+str(j))
                            for w in range(1,workstations+1):
                                selected_rows = merged_df.loc[merged_df['Workstation'] == w]
                                print()
                                print("Cycle of workstation "+str(w)+" is "+str(selected_rows['Time'].sum()))
                                print()
                                print(selected_rows)
                                print()
            if len(feasible_sol)==0:
                print('There is no feasible solutions.')
            else:
                print()
                print(feasible_sol)
            sys.stdout = sys.__stdout__
            
    def test_solution_print_console(self, outputs, workstations):
        print()
        print('The solution is tested')
        print()
        outputs = outputs.drop(columns=[col for col in outputs.columns if 'slack' in col or 'y' in col or 'Unnamed' in col])
        feasible_sol = []
        for j in range(len(outputs)):
            selected_response = outputs.loc[j]
            task_work_var = selected_response[selected_response == 1].index.tolist()
            if len(task_work_var) == len(self.tasks):
                sol_data = pd.DataFrame([item.split('_') for item in task_work_var], columns=['Task', 'Workstation'])
                sol_data['Task'] = sol_data['Task'].str[1:]
                sol_data['Workstation'] = sol_data['Workstation'].str[1:]
                sol_data = sol_data.astype(int)
                any_tasks = sol_data['Task'].isin(self.tasks).any() and all(name in sol_data['Task'].values for name in self.tasks)
                any_workstation = sol_data['Workstation'].isin(list(range(1, workstations+1))).any() and all(name in sol_data['Workstation'].values for name in list(range(1, workstations+1)))
                if any_tasks and any_workstation:
                    flag = True
                    for p in self.precedence:
                        initial_task = p[0]
                        arriving_task = p[1]
                        initial_workstation = sol_data.loc[sol_data['Task'] == initial_task, 'Workstation'].iloc[0]
                        final_workstation = sol_data.loc[sol_data['Task'] == arriving_task, 'Workstation'].iloc[0]
                        if initial_workstation > final_workstation :
                            flag = False
                            break
                    if flag:
                        feasible_sol.append(j)
                        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
                        #print(sol_data)
                        df = pd.DataFrame({'Task': self.tasks,'Time': self.times})
                        merged_df = pd.merge(df, sol_data, on='Task')
                        print("Solution: "+str(j))
                        for w in range(1,workstations+1):
                            selected_rows = merged_df.loc[merged_df['Workstation'] == w]
                            print()
                            print("Cycle of workstation "+str(w)+" is "+str(selected_rows['Time'].sum()))
                            print()
                            print(selected_rows)
                            print()
        if len(feasible_sol)==0:
            print('There are no feasible solutions.')
        else:
            print()
            print(feasible_sol)
          
        
