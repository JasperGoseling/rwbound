__author__ = 'BaiX'

# Jasper: making a small change and trying a commit from PyCharm

# A complete linear program that solves the optimization Problem 3.
# Try to find the coefficients c from a flow problem.

from pyomo.environ import *
from pyomo.opt import SolverFactory

import math
import itertools
import logging
logging.basicConfig(level=logging.DEBUG)

# Create a solver
optimizer = SolverFactory('gurobi')



############################
#
#
# Filename of data file
#
############################

# statespacedatafile = "rw_joint_statespace.dat"
# datafile = "rw_joint_data.dat"
# transitionprobability_datafile = "rw_joint_transition.dat"

# statespacedatafile = "tandem2d_statespace.dat"
# datafile = "tandem2d_data.dat"
# transitionprobability_datafile = "tandem2d_transition.dat"

statespacedatafile = "tandem2dB_statespace.dat"
datafile = "tandem2dB_data.dat"
transitionprobability_datafile = "tandem2dB_transition.dat"

# statespacedatafile = "tandem3d_statespace.dat"
# datafile = "tandem3d_data.dat"
# transitionprobability_datafile = "tandem3d_transition.dat"

# statespacedatafile = "couple3d_statespace.dat"
# datafile = "couple3d_data.dat"
# transitionprobability_datafile = "couple3d_transition.dat"


############################
#
#
# Initialize an abstract model
#
############################

# Create a model
model = AbstractModel('RW')
#model = ConcreteModel('RW')





############################
#
#
# Definition of the statespace
#  and its partitions
#
############################

model.dimen = Param(within=PositiveIntegers)
model.dimset = RangeSet(0,model.dimen-1,1) # starts at 0 to be consistent with Python indexing that starts at 0
model.extdimset = RangeSet(-1,model.dimen-1,1) # starts at -1, will be used to define affine (constant term at -1) functions in dimen dimensions
model.cuts = Set(model.dimset,within=NonNegativeIntegers,ordered=True)


# We want to use the dimension of the state space as the dimension of some sets. It
# seems that we can only do this when model.dimen is constructed. Therefore, we read
# the value of the model.dimen from the data file and load it into model.dimen. In
# addition we read the partition data
statespacedata = DataPortal(model=model)
statespacedata.load(filename=statespacedatafile,param=model.dimen,set=model.cuts) # ,
model.load(statespacedata)

# for notational convenience we define dimen outside the model
dimen = model.dimen.value



# we need to define dimset outside model, since model.dimset is not yet constructed
dimset = range(dimen) # starts at 0

#display(model.extdimset)





###########################################################################
#
#
# Definition of components, transitions and other parts of the random walk
#
###########################################################################

comps_init = [vector for vector in itertools.product( *[model.cuts[i] for i in dimset] )]
# print("comps_init: ", comps_init)
model.comps = Set(dimen=dimen, initialize=comps_init)

transitions_initial = [vector for vector in itertools.product((-1,0,1), repeat=dimen)]
model.transitions = Set(dimen=dimen, initialize=transitions_initial)



# In this part, we construct the transitions depending on the components to make sure that on the boundaries there are
# no transitions going outside the state space.

# First we construct all directions (-1,0,1) * (-1,0,1). Then we define a function on component,
# if one entry of the component is 0, then the directions corresponding to that entry will be (0,1)
def direction_init(component):
    direction_comps = [(-1, 0, 1)] * dimen
    for dim in range(len(component)):
        if component[dim] == 0:
            direction_comps[dim] = (0,1)
    return direction_comps

# print(direction_init((2,0)))

# print(direction_init((0,0)))

def transitions_init(component):
    direction_comps = direction_init(component)
    result = [vec for vec in itertools.product(*[direction_comps[i] for i in model.dimset])]
    return result

# print( "Transitions in dimensions in component (0,0): ", transitions_init((0,0)) )

# print( "Transitions in dimensions in component (0,2): ", transitions_init((0,2)) )



# Combine components and their transitions and make this list a set
def comps_transitions_init(model):
    for component in model.comps:
        transitions = transitions_init(component)
        for transition in transitions:
            yield component + transition
    return

# transitions_init = [vector for vector in itertools.product((-1,0,1), repeat=dimen)]
model.comps_transitions = Set(dimen = 2*dimen, initialize=comps_transitions_init)



model.transitionprobability = Param(model.comps_transitions, default=0)
model.perturbedtransitionprobability = Param(model.comps_transitions, default=0)





################################################################################################
#
# Define the coordinates, and a set containing a component + its transition + its coordinates
#
################################################################################################

def init_partcoords(model):
    # From model.cuts we generate partition
    # partition is a list of tuples, one tuple for each component in the partition
    # a tuple that describes a component is itself a tuple, with one entry for each dimension
    # the tuple for each dimension has two entries. If we consider the entry that corresponds to dimension DIM
    #    the first entry is the index into model.cuts[DIM] corresponding to this component
    #    the second entry is the value of the cut
    #print(partition)
    partition = [v for v in itertools.product(*[ [v for v in enumerate(model.cuts[dim])] for dim in dimset ])]
    for component in partition:
        # debug
        #print('--')
        #print(component)

        # find all the dimensions in which this component is bounded. What we need is that the index into model.cuts[dim] ,ie
        # component[dim][0], is not the last entry in model.cuts[dim]. '-1' is needed since the index start with 0.
        boundeddims = [dim for dim in dimset if component[dim][0]<(len(model.cuts[dim])-1) ]

        # next, we consider the dimensions in which this bounded component has a width/thickness larger than one
        nonsingulardims = [dim for dim in boundeddims if component[dim][1]<model.cuts[dim].value[component[dim][0]+1]-1]

        # debug
        #print(boundeddims)
        #print(nonsingulardims)

        coordstmp = [[component[dim][1]] for dim in dimset] # a list of lists
        for dim in nonsingulardims:
           coordstmp[dim].append(model.cuts[dim].value[component[dim][0]+1]-1)

        for c in itertools.product(*[coordstmp[dim] for dim in dimset] ) :
            yield c

    return

model.partcoords = Set(dimen=dimen,initialize=init_partcoords)





##############################################################################################
#
# Define the fine components, coordinates of the fine components, function TtoC
#
##############################################################################################

# Define the fine partition. The fine partition in each dimension contains all the coordinates, plus the last coordinate+1
list = [model.cuts[dim].value for dim in model.dimset]
list_total = []
for lst in list:
    list_new = []
    for index in range(len(lst)):
        list_new.append(lst[index])
        if lst[index]-1 >= 0 and lst[index]-1 not in lst:
            list_new.append(lst[index]-1)
        if index == len(lst)-1:
            list_new.append(lst[index]+1)
    list_new.sort()
    list_total.append(list_new)

fine_comps_init = [vector for vector in itertools.product( *[list_total[i] for i in model.dimset] )]
#print("fine_comps_init: ", fine_comps_init)
model.fine_comps = Set(dimen=dimen, initialize=fine_comps_init)



# Define a function. This function returns in which component(coarse partition) we will end, if we start in a
# component(fine partition) and go to a specific direction.
def TtoC(component,direction):
    component_new_list = [component[dim] + direction[dim] for dim in model.dimset]
    component_new = tuple(component_new_list)
    if component_new in comps_init:
        return component_new
    else:
        component_final_list = []
        for dim in model.dimset:
            list = []
            for item in model.cuts[dim].value:
                if item <= component_new[dim]:
                    list.append(item)
            component_final_list.append(max(list))
        component_final = tuple(component_final_list)
        return component_final



# Test if the TtoC function works
#for component in fine_comps_init:
#    transitions = transitions_init(component)
#    for transition in transitions:
#        component_new = TtoC(component,transition)
#        print("Fine component:", component, "Transition:", transition, "Coarse component after transition:", component_new)





##########################################################################
#
# We find the parameter c through an optimization problem
#
##########################################################################

transitionprobabilitydata = DataPortal(model=model)
transitionprobabilitydata.load(filename=transitionprobability_datafile,param=model.transitionprobability)
transitionprobabilitydata.load(filename=transitionprobability_datafile,param=model.perturbedtransitionprobability)
model.load(transitionprobabilitydata)

#print(model.transitionprobability[(1,0,1,0)])



def init_q(model,*args):
    return model.perturbedtransitionprobability.__getitem__(args) - model.transitionprobability.__getitem__(args)

model.q = Param(model.comps_transitions, initialize=init_q)



# We define a set that is called ttype. In the ttype set, we include all the e_i transition and all the transitions that
# have different probabilities between original and perturbed random walks.
def init_type(model):
    list = []
    for component in model.comps:
        transitions = transitions_init(component)
        for transition in transitions:
            if model.q.__getitem__(tuple(component)+transition) != 0:
                if transition not in list:
                    list.append(transition)
                transition_opp = tuple([-transition[dim] for dim in model.dimset])
                if transition_opp not in list:
                    list.append(transition_opp)
    list.remove(tuple([0]*dimen))
    for dim in model.dimset:
        vector = [0]*dimen
        vector[dim] = 1
        if tuple(vector) not in list:
            list.append(tuple(vector))
            vector_opp = [-vector[dim] for dim in model.dimset]
            list.append(tuple(vector_opp))
    return list
model.ttype = Set(dimen=dimen, initialize=init_type)



def init_type_plus(model):
    list = []
    for transition in model.transitions:
        ind = [0]*dimen
        for dim in model.dimset:
            if  transition[dim] > 0:
                ind[dim] = 1
        if ind != [0]*dimen:
            list.append(transition)
    list_new = []
    for item in list:
        if item in model.ttype:
            list_new.append(item)
    for item in list_new:
        if tuple([-item[dim] for dim in model.dimset]) in list_new:
            list_new.remove(item)
    for dim in model.dimset:
        vector = [0]*dimen
        vector[dim] = 1
        if tuple(vector) not in list_new:
            list_new.append(tuple(vector))
    return list_new
model.type_plus = Set(dimen=dimen, initialize=init_type_plus)



def ind_trans_comps(direction, component):
    transitions = transitions_init(component)
    if tuple(direction) in transitions:
        return 1
    else:
        return 0



# define a set containing N_{k(n)} union N_{k(n+u)} + u, where u in model.ttype

def func_fine_comps_type_state(fine_component,u):
    component = TtoC(fine_component,tuple([0]*dimen))
    transitions = transitions_init(component)
    list_temp = []
    for d in transitions:
        list_temp.append(d)
    new_component = TtoC(fine_component,u)
    new_transitions = transitions_init(new_component)
    for new_transition in new_transitions:
        sum_transition = tuple([new_transition[dim]+u[dim] for dim in model.dimset])
        if sum_transition not in list_temp:
            list_temp.append(sum_transition)
    return list_temp

#print(func_fine_comps_u_d((0,0),(0,0)))

def init_fine_comps_type_plus_state(model):
    for fine_component in model.fine_comps:
        component = TtoC(fine_component,tuple([0]*dimen))
        transitions = transitions_init(component)
        for u in transitions:
            if u in model.type_plus:
                set_d = func_fine_comps_type_state(fine_component,u)
                for vector in set_d:
                    yield fine_component+u+vector
    return

model.fine_comps_type_plus_state = Set(dimen=3*dimen, initialize=init_fine_comps_type_plus_state)



#i=0
#for fine_component in fine_comps_init:
#    component = TtoC(fine_component,(0,0))
#    transitions = transitions_init(component)
#    for u in transitions:
#        set_d = func_fine_comps_u_d(fine_component,u)
#        for vector in set_d:
#            i += 1
#            print(i,":", tuple(fine_component)+u+vector)





##################################################################################
#
# Variables
#
##################################################################################

model.c = Var(model.fine_comps_type_plus_state, model.ttype, domain=NonNegativeReals, initialize=0)





###################################################################
#
# Constraints
#
###################################################################

def init_flow_constr_expr(model,*args):
    fine_component = [item for item in args[0:dimen]]
    u = [item for item in args[dimen:2*dimen]]
    w = [item for item in args[2*dimen:3*dimen]]
    component = TtoC(fine_component, tuple([0]*dimen))
    transitions = transitions_init(component)
    complete_states = func_fine_comps_type_state(fine_component,u)
    result = 0
    for d in complete_states:
        if tuple([w[dim]-d[dim] for dim in model.dimset]) in model.ttype:
            result += model.c.__getitem__(tuple(fine_component+u)+d+tuple([w[dim]-d[dim] for dim in model.dimset]))
    for d in complete_states:
        if tuple([d[dim]-w[dim] for dim in model.dimset]) in model.ttype:
            result += -model.c.__getitem__(tuple(fine_component+u+w+[d[dim]-w[dim] for dim in model.dimset]))
    if tuple(w) in transitions:
        result += model.transitionprobability.__getitem__(component+tuple(w))
    new_component = TtoC(fine_component,u)
    new_transitions = transitions_init(new_component)
    if tuple([w[dim]-u[dim] for dim in model.dimset]) in new_transitions:
        result += -model.transitionprobability.__getitem__(new_component+tuple([w[dim]-u[dim] for dim in model.dimset]))
    return result

model.flow_constr_expr = Expression(model.fine_comps_type_plus_state, rule=init_flow_constr_expr)



def init_flow_constr(model,*args):
    return model.flow_constr_expr.__getitem__(args) == 0

model.flow_constr = Constraint(model.fine_comps_type_plus_state, rule=init_flow_constr)





#####################################################################
#
# Objective function
#
#####################################################################

"""def init_obj_func(model):
    result = 0
    return result"""

def init_obj_func(model):
    result = 0
    for fine_component in model.fine_comps:
        component = TtoC(fine_component,tuple([0]*dimen))
        transitions = transitions_init(component)
        for u in transitions:
            if u in model.type_plus:
                complete_states = func_fine_comps_type_state(fine_component,u)
                for d in complete_states:
                    new_component = TtoC(fine_component,d)
                    new_transitions = transitions_init(new_component)
                    for v in new_transitions:
                        if v in model.ttype:
                            result += model.c.__getitem__(fine_component+u+d+v)
    return result

model.obj_func = Objective(rule=init_obj_func, sense=minimize)





# # ###########################
# #
# # Create instance
# #
# # ###########################
#
#
#
# instance = model.create(datafile)
#
#
#
#
# ############################
# #
# # Print
# #
# ############################
#
# print("\n Abstract model: \n")
# #model.pprint()
#
# print("\n INSTANCE: \n")
# #instance.pprint()
#
#
# ############################
# #
# # Solve
# #
# ############################
#
# #print("\n PREPROCESS. \n")
# #instance.preprocess()
#
# print("\n SOLVE: \n")
# results = optimizer.solve(instance)
# instance.load(results)
#
# print("\n SOLUTION: \n")
# # instance.pprint()
# print(results)
