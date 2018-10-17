__author__ = 'BaiX'

# A complete linear program that solves the optimization problem.
# Find the coefficients c from a flow problem.
# Consider the non-negative coefficients, and generalize the flow problem.
# The bias terms are generalized, i.e. we consider D_u^t(n)

# Add a new comment and commit it.

from pyomo.environ import *
from pyomo.opt import SolverFactory

import math
import itertools

# a new comment
# Create a solver
optimizer = SolverFactory('gurobi')


#########################################
#
# Input from users.
#
#########################################

# The datafile has everything
# The statespacedatafile contains initialization for parameter dimen and set cuts
# The transitionprobability_datafile has parameter transitionprobability and perturbedtransitionprobability

# statespacedatafile = "statespaceYanting.dat"
# datafile = "dataYanting.dat"
# transitionprobability_datafile = "transitionYanting.dat"

# statespacedatafile = "statespaceTandem2d.dat"
# datafile = "dataTandem2d.dat"
# transitionprobability_datafile = "transitionTandem2d.dat"

# statespacedatafile = "statespaceTandem2d_finebias.dat"
# datafile = "dataTandem2d_finebias.dat"
# transitionprobability_datafile = "transitionTandem2d_finebias.dat"

# statespacedatafile = "statespaceRWJointDeparture.dat"
# datafile = "dataRWJointDeparture.dat"
# transitionprobability_datafile = "transitionRWJointDeparture.dat"

statespacedatafile = "tandem3d_statespace.dat"
datafile = "tandem3d_data.dat"
transitionprobability_datafile = "tandem3d_transition.dat"

# statespacedatafile = "couple3d_statespace.dat"
# datafile = "couple3d_data.dat"
# transitionprobability_datafile = "couple3d_transition.dat"

# Here we input the parameters needed for the invariant measure of the perturbed random walk.
# The parameter pi_bar is a list of tuples. Every sublist contains a set of parameters for a product-form invariant measure.
# The whole list gives us an invariant measure, the form of which is a sum of geometric terms.
# pi_bar = [(1/4, 1/4)]
# pi_bar = [(2/5,2/3)]
# pi_bar = [(-3+math.sqrt(33))/6, (-3+math.sqrt(33))/6]
rho = 1/3
pi_bar = [(rho, rho, rho)]

# The weights of the geometric terms in the sum
weight = [1]

# Input the performance measure below
def F(component, dim):
    if component == tuple([0]*dimen) and dim == -1:
    # if dim == 0:
    # if dim == 1:
        return 1
    else:
        return 0



#####################################
#
#
# Initialize an abstract model
#
#####################################

# Create a model
model = AbstractModel('RW')
#model = ConcreteModel('RW')
flow = AbstractModel()

# import the flow problem from flow_problem. Solve the problem and load the results
from flow_problem_type import model as flow_model

# Create an instance for the flow problem
flow_instance = flow_model.create(datafile)
#flow_instance.pprint()

# Solve the flow problem
flow_results = optimizer.solve(flow_instance)
#print(flow_results)
flow_instance.load(flow_results)
#flow_instance.pprint()



############################
#
# Definition of the statespace
#  and its partitions.
# Load the data from all the datafiles
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

# We want Python to generate all the expressions for the constraints automatically.
# We divide all the expressions to be generated into two groups.
# One group is the directed expressions, for example A, B since they depend on the component and the transition.
# Another one is the undirected expressions, which only depend on the component.
# The bound set is used to generate constraints to bound the bias terms both from below and also from above.

model.dstring = Set(initialize=["A","B"])
model.ustring = Set(initialize=["Fbar","FF"])
model.bound = Set(initialize=["lower","upper"])



###########################################################################
#
#
# Definition of components, transitions and other parts of the random walk
#
###########################################################################

comps_init = [vector for vector in itertools.product( *[model.cuts[i] for i in dimset] )]
#print("comps_init: ", comps_init)
model.comps = Set(dimen=dimen, initialize=comps_init)

transitions_initial = [vector for vector in itertools.product((-1,0,1), repeat=dimen)]
#print(transitions_initial)
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

def init_FF(model,*args):
    component = tuple([item for item in args[0:dimen]])
    dim = args[dimen]
    return F(component,dim)
model.FF = Param(model.comps, model.extdimset, initialize=init_FF)

# Load the data for transition probabilities
transitionprobabilitydata = DataPortal(model=model)
transitionprobabilitydata.load(filename=transitionprobability_datafile,param=model.transitionprobability)
transitionprobabilitydata.load(filename=transitionprobability_datafile,param=model.perturbedtransitionprobability)
transitionprobabilitydata.load(filename=transitionprobability_datafile,param=model.FF)

model.load(transitionprobabilitydata)



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



# Define a function to find the bounded dimensions for a component
def boundeddim(component):
    list = []
    for dim in model.dimset:
        cut_dim = model.cuts[dim].value
        index = cut_dim.index(component[dim])
        if index < len(cut_dim) - 1:
            list.append(dim)
    return list

#print(boundeddim((5,2)))



# This function returns in which coarse component a state is
def coordtoC(state):
    if state in comps_init:
        return state
    else:
        component_final_list = []
        for dim in model.dimset:
            list = []
            for item in model.cuts[dim].value:
                if item <= state[dim]:
                    list.append(item)
            component_final_list.append(max(list))
        component_final = tuple(component_final_list)
        return component_final

#print(coordtoC((0,0)))



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



# Define the TtoC function. This function returns in which component(coarse partition) we will end, if we start in a
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



def init_fine_partcoords(model):
    partition = [v for v in itertools.product(*[ [v for v in enumerate(list_total[dim])] for dim in model.dimset ])]
    for component in partition:
        # debug
        #print('--')
        #print(component)

        # find all the dimensions in which this component is bounded. What we need is that the index into model.cuts[dim] ,ie
        # component[dim][0], is not the last entry in model.cuts[dim]. '-1' is needed since the index start with 0.
        boundeddims = [dim for dim in dimset if component[dim][0]<(len(list_total[dim])-1) ]

        # next, we consider the dimensions in which this bounded component has a width/thickness larger than one
        nonsingulardims = [dim for dim in boundeddims if component[dim][1]<list_total[dim][component[dim][0]+1]-1]

        # debug
        #print(boundeddims)
        #print(nonsingulardims)

        coordstmp = [[component[dim][1]] for dim in dimset] # a list of lists
        for dim in nonsingulardims:
           coordstmp[dim].append(list_total[dim][component[dim][0]+1]-1)

        for c in itertools.product(*[coordstmp[dim] for dim in dimset] ) :
            yield c

    return

model.fine_partcoords = Set(dimen=dimen,initialize=init_fine_partcoords)



def fine_boundeddim(component):
    list = []
    for dim in model.dimset:
        cut_dim = list_total[dim]
        index = cut_dim.index(component[dim])
        if index < len(cut_dim) - 1:
            list.append(dim)
    return list



def init_fine_comps_unbddim(model):
    for fine_component in model.fine_comps:
        bddim = fine_boundeddim(fine_component)
        for dim in model.dimset:
            if dim not in bddim:
                yield tuple([fine_component[d] for d in model.dimset] + [dim])
    return

model.fine_comps_unbddim = Set(dimen=dimen+1, initialize=init_fine_comps_unbddim)



# This function returns in which fine component a state is
def coordtoT(state):
    if state in fine_comps_init:
        return state
    else:
        component_final_list = []
        for dim in model.dimset:
            list = []
            for item in list_total[dim]:
                if item <= state[dim]:
                    list.append(item)
            component_final_list.append(max(list))
        component_final = tuple(component_final_list)
        return component_final

#print(coordtoT((5,1)))



##########################################################################
#
# We find the parameter c through an optimization problem
#
##########################################################################



#print(model.transitionprobability[(1,0,1,0)])



def init_q(model,*args):
    return model.perturbedtransitionprobability.__getitem__(args) - model.transitionprobability.__getitem__(args)

model.q = Param(model.comps_transitions, initialize=init_q)



# We only consider bias terms D^t_u(n) if in any component(even not the one n belongs to), we make perturbation on direction u.
# Therefore, we call all those perturbed directions types and define a set type below containing all the perturbed directions and their opposite directions.

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

model.type = Set(dimen=dimen, initialize=init_type)



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
        if item in model.type:
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



# define a set containing N_{k(n)} union N_{k(n+u)} + u, where u in model.type

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



# Define a new set. An element in the set contains a fine component and a type_plus(provided that it is in the transitions of the fine component)
# and a state that is in N_{k(n)} union N_{k(n+u)} + u.
# This set is used later to pass value to c coefficients.

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



# The difference from the previous set is that an element now contains coordinates of a fine component.
# This set is used later when we write the expressions for the constraints.

def init_fine_coords_type_plus_state(model):
    for fine_coord in model.fine_partcoords:
        fine_component = coordtoT(fine_coord)
        component = TtoC(fine_component,tuple([0]*dimen))
        transitions = transitions_init(component)
        for u in transitions:
            if u in model.type_plus:
                set_d = func_fine_comps_type_state(fine_component,u)
                for vector in set_d:
                    yield fine_coord+u+vector
    return

model.fine_coords_type_plus_state = Set(dimen=3*dimen, initialize=init_fine_coords_type_plus_state)



# Obtain the parameter c from the file flow_problem.py, which describes an optimization flow problem solving c_{k,u,w,v}
# We consider
# D_u^{t+1}(n) = F(n+u) - F(n) + \sum_{w \in N_{k(n)}} \sum_{v \in N_{k(n+w)}} c_{k(n),u,w,v} D_v^t(n+w)

def init_c(model,*args):
    fine_component = [item for item in args[0:dimen]]
    u = [item for item in args[dimen:2*dimen]]
    d = [item for item in args[2*dimen:3*dimen]]
    v = [item for item in args[3*dimen:4*dimen]]
    return float(flow_instance.c.__getitem__(tuple(fine_component+u+d+v)))

model.c = Param(model.fine_comps_type_plus_state, model.type, initialize=init_c)



###################################
#
#
# Declaration of all the variables
#
###################################

model.A = Var(model.comps, model.type_plus, model.extdimset, domain=NonNegativeReals)
model.B = Var(model.comps, model.type_plus, model.extdimset, domain=NonNegativeReals)
# model.G = Var(model.comps, model.extdimset, domain=NonNegativeReals)
model.Fbar = Var(model.comps, model.extdimset, domain=NonNegativeReals)



################### End of declaration for sets and parameters. ###########################################





###################################################################################
#
# Objective function
#
###################################################################################

# In this function, we try to calculate the normalization constant. The idea is that we take the sum over all the states
# in state space. Divide 1 by the sum is the normalization constant.
def init_const(model):
    sum_result = 0
    for i in range(len(pi_bar)):
        result = 0
        for component in model.comps:
            result_part = 1
            bddim = boundeddim(component)
            for dim in model.dimset:
                if dim in bddim:
                    start = component[dim]
                    index = model.cuts[dim].value.index(component[dim])
                    end = model.cuts[dim].value[index+1]-1
                    result_part = result_part*(pi_bar[i][dim]**(end+1)-pi_bar[i][dim]**start)/ (pi_bar[i][dim]-1)
                else:
                    start = component[dim]
                    result_part = result_part*pi_bar[i][dim]**start/(1-pi_bar[i][dim])
            result += result_part
        sum_result += weight[i]*result
    return 1/sum_result

model.const = Param(initialize=init_const)



def obj_func_expr(model):
    # The final result equals sum over all the product-form measure times their weights
    sum_result = 0
    for i in range(len(pi_bar)):
        # This result is the sum for one product-form measure
        result = 0
        # For each component, we sum up the performance measures of all the states in this component.
        for component in model.comps:
            # A component can be bounded in some dimensions while unbounded on the others.
            # Get the bounded dimensions of the component
            bddim = boundeddim(component)
            # Since the performance measure is component-wise linear, let us assume that for the current component,
            # the coefficient w.r.t dimension k is f_k(dimensions are {1,2,...,n}), \ie
            # F(n) = f_0 + f_1 n_1 + f_2 n_2 + ... f_k n_k
            # The prodct-form invariant measure(ignore the normalization constant) is rho_1^{n_1} rho_2^{n_2} ... rho_k^{n_k}
            # Thus we want to calculate \sum_{n_1 = start_1}^{end_1}...\sum_{n_k = start_k}^{end_k} (f_0 + f_1 n_1 + f_2 n_2 + ... f_k n_k)rho_1^{n_1} rho_2^{n_2} ... rho_k^{n_k}
            # Note that those ends can be infinity if the component is unbounded on some dimensions.
            # The sum above = f_0 \sum_{n_1 = start_1}^{end_1}...\sum_{n_k = start_k}^{end_k}rho_1^{n_1} rho_2^{n_2} ... rho_k^{n_k}
            # + f_1 \sum_{n_1 = start_1}^{end_1}...\sum_{n_k = start_k}^{end_k} n_1 rho_1^{n_1} rho_2^{n_2} ... rho_k^{n_k} + ...
            # In the following part, we first calculate the big junk of sum(called result_part) for each extended dimension,
            # and then add it up to the result
            for c_dim in model.extdimset:
                result_part = 1
                # Since the invariant measure is product-form and separable, we can first take the sum on every dimension and then take the product.
                for dim in model.dimset:
                    # If dim is equal to the extended dimension we consider currently, then we need to calculate
                    # \sum_{start_k}^{end_k} n_k rho^{n_k}
                    # Otherwise, we should calculate \sum_{start_k}^{end_k} rho^{n_k}
                    if dim == c_dim:
                        # Depending on whether end_k is infinity or not, the formula is different.
                        if dim in bddim:
                            # Find the start_k and end_k for the component
                            start = component[dim]
                            index = model.cuts[dim].value.index(component[dim])
                            end = model.cuts[dim].value[index+1]-1
                            # Formula for \sum_{start_k}^{end_k} n_k rho^{n_k}
                            result_part = result_part*( pi_bar[i][dim]*(pi_bar[i][dim]**start + end*pi_bar[i][dim]**(end+1) - (end+1)*pi_bar[i][dim]**end) - start*(pi_bar[i][dim]-1)*pi_bar[i][dim]**start )/ (pi_bar[i][dim]-1)**2
                        else:
                            start = component[dim]
                            # Formula for \sum_{start_k}^{\infty} n_k rho^{n_k}
                            result_part = result_part*(pi_bar[i][dim]**start*(start*(-pi_bar[i][dim])+start+pi_bar[i][dim]))/(pi_bar[i][dim]-1)**2
                    else:
                        if dim in bddim:
                            start = component[dim]
                            index = model.cuts[dim].value.index(component[dim])
                            end = model.cuts[dim].value[index+1]-1
                            # Formula for \sum_{start_k}^{end_k} rho^{n_k}
                            result_part = result_part*(pi_bar[i][dim]**(end+1)-pi_bar[i][dim]**start)/ (pi_bar[i][dim]-1)
                        else:
                            start = component[dim]
                            # Formula for \sum_{start_k}^{\infty} rho^{n_k}
                            result_part = result_part*pi_bar[i][dim]**start/(1-pi_bar[i][dim])
                result += (model.Fbar.__getitem__(component+tuple([c_dim])))*result_part
        sum_result += weight[i]*result
    # Multiply the summation by the normalization constant
    return model.const*sum_result

model.obj_func = Objective(rule=obj_func_expr,sense=maximize)





################################################################################
#
# First define expressions for directed expression (A_i(n), B_i(n)) and
# undirected expression (Fbar(n),F(n), G(n))
#
################################################################################

# Define that D_u^t(n) = D_u^t(n) if u is type_plus, and D_u^t(n) = D^t_{-u}(n+u) if u is not type_plus.
# Therefore, the expressions are only defined for u in type_plus set.
# Here we generate the expressions A_u(n) and B_u(n).

def init_directed_expr(model,*args):
    coord = [a for a in args[0:dimen]]
    u = [a for a in args[dimen:2*dimen]]
    name = args[2*dimen]
    component = coordtoC(tuple(coord))
    if tuple(u) in transitions_init(component):
        if tuple(u) in model.type_plus:
            constant_term = getattr(model,name)[component + tuple(u+[-1])]
            linear_coefficients = [getattr(model,name)[component + tuple(u+[dim])] for dim in dimset]
            linear_term = summation(linear_coefficients,coord,index=range(dimen))
            return constant_term + linear_term
        else:
            u_opp =[-u[dim] for dim in model.dimset]
            fine_component = coordtoT(coord)
            new_component = TtoC(fine_component,u)
            new_state = [coord[dim]+u[dim] for dim in model.dimset]
            constant_term = getattr(model,name)[new_component + tuple(u_opp+[-1])]
            linear_coefficients = [getattr(model,name)[new_component + tuple(u_opp+[dim])] for dim in dimset]
            linear_term = summation(linear_coefficients,new_state,index=range(dimen))
            return -constant_term - linear_term
    else:
        return 0

model.directed_expr = Expression(model.fine_partcoords, model.type, model.dstring, initialize=init_directed_expr)



# Generate expressions F(n), G(n), Fbar(n)

def init_expr(model,*args):
    coord = [item for item in args[0:dimen]]
    component = coordtoC(tuple(coord))
    name = args[dimen]
    constant_term = getattr(model,name)[component+tuple([-1])]
    linear_coefficients = [getattr(model,name)[component+tuple([dim])] for dim in model.dimset]
    linear_term = summation(linear_coefficients,coord,index=model.dimset)
    return constant_term + linear_term

model.expr = Expression(model.fine_partcoords, model.ustring, initialize=init_expr)





#####################################################################
#
# Define expressions for directed shifted expressions (A_i(n+u), B_i(n+u))
# and undirected shifted expressions (F(n+e_i)...)
#
#####################################################################

# Expression for A_v(n+d) and B_v(n+d)

def init_directed_shifted_expr(model,*args):
    coord = [item for item in args[0:dimen]]
    u = [item for item in args[dimen:2*dimen]]
    d = [item for item in args[2*dimen:3*dimen]]
    v = [item for item in args[3*dimen:4*dimen]]
    name = args[4*dimen]
    fine_component = coordtoT(tuple(coord))
    new_component_d = TtoC(fine_component,d)
    if tuple(v) in transitions_init(new_component_d):
        if tuple(v) in model.type_plus:
            new_state = [coord[dim]+d[dim] for dim in model.dimset]
            constant_term = getattr(model,name)[new_component_d+tuple(v+[-1])]
            linear_coefficients = [getattr(model,name)[new_component_d+tuple(v+[dim])] for dim in model.dimset]
            linear_term = summation(linear_coefficients, new_state, index=model.dimset)
            return constant_term + linear_term
        else:
            v_opp = [-v[dim] for dim in model.dimset]
            new_component_dpv = TtoC(fine_component, tuple([d[dim]+v[dim] for dim in model.dimset]))
            new_state = [coord[dim]+d[dim]+v[dim] for dim in model.dimset]
            constant_term = getattr(model,name)[new_component_dpv+tuple(v_opp+[-1])]
            linear_coefficients = [getattr(model,name)[new_component_dpv+tuple(v_opp+[dim])] for dim in model.dimset]
            linear_term = summation(linear_coefficients, new_state, index=model.dimset)
            return - constant_term - linear_term
    else:
        return 0

model.directed_shifted_expr = Expression(model.fine_coords_type_plus_state, model.type, model.dstring, initialize=init_directed_shifted_expr)


# Expressions for F(n+u), Fbar(n+u),G(n+u). Some of the expressions may be useless.

def init_shifted_expr(model,*args):
    coord = [item for item in args[0:dimen]]
    u = [item for item in args[dimen:2*dimen]]
    name = args[2*dimen]
    fine_component = coordtoT(tuple(coord))
    component = coordtoC(tuple(coord))
    if tuple(u) in transitions_init(component):
        new_component = TtoC(fine_component,u)
        new_component_list = [new_component[dim] for dim in model.dimset]
        new_state = [coord[dim]+u[dim] for dim in model.dimset]
        linear_coefficients = [getattr(model,name)[tuple(new_component_list+[dim])] for dim in model.dimset]
        linear_term = summation(linear_coefficients, new_state, index=model.dimset)
        constant_term = getattr(model,name)[tuple(new_component_list+[-1])]
        return linear_term + constant_term
    else:
        return 0

model.shifted_expr = Expression(model.fine_partcoords, model.type_plus, model.ustring, initialize=init_shifted_expr)





###########################################################################################################
#
# Expressions for the constraints
#
###########################################################################################################

# Expressions for bounding the bias terms in the bounded dimensions.
# We need to make sure that the inequalities hold at all the fine coordinates.
# The "lower" part is to bound the bias from below and the "upper" part to bound from above.

def init_bound_bias_bd(model,*args):
    coord = [item for item in args[0:dimen]]
    u = [item for item in args[dimen:2*dimen]]
    bound = args[2*dimen]
    fine_component = coordtoT(tuple(coord))
    component = coordtoC(tuple(coord))
    transitions = transitions_init(component)
    if tuple(u) in transitions:
        complete_transitions = func_fine_comps_type_state(fine_component,u)
        result = model.shifted_expr.__getitem__(tuple(coord+u+["FF"])) - model.expr.__getitem__(tuple(coord+["FF"]))
        if bound == "lower":
            for d in complete_transitions:
                new_component = TtoC(fine_component,d)
                for v in model.type:
                    if v in transitions_init(new_component):
                        if v in model.type_plus:
                            result += -model.c.__getitem__(fine_component+tuple(u)+d+v)*model.directed_shifted_expr.__getitem__(tuple(coord+u)+d+v+tuple(["A"]))
                        else:
                            result += model.c.__getitem__(fine_component+tuple(u)+d+v)*model.directed_shifted_expr.__getitem__(tuple(coord+u)+d+v+tuple(["B"]))
            return -result - model.directed_expr.__getitem__(tuple(coord+u+["A"]))
        if bound == "upper":
            for d in complete_transitions:
                new_component = TtoC(fine_component,d)
                for v in model.type:
                    if v in transitions_init(new_component):
                        if v in model.type_plus:
                            result += model.c.__getitem__(fine_component+tuple(u)+d+v)*model.directed_shifted_expr.__getitem__(tuple(coord+u)+d+v+tuple(["B"]))
                        else:
                            result += -model.c.__getitem__(fine_component+tuple(u)+d+v)*model.directed_shifted_expr.__getitem__(tuple(coord+u)+d+v+tuple(["A"]))
            return result - model.directed_expr.__getitem__(tuple(coord+u+["B"]))
    else:
        return 0

model.bound_bias_bd = Expression(model.fine_partcoords, model.type_plus, model.bound, initialize=init_bound_bias_bd)



# For a fine component, if it is unbounded in one dimension, then we need to make sure that the slope in this direction is <= 0
# Below is the expression for the slope.
# The "lower" part is to bound the bias from below and the "upper" part to bound from above.

def init_bound_bias(model,*args):
    u = [item for item in args[0:dimen]]
    fine_component = [item for item in args[dimen:2*dimen]]
    unbddim = args[2*dimen]
    bound = args[2*dimen+1]
    component = TtoC(fine_component,tuple([0]*dimen))
    transitions = transitions_init(component)
    if tuple(u) in transitions:
        new_component_u = TtoC(fine_component,u)
        complete_transitions = func_fine_comps_type_state(fine_component,u)
        result = model.FF.__getitem__(new_component_u+tuple([unbddim])) - model.FF.__getitem__(component+tuple([unbddim]))
        if bound == "lower":
            for d in complete_transitions:
                new_component = TtoC(fine_component,d)
                for v in model.type:
                    if v in transitions_init(new_component):
                        if v in model.type_plus:
                            result += -model.c.__getitem__(tuple(fine_component+u)+d+v)*model.A.__getitem__(new_component+v+tuple([unbddim]))
                        else:
                            v_opp = [-v[dim] for dim in model.dimset]
                            new_component_dpv = TtoC(fine_component,tuple(d[dim]+v[dim] for dim in model.dimset))
                            result += -model.c.__getitem__(tuple(fine_component+u)+d+v)*model.B.__getitem__(new_component_dpv+tuple(v_opp+[unbddim]))
            return -result - model.A.__getitem__(component+tuple(u+[unbddim]))
        if bound == "upper":
            result = model.FF.__getitem__(new_component_u+tuple([unbddim])) - model.FF.__getitem__(component+tuple([unbddim]))
            for d in complete_transitions:
                new_component = TtoC(fine_component,d)
                for v in model.type:
                    if v in transitions_init(new_component):
                        if v in model.type_plus:
                            result += model.c.__getitem__(tuple(fine_component+u)+d+v)*model.B.__getitem__(new_component+v+tuple([unbddim]))
                        else:
                            v_opp = [-v[dim] for dim in model.dimset]
                            new_component_dpv = TtoC(fine_component,tuple(d[dim]+v[dim] for dim in model.dimset))
                            result += model.c.__getitem__(tuple(fine_component+u)+d+v)*model.A.__getitem__(new_component_dpv+tuple(v_opp+[unbddim]))
            return result - model.B.__getitem__(component+tuple(u+[unbddim]))
    else:
         return 0

model.bound_bias = Expression(model.type_plus, model.fine_comps_unbddim, model.bound, initialize=init_bound_bias)



# Expressions for the error bounds constraints.
# Again, we have constraints for bounded dimensions and unbounded dimensions

def init_bound_error_bd(model,*args):
    coord = [item for item in args[0:dimen]]
    bound = args[dimen]
    component = coordtoC(tuple(coord))
    transitions = transitions_init(component)
    result = model.expr.__getitem__(tuple(coord+["Fbar"])) - model.expr.__getitem__(tuple(coord+["FF"]))
    if bound == "lower":
        return 0
    if bound == "upper":
        for u in transitions:
            if u in model.type:
                if model.q.__getitem__(tuple(component)+u) >= 0:
                    if u in model.type_plus:
                        result += model.q.__getitem__(tuple(component)+u)*model.directed_expr.__getitem__(tuple(coord)+u+tuple(["B"]))
                    else:
                        result += model.q.__getitem__(tuple(component)+u)*model.directed_expr.__getitem__(tuple(coord)+u+tuple(["A"]))
                else:
                    if u in model.type_plus:
                        result += -model.q.__getitem__(tuple(component)+u)*model.directed_expr.__getitem__(tuple(coord)+u+tuple(["A"]))
                    else:
                        result += -model.q.__getitem__(tuple(component)+u)*model.directed_expr.__getitem__(tuple(coord)+u+tuple(["B"]))
        return result

model.bound_error_bd = Expression(model.fine_partcoords, model.bound, initialize=init_bound_error_bd)



def init_bound_error(model,*args):
    fine_component = [item for item in args[0:dimen]]
    component = TtoC(fine_component,tuple([0]*dimen))
    transitions = transitions_init(component)
    unbddim = args[dimen]
    bound = args[dimen+1]
    result = model.Fbar.__getitem__(component+tuple([unbddim])) - model.FF.__getitem__(component+tuple([unbddim]))
    if bound == "lower":
        return 0
    if bound == "upper":
        for u in transitions:
            if u in model.type:
                if model.q.__getitem__(tuple(component)+u) >= 0:
                    if u in model.type_plus:
                        result += model.q.__getitem__(tuple(component)+u)*model.B.__getitem__(component+u+tuple([unbddim]))
                    else:
                        u_opp = tuple([-u[dim] for dim in model.dimset])
                        new_component = TtoC(fine_component,u)
                        result += model.q.__getitem__(tuple(component)+u)*model.A.__getitem__(new_component+u_opp+tuple([unbddim]))
                else:
                    if u in model.type_plus:
                        result += -model.q.__getitem__(tuple(component)+u)*model.A.__getitem__(component+u+tuple([unbddim]))
                    else:
                        u_opp = tuple([-u[dim] for dim in model.dimset])
                        new_component = TtoC(fine_component,u)
                        result += - model.q.__getitem__(tuple(component)+u)*model.B.__getitem__(new_component+u_opp+tuple([unbddim]))
        return result

model.bound_error = Expression(model.fine_comps_unbddim, model.bound, initialize=init_bound_error)





#########################################################################################
#
# Define all the constraints by making the expression <= 0
#
#########################################################################################

# Bounds on bias terms

def init_constr_bound_bias_bd(model,*args):
    return model.bound_bias_bd.__getitem__(args) <= 0

model.Constr_bound_bias_bd = Constraint(model.fine_partcoords, model.type_plus, model.bound, rule=init_constr_bound_bias_bd)



def init_constr_bound_bias(model,*args):
    return model.bound_bias.__getitem__(args) <= 0

model.Constr_bound_bias = Constraint(model.type_plus, model.fine_comps_unbddim, model.bound, rule=init_constr_bound_bias)



# Bounds on errors

def init_constr_bound_error_bd(model,*args):
    return model.bound_error_bd.__getitem__(args) <= 0

model.Constr_bound_error_bd = Constraint(model.fine_partcoords, model.bound, rule=init_constr_bound_error_bd)



def init_constr_bound_error(model,*args):
    return model.bound_error.__getitem__(args) <= 0

model.Constr_bound_error = Constraint(model.fine_comps_unbddim, model.bound, rule=init_constr_bound_error)





############################
#
# Create instance
#
############################



instance = model.create(datafile)




############################
#
# Print
#
############################

#print("\n Abstract model: \n")
#model.pprint()

#print("\n INSTANCE: \n")
#instance.pprint()


############################
#
# Solve
#
############################

#print("\n PREPROCESS. \n")
#instance.preprocess()

print("\n SOLVE: \n")
results = optimizer.solve(instance)
instance.load(results)

print("\n SOLUTION: \n")
instance.pprint()
print(results)
