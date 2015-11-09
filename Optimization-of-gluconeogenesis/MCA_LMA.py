# coding: utf-8
__author__ = 'Nonchalant Dave'
import pysces
import numpy as np
from lmfit import minimize, Parameters

# Instantiate model object
mod = pysces.model('kouril3')

# Variable declarations
enzyme_total = mod.protGAPDH+mod.protALDPase+mod.protPGK+mod.protTIM
all_parameters = mod.parameters
enzyme_names = []
Params = []
F6ps = []

# Extract enzyme names from the model.
for word in all_parameters:
    if word.startswith('prot'):
        enzyme_names.append(word)
print "The system will be optimized by altering the relative enzyme activities of the following enzymes: ",enzyme_names


#Run the pysces model simulation and return max F6P value
def simulate(p,enzyme_total):
	#parameter values are normalized to add up to total enzyme constraint
    mod.protALDPase = enzyme_total*p[0]/(p[0]+p[1]+p[2]+p[3])
    mod.protGAPDH = enzyme_total*p[1]/(p[0]+p[1]+p[2]+p[3])
    mod.protPGK = enzyme_total*p[2]/(p[0]+p[1]+p[2]+p[3])
    mod.protTIM = enzyme_total*p[3]/(p[0]+p[1]+p[2]+p[3])
                  
    #Run the PYSCES simulation for 200 seconds
    mod.doSim(end=200.0, points=200.0)
   
    f6p_all = mod.data_sim.getSimData('f6p')
    f6p_max = max([col[1] for col in f6p_all])
    return round(f6p_max, 5)

#Objective function to be minimised using LMA
def objFunction(pi):
    #The list comprehension below is used as an ultra efficient way of storing a list of parameters from the parameter object 'pi'
    param_list = [pi[i].value for i in ['protALDPase','protGAPDH','protPGK','protTIM']]
    etot = 0.04995
    res = simulate(param_list,etot)
    result = [res]
    
    return np.resize(np.array([result[0]**(-1)]), 4) #note that the inverse of the maximum is returned for minimization

'''Monte Carlo Algorithm'''
Generate 1000 sets of random parameters between 0 and 1 and run the simulation
for iteration in range(2000):
    p = np.random.uniform(0,1,size=4)
    p = enzyme_total*(p/(p[0]+p[1]+p[2]+p[3]))
    
    #Add the random parameter set 'p' to a permanent list of all random parameter sets
    Params.append(p)

    f6p = simulate(p,enzyme_total)
    #add the resulting F6P value to a permanent list of parameters
    F6ps.append(f6p)

#Sorting F6P values and rearranging parameters to match accordingly
F6ps = np.asarray(F6ps)
Params = np.asarray(Params)
indices = F6ps.argsort()
Params = Params[indices]
F6ps = F6ps[indices]

'''Levenberg-Marquardt Algorithm'''
#setting initial parameters to be used for LM simulation
pi = Parameters()
pi.add('protALDPase', value=Params[-1,0], min=0.0, max=1.0)
pi.add('protGAPDH', value=Params[-1,1], min=0.0, max=1.0)
pi.add('protPGK', value=Params[-1,2], min=0.0, max=1.0)
pi.add('protTIM', value=Params[-1,3], min=0.0, max=1.0)

#Run the LM algorithm using lmfit:
result = minimize(objFunction,pi)
print '\n',"The following is the set of enzyme concentrations which result in optimal production of F6P for the reconstituted system:",'\n',result.params
print '\n',"The above parameters result in a theoretical F6P concentration of:",1/(result.residual[0]),"ug/mL"

