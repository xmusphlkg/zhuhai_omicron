import covasim as cv
import numpy as np
import pandas as pd


def vaccinate_by_age(sim):
    child = cv.true(sim.people.age <= 18) # cv.true() returns indices of people matching this condition, i.e. people under 50
    adult = cv.true((sim.people.age > 18) * (sim.people.age < 65)) # Multiplication means "and" here
    old    = cv.true(sim.people.age >= 65)
    inds = sim.people.uid # Everyone in the population -- equivalent to np.arange(len(sim.people))
    vals = np.ones(len(sim.people)) # Create the array
    vals[child] = 0.15
    vals[adult] = 0.62
    vals[old] = 0.66
    output = dict(inds=inds, vals=vals)
    return output

vaccine = cv.simple_vaccine(days=0, rel_sus=0.65, rel_symp=0.05, subtarget=vaccinate_by_age)


pars = dict(
)


sim1 = cv.Sim(label='Baseline')
sim2 = cv.Sim(interventions=vaccine, label='With age-targeted vaccine')
msim = cv.parallel(sim1, sim2)
msim.plot()
