import math
import random
import numpy as np
import scipy.stats as stats
from MyEnum import States

RECRUIT_TIMESTEPS = 50 # duration
RECRUIT_INIT_TIME = 20 # start at t = 20
RECRUIT_N_NUCS = 50 # number of nucleosomes recruited to
RECRUIT_CONV_TO = States.M_STATE


# number of events in a timestep determines prob of nothing
# RANDOM NUMBER RIGHT NOW
PROB_NOTHING = 0.25


#NOTHING_HAPPENED = 0.05
#CR_U_to_A = 0.1
#CR_A_to_U = 0.7
#CR_U_to_M = 0.1
#CR_M_to_U = 0.8
CR_U_to_A = 0.5
CR_A_to_U = 0.5
CR_U_to_M = 0.5
CR_M_to_U = 0.5

POWER = 2 # arbitrary number

def get_rate(old, new):
    if old == new:
        if old == States.M_STATE:
            return CR_U_to_M
        
        if old == States.A_STATE:
            return CR_U_to_A

        if old == States.U_STATE:
            return 0
        # can't have a U - to - U
    else:
        if old == States.M_STATE:
            return CR_M_to_U

        elif old == States.A_STATE:
            return CR_A_to_U

        elif old == States.U_STATE:
            if new == States.M_STATE:
                return CR_U_to_M

            if new == States.A_STATE:
                return CR_U_to_A

SEED = 1

random.seed(SEED)
np.random.seed(SEED)

GRAY = (0.662745,0.662745,0.662745)
RED = (0.545098,0,0)
BLUE = (0,0,0.545098)

def truncated_power_law(power, limit_neg, limit_pos):
    prob_direction = random.random()
    neg_prob = limit_neg / (limit_neg + limit_pos)
    if prob_direction < neg_prob:
        #print("NEGATIVE DIRECTION")
        # negative direction
        # limit_neg must be a positive number of nucs in neg direction
        x = np.arange(1, limit_neg+1, dtype='float')
        pmf = 1/x**power
        pmf /= pmf.sum()
        #print("neg:",limit_neg)
        dist = stats.rv_discrete(values=(range(1, limit_neg+1), pmf))
        return 0 - dist.rvs(size=1)
    else:
        #print("POSITIVE DIRECTION")
        # positive direction
        x = np.arange(1, limit_pos+1, dtype='float')
        pmf = 1/x**power
        pmf /= pmf.sum()
        dist = stats.rv_discrete(values=(range(1, limit_pos+1), pmf))
        return dist.rvs(size=1)

# use like this:
#a, m = 2, 10
#d = truncated_power_law(a=a, m=m)
#
#N = 10**4
#sample = d.rvs(size=N)

def state_to_color(state):
    if state == States.A_STATE:
        return RED

    if state == States.U_STATE:
        return GRAY

    if state == States.M_STATE:
        return BLUE

def exponential(mean):
    return - mean * math.log(random.random())

