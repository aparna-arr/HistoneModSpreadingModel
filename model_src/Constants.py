import math
import random

SEED = 1

GRAY = (0.662745,0.662745,0.662745)
RED = (0.545098,0,0)
BLUE = (0,0,0.545098)

def state_to_color(state):
    if state == States.A_STATE:
        return RED

    if state == States.U_STATE:
        return GRAY

    if state == States.M_STATE:
        return BLUE

def exponential(mean):
    return - mean * math.log(random.random())

