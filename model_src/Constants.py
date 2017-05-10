import math
import random

SEED = 1

GRAY = (0.662745,0.662745,0.662745)
RED = (0.545098,0,0)
BLUE = (0,0,0.545098)

class MyEnum():
    vals = ()
    enum_list = ()

    FAILURE = -1

    @classmethod
    def get_values(cls):
        return cls.vals

    @classmethod
    def get_enums(cls):
        return cls.enum_list

    @classmethod
    def string_to_enum(cls, string):
        if string in cls.vals:
            return cls.enum_list[cls.vals.index(string)]
        else:
            return cls.FAILURE

    @classmethod
    def enum_to_string(cls, enum):
        if enum in cls.enum_list:
            return cls.vals[cls.enum_list.index(enum)]
        else:
            return cls.FAILURE

class ProbSpread(MyEnum):
    RANDOM, POWERLAW = range(2)
    vals = ("rand", "powerlaw")
    enum_list = (RANDOM, POWERLAW)

class States(MyEnum):
    INIT_STATE, M_STATE, U_STATE, A_STATE = range(4)
    vals = ("R", "M", "U", "A")
    enum_list = (INIT_STATE, M_STATE, U_STATE, A_STATE)

class Domain(MyEnum):
    NONE, EQUAL_DEFAULT, USER_SET, BLEED_NONE, BLEED_USER_SET = range(5)
    vals = ("none", "equal", "set", "nobleed", "bleed_set")
    enum_list = (NONE, EQUAL_DEFAULT, USER_SET, BLEED_NONE, BLEED_USER_SET)

class Divisions(MyEnum):
    NONE, EQUAL_DEFAULT, USER_SET = range(3)
    vals = ("none", "equal", "set")
    enum_list = (NONE, EQUAL_DEFAULT, USER_SET)

class ProbConv():
    EQUAL_DEFAULT, MOD = range(2)
    vals = ("equal", "mod")
    enum_list = (EQUAL_DEFAULT, MOD)


def state_to_color(state):
    if state == States.A_STATE:
        return RED

    if state == States.U_STATE:
        return GRAY

    if state == States.M_STATE:
        return BLUE

def state_to_string(state):
    if state == States.A_STATE:
        return "A"

    if state == States.U_STATE:
        return "U"

    if state == States.M_STATE:
        return "M"

    if state == States.INIT_STATE:
        return "R"

def string_to_state(string):
    if string == "A":
        return States.A_STATE

    if string == "U":
        return States.U_STATE

    if string == "M":
        return States.M_STATE

    if string == "R":
        return States.INIT_STATE

def exponential(mean):
    return - mean * math.log(random.random())

