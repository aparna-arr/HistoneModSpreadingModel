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
    NONE, EQUAL_DEFAULT, USER_SET = range(3)
    vals = ("none", "equal", "set")
    enum_list = (NONE, EQUAL_DEFAULT, USER_SET)

class DomainBleed(MyEnum):
    NONE, USER_SET = range(2)
    vals = ("none", "set")
    enum_list = (NONE, USER_SET)

class Divisions(MyEnum):
    NONE, EQUAL_DEFAULT, USER_SET = range(3)
    vals = ("none", "equal", "set")
    enum_list = (NONE, EQUAL_DEFAULT, USER_SET)

class ProbConv(MyEnum):
    EQUAL_DEFAULT, MOD = range(2)
    vals = ("equal", "mod")
    enum_list = (EQUAL_DEFAULT, MOD)
