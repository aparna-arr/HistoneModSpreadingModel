import getopt
import sys
import Constants
from MyEnum import ProbSpread, States, Domain, DomainBleed, Divisions, ProbConv

class InputError(Exception):
    def __init__(self, opt, arg, msg):
        self.opt = opt
        self.arg = arg
        self.msg = msg

class Help(Exception):
    pass

def usage():
    print("usage: python3 Main.py [OPTIONS]")
    print("\t-h, --help\n\t\tprint usage statement and exit")
    print("\t-n, --nucleosomes <INT>\n\t\tnumber of nucleosomes\n\t\t[default: 60]")
    print("\t-e, --events <INT>\n\t\tnumber of events\n\t\t[default: 10000]")
    print("\t-f, --Fval <FLOAT>\n\t\tnoise level\n\t\t[default: 1]")
    print("\t-i, --initstate <A, U, M, R>\n\t\tinitial state of nucleosomes. R is random.\n\t\t[default: U]")
    print("\t-d, --divisions\n\t\tinclude divisions in simulation\n\t\t[default: False]")
    print("\t-o, --outfile <STRING>\n\t\tprint to outfile instead of interactive\n\t\t[default: interactive mode]")

    print("Advanced Options")
    print("\t--prob-spread <rand, powerlaw>\n\t\tProbability distribution for spreading of modification\n\t\t[default: rand]")

    print("\t--domain <none, equal, set>\n\t\tadd static domains of either equal size with user-set number of domains or custom sizes\n\t\t[default: none]")
    print("\t--domain-equal <INT>\n\t\tadd static domains of the same size. Input number of domains.\n\t\t[default:2]")
    print("\t--domain-set <comma separated list of integers>\n\t\tadd static domains of different, user-set sizes. Input comma-sep domain sizes.\n\t\t[default:n_nucleosomes/2,n_nucleosomes/2]")

    print("\t--domainbleed <none, set>\n\t\twhether to allow domain bleedthrough\n\t\t[default: none]")
    print("\t--domainbleed-prob <FLOAT>\n\t\tprobability of domain bleed through\n\t\t[default: 0.05]")

    print("\t--divisions-num <INT>\n\t\tnumber of divisions attempts\n\t\t[default: 5]")

#    print("\t--prob-conv <equal,mod>\n\t\tWhether to change feedback-induced conversion rates.\n\t\t[default: equal]")
    print("\t--prob-conv-mod <4 comma seperated floats>\n\t\tfeedback-induced conversion rates.\n\t\tFormat: <U to M, A to U, U to A, M to U>\n\t\t[default: sequal rates]")

def test_int(string):
    num = int(string)

    return num

def test_float(string):
    num = float(string)

    return num

def test_state(string):
    state = ""
    if string in ("A", "U", "M", "R"):
        return string_to_state(string)
    else:
        raise ValueError

def test_emptystr(string):
    if string != "":
        return string
    else:
        raise ValueError

def test_enum(string, classname):
    if string in classname.get_values():
        return classname.string_to_enum(string)
    else:
        raise ValueError

def parse_input(opts):
    default_n = 60
    inputs = {
        'n':default_n,
        'e':10000,
        'f':1,
        'i':States.U_STATE,
        'd':Divisions.NONE,
        'o':"",
        'adv': {
            'prob_spread' : ProbSpread.RANDOM,
            'domain' : Domain.NONE,
            'domainbleed' : DomainBleed.NONE,
            'prob_conv' : ProbConv.EQUAL_DEFAULT
            },
        'data': {
            'divisions':5,
            'prob_conv':[1,1,1,1],
            'domains':default_n,
            'domainbleed':0
            }
            }

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            raise Help
        if opt in ("-n", "--nucleosomes"):
            try:
                inputs['n'] = test_int(arg)
            except ValueError:
                raise InputError(opt, arg, "requires int!")
        elif opt in ("-e", "--events"):
            try:
                inputs['e'] = test_int(arg)
            except ValueError:
                raise InputError(opt, arg, "requires int!")
        elif opt in ("-f", "--Fval"):
            try:
                inputs['f'] = test_float(arg)
            except ValueError:
                raise InputError(opt, arg, "requires float!")
        elif opt in ("-i", "--initstate"):
            try:
                inputs['i'] = test_state(arg)
            except ValueError:
                raise InputError(opt, arg, "requires \"A\", \"U\", \"M\", or \"R\"!")

        elif opt in ("-d", "--divisions"):
            inputs['d'] = Divisions.EQUAL_DEFAULT
        elif opt in ("-o", "--outfile"):
            try:
                inputs['o'] = test_emptystr(arg)
            except ValueError:
                raise InputError(opt, arg, "requires STRING argument!")
        elif opt == "--prob-spread":
            try:
                inputs['adv']['prob_spread'] = test_enum(arg, ProbSpread)
            except ValueError:
                raise InputError(opt, arg, "must be in [" + ", ".join(ProbSpread.get_values()) + "]")
        elif opt == "--domain":
            try:
                inputs['adv']['domain'] = test_enum(arg, Domain)
            except ValueError:
                raise InputError(opt, arg, "must be in [" + ", ".join(Domain.get_values()) + "]")
        elif opt == "--domainbleed":
            try:
                inputs['adv']['domainbleed'] = test_enum(arg, DomainBleed)
            except ValueError:
                raise InputError(opt, arg, "must be in [" + ", ".join(DomainBleed.get_values()) + "]")
        elif opt == "--domain-equal":
            try:
                inputs['adv']['domain'] = Domain.EQUAL_DEFAULT
                inputs['data']['domains'] = test_int(arg) 
            except ValueError:
                raise InputError(opt, arg, "requires int!")
        elif opt == "--domain-set":
            raise RuntimeError("Unimplemented option:[",opt,"]")
        elif opt == "--domainbleed-prob":
            try:
                inputs['adv']['domainbleed'] = DomainBleed.USER_SET
                inputs['data']['domainbleed'] = test_float(arg)
            except ValueError:
                raise InputError(opt, arg, "requires float!")
        elif opt == "--divisions-num":
            try:
                inputs['d'] = Divisions.USER_SET
                inputs['data']['divisions'] = test_int(arg)
            except ValueError:
                raise InputError(opt, arg, "requires int!")
        elif opt == "--prob-conv-mod":
            raise RuntimeError("Unimplemented option:[",opt,"]")
        else:
            raise InputError(opt, arg, "unrecognized opt!")

    return inputs

def display_inputs(inputs):
    print("==========## INPUTS ##==========")
    print("Nucleosomes:", inputs['n'])
    print("Events:", inputs['e'])
    print("F-value:", inputs['f'])
    print("Init State:", States.enum_to_string(inputs['i']))
    print("Divisions:", inputs['d'])
    print("Outfile:", inputs['o'])
    print("================================")


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hn:e:f:i:do:", [
            "help",
            "nucleosomes=", 
            "events=",
            "Fval=",
            "initstate=",
            "divisions",
            "outfile=",
            "prob-spread=",
            "domain=",
            "domain-equal=",
            "domain-set=",
            "domainbleed=",
            "domainbleed-prob=",
            "divisions-num=",
            "prob-conv-mod="
            ])
    except getopt.GetoptError as err:
        print(str(err), file=sys.stderr)  
        usage()
        sys.exit(2)

    try:
        inputs = parse_input(opts)
    except InputError as e:
        print(type(e).__name__ + ":", "Opt [", e.opt, "] arg [", e.arg, "]:", e.msg, file=sys.stderr)
        usage()
        sys.exit(2)
    except Help:
        usage()
        sys.exit(2)
    except RuntimeError as e:
        print(type(e).__name__ + ":", e.arg)
        sys.exit(2)

    display_inputs(inputs) 

main()
