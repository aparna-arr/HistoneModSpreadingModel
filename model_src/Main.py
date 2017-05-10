import getopt
import sys
import Constants

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
    print("\t--domain <<num domains>, <INT comma-sep list of number nucleosomes> >\n\t\tadd static domains of either equal size with user-set number of domains or custom sizes\n\t\t[default: none]")
    print("\t--div <INT>\n\t\tnumber of divisions attempts\n\t\t[default: 5]")
    print("\t--prob-conv <4 comma seperated floats>\n\t\tfeedback-induced conversion rates. Format: <U to M, A to U, U to A, M to U>\n\t\t[default: sequal rates]")

def test_int(string):
    num = int(string)

    return num

def test_float(string):
    num = float(string)

    return num

def test_state(string):
    state = ""
    if string in ("A", "U", "M", "R"):
        return Constants.string_to_state(string)
    else:
        raise ValueError

def test_emptystr(string):
    if string != "":
        return string
    else:
        raise ValueError

def parse_input(opts):
    inputs = {
        'n':60,
        'e':10000,
        'f':1,
        'i':Constants.States.U_STATE,
        'd':Constants.Divisions.NONE,
        'o':"",
        'adv': {
            'prob_spread' : Constants.ProbSpread.RANDOM,
            'domain' : Constants.Domain.NONE,
            'prob_conv' : Constants.ProbConv.EQUAL_DEFAULT
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
            inputs['d'] = Constants.Divisions.EQUAL_DEFAULT
        elif opt in ("-o", "--outfile"):
            try:
                inputs['o'] = test_emptystr(arg)
            except ValueError:
                raise InputError(opt, arg, "requires STRING argument!")
        elif opt == "--prob-spread":
            if arg in ("rand", "powerlaw"):
                inputs['adv']['prob_spread'] = str_to_ProbSpread(arg)



        else:
            raise InputError(opt, arg, "unrecognized opt!")

    return inputs

def display_inputs(inputs):
    print("==========## INPUTS ##==========")
    print("Nucleosomes:", inputs['n'])
    print("Events:", inputs['e'])
    print("F-value:", inputs['f'])
    print("Init State:", Constants.state_to_string(inputs['i']))
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
            "outfile="
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

    print(Constants.ProbSpread.string_to_enum("powerlaw"))
    display_inputs(inputs) 

main()
