class Event:
    def __init__(self, a, n):
        self.nuc_index = random.randint(0, n - 1)
        self.nuc_ptr = None
        self.time = -1
        self.recruit_from = -1
        self.rand_change = 0
        self.a = a
        self.n = n

        prob = random.random()

        if prob < self.a:
            self.feedback_event()
        else:
            self.random_event()
#    def feedback_event(self):
        

class Nucleosome:
    def __init__(self, init_state):
        if init_state == States.INIT_STATE:
            self.state = int(random.choice(
                [States.U_STATE, States.M_STATE, States.A_STATE]
                ))
        else:
            self.state = init_state

    def change_state(self, new):
        if new > self.state:
            self.state += 1
        elif new < self.state:
            self.state -= 1

class Chromatin:
    def __init__(self, input_dat):
        self.dat = input_dat
        self.nucleosomes = [ Nucleosome(input_dat['i']) for x in range(input_dat['n']) ]

    def feedback_event(self):


    def generate_event(self):
        prob = random.random()

        if prob < self.a:
            self.feedback_event()
        else:
            self.random_event()
        
