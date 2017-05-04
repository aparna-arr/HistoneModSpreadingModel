import sys
import random
import math
import operator

import matplotlib.pyplot as plt
import matplotlib.animation as animation

## CONSTANTS ##
SEED = 1
MEAN_TIME = 5
random.seed(SEED)

GRAY = (0.662745,0.662745,0.662745)
RED = (0.545098,0,0)
BLUE = (0,0,0.545098)
##

## FUNCTIONS ##

class Nucleosome:
    def __init__(self):
        self.state = 0

    def change_state(self,new_state):
        if new_state > self.state:
            self.state += 1
        elif new_state < self.state:
            self.state -= 1

class Event:
    def __init__(self):
        self.nucleosome = -1
        self.time = -1
        self.recruit = False
        self.recruit_from = -1
        self.rand_change = 0
       
    def generate(self, N_NUCLEOSOMES, A):
        self.time = exponential()
        self.nucleosome = random.randint(0,N_NUCLEOSOMES - 1)

        prob = random.random()
        
        if prob < A:
            # recruitment
            self.recruit = True
            self.recruit_from = random.randint(0,N_NUCLEOSOMES - 1)
        else:
            # random event
            next_state_prob = random.random()
            if next_state_prob < 1/3:
                self.rand_change = -1
            elif next_state_prob < 2/3:
                self.rand_change = 1
            else:
                self.rand_change = 0
        
class Chromatin:
    def __init__(self, N_NUCLEOSOMES, TOT_EVENTS, A):
        self.N_NUCLEOSOMES = N_NUCLEOSOMES
        self.TOT_EVENTS = TOT_EVENTS
        self.A = A
        self.nucleosomes = [Nucleosome() for x in range(self.N_NUCLEOSOMES)]
        self.events = []
        self.order = []

        self.total = {
            'M':0,
            'A':0,
            'U':self.N_NUCLEOSOMES
                }

        self.colors = [GRAY] * self.N_NUCLEOSOMES
        self.updated = []


    def animate_nucs(self):
        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.9,0.9])

        cols = math.ceil(math.sqrt(self.N_NUCLEOSOMES))
        rows = math.ceil(self.N_NUCLEOSOMES/cols)

        x_int = 1/cols
        y_int = 1/rows

        x_vals = []
        y_vals = []

        for i in range(cols):
            for j in range(rows):
                x_vals.append(i*x_int)
                y_vals.append(j*y_int)

        p = ax.scatter(x_vals, y_vals, facecolors = self.colors)

        def update_animation(i):
            self.run_event()

            if i % int(self.TOT_EVENTS/5) == 0:
                self.divide()

            p.set_facecolors(self.colors)
            

        an = animation.FuncAnimation(fig,update_animation,frames=self.TOT_EVENTS,interval=10,repeat=False)

        plt.show()

    def divide(self):
        for i in range(self.N_NUCLEOSOMES):
            if random.random() <= 0.5:
                prev = self.nucleosomes[i].state
                self.nucleosomes[i].state = 0
                self.update(prev,0)

                self.updated.append(i)
                self.colors[i] = GRAY


    def print_nucleosomes(self):
        n_string = ''

        for n in self.nucleosomes:
            if n.state == -1:
                n_string += 'M'
            elif n.state == 0:
                n_string += 'U'
            elif n.state == 1:
                n_string += 'A'

        print (n_string)

    def run(self):
        fp = open("results.txt", "w")
        fp.write('Event\tM\tU\tA\n')

        divisions = []
        event_count = 0
        while len(self.order) > 0:
            self.run_event()

            fp.write( str(event_count) + '\t' + str(self.total['M']/self.N_NUCLEOSOMES) + '\t' + str(self.total['U']/self.N_NUCLEOSOMES) + '\t' + str(self.total['A']/self.N_NUCLEOSOMES) + '\n')

            if event_count % int(self.TOT_EVENTS/5) == 0:
                self.divide()
                divisions.append(str(event_count))

            self.update_animation()

            event_count += 1

        fp.close()
        print("========== SUMMARY ==========")
        print("M:",self.total['M'], "U:",self.total['U'], "A:",self.total['A'])
        print(', '.join(divisions))

    def generate_events(self):
        index = 0
        tmp_order = []

        while len(self.events) < self.TOT_EVENTS:
            self.events.append(Event())
            self.events[index].generate(self.N_NUCLEOSOMES, self.A)
            curr_time = self.events[index].time
    
            tmp_order.append({'index' : index, 'time' : curr_time})

            index += 1

        self.order = [x['index'] for x in sorted(tmp_order, key=operator.itemgetter('time'), reverse=True)]


    def run_event(self):
        event_index = self.order.pop()
        curr_time = self.events[event_index].time
        curr_nucleosome = self.events[event_index].nucleosome
        
        state = 0

        if self.events[event_index].recruit:
            next_n = self.events[event_index].recruit_from
            state = self.nucleosomes[next_n].state
#            print("Recruit is TRUE, next_n:",next_n)
        else:
            state = self.events[event_index].rand_change
#            print("Recruit is FALSE")

        old_state = self.nucleosomes[curr_nucleosome].state
        self.nucleosomes[curr_nucleosome].change_state(state)
        new_state = self.nucleosomes[curr_nucleosome].state

            

#        print("old state:",old_state,"new_state:",new_state, "change_state:", state)

        self.updated.append(curr_nucleosome)
        if new_state == 1:
            self.colors[curr_nucleosome] = RED
        elif new_state == 0:
            self.colors[curr_nucleosome] = GRAY
        elif new_state == -1:
            self.colors[curr_nucleosome] = BLUE

        self.update(old_state, new_state)

        print("Time:", curr_time, "event index:", event_index, "curr_n", curr_nucleosome)
        print("M:",self.total['M'], "U:",self.total['U'], "A:",self.total['A'])
        self.print_nucleosomes()
    
    def update(self, old, new):
        if old == new:
            return
        
        if old == -1:
            self.total['M'] -= 1
        elif old == 0:
            self.total['U'] -= 1
        elif old == 1:
            self.total['A'] -= 1

        if new == -1:
            self.total['M'] += 1
        elif new == 0:
            self.total['U'] += 1
        elif new == 1:
            self.total['A'] += 1

def exponential( ):
    return - MEAN_TIME * math.log(random.random())

def main( ):

    if len(sys.argv) < 4:
        print("usage: <n_nucleosomes> <n_events> <F>", file=sys.stderr)
        sys.exit(1)

    N_NUCLEOSOMES = int(sys.argv[1])
    TOT_EVENTS = int(sys.argv[2])
    f = float(sys.argv[3])
    A = f/(f+1)

    print("INPUTS: N_NUCLEOSOMES:",N_NUCLEOSOMES,"TOT_EVENTS:",TOT_EVENTS,"A:",A)

    chromatin = Chromatin(N_NUCLEOSOMES, TOT_EVENTS, A)
    chromatin.generate_events()
#    chromatin.run()
    chromatin.animate_nucs()

main()

