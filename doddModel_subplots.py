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
    def __init__(self, init_states):
        #self.state = 0

        if -1 <= init_states <= 1:
            self.state = init_states
        else:
            self.state = int(random.randint(-1,1))

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
    def __init__(self, N_NUCLEOSOMES, TOT_EVENTS, A, F, init_nuc_states, divisions, outfile):
        self.N_NUCLEOSOMES = N_NUCLEOSOMES
        self.TOT_EVENTS = TOT_EVENTS
        self.A = A
        self.F = F
        self.nucleosomes = [Nucleosome(init_nuc_states) for x in range(self.N_NUCLEOSOMES)]
        self.events = []
        self.order = []

        self.total = {'M':0, 'A':0, 'U':0}
        self.colors = []
        self.outfile = outfile

        self.divisionbool = divisions
        for n in self.nucleosomes:
            #print(n.state)
            if n.state == -1:
                self.total['M'] += 1
                self.colors.append(BLUE)
            elif n.state == 0:
                self.total['U'] += 1
                self.colors.append(GRAY)
            elif n.state == 1:
                self.total['A'] += 1
                self.colors.append(RED)


#        self.total = {
#            'M':0,
#            'A':0,
#            'U':self.N_NUCLEOSOMES
#                }

        self.updated = []


    def animate_nucs(self):
        fig = plt.figure()
#        ax = fig.add_axes([0.1,0.1,0.9,0.9])

        ax1 = fig.add_subplot(2,2,1)
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_title("N = " + str(self.N_NUCLEOSOMES) + " F = " + str(self.F))

        ax2 = fig.add_subplot(2,2,2)
        ax2.set_ylim([-5,105])
        ax2.set_xlim([-self.TOT_EVENTS/40,self.TOT_EVENTS + self.TOT_EVENTS/40])
        ax2.set_xlabel("Events")
        ax2.set_ylabel("% Nucleosomes")


        fig.tight_layout()

        cols = math.ceil(math.sqrt(self.N_NUCLEOSOMES))
        rows = math.ceil(self.N_NUCLEOSOMES/cols)

        x_int = 1/cols
        y_int = 1/rows

        x_vals = []
        y_vals = []

        count = 0
        for i in range(rows):
            for j in range(cols):
                if count >= 60:
                    break

                x_vals.append(j*x_int)
                y_vals.append(i*y_int)

                count += 1

        p = ax1.scatter(x_vals, y_vals, facecolors = self.colors)

        ### lines subplot ###
        lines = []
        lobj_1 = ax2.plot([],[],lw=2,color="red", label="A")[0]
        lobj_2 = ax2.plot([],[],lw=2,color="blue", label="M")[0]
        lines.append(lobj_1)
        lines.append(lobj_2)

        ax2.legend(loc="upper right")

        for line in lines:
            line.set_data([],[])


        line1_x = []
        line1_y = []

        line2_x = []
        line2_y = []

        curr = -1

        ## NEED THIS TRUST ME
        def init_an():
            return p, (*lines)

        def update_animation(i):
            # because the thing tries to loop EVEN WHEN I SET REPEAT = FALSE
            if i == 0 and len(self.order) == 0:
                return p,(*lines) 

            #print("ON EVENT:",i)

            self.run_event()

            if self.divisionbool and i % int(self.TOT_EVENTS/5) == 0 and i != 0:
                self.divide()
                lobj_3 = ax2.plot([],[],lw=1,ls="dotted",color="black")[0]
                lines.append(lobj_3)
                lines[len(lines)-1].set_data([i,i], [-5,105])

            line1_x.append(i)
            line1_y.append(self.total['A']/self.N_NUCLEOSOMES*100)

            line2_x.append(i)
            line2_y.append(self.total['M']/self.N_NUCLEOSOMES*100)

            lines[0].set_data(line1_x, line1_y)
            lines[1].set_data(line2_x, line2_y)

            p.set_facecolors(self.colors)
            
            return p, (*lines)


        an = animation.FuncAnimation(fig,update_animation,init_func=init_an,frames=self.TOT_EVENTS,interval=1,repeat=False, blit=True)

        
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=200, metadata=dict(artist='Me'), bitrate=1800)

        #plt.show()
        an.save(self.outfile, writer=writer)
        

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
        #print ("LEN ORDER IS:",len(self.order))
        event_index = self.order.pop()
        curr_time = self.events[event_index].time
        curr_nucleosome = self.events[event_index].nucleosome
        
        state = 0

        if self.events[event_index].recruit:
            next_n = self.events[event_index].recruit_from
            state = self.nucleosomes[next_n].state
            if state == 0:
                state = self.nucleosomes[curr_nucleosome].state
                # no change if the recruit-from nucleosome is unmodified
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

        #print("Time:", curr_time, "event index:", event_index, "curr_n", curr_nucleosome)
        #print("M:",self.total['M'], "U:",self.total['U'], "A:",self.total['A'])
        #self.print_nucleosomes()
    
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

    if len(sys.argv) < 7:
        print("usage: <n_nucleosomes> <n_events> <F> <init state: -2: random, -1: M, 0: U, 1: A> <divisions y/n: 0:n 1:y> <outfilename>", file=sys.stderr)
        sys.exit(1)

    N_NUCLEOSOMES = int(sys.argv[1])
    TOT_EVENTS = int(sys.argv[2])
    f = float(sys.argv[3])
    A = f/(f+1)
    initstate = int(sys.argv[4])

    out = sys.argv[6]

    divisions = False
    if int(sys.argv[5]) == 0:
        divisions = False
    elif int(sys.argv[5]) == 1:
        divisions = True

    print("INPUTS: N_NUCLEOSOMES:",N_NUCLEOSOMES,"TOT_EVENTS:",TOT_EVENTS,"A:",A, "init_state:",initstate, "divisions:",divisions)

    chromatin = Chromatin(N_NUCLEOSOMES, TOT_EVENTS, A, f, initstate, divisions, out)
    chromatin.generate_events()
#    chromatin.run()
    chromatin.animate_nucs()

main()

