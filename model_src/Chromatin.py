import Constants
import random
import numpy as np
import operator
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
from MyEnum import ProbSpread, States, Domain, DomainBleed, Divisions, ProbConv

#class Event:
#    def __init__(self, a, n):
#        self.nuc_index = random.randint(0, n - 1)
#        self.nuc_ptr = None
#        self.time = -1
#        self.recruit_from = -1
#        self.rand_change = 0
#        self.a = a
#        self.n = n
#
#        prob = random.random()
#
#        if prob < self.a:
#            self.feedback_event()
#        else:
#            self.random_event()
##    def feedback_event(self):
        

class Nucleosome:
    def __init__(self, init_state, left, right):
        if init_state == States.INIT_STATE:
            self.state = int(random.choice(
                [States.U_STATE, States.M_STATE, States.A_STATE]
                ))
        else:
            self.state = init_state

        # domain limits
        self.left_limit = left
        self.right_limit = right

    def change_state(self, new):
        if new > self.state:
            self.state += 1
        elif new < self.state:
            self.state -= 1

class Chromatin:
    def __init__(self, input_dat):
        self.dat = input_dat
        self.events = []
#        self.order = []
        self.totals = {States.M_STATE:0, States.A_STATE:0, States.U_STATE:0}
        self.colors = []

        self.nucleosomes = []

        self.TIME = 0
        self.time_array = []

        if input_dat['adv']['domain'] == Domain.NONE:
            self.nucleosomes = [ Nucleosome(input_dat['i'], 0, input_dat['n']-1) for x in range(input_dat['n']) ]
        elif input_dat['adv']['domain'] == Domain.EQUAL_DEFAULT:
            domain_sizes = math.floor(input_dat['n'] / input_dat['data']['domains'])
            #print ("domain sizes:", domain_sizes)
            # last domain gets the leftover
            #domain_leftover = input_dat['n'] % input_dat['data']['domains']

            num_d = input_dat['data']['domains']
            print("num_d:", num_d)

            for i in range(input_dat['n']):
                low_index = math.floor(i / domain_sizes) * domain_sizes
                high_index = low_index + domain_sizes - 1

                if high_index >= input_dat['n']:
                    high_index = input_dat['n'] - 1

                self.nucleosomes.append(Nucleosome(input_dat['i'], low_index, high_index))
#            self.nucleosomes = [ Nucleosome(
#                input_dat['i'], 
#                math.floor(i / num_d) * domain_sizes, 
#                math.floor(i / num_d) * domain_sizes + \
#                        math.floor(math.ceil(i / num_d) / input_dat['n']) * domain_leftover + \
#                        (1 - math.floor(math.ceil(i / num_d))) * math.ceil(i / num_d) * domain_sizes - 1
#                ) for i in range(input_dat['n']) ]
        
        for nuc in self.nucleosomes:
            self.colors.append(Constants.state_to_color(nuc.state))
            self.totals[nuc.state] += 1

    ## EVENTS ##
    def pick_random_nuc(self, curr_nuc):

        left = self.nucleosomes[curr_nuc].left_limit 
        right = self.nucleosomes[curr_nuc].right_limit

        if self.dat['adv']['domainbleed'] != DomainBleed.NONE:
            prob = random.random()

            if prob < self.dat['data']['domainbleed']:
                # NOTE only allowing bleedthrough into adjacent domains right now
                if left > 0:
                    left = self.nucleosomes[curr_nuc - 1].left_limit
                if right < self.dat['n'] - 1:
                    left = self.nucleosomes[curr_nuc + 1].right_limit

        if self.dat['adv']['prob_spread'] == ProbSpread.RANDOM :
            #print("left:", left, "right:",right, "nuc", curr_nuc)
            seq = list(range(left,curr_nuc)) + list(range(curr_nuc+1, right))
            return random.choice(seq)
        elif self.dat['adv']['prob_spread'] == ProbSpread.POWERLAW :
            left = curr_nuc - self.nucleosomes[curr_nuc].left_limit  
            right = self.nucleosomes[curr_nuc].right_limit - curr_nuc

            #dist = Constants.truncated_power_law(1, left, right)
            #result = dist.rvs(size=1) + curr_nuc
            result = Constants.truncated_power_law(1, left, right)
            #print("limits: left", left, "right:", right, "result",result)
            return int(result + curr_nuc)


    def feedback_event(self, curr_nuc):
        return self.pick_random_nuc(curr_nuc)

#    def random_event(self):
#        return -1 


    def generate_event(self):
        curr_nuc = random.randint(0, self.dat['n'] - 1)
        prob = random.random()
        a = self.dat['f']/(self.dat['f'] + 1)

        next_event = -1 # random event

        if prob < a:
            next_event = self.feedback_event(curr_nuc)

        self.events.append({'nuc' : curr_nuc, 'event' : next_event})
        
        #return len(self.events) - 1

    def generate_all(self):
#        tmp_order = []
        
        while len(self.events) < self.dat['e']:
            events_index = self.generate_event()

#            beg_time = np.random.exponential(1/Constants.AVG_CR_TIME)
#            tmp_order.append({'index':events_index, 'time':beg_time})

#        self.order = [ x for x in sorted(tmp_order, key=operator.itemgetter('time'), reverse=True) ]

    def run_event(self, t):
#        event_time_dict = self.order.pop()
#        index = event_time_dict['index']
#        beg_time = event_time_dict['time']

#        event_info = self.events[index]
        event_info = self.events.pop(0)
        curr_nuc = event_info['nuc']
        event = event_info['event']

        
#        t_next = np.random.exponential(Constants.AVG_CR_TIME)

        t_next = -1

#        print("========")
#        print("T_NEXT", t_next)
#        print("T", t)
#        curr_t = t + t_next
#        print("CURR_T", curr_t)
#        print("========")

        old_state = self.nucleosomes[curr_nuc].state
        new_state = old_state

        if event == -1:
            # random event
            possible_states = {States.A_STATE, States.U_STATE, States.M_STATE}
            options = list(possible_states.difference({old_state}))

            prob = random.random()
#            if prob < 1/3:
            if prob < 1/2:
                new_state = options[0]
                t_next = np.random.exponential(Constants.get_rate(old_state, new_state))
#            elif prob < 2/3:
            else:
                new_state = options[1]
                t_next = np.random.exponential(Constants.get_rate(old_state, new_state))
        else:
            recruit_from = event
            rec_state = self.nucleosomes[recruit_from].state

            if rec_state == States.U_STATE:
                # want the rate for a no-attempt
#                t_next = np.random.exponential(Constants.NOTHING_HAPPENED)
                t_next = 0
                new_state = old_state
            else:
                new_state = rec_state
                t_next = np.random.exponential(Constants.get_rate(old_state, new_state))


        if t_next < 0:
            print("ERROR t_next < 0!!!")

        curr_t = t + t_next
        
        #duration_time = Constants.AVG_CR_TIME/10 # JUST FOR NOW
        #if duration_time + beg_time > curr_time
            
        if new_state == old_state :
            self.generate_event()
#            print("PASSING CURR_T AS T", curr_t)
            return self.run_event(curr_t)
#            print("FROM IF: returning ", curr_time)
#            return curr_time
        else:
            self.nucleosomes[curr_nuc].change_state(new_state)
            new_assigned = self.nucleosomes[curr_nuc].state
            self.colors[curr_nuc] = Constants.state_to_color(new_assigned)
            self.update(old_state, new_assigned)
#            print("FROM ELSE: returning ", curr_t)
            return curr_t
            #        if new_state != old_state and new_state != States.U_STATE:
            
    
    ## END EVENTS ##
    def divide(self):
        for i in range(self.dat['n']):
            if random.random() <= 0.5:
                prev = self.nucleosomes[i].state
                self.nucleosomes[i].state = States.U_STATE
                self.update(prev, States.U_STATE)
                self.colors[i] = Constants.state_to_color(States.U_STATE)

    def update(self, old, new):
        if old == new:
            return

        self.totals[old] -= 1
        self.totals[new] += 1

    ## Animation ##
    def animate(self):
        # base figure
        fig = plt.figure()
        # axes
        # figure 1: nucleosomes
        ax1 = fig.add_subplot(2,2,1)
        ax1.set_xticks([])
        ax1.set_yticks([])

        titlestr = "N = " + str(self.dat['n']) + " F = " + str(self.dat['f'])
        ax1.set_title(titlestr)

        # figure 2: proportion vs events
        ax2 = fig.add_subplot(2,2,2)
        ax2.set_ylim([-5,105])
        ax2.set_xlim(self.dat['e']/40, self.dat['e'] + self.dat['e']/40)
        ax2.set_xlabel("Events")
        ax2.set_ylabel("% Nucleosomes")

        ax3 = fig.add_subplot(2,2,3)
        # adjust spacing between figures
        fig.tight_layout()

        # calculate size of square for nucleosomes in fig 1
        cols = math.ceil(math.sqrt(self.dat['n']))
        rows = math.ceil(self.dat['n']/cols)

        x_int = 1/cols
        y_int = 1/rows

        x_vals = []
        y_vals = []

        count = 0
        for i in range(rows):
            for j in range(cols):
                if count >= self.dat['n']:
                    break
                
                x_vals.append(j * x_int)
                y_vals.append(i * y_int)
                count += 1

        # initialize fig1 as a scatterplot
        scat = ax1.scatter(x_vals, y_vals, facecolors = self.colors)

        # initialize fig2 as a lineplot
        lines = []
        lines.append(ax2.plot([], [], lw = 2, color = "red", label = "A")[0])
        lines.append(ax2.plot([], [], lw = 2, color = "blue", label = "M")[0])
        ax2.legend(loc = "upper right")

        for line in lines:
            line.set_data([],[])

        linex = []
        liney = [[], []]

        def init_an():
            '''initialize animation: this is needed'''
            return scat, (*lines)

        def update_an(i):
            '''update function for animation'''
#            if i == 0 or len(self.order) == 0:
            if i == 0 or len(self.events) == 0:
                return scat, (*lines)

#            print("PASSING self.TIME as t:", self.TIME)
            self.TIME = self.run_event(self.TIME)
            self.time_array.append(self.TIME)
            #print("Frame:",i)
            #print("Time is:", self.TIME)
#            print("After event, self.TIME is", self.TIME)

            if self.dat['d'] != Divisions.NONE and i != 0:
                div = self.dat['data']['divisions']
                if i % int(self.dat['e'] / div) == 0:
                    self.divide()
                    lines.append(ax2.plot([], [], lw = 1, ls = "dotted", color = "black")[0])
                    lines[len(lines) - 1].set_data([i, i], [-5, 105])

            linex.append(i)
            liney[0].append(self.totals[States.A_STATE] / self.dat['n'] * 100)
            liney[1].append(self.totals[States.M_STATE] / self.dat['n'] * 100)

            lines[0].set_data(linex, liney[0])
            lines[1].set_data(linex, liney[1])

            scat.set_facecolors(self.colors)

            return scat, (*lines)

        anim = animation.FuncAnimation(fig, update_an, init_func = init_an, frames = self.dat['e'], interval = 1, repeat = False, blit=True)

        if self.dat['o'] != "":
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps = 200, metadata = dict(artist = 'Me'), bitrate = 1800)
            anim.save(self.dat['o'] + ".mp4", writer = writer)
            #print(self.time_array)
            fp = open(self.dat['o'] + ".times", "w")
            string = '\n'.join(str(x) for x in self.time_array)
            fp.write(string)
            fp.close()
        else:
            plt.show()

