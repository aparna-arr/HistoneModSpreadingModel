import Constants
import random
import numpy as np
import operator
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.stats as stats
import math
from MyEnum import ProbSpread, States, Domain, DomainBleed, Divisions, ProbConv

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
        self.totals = {
                States.M_STATE:0, 
                States.A_STATE:0, 
                States.U_STATE:0
                }

        self.colors = []
        self.nucleosomes = []

        self.TIME = 0
        self.time_array = []

        self.prob_mat = np.zeros(shape=(input_dat['n'],input_dat['n']))
        self.M_mat = np.zeros(shape=(input_dat['n']))
        self.A_mat = np.zeros(shape=(input_dat['n']))

        # run initialization functions
        self.init_nucs(input_dat['adv']['domain'], input_dat['n'], input_dat['i'], input_dat['data']['domains'])

        self.init_prob_mat(input_dat['n'], input_dat['adv']['domainbleed'], input_dat['data']['domainbleed'])

        self.init_colors_and_state_mats(input_dat['n'])


    ## init functions ##
    def init_colors_and_state_mats(self, n):
        for i in range(n):
            curr_state = self.nucleosomes[i].state
            self.colors.append(Constants.state_to_color(curr_state))
            self.totals[curr_state] += 1

            if curr_state == States.M_STATE:
                self.M_mat[i] = 1
            if curr_state == States.A_STATE:
                self.A_mat[i] = 1

    def init_prob_mat(self, n_nucs, db_enum, db_val):
        for col in range(n_nucs):
            left = self.nucleosomes[col].left_limit
            right = self.nucleosomes[col].right_limit
            
            for row in range(n_nucs):
                if col == row:
                    continue
                
                prob = self.calc_prob(col, row, right, left, db_enum, db_val, n_nucs)

                if prob == -1:
                    continue
                elif prob == -2:
                    break

                self.prob_mat[col,row] = prob

    def init_nucs(self, domain_enum, n, initstate, num_domains):
        if domain_enum == Domain.NONE:
            self.nucleosomes = [ 
                    Nucleosome(initstate, 0, n) for x in range(n) 
                    ]
        elif domain_enum == Domain.EQUAL_DEFAULT:
            domain_sizes = math.floor(n / num_domains)
            num_d = num_domains

            for i in range(n):
                low_index = math.floor(i / domain_sizes) * domain_sizes
                high_index = low_index + domain_sizes - 1

                if high_index >= n :
                    high_index = n - 1

                self.nucleosomes.append(
                        Nucleosome(initstate, low_index, high_index)
                        )


    def calc_prob(self, col, row, right, left, domainbleed_enum, domainbleed_val, n):
        if domainbleed_enum == DomainBleed.NONE:
            if row < left:
                return -1 # continue
            elif row > right:
                return -2 # break
            
            if row < col:
                return stats.powerlaw.ppf( (col - row) / (col - left), Constants.POWER)
            else:
                return stats.powerlaw.ppf( (row - col) / (right - col), Constants.POWER)
            
        else:
            extreme_left = left
            extreme_right = right

            if left > 0:
                extreme_left = self.nucleosomes[left - 1].left_limit

            if right < n - 1:
                extreme_right = self.nucleosomes[right + 1].right_limit

            if row < extreme_left:
                return -1 # continue

            if row > extreme_right:
                return -2 # break

            prob = 0

            if row < left:
                prob = domainbleed_val * stats.powerlaw.ppf( 
                        (col - row) / (col - extreme_left), Constants.POWER
                        )
            elif row > right:
                prob = domainbleed_val * stats.powerlaw.ppf(
                        (row - col) / (extreme_right - col), Constants.POWER
                        )
            elif row < col:
                prob = stats.powerlaw.ppf( 
                        (col - row) / (col - left), Constants.POWER
                        )
            elif row > col:
                prob = stats.powerlaw.ppf( 
                        (row - col) / (right - col), Constants.POWER
                        )

            return prob

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
                    right = self.nucleosomes[curr_nuc + 1].right_limit

        if self.dat['adv']['prob_spread'] == ProbSpread.RANDOM :
            seq = list(range(left,curr_nuc)) + list(range(curr_nuc+1, right))
            return random.choice(seq)
        elif self.dat['adv']['prob_spread'] == ProbSpread.POWERLAW :
            left = curr_nuc - self.nucleosomes[curr_nuc].left_limit  
            right = self.nucleosomes[curr_nuc].right_limit - curr_nuc
            result = Constants.truncated_power_law(1, left, right)
            return int(result + curr_nuc)


    def feedback_event(self, curr_nuc):
        return self.pick_random_nuc(curr_nuc)

    def generate_event(self):
        curr_nuc = random.randint(0, self.dat['n'] - 1)
        prob = random.random()
        a = self.dat['f']/(self.dat['f'] + 1)

        next_event = -1 # random event

        if prob < a:
            next_event = self.feedback_event(curr_nuc)

        self.events.append({'nuc' : curr_nuc, 'event' : next_event})

    def generate_all(self):
        while len(self.events) < self.dat['e']:
            events_index = self.generate_event()

    def run_event(self, t):
        event_info = self.events.pop(0)
        curr_nuc = event_info['nuc']
        event = event_info['event']

        t_next = -1
        old_state = self.nucleosomes[curr_nuc].state
        new_state = old_state

        if event == -1:
            # random event
            possible_states = {States.A_STATE, States.U_STATE, States.M_STATE}
            options = list(possible_states.difference({old_state}))

            prob = random.random()

            if prob < 1/2:
                new_state = options[0]
                t_next = np.random.exponential(Constants.get_rate(old_state, new_state))
            else:
                new_state = options[1]
                t_next = np.random.exponential(Constants.get_rate(old_state, new_state))
        else:
            recruit_from = event
            rec_state = self.nucleosomes[recruit_from].state

            if rec_state == States.U_STATE:
                # want the rate for a no-attempt
                t_next = 0
                new_state = old_state
            else:
                new_state = rec_state
                t_next = np.random.exponential(Constants.get_rate(old_state, new_state))


        if t_next < 0:
            print("ERROR t_next < 0!!!")

        curr_t = t + t_next
            
        if new_state == old_state :
            self.generate_event()
            return self.run_event(curr_t)
        else:
            self.nucleosomes[curr_nuc].change_state(new_state)
            new_assigned = self.nucleosomes[curr_nuc].state
            self.colors[curr_nuc] = Constants.state_to_color(new_assigned)
            self.update(old_state, new_assigned)
            return curr_t
            
    def handle_timers(self, index, old, new, timers):
        t_next = int(np.random.exponential(Constants.get_rate(old, new)))

        if t_next > 0:
            timers[index] = {
                'timer' : t_next,
                'old' : old,
                'new' : new
                    }
        else:
            self.update(old, new, index)

        
    ## END EVENTS ##
    ## Timestep simulation
    def timesim(self, n_nucs):
        TOT_EVENTS = 10000
        EVENTS_PER_TIMESTEP = 10
        TOT_TIMESTEPS = EVENTS_PER_TIMESTEP * TOT_EVENTS
        prob_event = EVENTS_PER_TIMESTEP / n_nucs
    
        timers = {}
        nuc_index_seq = [ x for x in range(n_nucs) ]
        curr_limit = n_nucs
        #seq_nuc_indicies = [ x for x in range(n_nucs) ]
        for t in range(TOT_TIMESTEPS):
            # handle timers first
            delete_timers = []
            for nuc_index in timers:
                if timers[nuc_index]['timer'] == 0:
                    # add back to pool & update
                    nuc_index_seq.append(nuc_index)
                    curr_limit += 1

                    old = timers[nuc_index]['old']
                    new = timers[nuc_index]['new']
                    
                    self.update(old, new, nuc_index)
                    delete_timers.append(nuc_index)
                else:
                    timers[nuc_index]['timer'] -= 1

            for timer in delete_timers:
                del timers[timer]

            num_events = int(np.random.poisson(EVENTS_PER_TIMESTEP))
            nucs_w_event = rand.choice(seq_nuc_indicies, num_events)

            a = self.dat['f']/(self.dat['f'] + 1)
            num_rand_events = np.random.poisson(EVENTS_PER_TIMESTEP * a)
            nucs_w_rand_event = rand.choice(nucs_w_event, num_rand_events)
            nucs_w_feedback_event = list( set(nucs_w_event) - set(nucs_w_rand_event) )

            for nuc in nucs_w_rand_event:
                old = self.nucleosomes[nuc].state

                if old == States.U_STATE:
                    if rand.random() < 0.5:
                        self.nucleosomes[nuc].state = States.A_STATE
                        self.update(old, States.A_STATE, nuc)
                    else:
                        self.nucleosomes[nuc].state = States.M_STATE
                        self.update(old, States.M_STATE, nuc)
                else:
                    self.nucleosomes[nuc].state = States.U_STATE
                    self.update(old, States.U_STATE, nuc)
            
            ## VERY IMPORTANT that M_mat and A_mat are np ARRAYS, not matricies! WE DO NOT WANT DOT PRODUCT
            prob_conv_M = self.prob_mat[nucs_w_feedback_event,:] * self.M_mat
            
            prob_conv_A = self.prob_mat[nucs_w_feedback_event,:] * self.A_mat

            tot_prob_per_nuc_M = np.sum(prob_conv_M, axis = 1)
            tot_prob_per_nuc_A = np.sum(prob_conv_A, axis = 1)

            for nuc in range(len(nucs_w_feedback_event)):
                curr_state = self.nucleosomes[nucs_w_feedback_event[nuc]].state
                if curr_state == States.M_STATE:
                    if rand.random() < tot_prob_per_nuc_A[nuc] :
                        self.nucleosomes[nucs_w_feedback_event[nuc]].state = States.U_STATE
                        update(curr_state, States.U_STATE)
                elif curr_state == States.A_STATE:
                    if rand.random() < tot_prob_per_nuc_M[nuc] :
                        self.nucleosomes[nucs_w_feedback_event[nuc]].state = States.U_STATE
                        update(curr_state, States.U_STATE)
                else: # U state
                    if tot_prob_per_nuc_A[nuc] != 0 and \
                            tot_prob_per_nuc_M[nuc] != 0:

                        A_prob = tot_prob_per_nuc_A[nuc]
                        M_prob = tot_prob_per_nuc_M[nuc]
                        
                        added_prob = A_prob + M_prob

                        if added_prob > 1:
                            # normalize to 1
                            scaling = 1 / added_prob 
                            A_prob *= scaling
                            M_prob *= scaling

                        cumsum = np.cumsum([1 - A_prob - M_prob, A_prob, M_prob])
                        int_sums = cumsum < rand.random()
                        index = np.sum(int_sums.astype(int))

                        if index == 1:
                            self.nucleosomes[nucs_w_feedback_event[nuc]].state = States.A_STATE
                        elif index == 2:
                            self.nucleosomes[nucs_w_feedback_event[nuc]].state = States.M_STATE
                            update(curr_state, States.M_STATE)
                        # else nothing


    ##
    def keep_index_sort(self, ar):
        key_ar = [(ar[i],i) for i in range(len(ar))]

        sorted_ar = sorted(key_ar, key = operator.itemgetter(0))
        index_ar = []
        val_ar = []
        for val, index in sorted_ar:
            val_ar.append(val)
            index_ar.append(index)

        return (val_ar, index_ar)


    def divide(self):
        for i in range(self.dat['n']):
            if random.random() <= 0.5:
                prev = self.nucleosomes[i].state
                self.nucleosomes[i].state = States.U_STATE
                self.update(prev, States.U_STATE)
                self.colors[i] = Constants.state_to_color(States.U_STATE)

    def update(self, old, new, i):
        if old == new:
            return

        self.totals[old] -= 1
        self.totals[new] += 1
        self.colors[i] = Constants.state_to_color(new)
    
        self.nucleosomes[i].state = new

        if old == States.M_STATE:
            self.M_mat[i] = 0
        elif old == States.A_STATE:
            self.A_mat[i] = 0

        if new == States.M_STATE:
            self.M_mat[i] = 1
        elif new == States.A_STATE:
            self.A_mat[i] = 1

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

