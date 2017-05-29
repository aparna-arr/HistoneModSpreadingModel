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


    ## helper functions ##

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

    def handle_timers(self, index, old, new, timers, nuc_seq, map_seq, lim):
        t_next = int(np.random.exponential(Constants.get_rate(old, new)))

        if t_next > 0:
            timers[index] = {
                'timer' : t_next,
                'old' : old,
                'new' : new
                    }
            return self.fake_del(map_seq, nuc_seq, lim, index)

        else:
            self.update(old, new, index)
            return lim

    def fake_del(self, map_seq, nuc_seq, lim, del_elem):
        if lim >= 0:
            lim -= 1

            tmp = nuc_seq[lim]

            del_elem_index = map_seq[del_elem]
            nuc_seq[lim] = nuc_seq[del_elem_index]
            map_seq[del_elem] = lim

            nuc_seq[del_elem_index] = tmp
            map_seq[tmp] = del_elem_index

            return lim
        else:
            print("ERROR trying to delete but everything is deleted!")

    def fake_readd(self, map_seq, nuc_seq, lim, add_elem):
        if lim < len(nuc_seq):
            add_elem_index = map_seq[add_elem]
            tmp = nuc_seq[lim]

            nuc_seq[lim] = nuc_seq[add_elem_index]
            map_seq[add_elem] = lim

            nuc_seq[add_elem_index] = tmp
            map_seq[tmp] = add_elem_index

            return lim + 1
        else:
            print("ERROR trying to re-add but nothing deleted!")

    def print_nucs(self, fp):
        for n in self.nucleosomes:
            fp.write(States.enum_to_string(n.state))

        fp.write("\n")
        
    ## Timestep simulation
    def timesim(self, n_nucs):
        EVENTS_PER_TIMESTEP = 10
        TOT_TIMESTEPS = self.dat['t'] 
        prob_event = EVENTS_PER_TIMESTEP / n_nucs
    
        timers = {}
        map_to_seq = [ x for x in range(n_nucs) ]
        nuc_index_seq = [ x for x in range(n_nucs) ]
        lim = n_nucs
        fp = open(self.dat['o'] + ".txt", "w")

        for t in range(TOT_TIMESTEPS):
            print("timestep:",t)
            # handle timers first
            self.print_nucs(fp)

            delete_timers = []
            for nuc_index in timers:
                if timers[nuc_index]['timer'] == 0:
                    # add back to pool & update
                    lim = self.fake_readd(map_to_seq, nuc_index_seq, lim, nuc_index)

                    old = timers[nuc_index]['old']
                    new = timers[nuc_index]['new']
                    
                    self.update(old, new, nuc_index)
                    delete_timers.append(nuc_index)
                else:
                    timers[nuc_index]['timer'] -= 1

            for timer in delete_timers:
                del timers[timer]


            if self.dat['d'] != Divisions.NONE and t != 0:
                div_rate = self.dat['data']['divisions']

                if t % div_rate == 0:
                    #print("DIVIDING")
                    num_nucs_replaced = int(np.random.poisson(lim / 2))

                    if num_nucs_replaced > lim:
                        num_nucs_replaced = lim

                    nucs_replaced = random.sample(nuc_index_seq[:lim], num_nucs_replaced)
                        #lim = self.handle_timers(nuc, old, States.A_STATE, timers, nuc_index_seq, map_to_seq, lim)

                    for nuc in nucs_replaced:
                        old = self.nucleosomes[nuc].state
                        #lim = self.handle_timers(nuc, old, States.U_STATE, timers, nuc_index_seq, map_to_seq, lim)

                        self.update(old, new, nuc)
                    continue

            if t >= Constants.RECRUIT_INIT_TIME and \
                    t <= (Constants.RECRUIT_INIT_TIME + 
                            Constants.RECRUIT_TIMESTEPS):
                start_nuc = int(self.dat['n']/2 - Constants.RECRUIT_N_NUCS/2)
                end_nuc = start_nuc + Constants.RECRUIT_N_NUCS
                #print("start_nuc, end_nuc:", start_nuc, end_nuc)
                for i in range(start_nuc, end_nuc):
                    nuc_seq_i = map_to_seq[i]
                    if nuc_seq_i < lim:
                        lim = self.handle_timers(i, self.nucleosomes[i].state, Constants.RECRUIT_CONV_TO, timers, nuc_index_seq, map_to_seq, lim)
                
            num_events = int(np.random.poisson(EVENTS_PER_TIMESTEP * (lim / n_nucs)))

            if num_events >= lim:
                num_events = lim 

            nucs_w_event = random.sample(nuc_index_seq[:lim], num_events)

            a = 1/(self.dat['f'] + 1)
            
            num_rand_events = np.random.poisson(EVENTS_PER_TIMESTEP * (lim / n_nucs) * a)


            #print("num events total", len(nucs_w_event))
            #print("a", a, "adj poisson center:", EVENTS_PER_TIMESTEP * (lim / n_nucs) * a, "num rand events:", num_rand_events)
    
            if num_rand_events > len(nucs_w_event):
                num_rand_events = len(nucs_w_event)

            nucs_w_rand_event = random.sample(nucs_w_event, num_rand_events)
            nucs_w_feedback_event = list( set(nucs_w_event) - set(nucs_w_rand_event) )

            #print("num feedback:", len(nucs_w_feedback_event))

            for nuc in nucs_w_rand_event:
                old = self.nucleosomes[nuc].state

                if old == States.U_STATE:
                    if random.random() < 0.5:
                        lim = self.handle_timers(nuc, old, States.A_STATE, timers, nuc_index_seq, map_to_seq, lim)
                    else:
                        lim = self.handle_timers(nuc, old, States.M_STATE, timers, nuc_index_seq, map_to_seq, lim)
                else:
                    lim = self.handle_timers(nuc, old, States.U_STATE, timers, nuc_index_seq, map_to_seq, lim)
            
            ## VERY IMPORTANT that M_mat and A_mat are np ARRAYS, not matricies! WE DO NOT WANT DOT PRODUCT
            prob_conv_M = self.prob_mat[nucs_w_feedback_event,:] * self.M_mat
            
            prob_conv_A = self.prob_mat[nucs_w_feedback_event,:] * self.A_mat

            tot_prob_per_nuc_M = np.sum(prob_conv_M, axis = 1)
            tot_prob_per_nuc_A = np.sum(prob_conv_A, axis = 1)

            for nuc in range(len(nucs_w_feedback_event)):
                curr_state = self.nucleosomes[nucs_w_feedback_event[nuc]].state
                if curr_state == States.M_STATE:
                    if random.random() < tot_prob_per_nuc_A[nuc] :
                        lim = self.handle_timers(nucs_w_feedback_event[nuc], curr_state, States.U_STATE, timers, nuc_index_seq, map_to_seq, lim)
                elif curr_state == States.A_STATE:
                    if random.random() < tot_prob_per_nuc_M[nuc] :
                        lim = self.handle_timers(nucs_w_feedback_event[nuc], curr_state, States.U_STATE, timers, nuc_index_seq, map_to_seq, lim)
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
                        int_sums = cumsum < random.random()
                        index = np.sum(int_sums.astype(int))

                        if index == 1:
                            lim = self.handle_timers(nucs_w_feedback_event[nuc], curr_state, States.A_STATE, timers, nuc_index_seq, map_to_seq, lim)
                        elif index == 2:
                            lim = self.handle_timers(nucs_w_feedback_event[nuc], curr_state, States.M_STATE, timers, nuc_index_seq, map_to_seq, lim)
                        # else nothing

        fp.close()
    ##
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

