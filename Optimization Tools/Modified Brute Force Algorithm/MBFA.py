# Jackson Morgan
# MBFA for Knapsack Problem

import random

class Knapsack:
    def __init__(self, values, weights, max_weight):
        self.values = values
        self.weights = weights
        self.max_weight = max_weight

class candidate:
    def __init__(self, knapsack):

        # initialize counter variables and indexing lists:
        tot_weight = 0
        tot_value = 0
        used_items = []
        w = []
        v = []
        
        # build candidate solution:
        while tot_weight < knapsack.max_weight:
            i = random.randint(0, len(knapsack.weights)-1) # random number to index lists
            if i not in used_items: # check that item has not already been packed
                used_items.append(i) # add index to list of used items
                tot_weight = tot_weight + knapsack.weights[i] # add weight to total
                tot_value = tot_value + knapsack.values[i] # add value to total
                w.append(knapsack.weights[i]) # add weight of this item to list of items in this combo
                v.append(knapsack.values[i]) # add value of this item to list of items in this combo
        
        # remove last item from totals so that they accurately portray a solution <= max weight
        self.tot_w = tot_weight - knapsack.weights[i]
        self.tot_v = tot_value - knapsack.values[i]
        self.w = w[:-1]
        self.v = v[:-1]