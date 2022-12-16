import random
import bisect
import copy
import sys
from typing import Union

import numpy as np

LOW_PROB_THRESH = 1e-12


class DiscreteDistribution:
    def __init__(self, weights, values, degenerate_val=None, method='bisect'):

        # some sanity checking
        if not len(weights) or not len(values):
            print('\nWarning: weight or value vector given to DiscreteDistribution() are 0-length.\n')
            values = [0]
            weights = [0]

        self.method = method
        sum_weight = float(sum(weights))

        # if probability of all input events is 0, consider it degenerate and always return the first value
        if sum_weight < LOW_PROB_THRESH:
            self.degenerate = values[0]
        else:
            self.weights = [n / sum_weight for n in weights]
            self.values = values
            if len(self.values) != len(self.weights):
                print('\nError: length and weights and values vectors must be the same.\n')
                exit(1)
            self.degenerate = degenerate_val

            if self.method == 'alias':
                len_weights = len(self.weights)
                prob_vector = np.zeros(len_weights)
                count_vector = np.zeros(len_weights, dtype=np.int)
                smaller = []
                larger = []
                for kk, prob in enumerate(self.weights):
                    prob_vector[kk] = len_weights * prob
                    if prob_vector[kk] < 1.0:
                        smaller.append(kk)
                    else:
                        larger.append(kk)
                while len(smaller) > 0 and len(larger) > 0:
                    small = smaller.pop()
                    large = larger.pop()
                    count_vector[small] = large
                    prob_vector[large] = (prob_vector[large] + prob_vector[small]) - 1.0
                    if prob_vector[large] < 1.0:
                        smaller.append(large)
                    else:
                        larger.append(large)

                self.a1 = len(count_vector) - 1
                self.a2 = count_vector.tolist()
                self.a3 = prob_vector.tolist()

            elif self.method == 'bisect':
                self.cum_prob = np.cumsum(self.weights).tolist()[:-1]
                self.cum_prob.insert(0, 0.)

            else:
                print("\nUnknown discreet distribution method.\n")

    def __str__(self):
        return str(self.weights) + ' ' + str(self.values) + ' ' + self.method

    def sample(self) -> Union[int, float]:
        """
        This is one of the slowest parts of the code. Or it just gets hit the most times. Will need
        to investigate at some point.
        :return: Since this function is selecting an item from a list, and the list could theoretically be anything,
        then in a broad sense this function returns a list item or a generic object. But I'm fairly confident that most
        of these uses will be lists of ints or floats, but will investigate further
        """

        if self.degenerate is not None:
            return self.degenerate

        else:

            if self.method == 'alias':
                random1 = random.randint(0, self.a1)
                random2 = random.random()
                if random2 < self.a3[random1]:
                    return self.values[random1]
                else:
                    return self.values[self.a2[random1]]

            elif self.method == 'bisect':
                r = random.random()
                return self.values[bisect.bisect(self.cum_prob, r) - 1]


# takes k_range, lambda, [0,1,2,..], returns a DiscreteDistribution object
# with the corresponding to a poisson distribution

def poisson_list(k_range, input_lambda):
    min_weight = 1e-12
    if input_lambda < min_weight:
        return DiscreteDistribution([1], [0], degenerate_val=0)
    log_factorial_list = [0.0]
    for k in k_range[1:]:
        log_factorial_list.append(np.log(float(k)) + log_factorial_list[k - 1])
    w_range = [np.exp(k * np.log(input_lambda) - input_lambda - log_factorial_list[k]) for k in k_range]
    w_range = [n for n in w_range if n >= min_weight]
    if len(w_range) <= 1:
        return DiscreteDistribution([1], [0], degenerate_val=0)
    return DiscreteDistribution(w_range, k_range[:len(w_range)])
