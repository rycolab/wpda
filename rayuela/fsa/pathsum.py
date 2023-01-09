from typing import Dict, Type
from collections import defaultdict
from functools import reduce
from frozendict import frozendict

import numpy as np
from numpy import linalg as LA

class Strategy:
    VITERBI = 1
    BELLMANFORD = 2
    DIJKSTRA = 3
    LEHMANN = 4
    JOHNSON = 5
    FIXPOINT = 6
    DECOMPOSED_LEHMANN = 7
    VITERBI_FAILURE_ARCS = 8
    VITERBI_FAILURE_ARCS_CRF = 9

class Pathsum:
    def __init__(self, fsa):

        # basic FSA stuff
        self.fsa = fsa
        self.R = fsa.R
        self.N = self.fsa.num_states

        # state dictionary
        self.I = {}
        for n, q in enumerate(self.fsa.Q):
            self.I[q] = n

        # lift into the semiring
        self.W = self.lift()