from frozendict import frozendict

import numpy as np
from numpy import linalg as LA


class Strategy:
    LEHMANN = 4


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

    def _convert(self):
        mat = np.zeros((self.N, self.N))
        for n in range(self.N):
            for m in range(self.N):
                mat[n, m] = self.W[n, m].score
        return mat

    def max_eval(self):
        # computes the largest eigenvalue
        mat = self._convert()
        if len(mat) == 0:
            return 0.0
        vals = []
        for val in LA.eigvals(mat):
            vals.append(np.abs(val))
        return np.max(vals)

    def lift(self):
        """creates the weight matrix from the automaton"""
        W = self.R.zeros(self.N, self.N)
        for p in self.fsa.Q:
            for a, q, w in self.fsa.arcs(p):
                W[self.I[p], self.I[q]] += w
        return W

    def pathsum(self, strategy):

        if strategy == Strategy.LEHMANN:
            return self.lehmann_pathsum()

        else:
            raise NotImplementedError


    def _lehmann(self, zero=True):
        """
        Lehmann's (1977) algorithm.
        """

        # initialization
        V = self.W.copy()
        U = self.W.copy()

        # basic iteration
        for j in range(self.N):
            V, U = U, V
            V = self.R.zeros(self.N, self.N)
            for i in range(self.N):
                for k in range(self.N):
                    # i ➙ j ⇝ j ➙ k
                    V[i, k] = U[i, k] + U[i, j] * U[j, j].star() * U[j, k]

        # post-processing (paths of length zero)
        if zero:
            for i in range(self.N):
                V[i, i] += self.R.one

        return V

    def lehmann(self, zero=True):
        # TODO: check we if we can't do away with this method.

        V = self._lehmann(zero=zero)

        W = {}
        for p in self.fsa.Q:
            for q in self.fsa.Q:
                if p in self.I and q in self.I:
                    W[p, q] = V[self.I[p], self.I[q]]
                elif p == q and zero:
                    W[p, q] = self.R.one
                else:
                    W[p, q] = self.R.zero

        return frozendict(W)



