from collections import defaultdict as dd
from itertools import product
from frozendict import frozendict
from functools import reduce

from rayuela.base.semiring import Boolean
from rayuela.base.symbol import Sym, ε

from rayuela.fsa.state import State

from rayuela.cfg.cfg import CFG
from rayuela.cfg.nonterminal import NT, S, bottom
from rayuela.cfg.production import Production

class PDA:

    def __init__(self, R=Boolean):

        # DEFINITION
        # A weighted pushdown automaton is a 5-tuple <R, Σ, Q, δ, λ, ρ> where
        # • R is a semiring;
        # • Σ is an alphabet of symbols;
        # • Q is a finite set of states;
        # • Γ is an alphabet of stack symbols;
        # • δ is a finite relation Q × Σ × Q × Γ* × Γ* → R;
        # • λ is an initial weight function;
        # • ρ is a final weight function.

        # NOTATION CONVENTIONS
        # • single states (elements of Q) are denoted q
        # • multiple states not in sequence are denoted, p, q, r, ...
        # • multiple states in sequence are denoted i, j, k, ...
        # • symbols (elements of Σ) are denoted lowercase a, b, c, ...
        # • single weights (elements of R) are denoted w
        # • multiple weights (elements of R) are denoted u, v, w, ...

        # semiring
        self.R = R

        # finite set of states
        self.Q = set([])

        # alphabet of symbols
        self.Sigma = set([])

        # alphabet of stack symbols
        self.Gamma = set([])

        # unique ending stack symbol
        self.S = S

        # transition function : Q × Σ × Q × Γ × Γ* → R
        self.δ = dd(lambda : dd(lambda : dd(lambda : dd(lambda : dd(lambda : self.R.zero)))))

        # initial weight function
        self.λ = R.chart()

        # final weight function
        self.ρ = R.chart()

    def set_I(self, q, w):
        self.add_state(q)
        self.λ[q] = w

    def set_F(self, q, w):
        self.add_state(q)
        self.ρ[q] = w

    def add_I(self, q, w):
        self.add_state(q)
        self.λ[q] += w

    def add_F(self, q, w):
        self.add_state(q)
        self.ρ[q] += w

    def freeze(self):
        self.Sigma = frozenset(self.Sigma)
        self.Gamma = set(self.Gamma)
        self.Q = frozenset(self.Q)
        self.δ = frozendict(self.δ)
        self.λ = frozendict(self.λ)
        self.ρ = frozendict(self.ρ)

    @property
    def I(self):
        for q, w in self.λ.items():
            if w != self.R.zero:
                yield q, w

    @property
    def F(self):
        for q, w in self.ρ.items():
            if w != self.R.zero:
                yield q, w

    @property
    def arcs(self):
        for i in self.δ:
            for a in self.δ[i]:
                for j in self.δ[i][a]:
                    for head in self.δ[i][a][j]:
                        for body, w in self.δ[i][a][j][head].items():
                            if w != self.R.zero:
                                yield (i, a, j, head, body), w

    def in_cnf(self, strategy="bottom-up"):
        """Returns true if the PDA is in CNF normal form (as defined in EMNLP 2022
        submission `Algorithms for Weighted Pushdown Automata`)"""

        if strategy == "bottom-up":
            i = [p for p, _ in self.I]
            f = [p for p, _ in self.F]

            assert len(i) == 1, "A bottom-up WPDA can have a single initial state"
            assert len(f) == 1, "A bottom-up WPDA can have a single final state"

            for (_, a, _, head, body), _ in self.arcs:
                # push exactly 1
                if len(head) != 1:
                    return False
                if a == ε:
                    # pop exactly 2 for non-scanning transitions
                    if body is None or len(body) != 2:
                        return False
                else:
                    # pop at most 2 for scanning transitions
                    if body is not None and len(body) > 2:
                        return False
            return True
        elif strategy == "top-down":
            i = [p for p, _ in self.I]
            f = [p for p, _ in self.F]

            assert len(i) == 1, "A top-down WPDA can have a single initial state"
            assert len(f) == 1, "A top-down WPDA can have a single final state"

            for (_, a, _, head, body), _ in self.arcs:
                # pop exactly 1
                if body is None or len(body) != 1:
                    return False
                if a == ε:
                    # push exactly 2 for non-scanning transitions
                    if head is None or len(head) != 2:
                        return False
                else:
                    # push at most 2 for scanning transitions
                    if head is not None and len(head) > 2:
                        return False
            return True
        else:
            raise NotImplementedError

    def add_state(self, q):
        self.Q.add(q)

    def add_states(self, Q):
        for q in Q:
            self.add_state(q)

    def add(self, w, i, a, j, head, *body):
        assert isinstance(head, tuple)
        assert isinstance(body, tuple)
        self.add_states([i, j])
        self.Sigma.add(a)

        for X in head: self.Gamma.add(X)
        for X in body: self.Gamma.add(X)

        if len(body) > 0:
            for elem in body:
                self.Gamma.add(elem)
            self.δ[i][a][j][head][tuple(body)] += w
        else:
            self.δ[i][a][j][head][()] += w

    def spawn(self):
        return PDA(self.R)

    def remove_nullary(self, strategy="bottom-up"):
        from rayuela.pda.transformer import Transformer
        return Transformer.remove_nullary(self, strategy=strategy)

    def remove_unary(self, strategy="bottom-up"):
        from rayuela.pda.transformer import Transformer
        return Transformer.remove_unary(self, strategy=strategy)

    def binarize(self, strategy="bottom-up"):
        from rayuela.pda.transformer import Transformer
        return Transformer.binarize(self, strategy=strategy)

    def new_state(self):
        counter = len(self.Q)

        q = State(f"@{counter}")
        while q in self.Q:
            counter += 1
            q = State(f"@{counter}")

        return q

    def to_cnf(self, strategy="bottom-up"):
        """
        Converts a WPDA into a normal form
        analogous to CNF
        """
        npda = self.binarize(strategy=strategy)
        npda = npda.remove_nullary(strategy=strategy)
        npda = npda.remove_unary(strategy=strategy)
        return npda


    def __str__(self):
        return "\n".join([(f"{i} -- {a} --> {j}\t["+", ".join(map(str, body))+"]\t["+", ".join(map(str, head))+f"]\t{w}") for (i, a, j, head, body), w in self.arcs])

    def get_len(self):
        pda_str = "\n".join([(f"{i} -- {a} --> {j}\t["+", ".join(map(str, body))+"]\t["+", ".join(map(str, head))+f"]\t{w}") for (i, a, j, head, body), w in self.arcs])
        return len(pda_str.split("\n"))