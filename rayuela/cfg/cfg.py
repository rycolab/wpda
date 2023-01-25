import copy
from frozendict import frozendict

from rayuela.base.semiring import Semiring, Boolean
from rayuela.base.symbol import Sym, ε

from rayuela.fsa.state import State
from rayuela.fsa.fsa import FSA

from rayuela.cfg.nonterminal import NT, S
from rayuela.cfg.production import Production
from rayuela.cfg.exceptions import InvalidProduction


class CFG:

    def __init__(self, R=Boolean):

        # A weighted context-free grammar is a 5-tuple <R, Σ, V, P, S> where
        # • R is a semiring;
        # • Σ is an alphabet of terminal symbols;
        # • V is an alphabet of non-terminal symbols;
        # • P is a finite relation V × (Σ ∪ V)* × R;
        # • S ∈ V is a distinguished started symbol.

        # semiring
        self.R = R

        # alphabet
        self.Sigma = set([])

        # non-terminals
        self.V = set([S])

        # productions
        self._P = self.R.chart()

        # unique start symbol
        self.S = S

        # unary FSA
        self.unary_fsa = None

    @property
    def terminal(self):
        for p, w in self.P:
            (head, body) = p
            if len(body) == 1 and isinstance(body[0], Sym):
                yield p, w

    @property
    def unary(self):
        for p, w in self.P:
            (head, body) = p
            if len(body) == 1 and isinstance(body[0], NT):
                yield p, w

    @property
    def binary(self):
        for p, w in self.P:
            (head, body) = p
            if len(body) == 2 and isinstance(body[0], NT) \
                    and isinstance(body[1], NT):
                yield p, w

    @property
    def size(self):
        size = 0
        for (_, body), _ in self.P:
            size += len(body) + 1
        return size

    @property
    def num_rules(self):
        return len(self._P)

    def w(self, p):
        return self._P[p]

    def spawn(self):
        return CFG(R=self.R)

    def make_unary_fsa(self):
        one = self.R.one
        fsa = FSA(R=self.R)

        # add a state for every non-terminal
        for X in self.V:
            fsa.add_state(State(X))

        # add arcs between every pair of unary rules
        for (head, body), w in self.unary:
            fsa.add_arc(State(body[0]), ε, State(head), w)

        # add initial and final weight one for every state
        for q in list(fsa.Q):
            fsa.set_I(q, one)
            fsa.set_F(q, one)

        self.unary_fsa = fsa

    def eps_partition(self):
        """ makes a new grammar can only generate epsilons """
        ecfg = self.spawn()
        ncfg = self.spawn()

        def has_terminal(body):
            for elem in body:
                if isinstance(elem, Sym) and elem != ε:
                    return True
            return False

        for p, w in self.P:
            head, body = p
            if has_terminal(body):
                ncfg.add(w, head, *body)
            elif len(body) == 1 and body[0] == ε:
                ecfg.add(w, head, *body)
            else:
                ncfg.add(w, head, *body)
                ecfg.add(w, head, *body)

        return ecfg, ncfg

    @property
    def P(self):
        for p, w in self._P.items():
            yield p, w

    def P_byhead(self, X, unary=True):
        for p, w in self._P.items():
            if X == p.head:
                if not unary and len(p.body) == 1 and isinstance(p.body[0], NT):
                    continue
                yield p, w

    def add(self, w, head, *body):
        if not isinstance(head, NT):
            raise InvalidProduction
        if not isinstance(w, Semiring):
            w = self.R(w)

        self.V.add(head)

        for elem in body:
            if isinstance(elem, NT):
                self.V.add(elem)
            elif isinstance(elem, Sym) and elem != ε:
                self.Sigma.add(elem)
            elif elem != ε:
                raise InvalidProduction

        self._P[Production(head, body)] += w

    def get_productions(self):
        return self._P

    def freeze(self):
        self.Sigma = frozenset(self.Sigma)
        self.V = frozenset(self.V)
        self._P = frozendict(self._P)

    def copy(self):
        return copy.deepcopy(self)


    def __str__(self):
        return "\n".join(f"{p}\t{w}" for (p, w) in
                         sorted(self.P, key=lambda x: (len(str(x[0].head)), str(x[0].head), len(str(x[0])))))
    # return "\n".join(f"{p}" for (p, w) in sorted(self.P, key=lambda x: (len(str(x[0].head)), str(x[0].head), len(str(x[0])))))

