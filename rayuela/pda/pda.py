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


class Run:

    def __init__(self, pda, length):
        self.pda = pda
        self.index = 0
        self.σ = tuple()
        self.inputs = list(range(1,length+1))
        self.actions = None

    def moves(self):
        raise NotImplementedError

    def step(self, i, a, j, gg, g):
        assert self.pda.δ[i][a][j][gg][ggʼ]
        self.σ = self.σ[:len(self.σ)-len(gg)]
        self.σ += tuple([g])
        self.index += 1

    def get_options(self, buffer, state, stack):

        if not buffer:
            b = None
        else:
            b = buffer[0]

        options = []
        for (i, a, j, head, body), w in self.pda.arcs:
            if str(state) == str(i) \
                and (str(b) == str(a) or a == ε) \
                and (str(list(reversed(stack[:len(body)]))) == str(list(body)) or len(body) == 0):
                options.append((i, a, j, head, body))

        return options

    def update_stack(self, stack, head, body):
        if not head:
            # pop from stack
            n_stack = stack[1:]
        else:
            if not body:
                # push on stack
                n_stack = stack[0:]
            else:
                # pop from stack
                n_stack = stack[len(body):]
            for b in reversed(head):
                n_stack.insert(0, b)

        return n_stack

    def find_path(self, state=0, stack=[], buffer=[], actions=[]):
        if stack == [S] and not buffer:
            self.actions = actions
        else:
            options = self.get_options(buffer, state, stack)
            if options:
                for (i, a, j, head, body) in options:
                    if a != ε:
                        n_buffer = buffer[1:]
                    else:
                        n_buffer = buffer[0:]

                    n_stack = self.update_stack(stack, head, body)
                    n_actions = actions.copy()
                    self.find_path(j, n_stack, n_buffer, n_actions+[(i, a, j, head, body)])

    def __merge(self, head, body):
        if body:
            return self.__dep(self.pda.labels[Production(head[0], body)])
        else:
            return ""

    def __dep(self, label):
        if label is None:
            return ""
        else:
            (n1, n2, direction) = label
            if direction == "r":
                return f"{n2}→{n1}"
            elif direction == "l":
                return f"{n1}←{n2}"
            else:
                raise AssertionError("Invalid Direction")

    def oracle(self):
        return [" ".join(map(str, head)) + " → " + " ".join(map(str, body)) if body else " ".join(map(str, head)) + " → " + str(a) for (_, a, _, head, body) in self.actions]

    def __str__(self):
        # return "\n".join([(f"{i} -- {a} --> {j}\t[" + ", ".join(map(str, body)) + "]\t[" + ", ".join(
        # 	map(str, head)) + f"]") for (i, a, j, head, body) in self.actions])

        # return "\n".join([" ".join(map(str, head)) + " → " + " ".join(map(str, body)) if body else " ".join(map(str, head)) + " → " + str(a) for (_, a, _, head, body) in self.actions])

        return "\n".join(
            [(f"{i} -- {a} --> {j}\t[" + ", ".join(map(str, body)) + "]\t[" + ", ".join(
                map(str, head)) + f"]\t{self.__merge(head, body)}")
             for (i, a, j, head, body) in self.actions])

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

    def trim(self):
        from rayuela.pda.transformer import Transformer
        return Transformer.trim(self)

    def push(self):
        from rayuela.pda.transformer import Transformer
        return Transformer.push(self)

    def pathsum(self):
        from rayuela.pda.pathsum import Pathsum
        return Pathsum(self).sum()

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

    def to_push(self):
        """
        Convert th PDA to a normal form that pushes *at most* 1
        stack symbol at every time step and which terminates
        when we deposit an S on the stack.
        """

        zero, one = self.R.zero, self.R.one
        npda = self.spawn()
        # TODO: make this a fresh symbol
        Sʼ = NT("Sʼ")

        # head = stack elements that are pushed
        # body = stack elements that are popped
        for (i, a, k, head, body), w in self.arcs:
            # remove strip S
            #head = tuple([Sʼ if x == S else x for x in head])
            #body = tuple([Sʼ if x == S else x for x in body])

            if len(head) <= 1:
                npda.add(w, i, a, k, head, *body)
            else:
                npda.add(w, i, a, self.new_state(), (head[0],), *body)
                rest = list(head[1:])
                while rest:
                    X = rest.pop()
                    if len(rest) == 0:
                        npda.add(one, self.new_state(), a, k, (X,))
                    else:
                        npda.add(one, self.new_state(), a, self.new_state())

        for p, w in self.I:
            npda.set_I(p, w)
        for p, w in self.F:
            npda.set_F(p, w)

        # Add new S
        #for p, w in self.F:
        #	npda.add(one, p, ε, p, (S,))


        return npda

    def to_pop(self):
        """
        Convert the PDA to a normal form that pops *at most* 1
        stack symbol at every time step.
        """

        # TODO: add pop normal form

        zero, one = self.R.zero, self.R.one
        npda = self.spawn()

        # head = stack elements that are pushed
        # body = stack elements that are popped
        for (i, a, k, head, body), w in self.arcs:
            if len(body) <= 1:
                npda.add(w, i, a, k, head, *body)
            else:
                npda.add(w, i, a, self.new_state(), head, body[0])
                rest = list(body[1:])
                while rest:
                    X = rest.pop()
                    if len(rest) == 0:
                        npda.add(one, self.new_state(), a, k, (), (X,))
                    else:
                        npda.add(one, self.new_state(), a, self.new_state(), (), (X,))

        for p, w in self.I:
            npda.set_I(p, w)
        for p, w in self.F:
            npda.set_F(p, w)

        return npda

    def to_cfg(self):
        return self.to_push()._to_cfg_push()

    def _to_cfg_push(self):
        # TODO: The assumptions are unchecked:
        # • We assume that we cannot both scan and pop (from the stack)
        # • We assume that we cannot push more than on element
        # Add a check to ensure this is enforced


        # Timmy says:
        # Compute transitive closure of the PDA *ignoring* the stack
        # this results in a FSA which "overestimates" pairwise reachability.
        # It's a tighter estimate than the cross product, used below.

        ncfg = CFG(R=self.R)
        recovery = {}

        for j, k, a in product(self.Q, self.Q, self.Sigma.union([ε])):
            for ohead in self.δ[j][a][k]:
                for obody, w in self.δ[j][a][k][ohead].items():
                    head = tuple([bottom if x == S else x for x in ohead])
                    body = tuple([bottom if x == S else x for x in obody])

                    # pushes zero
                    if len(head) == 0:

                        # pops zero
                        if len(body) == 0:
                            for i, X in product(self.Q, self.Gamma):
                                X = bottom if X == S else X
                                nbody = [ NT((i, X, j)) ]
                                if a != ε: nbody = nbody + [a]
                                ncfg.add(w, NT((i, X, k)), *nbody)

                                # recovery
                                p = Production(NT((i, X, k)), tuple(nbody))
                                recovery[p] = (j, a, k, ohead, obody)

                        # pops at least one
                        else:
                            for span in product(self.Q, repeat=len(body)):
                                nbody = [ NT((n, Y, m)) for Y, (n, m) in zip(body, zip(span, span[1:]+(j,))) ]
                                if a != ε: nbody += [a]
                                for i, X in product(self.Q, self.Gamma):
                                    X = bottom if X == S else X
                                    ncfg.add(w, NT((i, X, k)), NT((i, X, span[0])), *nbody)

                                    # recovery
                                    p = Production(NT((i, X, k)), tuple([NT((i, X, span[0]))]+nbody))
                                    recovery[p] = (j, a, k, ohead, obody)

                    # pushes one
                    elif len(head) == 1:

                        # pops zeros
                        if len(body) == 0:
                            ncfg.add(w, NT((j, head[0], k)), a)

                            # recovery
                            p = Production(NT((j, head[0], k)), (a,))
                            recovery[p] = (j, a, k, ohead, obody)

                        # pops at least one
                        else:
                            for span in product(self.Q, repeat=len(body)):
                                nbody = [ NT((n, Y, m)) for Y, (n, m) in zip(body, zip(span, span[1:]+(j,))) ]
                                if a != ε: nbody += [a]
                                ncfg.add(w, NT((span[0], head[0], k)), *nbody)

                                # recovery
                                p = Production(NT((span[0], head[0], k)), tuple(nbody))
                                recovery[p] = (j, a, k, ohead, obody)

        # wrap-up (create start state)
        # TODO: make this work for non-commutative semirings
        heads = set([])
        for p, q in product(self.Q, self.Q):
            heads.add(NT((p, bottom, q)))

        for p, w in self.I:
            for q, v in self.F:
                ncfg.add(w*v, S, NT((p, bottom, q)))

        # TODO: check if these are necessary
        # ncfg = ncfg.trim()
        # ncfg.make_unary_fsa()
        return ncfg, recovery

    def _to_cfg_pop(self):

        ncfg = CFG(R=self.R)

        for j, k, a in product(self.Q, self.Q, self.Sigma.union([ε])):
            for head in self.δ[j][a][k]:
                for body, w in self.δ[j][a][k][head].items():

                    # pops zero
                    if len(body) == 0:

                        # push zero
                        if len(head) == 0:
                            for i, X in product(self.Q, self.Gamma):
                                nhead = [ NT((i, X, j)) ]
                                if a != ε: nhead = [a] + nhead
                                ncfg.add(w, NT((i, X, k)), *nhead)

                        # push at least one
                        else:
                            for span in product(self.Q, repeat=len(head)):
                                nhead = [ NT((n, Y, m)) for Y, (n, m) in zip(head, zip((k,)+span, span)) ]
                                if a != ε: nhead = [a] + nhead
                                for i, X in product(self.Q, self.Gamma):
                                    ncfg.add(w, NT((i, X, span[-1])), NT((i, X, j)), *nhead)

                    # pops one
                    elif len(body) == 1:

                        # pops zeros
                        if len(head) == 0:
                            ncfg.add(w, NT((j, body[0], k)), a)

                        # pushes at least one
                        else:
                            for span in product(self.Q, repeat=len(head)):
                                nhead = [ NT((n, Y, m)) for Y, (n, m) in zip(head, zip((k,)+span, span)) ]
                                if a != ε: nhead = [a] + nhead
                                ncfg.add(w, NT((j, body[0], span[-1])), *nhead)

        # wrap-up (create start state)
        # TODO: make this work for non-commutative semirings
        for p, w in self.I:
            for q, v in self.F:
                ncfg.add(w*v, S, NT((p, S, p)))


        ncfg = ncfg.trim()
        ncfg.make_unary_fsa()
        return ncfg

    def __str__(self):
        return "\n".join([(f"{i} -- {a} --> {j}\t["+", ".join(map(str, body))+"]\t["+", ".join(map(str, head))+f"]\t{w}") for (i, a, j, head, body), w in self.arcs])

    def get_len(self):
        pda_str = "\n".join([(f"{i} -- {a} --> {j}\t["+", ".join(map(str, body))+"]\t["+", ".join(map(str, head))+f"]\t{w}") for (i, a, j, head, body), w in self.arcs])
        return len(pda_str.split("\n"))