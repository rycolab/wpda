from collections import defaultdict as dd

from rayuela.base.symbol import ε

from rayuela.cfg.nonterminal import Slash, NT, S

from rayuela.pda.allsum import Allsum

from rayuela.fsa.state import State
from rayuela.fsa.fsa import FSA
from rayuela.fsa.pathsum import Pathsum

from itertools import product

class Transformer:


    def __init__(self):
        pass


    def binarize(pda, strategy="bottom-up"):
        """Creates an equivalent WPDA whose transitions are k-pop, 1-push, k <= 2,
        for the bottom-up case and 1-pop, k-push, k <= 2, for the top-down case"""
        if strategy == "bottom-up":
            npda = pda.spawn()
            # add all states initially so that we can generate
            # new states when adding transitions
            npda.add_states(pda.Q)

            for (p, a, q, head, body), w in pda.arcs:
                if body == None or len(body) <= 2:
                    npda.add(w, p, a, q, head, *body)
                else:
                    r = npda.new_state()
                    npda.add(w, p, a, r, (body[1], ), *(body[0], body[1]))
                    for i in range(1, len(body) - 1):
                        t = npda.new_state()
                        npda.add(npda.R.one, r, ε, t, (body[i + 1],), *(body[i], body[i + 1]))
                        r = t
                    npda.add(npda.R.one, r, ε, q, (head[0],), *(body[-2], body[-1]))

            return npda

        elif strategy == "top-down":
            raise NotImplementedError
        else:
            raise NotImplementedError

    def _unary_eps_partition(pda, strategy="bottom-up"):
        """Constructs an FSA E with states of the form (p, X) ∈ Q, Γ
        and transitions (p, X) -- ε/w --> (q, Y) from a WPDA with unary
        transitions, that is, transitions of the form p, X -- ε/w --> q, Y
        """
        E = FSA(R=pda.R)

        if strategy == "bottom-up":
            for (p, a, q, head, body), w in pda.arcs:
                if a == ε and body is not None and len(body) == 1 and len(head) == 1:
                    E.add_arc(i=(p, body[0]), a=ε, j=(q, head[0]), w=w)
        elif strategy == "top-down":
            for (p, a, q, head, body), w in pda.arcs:
                #TODO: check this again
                if a == ε and body is not None and len(body) == 1 and len(head) == 1:
                    E.add_arc(i=(q, head[0]), a=ε, j=(p, body[0]), w=w)
        else:
            raise NotImplementedError

        return E

    def remove_unary(pda, strategy="bottom-up"):
        """
        Constructs an equivalent WPDA by removing unary transitions,
        that is, transitions of the form p, X -- ε/w --> q, Y
        """
        npda = pda.spawn()
        E = Transformer._unary_eps_partition(pda, strategy=strategy)
        U = Pathsum(E).lehmann(zero=False)

        for q, w in pda.I:
            npda.set_I(q, w)

        for q, w in pda.F:
            npda.set_I(q, w)

        if strategy == "bottom-up":
            for (p, a, q, head, body), w in pda.arcs:
                if a != ε:
                    X = head[0]
                    for r in pda.Q:
                        for Y in pda.Gamma:
                            if (State((q, X)), State((r, Y))) in U:
                                npda.add(w * U[State((q, X)), State((r, Y))], p, a, r, (Y,), *body)
                elif body is not None and len(body) != 1:
                    npda.add(w, p, a, q, head, *body)
        #TODO: check this again
        elif strategy == "top-down":
            for (p, a, q, head, body), w in pda.arcs:
                if a != ε:
                    X = body[0]
                    for r in pda.Q:
                        for Y in pda.Gamma:
                            if (State((p, X)), State((r, Y))) in U:
                                npda.add(U[State((p, X)), State((r, Y))] * w, r, a, q, head, *(Y,))
                elif len(head) != 1:
                    npda.add(w, p, a, q, head, *body)
        else:
            raise NotImplementedError

        return npda

    def _precomputation(pda, strategy="bottom-up"):
        """Computes the weights of all non-scanning push computations"""
        eps_pda = pda.spawn()
        for (p, a, q, head, body), w in pda.arcs:
            if a == ε:
                eps_pda.add(w, p, a, q, head, *body)
        U = Allsum(eps_pda).table(strategy=strategy)

        V = dd(lambda: dd(lambda: dd(lambda: dd(lambda: pda.R.zero))))
        for p, r, s in product(pda.Q, repeat=3):
            for Y, Z in product(pda.Gamma, repeat=2):
                V[p][Y][Z][s] += U[p][Y][r] * U[r][Z][s]
        return U, V

    def remove_nullary(pda, strategy="bottom-up"):
        """Constructs an equivalent WPDA by removing nullary transitions,
            that is, transitions of the form p, ε -- ε/w --> q, X"""

        if strategy == "bottom-up":
            npda = pda.spawn()
            one = pda.R.one

            U, V = Transformer._precomputation(pda, strategy=strategy)
            npda.S = Slash(False, pda.S)

            # set the same initial & final weights
            # TODO: this is a bit off
            for q, w in pda.I:
                npda.set_I(q, w)

            for q, w in pda.F:
                npda.set_F(q, w)

            for (p, a, q, head, body), w in pda.arcs:
                # PDA must be binarized
                assert body == None or len(body) <= 2
                repeat = 0 if body == None else len(body)
                for eps_vals in product((True, False), repeat=repeat):
                    eps = all(eps_vals) and (a == ε)
                    nbody = tuple([Slash(eps_val, nt) for eps_val, nt in zip(eps_vals, body)])
                    # remove all rules with X_ε on the RHS
                    if not eps:
                        # p, Y_ε -- a/w --> q, X_¬ε
                        if eps_vals == (True,):
                            for s, t in product(pda.Q, repeat=2):
                                nhead = (Slash((s, t), head[0]),)
                                npda.add(U[t][body[0]][p] * w, s, a, q, nhead)
                        # p, Y_ε, Z_ε -- a/w --> q, X_¬ε
                        elif eps_vals == (True, True):
                            for s, t in product(pda.Q, repeat=2):
                                nhead = (Slash((s, t), head[0]),)
                                npda.add(V[t][body[0]][body[1]][p] * w, s, a, q, nhead)
                        # p, Y_¬ε, Z_ε -- a/w --> q, X_¬ε
                        elif eps_vals == (False, True):
                            nhead = (Slash(False, head[0]),)
                            nbody = (Slash(False, body[0]),)
                            for t in pda.Q:
                                npda.add(U[t][body[1]][p] * w, t, a, q, nhead, *nbody)
                        # p, Y_ε, Z_¬ε -- a/w --> q, X_¬ε
                        elif eps_vals == (True, False):
                            for r, s, t in product(pda.Q, repeat=3):
                                nhead = (Slash((r, s), head[0]),)
                                nbody = (Slash((r, t), body[1]),)
                                npda.add(U[s][body[0]][t] * w, p, a, q, nhead, *nbody)
                        elif eps_vals == ():
                            for s in pda.Q:
                                nhead = (Slash((s, p), head[0]),)
                                npda.add(w, s, a, q, nhead, *nbody)
                        else:
                            nhead = (Slash(False, head[0]),)
                            npda.add(w, p, a, q, nhead, *nbody)

            # q, X_pp -- a/w --> q, X_¬ε
            for p, q in product(pda.Q, repeat=2):
                for X in pda.Gamma:
                    nhead = (Slash(False, X),)
                    nbody = (Slash((p, p), X),)
                    npda.add(one, q, ε, q, nhead, *nbody)
            return npda

        elif strategy == "top-down":
            npda = pda.spawn()
            raise NotImplementedError

        else:
            raise NotImplementedError

if __name__ == '__main__':

    from rayuela.base.symbol import Sym
    from rayuela.base.semiring import Real
    from rayuela.cfg.nonterminal import NT, S
    from rayuela.pda.pda import PDA
    from rayuela.pda.parser import Parser

    pda = PDA(R=Real)

    pda.set_I(State('0'), Real.one)
    pda.set_F(State('3'), Real.one)

    pda.add(Real(0.18), State('0'), Sym("a"), State('1'), (NT("Y"),))
    pda.add(Real(0.23), State('1'), ε, State('2'), (NT("X"),), *(NT("Y"),))
    pda.add(Real(0.23), State('2'), Sym("a"), State('3'), (S,), *(NT("X"),))

    print(pda)

    parser = Parser(pda)
    print(parser.parse("aa", strategy="bottom-up"))


