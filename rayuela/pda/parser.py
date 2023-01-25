from itertools import product
from collections import defaultdict as dd

from rayuela.base.symbol import ε

from rayuela.pda.allsum import Allsum

class Parser:

    def __init__(self, pda):
        self.pda = pda
        self.R = self.pda.R

    def _bottom_up_parsing(self, input):
        """Implements the bottom-up parsing algorithm from EMNLP 2022
        submission `Algorithms for Weighted Pushdown Automata`"""
        # WPDA must be in normal form
        assert self.pda.in_cnf(strategy="bottom-up")
        print(self.pda.in_cnf(strategy="bottom-up"))
        N = len(input)

        # initialization
        β = dd(lambda: dd(lambda: dd(lambda: dd(lambda: dd(lambda: self.R.zero)))))

        #0-pop, 1-push rule
        for i in range(N):
            for (p, a, q, head, body), w in self.pda.arcs:
                if (body == None or len(body) == 0) and a.sym == input[i]:
                    β[i][p][head[0]][i + 1][q] += w

        for span in range(2, N + 1):
            for i in range(N - span + 1):
                j = i + span
                # 1-pop, 1-push
                for p in self.pda.Q:
                    for (r, a, q, head, body), w in self.pda.arcs:
                        if body is not None:
                            if len(body) == 1 and a.sym == input[j-1]:
                                β[i][p][head[0]][j][q] += β[i][p][body[0]][j - 1][r] * w
                # 2-pop, 1-push
                for p, r in product(self.pda.Q, repeat=2):
                    for (s, a, q, head, body), w in self.pda.arcs:
                        if body is not None:
                            if len(body) == 2:
                                # a = s_j
                                if a.sym == input[j - 1]:
                                    for k in range(i + 1, j - 1):
                                        β[i][p][head[0]][j][q] += β[i][p][body[0]][k][r] * β[k][r][body[1]][j - 1][s] * w
                                elif a == ε:
                                    for k in range(i + 1, j):
                                        β[i][p][head[0]][j][q] += β[i][p][body[0]][k][r] * β[k][r][body[1]][j][s] * w

        # PDA must have a single initial and final state
        i = [p for p, _ in self.pda.I]
        f = [p for p, _ in self.pda.F]
        assert len(i) == 1
        assert len(f) == 1

        return β[0][i[0]][self.pda.S][N][f[0]]

    def _top_down_parsing(self, input):
        """Implements the top-down parsing algorithm from EMNLP 2022
        submission `Algorithms for Weighted Pushdown Automata`"""
        # WPDA must be in normal form
        assert self.pda.in_cnf(strategy="top-down")

        N = len(input)

        # initialization
        β = dd(lambda: dd(lambda: dd(lambda: dd(lambda: dd(lambda: self.R.zero)))))

        #1-pop, 0-push rule
        for i in range(1, N + 1):
            for (p, a, q, head, body), w in self.pda.arcs:
                if len(head) == 0 and a.sym == input[i]:
                    β[i - 1][p][head[0]][i][q] += w

        for span in range(2, N + 1):
            for j in range(N - span + 2):
                i = j - span
                # 1-pop, 1-push
                for q in self.pda.Q:
                    for (p, a, r, head, body), w in self.pda.arcs:
                        if len(body) == 1 and a.sym == input[i]:
                            β[i][p][body[0]][j][q] += β[i + 1][r][head[0]][j][q] * w
                # # 2-pop, 1-push
                for s, q in product(self.pda.Q, repeat=2):
                    for (p, a, r, head, body), w in self.pda.arcs:
                        if len(head) == 2:
                            if a.sym == input[j]:
                                for k in range(i + 1, j - 1):
                                    β[i][p][body[0]][j][q] += β[i + 1][r][head[0]][k][s] * β[k][s][head[1]][j][q] * w
                            elif a == ε:
                                for k in range(i, j - 1):
                                    β[i][p][body[0]][j][q] += β[i][r][head[0]][k][s] * β[k][s][head[1]][j][q] * w

        # PDA must have a single initial and final state
        i = [p for p, _ in self.pda.I]
        f = [p for p, _ in self.pda.F]
        assert len(i) == 1
        assert len(f) == 1

        return β[0, i[0], self.pda.S, N, f[0]]

    def _bottom_up_parsing_v2(self, input):
        """Implements the bottom-up parsing algorithm from EMNLP 2022
        submission `Algorithms for Weighted Pushdown Automata`"""
        npda, eps_pda = self.pda.spawn(), self.pda.spawn()

        for q, w in self.pda.I:
            npda.set_I(q, w)
        for q, w in self.pda.F:
            npda.set_F(q, w)

        for (p, a, q, head, body), w in self.pda.arcs:
            if a == ε:
                eps_pda.add(w, p, a, q, head, *body)
            else:
                npda.add(w, p, a, q, head, *body)

        N = len(input)
        # initialization
        β = dd(lambda: dd(lambda: dd(lambda: dd(lambda: dd(lambda: self.R.zero)))))

        # Initialize 0 len spans to weights of non-scanning push computations
        U = Allsum(eps_pda).table(strategy="bottom-up")

        for i in range(N + 1):
            for p, q in product(self.pda.Q, repeat=2):
                for X in self.pda.Gamma:
                    β[i][p][X][i][q] += U[p][X][q]

        #0-pop, 1-push rule
        for i in range(N):
            for (p, a, q, head, body), w in npda.arcs:
                if (body == None or len(body) == 0) and a.sym == input[i]:
                    β[i][p][head[0]][i + 1][q] += w

        for span in range(1, N + 1):
            for i in range(N - span + 1):
                j = i + span
                # 1-pop, 1-push
                for p in self.pda.Q:
                    for (r, a, q, head, body), w in npda.arcs:
                        if body is not None:
                            if len(body) == 1 and a.sym == input[j-1]:
                                β[i][p][head[0]][j][q] += β[i][p][body[0]][j - 1][r] * w
                # 2-pop, 1-push
                for p, r in product(self.pda.Q, repeat=2):
                    for (s, a, q, head, body), w in npda.arcs:
                        if body is not None:
                            if len(body) == 2:
                                # a = s_j
                                if a.sym == input[j - 1]:
                                    for k in range(i, j):
                                        β[i][p][head[0]][j][q] += β[i][p][body[0]][k][r] * β[k][r][body[1]][j - 1][s] * w

        # PDA must have a single initial and final state
        i = [p for p, _ in self.pda.I]
        f = [p for p, _ in self.pda.F]
        assert len(i) == 1
        assert len(f) == 1

        return β[0][i[0]][self.pda.S][N][f[0]]

    def parse(self, input, strategy="bottom-up"):
        if strategy == "bottom-up":
            return self._bottom_up_parsing_v2(input)
        elif strategy == "top-down":
            return self._top_down_parsing(input)
        else:
            raise NotImplementedError

