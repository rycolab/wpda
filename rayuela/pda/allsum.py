from collections import defaultdict as dd
from itertools import product

from rayuela.base.semiring import Real, Rational

class Allsum:

	def __init__(self, pda):
		self.pda = pda
		self.R = pda.R

	def sum(self, strategy="bottom-up"):
		i = [p for p, _ in self.pda.I]
		f = [p for p, _ in self.pda.F]
		assert len(i) == 1
		assert len(f) == 1
		table = self.table(strategy)
		return table[i[0]][self.pda.S][f[0]]

	def table(self, strategy="bottom-up"):
		if strategy == "bottom-up":
			return self.forwardchain()
		elif strategy == "top-down":
			raise NotImplementedError
		else:
			raise NotImplementedError

	def _judge_of_the_change(self, U, V, tol):

		if self.pda.R is Real or self.pda.R is Rational:
			total = 0.0
			for p, q in product(self.pda.Q, repeat=2):
				for X in self.pda.Sigma:
					val1, val2 = U[p][X][q], V[p][X][q]
					total += abs(float(val1) - float(val2))
				if total < tol:
					return True
				return False
		elif self.R.idempotent:
			for p, q in product(self.pda.Q, repeat=2):
				for X in self.pda.Gamma:
					if U[p][X][q] != V[p][X][q]:
						return False
			return True
		else:
			raise NotImplementedError

	def _bottom_up_step(self, V):
		U = dd(lambda : dd(lambda : dd(lambda : self.R.zero)))

		for (p, a, q, head, body), w in self.pda.arcs:
			if body is not None and len(body) > 0:
				for states in product(self.pda.Q, repeat=len(body)):
					update = w
					for i in range(len(states) - 1):
						update *= V[states[i]][body[i]][states[i+1]]
					update *= V[states[-1]][body[-1]][p]

					U[states[0]][head[0]][q] += update
			else:
				update = w
				U[p][head[0]][q] += update

		return U

	def forwardchain(self, tol=1e-100, timeout=1000):
		V = dd(lambda : dd(lambda : dd(lambda : self.R.zero)))

		counter = 0
		while counter < timeout:

			U = self._bottom_up_step(V)
			if self._judge_of_the_change(U, V, tol):
				return V
			V = U
			counter += 1

		return V