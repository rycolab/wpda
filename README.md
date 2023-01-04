# Algorithms for Weighted Pushdown Automata

Code accompanying the EMNLP 2022 publication "[Algorithms for Weighted Pushdown Automata](https://arxiv.org/abs/2210.06884)".

---
## Installation
```bash
$ git clone https://github.com/rycolab/wpda.git
$ cd wpda
$ pip install -e .
```
---
## Example
```python
from rayuela.base.symbol import Sym
from rayuela.base.semiring import Real
from rayuela.fsa.state import State
from rayuela.cfg.nonterminal import NT, S
```
### Random WPDA
Create a random WPDA:
```python
from rayuela.pda.random_pda import random_pda

Sigma = [Sym("a"), Sym("b"), Sym("c")]
Gamma = [S, NT("X"), NT("Y"), NT("Z")]
Q = [State('0'), State('1'), State('2')]

pda = random_pda(Sigma=Sigma, Gamma=Gamma, Q=Q, R=Real, num_transitions=3)
print(pda)
```
Output:
```bash
1 -- a --> 0	[Y, S]	[X, S, X, Z]	0.036
2 -- b --> 0	[X]	[S]	0.135
2 -- a --> 0	[S, Z, X, Y]	[]	0.092
```
### Stringsums
```python
from rayuela.pda.pda import PDA
from rayuela.pda.parser import Parser

pda = PDA(R=Real)
# initial and final states
pda.set_I(State('0'), Real.one)
pda.set_F(State('2'), Real.one)
# add transitions
pda.add(Real(0.18), State('0'), Sym("a"), State('1'), (NT("X"),))
pda.add(Real(0.23), State('1'), Sym("a"), State('2'), (S,), *(NT("X"),))

print(pda)
# compute stringsum
parser = Parser(pda)
print(parser.parse("aa", strategy="bottom-up"))
```
Output:
```bash
0 -- a --> 1	[]	[X]	0.18
1 -- a --> 2	[X]	[S]	0.23

0.0414
```
---
## Cite

Alexandra Butoi, Brian DuSell, Tim Vieira, Ryan Cotterell, and David Chiang. 2022. [Algorithms for weighted pushdown
automata](https://arxiv.org/abs/2210.06884). In _Proceedings of the 2022 Conference on Empirical Methods in
Natural Language Processing_. Association for Computational Linguistics.
```
@inproceedings{butoi+al.emnlp22,
     abstract = {Weighted pushdown automata (WPDAs) are at the core of many natural language processing tasks, like syntax-based statistical machine translation and transition-based dependency parsing. As most existing dynamic programming algorithms are designed for context-free grammars (CFGs), algorithms for PDAs often resort to a PDA-to-CFG conversion. In this paper, we develop novel algorithms that operate directly on WPDAs. Our algorithms are inspired by Lang's algorithm, but use a more general definition of pushdown automaton and either reduce the space requirements by a factor of |Γ|(the size of the stack alphabet) or reduce the runtime by a factor of more than |Q| (the number of states). When run on the same class of PDAs as Lang's algorithm, our algorithm is both more space-efficient by a factor of |Γ| and more time-efficient by a factor of |Q|⋅|Γ|.},
     address = {Abu Dhabi, United Arab Emirates},
     arxiv = {https://arxiv.org/abs/2210.06884},
     author = {Butoi, Alexandra and 
    DuSell, Brian and 
    Vieira, Tim and 
    Cotterell, Ryan and 
    Chiang, David},
     booktitle = {Proceedings of the 2022 Conference on Empirical Methods in Natural Language Processing},
     code = {https://github.com/rycolab/wpda},
     month = {December},
     publisher = {Association for Computational Linguistics},
     title = {Algorithms for Weighted Pushdown Automata},
     url = {https://arxiv.org/abs/2210.06884},
     venue = {EMNLP},
     year = {2022}
}
```
---
## Contact
To ask questions or report problems, please open an [issue](https://github.com/rycolab/wpda/issues).


