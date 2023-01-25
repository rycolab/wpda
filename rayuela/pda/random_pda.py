import random

from rayuela.base.symbol import ε
from rayuela.base.misc import _random_weight as rw

from rayuela.pda.pda import PDA


def random_pda(Sigma, Gamma, Q, R, num_transitions=5):
    pda = PDA(R=R)
    Sigma = set(Sigma).union({ε})

    for _ in range(num_transitions):
        p, q = random.choice(Q), random.choice(Q)
        head, body = [], []
        a = random.choice(list(Sigma))
        for _ in range(random.randint(0, 5)):
            head.append(random.choice(Gamma))
        for _ in range(random.randint(0, 5)):
            body.append(random.choice(Gamma))
        pda.add(rw(R), p, a, q, tuple(head), *tuple(body))

    return pda


def random_bottom_up(Sigma, Gamma, Q, R, num_transitions=5):
    pda = PDA(R=R)
    Sigma = set(Sigma).union({ε})

    i, f = random.choice(Q), random.choice(Q)
    pda.set_I(i, R.one)
    pda.set_F(f, R.one)

    for _ in range(num_transitions):
        p, q = random.choice(Q), random.choice(Q)
        head, body = random.choice(Gamma), []
        a = random.choice(list(Sigma))

        for _ in range(random.randint(0, 5)):
            body.append(random.choice(Gamma))
        pda.add(rw(R), p, a, q, (head, ), *tuple(body))

    return pda


def random_top_down(Sigma, Gamma, Q, R, num_transitions=5):
    pda = PDA(R=R)
    Sigma = set(Sigma).union({ε})

    i, f = random.choice(Q), random.choice(Q)
    pda.set_I(i, R.one)
    pda.set_F(f, R.one)

    for _ in range(num_transitions):
        p, q = random.choice(Q), random.choice(Q)
        head, body = [], random.choice(Gamma)
        a = random.choice(list(Sigma))

        for _ in range(random.randint(0, 5)):
            head.append(random.choice(Gamma))
        pda.add(rw(R), p, a, q, tuple(head), *(body, ))

    return pda


def random_bottom_up_cnf(Sigma, Gamma, Q, R, num_transitions=5):
    pda = PDA(R=R)
    Sigma = set(Sigma).union({ε})

    i, f = random.choice(Q), random.choice(Q)
    pda.set_I(i, R.one)
    pda.set_F(f, R.one)

    for _ in range(num_transitions):
        p, q = random.choice(Q), random.choice(Q)
        head, body = random.choice(Gamma), []
        a = random.choice(list(Sigma))
        if a == ε:
            body = (random.choice(Gamma), random.choice(Gamma))
        else:
            for _ in range(random.randint(0, 2)):
                body.append(random.choice(Gamma))

        pda.add(rw(R), p, a, q, (head,), *tuple(body))

    return pda


def random_top_down_cnf(Sigma, Gamma, Q, R, num_transitions=5):
    pda = PDA(R=R)
    Sigma = set(Sigma).union({ε})

    i, f = random.choice(Q), random.choice(Q)
    pda.set_I(i, R.one)
    pda.set_F(f, R.one)

    for _ in range(num_transitions):
        p, q = random.choice(Q), random.choice(Q)
        head, body = [], random.choice(Gamma)
        a = random.choice(list(Sigma))
        if a == ε:
            head = (random.choice(Gamma), random.choice(Gamma))
        else:
            for _ in range(random.randint(0, 2)):
                head.append(random.choice(Gamma))

        pda.add(rw(R), p, a, q, tuple(head), *(body, ))

    return pda