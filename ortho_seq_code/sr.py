import numpy as np

import collections
import functools


class memoized(object):
    '''Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    '''

    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)
        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value

    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__

    def __get__(self, obj, objtype):
        '''Support instance methods.'''
        return functools.partial(self.__call__, obj)


@memoized
def recursive_array_builder(ord_one, ord_two, d):
    if ord_one - ord_two == 0:
        return 0
    else:
        return [
            recursive_array_builder(ord_one - 1, ord_two, d) for i in range(d)]


def for_loop_builder(
        ip, depth, ord_one, ord_two, d, index_list, ar_one, ar_two):
    if depth == 0:
        command, sec_comm = command_builder(ord_one, ord_two, index_list)
        result = eval(sec_comm)

        if ord_one == ord_two:
            ip += result
        else:
            exec(command + "result")

        return ip
    else:
        for a in range(d):
            index_list[ord_one - depth] = a
            ip = for_loop_builder(
                ip, depth - 1, ord_one, ord_two, d, index_list, ar_one, ar_two)

        return ip


def command_builder(ord_one, ord_two, index_list):
    new_command = "ip"
    for i in range(ord_one - ord_two):
        # index_list contains the current
        # position in the loop for each i variable
        new_command += "[" + str(index_list[i]) + "]"

    first_com = new_command

    new_command += " = " + first_com + " + "  # can't do += when using exec

    # ar_one and ar_two are like a and b in inner_general
    second_command = "ar_one"

    for i in range(ord_one):
        second_command += "[" + str(index_list[i]) + "]"

    second_command += "*ar_two"

    for j in range(ord_two):
        second_command += "[" + str(index_list[ord_one - ord_two + j]) + "]"  #

    return new_command, second_command


def inner_general(a, b):
    if type(a) == list:
        a = np.array(a)

    if type(b) == list:
        b = np.array(b)

    ord_one = len(a.shape)
    ord_two = len(b.shape)
    a = a.tolist()  # tolist turns np.array to a normal python list
    b = b.tolist()
    d = len(a)
    # print(d)
    # print(ord_one)
    ip = np.array(recursive_array_builder(ord_one, ord_two, d))
    ip = ip.tolist()
    ip = for_loop_builder(
        ip, ord_one, ord_one, ord_two, d, [0] * ord_one, a, b)
    return ip


@memoized
def outer_array_builder(order, d):
    if order == 0:
        return 0
    else:
        return [outer_array_builder(order - 1, d) for i in range(d)]


def outer_loop_builder(
        op, depth, ord_one, ord_two, d, index_list, ar_one, ar_two):
    if depth == 0:
        first_com = outer_command_builder(ord_one, ord_two, index_list)
        exec(first_com)

        return op
    else:
        for a in range(d):
            index_list[ord_one - depth] = a
            op = outer_loop_builder(
                op, depth - 1, ord_one, ord_two, d, index_list, ar_one, ar_two)

        return op


def outer_command_builder(ord_one, ord_two, index_list):
    new_command = "op"
    for i in range(ord_one + ord_two):
        # index_list contains the current position in loop for each i variable
        new_command += "[" + str(index_list[i]) + "]"

    new_command += " = "
    # ar_one and ar_two are like a and b in inner_general
    second_command = "ar_one"
    for i in range(ord_one):
        second_command += "[" + str(index_list[i]) + "]"

    second_command += "*ar_two"

    for j in range(ord_two):
        second_command += "[" + str(index_list[ord_one + j]) + "]"  #

    return new_command + second_command


def outer_general(A, B):
    if type(A) == list():
        A = np.array(A)

    if type(B) == list():
        B = np.array(B)

    ord_one = len(A.shape)
    ord_two = len(B.shape)
    d = len(A)
    A = A.tolist()
    B = B.tolist()
    order = ord_one + ord_two
    op = np.array(outer_array_builder(order, d))
    op = op.tolist()
    op = outer_loop_builder(op, order, ord_one, ord_two, d, [0] * order, A, B)
    return np.array(op)
