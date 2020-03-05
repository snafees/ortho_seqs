import numpy as np
from numpy import matrix
from numpy import linalg
from numpy import *
from numpy.linalg import *
import subprocess

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
        return [recursive_array_builder(ord_one - 1, ord_two, d) for i in range(d)]
        # uses the concept of list comprehension]


def for_loop_builder(ip, depth, ord_one, ord_two, d, index_list, ar_one, ar_two):
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
            ip = for_loop_builder(ip, depth - 1, ord_one, ord_two, d, index_list, ar_one, ar_two)

        return ip


# def for_loop_builder(ip, depth, max_depth, ord_one, ord_two, d, index_list, ar_one, ar_two):
#   if depth == 1: # If on last loop...
#      for a in range(d): # The loop
#         index_list[max_depth - depth] = a # Set the value for index_list. Index list tracks all values for every iteration of the loop.
#         first_com, second_com = command_builder(ord_one, ord_two, index_list)
#         result = eval(second_com) # Evaluate the multiplication part of the expression.
#         exec(first_com + "result") # Assign the result.
#
#      return ip
#   else: # Otherwise, build another loop level
#      for a in range(d): # Loop level
#         index_list[max_depth - depth] = a
#         ip = for_loop_builder(ip, depth - 1, max_depth, ord_one, ord_two, d, index_list, ar_one, ar_two) # Calls for the next level of the loop series.
#
#      return ip



def command_builder(ord_one, ord_two, index_list):
    new_command = "ip"
    for i in range(ord_one - ord_two):
        new_command += "[" + str(
                index_list[i]) + "]"  # index_list contains the current position in the loop for each i variable

    first_com = new_command

    new_command += " = " + first_com + " + "  # can't do += when using exec

    second_command = "ar_one"  # ar_one and ar_two are like a and b in inner_general

    for i in range(ord_one):
        second_command += "[" + str(index_list[i]) + "]"

    second_command += "*ar_two"

    for j in range(ord_two):
        second_command += "[" + str(index_list[ord_one - ord_two + j]) + "]"  #

    return new_command, second_command


def inner_general(a, b):
    if type(a) == type(list()):
        a = np.array(a)

    if type(b) == type(list()):
        b = np.array(b)

    ord_one = len(a.shape)
    ord_two = len(b.shape)
    a = a.tolist()  # tolist turns np.array to a normal python list
    b = b.tolist()
    d = len(a)
    # print(d)
    # print(ord_one)
    ip = array(recursive_array_builder(ord_one, ord_two, d))
    ip = ip.tolist()
    ip = for_loop_builder(ip, ord_one, ord_one, ord_two, d, [0] * ord_one, a, b)
    return ip


#



@memoized
def outer_array_builder(order, d):
    if order == 0:
        return 0
    else:
        return [outer_array_builder(order - 1, d) for i in range(d)]


def outer_loop_builder(op, depth, ord_one, ord_two, d, index_list, ar_one, ar_two):
    if depth == 0:
        first_com = outer_command_builder(ord_one, ord_two, index_list)
        exec(first_com)

        return op
    else:
        for a in range(d):
            index_list[ord_one - depth] = a
            op = outer_loop_builder(op, depth - 1, ord_one, ord_two, d, index_list, ar_one, ar_two)

        return op


def outer_command_builder(ord_one, ord_two, index_list):
    new_command = "op"
    for i in range(ord_one + ord_two):
        new_command += "[" + str(
                index_list[i]) + "]"  # index_list contains the current position in the loop for each i variable

    new_command += " = "

    second_command = "ar_one"  # ar_one and ar_two are like a and b in inner_general

    for i in range(ord_one):
        second_command += "[" + str(index_list[i]) + "]"

    second_command += "*ar_two"

    for j in range(ord_two):
        second_command += "[" + str(index_list[ord_one + j]) + "]"  #

    return new_command + second_command


def outer_general(A, B):
    if type(A) == type(list()):
        A = np.array(A)

    if type(B) == type(list()):
        B = np.array(B)

    ord_one = len(A.shape)
    ord_two = len(B.shape)
    d = len(A)
    A = A.tolist()
    B = B.tolist()
    order = ord_one + ord_two
    op = array(outer_array_builder(order, d))
    op = op.tolist()
    op = outer_loop_builder(op, order, ord_one, ord_two, d, [0] * order, A, B)
    return np.array(op)


    # def inner11(v1,v2):
    # 	return inner_general(v1, v2)

    # def inner11(v1,v2):
    #     d = len(v1)
    #     ip = 0
    #     for i in range(d):
    #         ip += v1[i]*v2[i]
    #     return ip
    #
    # def inner21(m,v):
    # 	d = len(m)
    # 	ip = array([0.0 for k in range(d)])
    # 	print("checking inner 21: ")
    # 	print(ip)
    # 	for i in range(d):
    # 		for j in range(d):
    # 			ip[i] += m[i][j]* v[j]
    # 	print("checking inner part 2 21: ")
    # 	print(ip)
    # 	return ip
    #
    # def inner22(m1,m2):
    # 	d = len(m1)
    # 	ip = 0
    # 	for i in range(d):
    # 		for j in range(d):
    # 			ip += m1[i][j]* m2[i][j]
    # 	return ip
    #
    # def inner31(t3,t1):
    # 	d = len(t3)
    # 	ip = array([[0.0 for i in range(d)] for j in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				ip[i1][i2] += t3[i1][i2][i3] * t1[i3]
    # 	return ip
    #
    # def inner32(t3,t2):
    # 	d = len(t3)
    # 	ip = array([0.0 for i in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				ip[i1] += t3[i1][i2][i3] * t2[i2][i3]
    # 	return ip
    #
    # def inner33(t3a,t3b):
    # 	d = len(t3a)
    # 	ip = 0
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				ip += t3a[i1][i2][i3] * t3b[i1][i2][i3]
    # 	return ip
    #
    # def inner41(t4,t1):
    # 	d = len(t4)
    # 	ip = array([[[0.0 for i in range(d)] for j in range(d)] for k in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					ip[i1][i2][i3] += t4[i1][i2][i3][i4] * t1[i4]
    # 	return ip
    #
    # def inner42(t4,t2):
    # 	d = len(t4)
    # 	ip = array([[0.0 for i in range(d)] for j in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					ip[i1][i2] += t4[i1][i2][i3][i4] * t2[i3][i4]
    # 	return ip
    #
    # def inner43(t4,t3):
    # 	d = len(t4)
    # 	ip = array([0.0 for i in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					ip[i1] += t4[i1][i2][i3][i4] * t3[i2][i3][i4]
    # 	return ip
    #
    # def inner44(t4a,t4b):
    # 	d = len(t4a)
    # 	ip = 0
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					ip += t4a[i1][i2][i3][i4] * t4b[i1][i2][i3][i4]
    # 	return ip
    #
    # def inner51(t5,t1):
    # 	d = len(t5)
    # 	ip = array([[[[0.0 for i in range(d)] for j in range(d)] for k in range(d)] for l in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						ip[i1][i2][i3][i4] += t5[i1][i2][i3][i4][i5] * t1[i5]
    # 	return ip
    #
    # def inner52(t5,t2):
    # 	d = len(t5)
    # 	ip = array([[[0.0 for i in range(d)] for j in range(d)] for k in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						ip[i1][i2][i3] += t5[i1][i2][i3][i4][i5] * t2[i4][i5]
    # 	return ip
    #
    # def inner53(t5,t3):
    # 	d = len(t5)
    # 	ip = array([[0.0 for i in range(d)] for j in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						ip[i1][i2] += t5[i1][i2][i3][i4][i5] * t3[i3][i4][i5]
    # 	return ip
    #
    # def inner54(t5,t4):
    # 	d = len(t5)
    # 	ip = array([0.0 for i in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						ip[i1] += t5[i1][i2][i3][i4][i5] * t4[i2][i3][i4][i5]
    # 	return ip
    #
    # def inner55(t5a,t5b):
    # 	d = len(t5a)
    # 	ip = 0
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						ip += t5a[i1][i2][i3][i4][i5] * t5b[i1][i2][i3][i4][i5]
    # 	return ip
    #
    # def inner61(t6,t1):
    # 	d = len(t6)
    # 	ip = array([[[[[0.0 for i in range(d)] for j in range(d)] for k in range(d)] for l in range(d)] for m in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						for i6 in range(d):
    # 							ip[i1][i2][i3][i4][i5] += t6[i1][i2][i3][i4][i5][i6] * t1[i6]
    # 	return ip
    #
    # def inner62(t6,t2):
    # 	d = len(t6)
    # 	ip = array([[[[0.0 for i in range(d)] for j in range(d)] for k in range(d)] for l in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						for i6 in range(d):
    # 							ip[i1][i2][i3][i4] += t6[i1][i2][i3][i4][i5][i6] * t2[i5][i6]
    # 	return ip
    #
    # def inner63(t6,t3):
    # 	d = len(t6)
    # 	ip = array([[[0.0 for i in range(d)] for j in range(d)] for k in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						for i6 in range(d):
    # 							ip[i1][i2][i3] += t6[i1][i2][i3][i4][i5][i6] * t3[i4][i5][i6]
    # 	return ip
    #
    # def inner64(t6,t4):
    # 	d = len(t6)
    # 	ip = array([[0.0 for i in range(d)] for j in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						for i6 in range(d):
    # 							ip[i1][i2] += t6[i1][i2][i3][i4][i5][i6] * t4[i3][i4][i5][i6]
    # 	return ip
    #
    # def inner65(t6,t5):
    # 	d = len(t6)
    # 	ip = array([0.0 for i in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						for i6 in range(d):
    # 							ip[i1] += t6[i1][i2][i3][i4][i5][i6] * t5[i2][i3][i4][i5][i6]
    # 	return ip
    #
    # def inner66(t6a,t6b):
    # 	d = len(t6a)
    # 	ip = 0
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						for i6 in range(d):
    # 							ip += t6a[i1][i2][i3][i4][i5][i6] * t6b[i1][i2][i3][i4][i5][i6]
    # 	return ip
    #
    # def outer11(a,b):
    # 	d = len(a)
    # 	op = array([[0.0 for i in range(d)] for j in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			op[i1][i2] = a[i1]*b[i2]
    # 	return op
    #
    #
    # def outer21(a,b):
    # 	d = len(a)
    # 	op = array([[[0.0 for i in range(d)] for j in range(d)] for k in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				op[i1][i2][i3] = a[i1][i2]*b[i3]
    # 	return op
    #
    # def outer12(a,b):
    # 	d = len(a)
    # 	op = array([[[0.0 for i in range(d)] for j in range(d)] for k in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				op[i1][i2][i3] = a[i1]*b[i2][i3]
    # 	return op
    #
    # def outer22(a,b):
    # 	d = len(a)
    # 	op = array([[[[0.0 for i in range(d)] for j in range(d)] for k in range(d)] for l in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					op[i1][i2][i3][i4] = a[i1][i2]*b[i3][i4]
    # 	return op
    #
    # def outer31(a,b):
    # 	d = len(a)
    # 	op = array([[[[0.0 for i in range(d)] for j in range(d)] for k in range(d)] for l in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					op[i1][i2][i3][i4] = a[i1][i2][i3]*b[i4]
    # 	return op
    #
    # def outer32(a,b):
    # 	d = len(a)
    # 	op = array([[[[[0.0 for i in range(d)] for j in range(d)] for k in range(d)] for l in range(d)]  for m in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						op[i1][i2][i3][i4][i5] = a[i1][i2][i3]*b[i4][i5]
    # 	return op
    #
    # def outer33(a,b):
    # 	d = len(a)
    # 	op = array([[[[[[0.0 for i in range(d)] for j in range(d)] for k in range(d)] for l in range(d)]  for m in range(d)] for n in range(d)])
    # 	for i1 in range(d):
    # 		for i2 in range(d):
    # 			for i3 in range(d):
    # 				for i4 in range(d):
    # 					for i5 in range(d):
    # 						for i6 in range(d):
    # 							op[i1][i2][i3][i4][i5][i6] = a[i1][i2][i3]*b[i4][i5][i6]
    # 	return op
