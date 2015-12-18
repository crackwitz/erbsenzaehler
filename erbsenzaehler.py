# -*- coding: utf-8 -*-
from __future__ import division
import sys
import serial
import math
import pprint; pp = pprint.pprint


def argmin(arg, key=None):
	if not arg:
		return None

	if key is None:
		keyfn = lambda k: arg[k]
	else:
		keyfn = lambda k: key(arg[k])

	if isinstance(arg, dict):
		return min(arg, key=keyfn)
	else:
		return min(xrange(len(arg)), key=keyfn)

def mean(seq):
	lst = list(seq)
	if not lst:
		return None
	else:
		return sum(lst) / len(lst)

def ValueStream(fd, scale=None):
	while True:
		try:
			line = fd.next()
		except StopIteration:
			continue

		try:
			value = int(line)
			if scale is not None:
				value *= scale
			yield value


		except ValueError, e:
			continue

class StableValue(object):
	def __init__(self, range=None, threshold=None, history=4):
		self.history = history
		self.threshold = threshold
		self.range = range
		assert not (range is None and threshold is None)
		self.clear()

	def clear(self):
		self.values = (None,) * self.history
	
	def update(self, newval):
		self.values = self.values[1:] + (newval,)
		return self.stable

	def is_stable(self):
		if any(x is None for x in self.values):
			return False

		u = min(self.values)
		v = max(self.values)
		if self.range:
			return (v-u) <= self.range
		elif self.threshold:
			return (v-u) < self.threshold * (u+v) * 0.5

	stable = property(is_stable)

	def __float__(self):
		return mean(v for v in self.values if v is not None) or 0.0

class AverageValue(object):
	def __init__(self, alpha=0.05, fmt="f", init=None):
		self.fmt = fmt
		self.alpha = alpha
		self.clear()
		if init is not None:
			self.update(init)

	def is_valid(self):
		return (self.value is not None)

	def clear(self):
		self.value = None
		self.mdev = 0

	def update(self, newval):
		if self.is_valid():
			dev = newval - self.value
			self.value += dev * self.alpha
			self.mdev += (abs(dev) - self.mdev) * self.alpha
		else:
			self.value = newval
			self.mdev = 0

	def __float__(self):
		return float(self.value or 0.0)

	def __str__(self):
		res = self.fmt.format(float(self), self.mdev or 0)
		
		if not self.is_valid():
			res = res.replace('0', '-')

		return res

def DeltaStream(input, tare=0.1):
	zero = AverageValue()
	increment = StableValue(range=tare)

	# init
	for newval in input:
		if increment.update(newval):
			zero.update(float(increment))
			increment.clear()
			break

	for newval in input:
		delta = newval - float(zero)

		if abs(delta) <= tare:
			zero.update(newval)
			increment.clear()
			continue

		if not increment.update(delta):
			# not stable
			continue

		# stable. reset zero and return averaged delta.
		zero.clear()
		zero.update(newval)
		yield float(increment)

class Cluster(object):
	# http://www.johndcook.com/blog/standard_deviation/
	def __init__(self, mdev0=None, values=()):
		self.count = 0
		self.mean = 0.0
		self.mdev = mdev0 or 0.0
		for v in values:
			self.add(v)

	def add(self, value, count=1):
		dev = value - self.mean

		self.mean = (self.mean * self.count + value * count) / (self.count + count)
		self.count += count

		if self.count > 1:
			self.mdev += (abs(dev) - self.mdev) * count / self.count

	def remove(self, count=1):
		self.count -= count

	def __cmp__(self, other):
		return (self.mean > other.mean) - (self.mean < other.mean)

	def __repr__(self):
		return "<Cluster {0:.2f} ~{1:.2f} n={2:d}>".format(self.mean, self.mdev, self.count)

class Mixture(object):
	def __init__(self, mdev0, th_merge=3):
		self.mdev0 = mdev0
		self.th_merge = th_merge
		self.clusters = set()

	def evaluate(self, value):
		err = (lambda x:
			abs(x - round(x)))
		score = (lambda c:
			err(abs(value) / c.mean) * c.mean / c.mdev)

		scores = {
			c: score(c)
			for c in self.clusters
			if (value > 0 and c.count >= 0) or (value < 0 and c.count > 0)
		}

		best = argmin(scores)
		if best in scores:
			count = int(round(abs(value) / best.mean))
		else:
			count = None

		return best, count, scores.get(best, None)

	def add(self, value):
		(best, count, sigmas) = self.evaluate(value)
		does_fit = (best is not None) and (sigmas <= self.th_merge)

		if value > 0:
			# check if it fits anywhere
			if does_fit:
				best.add(value/count, count=count)
			else:
				self.clusters.add(Cluster(values=[value], mdev0=self.mdev0(value)))

		elif value < 0:
			if does_fit:
				best.remove(count=count)
				if best.count == 0:
					self.clusters.remove(best)
			else:
				self.clusters = set()
			return


	def __repr__(self):
		return "<Mixture {0!r}>".format(sorted(self.clusters))


if __name__ == '__main__':
	comport = sys.argv[1]

	ser = serial.Serial(comport, baudrate=9600, timeout=0.2)
	assert ser.isOpen()

	scale = 720.0e-6 # calibration parameter
	deltas = DeltaStream(ValueStream(ser, scale=scale), tare=0.1)
	mixture = Mixture(mdev0=(lambda value: 0.02 * value), th_merge=3)

	for increment in deltas:
		print "{0:+6.2f}".format(increment)
		#if increment > 0:
		mixture.add(increment)

		print mixture
