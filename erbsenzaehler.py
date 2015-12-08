# -*- coding: utf-8 -*-
from __future__ import division
import sys
import serial


def argmin(dictionary, key=None):
	if not dictionary:
		return None

	if key is None:
		keyfn = lambda k: dictionary[k]
	else:
		keyfn = lambda k: key(dictionary[k])

	return min(dictionary, key=keyfn)

class SettlingValue(object):
	def __init__(self, range=None, threshold=None, history=5):
		self.threshold = threshold
		self.range = range
		assert not (range is None and threshold is None)
		self.values = (None,) * history
	
	def update(self, newval):
		self.values = self.values[1:] + (newval,)

	def is_settled(self):
		if not all(x is not None for x in self.values):
			return False

		#nn = [x for x in self.values if x is not None]
		nn = self.values
		u = min(nn)
		v = max(nn)
		if self.range:
			return (v-u) < self.range
		elif self.threshold:
			return (v-u) < self.threshold * (u+v) * 0.5

	settled = property(is_settled)

	def __float__(self):
		return float(self.values[-1] or 0.0)

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

class ValueStream(object):
	def __init__(self, fd):
		self.fd = fd

	def __iter__(self):
		return self

	def next(self):
		while True:
			line = self.fd.readline()

			try:
				return int(line)
			except ValueError, e:
				continue

class MixedBag(object):
	def __init__(self, new_factor=0.1):
		self.new_factor = new_factor
		self.clear()

	def clear(self):
		self.kinds = 0
		self.counts = []
		self.weights = []

	def add_kind(self, weight=None, count=0):
		self.counts.append(count)
		self.weights.append(AverageValue(init=weight, alpha=1, fmt="{0:+7.3f} g"))
		self.kinds += 1

	def remove_kind(self, kind):
		del self.counts[kind]
		del self.weights[kind]
		self.kinds -= 1

	def __iter__(self):
		for k in xrange(self.kinds):
			yield float(self.weights[k]), self.counts[k]

	def total(self, *keys):
		if not keys:
			return sum( w*c for k,(w,c) in enumerate(self))
		else:
			return sum( w*c for k,(w,c) in enumerate(self) if k in keys)

	def update(self, new_total):
		# categorize, or init new category
		old_total = self.total()
		update = new_total - old_total

		multiples = {
			k: update / float(self.weights[k])
			for k in xrange(self.kinds)
		}

		# best match in terms of being a full multiple
		closest = argmin(multiples, key=lambda m: abs(m - round(m)))
		if closest is None:
			self.add_kind(update, 1)
			return

		# update all equally
		if all(round(multiples[k]) == 0 and abs(multiples[k]) <= self.new_factor for k in multiples):
			factor = new_total / old_total
			for w in self.weights:
				w.update(float(w) * factor)
			return

		dcount = int(round(multiples[closest]))
		error = multiples[closest] - dcount

		if abs(error) > self.new_factor:
			if update > 0:
				self.add_kind(update, 1)
				return

		if self.counts[closest] + dcount >= 0:
			# compute what must be the remainder
			kindsum = update + self.total(closest)
			self.counts[closest] += dcount
			if self.counts[closest] == 0:
				self.remove_kind(closest)
			elif self.counts[closest] > 0:
				self.weights[closest].update(kindsum / self.counts[closest])


if __name__ == '__main__':
	scale = 720.0e-6 # calibration parameter

	current = SettlingValue(range=0.1, history=4) # grams

	zero = AverageValue(fmt=u"{0:9.1f}") # in counts
	tarethreshold = 0.1 # grams

	# usable weight, neq 0
	weight = AverageValue(alpha=0.1, fmt=u"{0:+7.2f} g") # grams

	bag = MixedBag()

	# obsolete
	singleweight = AverageValue(alpha=0.2, fmt="{0:+7.3f} g") # grams
	count = 0

	comport = sys.argv[1]

	ser = serial.Serial(comport, baudrate=9600, timeout=0.2)
	assert ser.isOpen()

	for rawvalue in ValueStream(ser):
		print u"raw: {0:8.0f}, zero {1:s}, current {2:+7.2f} g, weight {3:s} -> {4:s}".format(
			rawvalue,
			zero,
			float(current),
			weight,
			", ".join("{1}x {0:.2f}".format(w,c) for w,c in bag)
		)

		currentvalue = (rawvalue - float(zero)) * scale
		current.update(currentvalue)

		if not current.is_settled():
			weight.clear()
			continue

		if not zero.is_valid() or abs(float(current)) < tarethreshold:
			zero.update(rawvalue)
			weight.clear()
			bag.clear()
			count = 0
			continue

		weight.update(float(current))

		#if not weight.is_settled():
		#	continue

		if float(weight) < 0:
			continue

		bag.update(float(weight))

