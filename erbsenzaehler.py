# -*- coding: utf-8 -*-
from __future__ import division
import sys
import serial


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
	# this should be doing a mixture of gaussians
	# this should keep track of all increments to recompute sigmas
	# this should recluster on demand or when needed

	# TODO: need a new class to emit increments and track/calibrate otherwise

	def __init__(self, new_abs=0.1, new_rel=0.1):
		# exactly one candidate category should score within relative range for a count
		# for best candidate, update / round(count) should match exactly one weight, its own
		# thresholds to trigger a new category (both need to be overstepped)
		self.th_abs = new_abs
		self.th_rel = new_rel

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
		# check if it fits into exactly one category
		# if not, could it be a new category?
		# otherwise, do nothing


		# categorize, or init new category
		old_total = self.total()
		update = new_total - old_total

		# update all equally
		if abs(update) < self.th_abs:
			factor = new_total / old_total
			for w in self.weights:
				w.update(float(w) * factor)
			return

		multiples = [ update / float(w) for w in self.weights ]
		rel_errors = [abs(m - round(m)) for m in multiples ] # fractional error to nearest multiple
		abs_errors = [ e * float(w) for e,w in zip(rel_errors, self.weights) ] # how much that really is

		if sum(e < self.th_rel for e in rel_errors) > 1:
			# ambiguous. can't work with that.
			return

		# best (only) match in terms of being a full multiple
		closest = argmin(abs_errors)
		if closest is None:
			self.add_kind(update, 1)
			return

		dcount = int(round(multiples[closest]))
		error_rel = rel_errors[closest]
		error_abs = abs_errors[closest]

		if error_abs > self.th_abs:
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

	zero = AverageValue(alpha=0.02, fmt=u"{0:9.1f}") # in counts
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

