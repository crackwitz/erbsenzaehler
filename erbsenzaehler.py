from __future__ import division
import sys
import serial


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
	def __init__(self, alpha=0.05, fmt="f"):
		self.fmt = fmt
		self.alpha = alpha
		self.clear()

	def is_valid(self):
		return (self.value is not None)

	def clear(self):
		self.value = None

	def update(self, newval):
		if self.is_valid():
			self.value += (newval - self.value) * self.alpha
		else:
			self.value = newval

	def __float__(self):
		return float(self.value or 0.0)

	def __str__(self):
		res = self.fmt.format(float(self))
		
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

scale = 720.0e-6

current = SettlingValue(range=0.05, history=4) # grams

zero = AverageValue(fmt="{0:9.1f}") # in counts
tarethreshold = 0.1 # grams

weight = AverageValue(alpha=0.1, fmt="{0:+7.2f} g") # grams

singleweight = AverageValue(alpha=0.2, fmt="{0:+7.3f} g") # grams
count = 0

if __name__ == '__main__':
	comport = sys.argv[1]

	ser = serial.Serial(comport, baudrate=9600, timeout=0.2)
	assert ser.isOpen()

	for rawvalue in ValueStream(ser):
		print "raw: {0:8.0f}, zero {1:s}, current {2:+7.2f} g, weight {3:s}, {4:6.2f} x {5:s}".format(
			rawvalue,
			zero,
			float(current),
			weight,
			count,
			singleweight
		)

		currentvalue = (rawvalue - float(zero)) * scale
		current.update(currentvalue)

		if not current.is_settled():
			weight.clear()
			continue

		if not zero.is_valid() or abs(float(current)) < tarethreshold:
			zero.update(rawvalue)
			weight.clear()
			singleweight.clear()
			count = 0
			continue

		weight.update(float(current))

		#if not weight.is_settled():
		#	continue

		if float(weight) < 0:
			continue

		if not singleweight.is_valid():
			count = 1.0
		else:
			count = float(weight) / float(singleweight)

		rcount = round(count)

		if rcount > 0:
			singleweight.update(float(weight) / rcount)

		# if abs(currentvalue) < tarethreshold:
		# 	zero.update(raw)
		# 	weight.clear()
		# 	continue

		# weight.update(currentvalue)

		# if not weight.is_settled():
		# 	continue

		# if singleweight.is_valid():
		# 	count = weight / singleweight
		# 	singleweight.update(weight / round(count))
		# else:
		# 	singleweight.update(weight)

