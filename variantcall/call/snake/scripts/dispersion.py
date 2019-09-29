#!/usr/bin/python
import math

# copy-pasted a statistics calculator
class Dispersion(object):
    """>>> a = Dispersion()
       >>> a.update(2)
       >>> a.S
       0.0
       >>> a.mean
       2.0
       >>> a.update(3)
       >>> a.mean
       2.5
       >>> import math
       >>> a.stddev() == math.sqrt(2)/2
       True
    """
    def __init__(self):
        self.n                = 0
        self.max     = float("-infinity")
        self.min     = float("+infinity")
        self.mean    = 0.0
        self.S  = 0.0
        self.sum     = 0.0

    def variance(self):
        if self.n == 1: return 0.0
        return self.S / (self.n - 1)

    def stddev(self):
        return math.sqrt(self.variance())

    def update(self, val):
        self.n += 1
        self.sum += val

        if val > self.max:
            self.max = val
        if val < self.min:
            self.min = val

        old_mean = self.mean

        self.mean = self.sum / self.n

        # TAOCP Vol 2 4.2.2
        # S(1) := 0; S(k) := S(k-1) + (x(k) - M(k-1)) * (x(k) - M(k))
        # sigma(k) = sqrt(S(k)/k-1)
        self.S = self.S + (val - old_mean) * (val - self.mean)

    def __str__(self):
        return "min=%(min)f max=%(max)f n=%(n)d avg=%(mean)f sum=%(sum)f std=%(std)f" % {
            "min":self.min, "max":self.max, "n":self.n, "sum":self.sum, "mean":self.mean, "std":self.stddev()
            }

if __name__ == "__main__":
    import sys

    spread = Dispersion()
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        num = float(line)
        spread.update(num)
    print("%s" % (spread,))
