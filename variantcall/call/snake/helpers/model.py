
class Contig(object):
    # pos starts at 1
    # end is inclusive in range
    # length is length of containing chromosome (same for consecutive contigs within one chromosome)
    __slots__ = ("chrom", "pos", "end", "length")
    def __init__(self, chrom, pos, end, length):
        self.chrom = chrom
        self.pos = pos
        self.end = end
        self.length = length

    def to_interval(self):
        return "%(chrom)s:%(pos)09d-%(end)09d" % {
            "chrom": self.chrom,
            "pos": self.pos,
            "end": self.end
        }

    @property
    def size(self):
        return self.end - self.pos + 1

    def to_file(self):
        return "%(chrom)s %(pos)d %(end)d %(len)d" % {
            "chrom": self.chrom,
            "pos": self.pos,
            "end": self.end,
            "len": self.length
        }

    def to_gatk_interval(self):
        # gatk interval ranges are end-inclusive
        return "%s:%09d-%09d" % (self.chrom, self.pos, self.end)

    def to_wildcard(self):
        """returns a string usable in a filename"""
        return "%s-%09d-%09d" % (self.chrom, self.pos, self.end)

    @classmethod
    def from_wildcard(klass, chr_start_end):
        """
        transforms Chr01-000000001-100000000 into  -> Contig("Chr01", 1, 100000000, 0)
        """
        toks = chr_start_end.split("-")
        try:
            return klass(toks[0], int(toks[1], 10), int(toks[2], 10), 0)
        except (ValueError, IndexError):
            raise Exception("wildcard %s is not a valid contig specification" % (chr_start_end,))

    def same_start_end(self, other):
        """compare this contig with the other. return True iff they both
           are on the same chromosome, and have same start and end pos.
        """
        if self.chrom != other.chrom:
            return False
        if self.pos != other.pos:
            return False
        if self.end != other.end:
            return False
        return True

    @staticmethod
    def contig_window(contig_slice, a, b):
        """generate contigs from the array contig_slice, starting from a ending at b,
           inclusively.

           yields (i, contig)
        """
        start_i = -1
        end_i = -1

        for contig_i, contig in enumerate(contig_slice):
            if start_i < 0:
                if not contig.same_start_end(a):
                    # skip until we find start
                    continue
                # found start
                start_i = contig_i

            yield (contig_i, contig)

            if not contig.same_start_end(b):
                continue
            # found end
            end_i = contig_i
            break
        else:
            if start_i < 0:
                raise Exception("Cannot find start contig %s in slice" % (a.to_wildcard(),))
            if end_i < 0:
                raise Exception("Cannot find end contig %s in slice" % (a.to_wildcard(),))

