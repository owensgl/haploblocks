#!/usr/bin/env python

"""
  A dirty script that combines information from two complementary Tab-separated value
  files which have metadata about samples in the experiment.

  outputs metadata in yaml form:

  samples:
    samplename: {k:v, k:v}
    ...
"""
  

import sys
import yaml
import re



POPULATIONS="samples/populations.txt"
SPECIES="samples/WGS2018.txt"

all_samples = {}

with open(SPECIES, "r") as fd:
    for lineno, line in enumerate(fd):
        if lineno == 0:
            continue
        line = line.strip()
        bam, samplename, species, population, is_wild, project = line.split('\t')
        all_samples[samplename] = {
            "name": samplename,
            "species": species,
            "population": population,
            "is_wild": is_wild,
            "project": project
        }


sample_re = re.compile(r"([^0-9]+)([0-9]+)-([0-9]+)")
with open(POPULATIONS, "r") as fd:
    for lineno, line in enumerate(fd):
        if lineno == 0: continue
        line = line.strip()
        individualdash, population = line.split('\t')
        m = sample_re.match(individualdash)

        if m is None:
            samplename = individualdash
            info = all_samples.get(samplename, {"name": samplename})
            info["population"] = population
        else:
            prefix, start, end = m.groups()
            padlength = len(start)
            start = int(start, 10)
            end = int(end, 10)

            for i in range(start, end + 1):
                num_s = str(i)
                samplename = prefix + ("0" * (padlength - len(num_s))) + num_s
                info = all_samples.get(samplename, {"name": samplename})
                info["population"] = population

for samplename in all_samples:
    sample = all_samples[samplename]
    sample.setdefault('project',    None)
    sample.setdefault('is_wild',    None)
    sample.setdefault('population', None)
    sample.setdefault('species',    None)

obj = {"samples": all_samples}
yaml.dump(obj, sys.stdout)

