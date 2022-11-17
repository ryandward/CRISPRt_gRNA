#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import collections
import copy
import itertools
import logging
import os.path
import re
import subprocess
import sys

from Bio import SeqIO


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class Error(Exception):
  pass

class SampleError(Error):
  pass


DNA_PAIRINGS = str.maketrans('atcgATCG', 'tagcTAGC')

def revcomp(x):
  return x.translate(DNA_PAIRINGS)[::-1]


def bool_from_rev(x):
  if type(x) is bool:
    return x
  if x == 'rev':
    return True
  elif x == 'fwd':
    return False
  else:
    logging.fatal('Unrecognized bool signifier {x}'.format(**vars()))


def bool_from_sense(x):
  if type(x) is bool:
    return x
  if x == 'sense':
    return True
  elif x == 'anti':
    return False
  else:
    logging.fatal('Unrecognized bool signifier {x}'.format(**vars()))


def none_or_bool(x):
  if x == 'None':
    return None
  if type(x) is bool:
    return x
  elif x == 'True':
    return True
  elif x == 'False':
    return False
  else:
    logging.fatal('Unrecognized bool signifier {x}'.format(**vars()))


def none_or_int(x):
  if x == 'None':
    return None
  if type(x) is int:
    return x
  else:
    return int(x)


def none_or_str(x):
  if x == 'None':
    return None
  if type(x) is str:
    return x
  else:
    return str(x)


class grna_target(object):
  def __init__(self, target, pam, chrom, start, end, reverse):
    self.gene = None
    self.offset = None
    self.target = str(target)
    self.pam = str(pam)
    self.chrom = str(chrom)
    self.start = int(start)
    self.end = int(end)
    self.reverse = bool_from_rev(reverse)
    self.sense_strand = None
    self.specificity = 0

  @classmethod
  def from_tsv(cls, tsv, sep='\t'):
    """Alternate factory constructor from serialized grna_target string.

    Intended for reading back in the result of grna_target.__str__()
    """
    (gene,
     offset,
     target,
     pam,
     chrom,
     start,
     end,
     reverse,
     sense_strand,
     specificity) = tsv.strip().split(sep)
    t = grna_target(target, pam, chrom, start, end, reverse)
    t.gene = none_or_str(gene)
    t.offset = none_or_int(offset)
    t.sense_strand = bool_from_sense(sense_strand)
    t.specificity = none_or_int(specificity)
    return t

  @classmethod
  def header(cls, sep='\t'):
    return sep.join([
      'gene',
      'offset',
      'target',
      'pam',
      'chrom',
      'start',
      'end',
      'repldir',
      'transdir',
      'specificity'])

  def __str__(self, sep='\t'):
    if self.reverse:
      rev_str = 'rev'
    else:
      rev_str = 'fwd'
    if self.sense_strand:
      sst_str = 'sense'
    else:
      sst_str = 'anti'
    return sep.join([str(x) for x in [
      self.gene,
      self.offset,
      self.target,
      self.pam,
      self.chrom,
      self.start,
      self.end,
      rev_str,
      sst_str,
      self.specificity]])

  def id_str(self, sep=';'):
    return sep.join([str(x) for x in [
      self.target,
      self.pam,
      self.chrom,
      self.start,
      self.reverse]])

  def sequence_with_pam(self):
    """DNA sequence with trailing PAM in place."""
    return self.target + self.pam

  def pam_with_sequence(self):
    """PAM with DNA sequence trailing in place."""
    return self.pam + self.target
