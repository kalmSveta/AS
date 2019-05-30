from numpy import argmin, unravel_index, full, empty, load, set_printoptions, argwhere, argsort
from math import ceil
import binascii, itertools, sys, getopt
import time
from functools import partial
from re import sub
import pandas as pd

inf = float('inf')

# Dictionary for nts (used in 1x1, 2x1, 2x2 loops in last 2 dims)
Dic_nt = {'@': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}
# Dictionary of basepairs (used in 1x1, 2x1, 2x2 loops in first 2 dims)
Dic_bp = {'CG': 0, 'GC': 1, 'GT': 2, 'TG': 3, 'AT': 4, 'TA': 5}
# Adding for long bulges
TerminalAU = 50
# Adding for assymetric loops
ninio = 100
# Internal loop
internal_loop = [inf, inf, 100., 100., 110., 200., 200., 210., 230., 240.]

stacking_matrix = load("../lib/stacking_matrix.npy")
bulge_list = load("../lib/bulge_list.npy")
intl11_matrix = load("../lib/intl11_matrix.npy")
intl12_matrix = load("../lib/intl12_matrix.npy")
intl22_matrix = load("../lib/intl22_matrix.npy")
intl13_matrix = load("../lib/intl13_matrix.npy")
intl23_matrix = load("../lib/intl23_matrix.npy")
intl33_matrix = load("../lib/intl33_matrix.npy")


def SplitStructure(structure):
    structure1 = structure[: structure.find(')')]
    structure2 = structure[structure.find(')'):]
    return (structure1, structure2)


def AddDelitions(structure, seq, seq_compl, mut_coord, mut_seq):
    part1, part2 = SplitStructure(structure)
    i = 0
    while i != max(len(part1), len(part2)):
        if part1[len(part1) - i - 1] == '.' and part2[i] != '.':
            if mut_seq == 'seq_compl' and mut_coord >= i:
                mut_coord += 1
            part2 = part2[: i] + '-' + part2[i:]
            seq_compl = seq_compl[: i] + '-' + seq_compl[i:]
        elif part2[i] == '.' and part1[len(part1) - i - 1] != '.':
            if mut_seq == 'seq' and mut_coord > len(part1) - i - 1:
                mut_coord += 1
            part1 = part1[: len(part1) - i] + '-' + part1[len(part1) - i:]
            seq = seq[: len(seq) - i] + '-' + seq[len(seq) - i:]
        i += 1
    return (part1, part2, seq, seq_compl, mut_coord)


def Mutate(seq, seq_compl, mut_coord, mut_seq, mut_to, part1, part2):
    if mut_coord > len(seq):
        print('Mut coordsinate is too big!')
        return (0, 0, 0, 0)
    else:
        seq = list(seq)
        seq_compl = list(seq_compl)
        part1 = list(part1)
        part2 = list(part2)
        if mut_seq == 'seq':  # mutation in seq
            if part1[mut_coord] == '(':  # affects stacking pair
                if (seq_compl[len(part2) - mut_coord - 1] == 'T' and (
                        seq[mut_coord] == 'A' and mut_to == 'G' or seq[mut_coord] == 'G' and mut_to == 'A')) or (
                        seq_compl[len(part2) - mut_coord - 1] == 'G' and (
                        seq[mut_coord] == 'T' and mut_to == 'C' or seq[
                    mut_coord] == 'C' and mut_to == 'T')):  # doesn't disrupt it
                    seq[mut_coord] = mut_to
                else:  # disrupts it
                    seq[mut_coord] = mut_to
                    part1[mut_coord] = '.'
                    part2[len(part2) - mut_coord - 1] = '.'
            else:  # affects bulge or loop
                if part1[mut_coord] == '.' and part2[len(part2) - mut_coord - 1] == '.':
                    # affects loop  1x1 or 2x2
                    pair = mut_to + seq_compl[len(part2) - mut_coord - 1]
                    if pair in Dic_bp:  # restores stacking
                        seq[mut_coord] = mut_to
                        part1[mut_coord] = '('
                        part2[len(part2) - mut_coord - 1] = ')'
                    elif part1[mut_coord - 1] == '-':
                        # affects loop1x2
                        pair = mut_to + seq_compl[len(part2) - mut_coord]
                        if pair in Dic_bp:  # restores stacking
                            seq[mut_coord] = mut_to
                            part1[mut_coord] = '('
                            part2[len(part2) - mut_coord] = ')'
                        else:  # doesn't restore stacking
                            seq[mut_coord] = mut_to
                    else:  # doesn't restore stacking
                        seq[mut_coord] = mut_to
                elif part1[mut_coord] == '.' and part2[len(part2) - mut_coord - 1] == '-':
                    if part2[len(part2) - mut_coord - 2] == '.':
                        # affects loop 2x1
                        pair = mut_to + seq_compl[len(part2) - mut_coord - 2]
                        if pair in Dic_bp:  # restores stacking
                            seq[mut_coord] = mut_to
                            part1[mut_coord] = '('
                            part2[len(part2) - mut_coord - 2] = ')'
                        else:  # doesn't restore stacking
                            seq[mut_coord] = mut_to
                    else:  # bulge
                        print('bulge')
                        seq[mut_coord] = mut_to
        elif mut_seq == 'seq_compl':
            if part2[mut_coord] == ')':
                if (seq[len(part1) - mut_coord - 1] == 'T' and (
                        seq_compl[mut_coord] == 'A' and mut_to == 'G' or seq_compl[
                    mut_coord] == 'G' and mut_to == 'A')) or (seq[len(part1) - mut_coord - 1] == 'G' and (
                        seq_compl[mut_coord] == 'C' and mut_to == 'T' or seq_compl[
                    mut_coord] == 'T' and mut_to == 'C')):

                    seq_compl[mut_coord] = mut_to
                else:

                    seq_compl[mut_coord] = mut_to
                    part2[mut_coord] = '.'
                    part1[len(part1) - mut_coord - 1] = '.'
            else:
                if part1[len(part1) - mut_coord - 1] == '.' and part2[mut_coord] == '.':
                    pair = mut_to + seq[len(part1) - mut_coord - 1]
                    if pair in Dic_bp:
                        seq_compl[mut_coord] = mut_to
                        part1[len(part1) - mut_coord - 1] = '('
                        part2[mut_coord] = ')'
                    elif part2[mut_coord + 1] == '-':
                        # affects loop 2x1
                        pair = mut_to + seq[len(part1) - mut_coord - 2]
                        if pair in Dic_bp:  # restores stacking
                            seq_compl[mut_coord] = mut_to
                            part1[len(part1) - mut_coord - 2] = '('
                            part2[mut_coord] = ')'
                        else:  # doesn't restore stacking
                            seq_compl[mut_coord] = mut_to
                    else:  # doesn't restore stacking
                        seq_compl[mut_coord] = mut_to
                elif part1[len(part1) - mut_coord - 1] == '-' and part2[mut_coord] == '.':
                    if part1[len(part1) - mut_coord] == '.':
                        # loop 1x2
                        pair = mut_to + seq[len(part1) - mut_coord]
                        if pair in Dic_bp:  # restores stacking
                            seq_compl[mut_coord] = mut_to
                            part1[len(part1) - mut_coord] = '('
                            part2[mut_coord] = ')'
                        else:  # doesn't restore stacking
                            seq_compl[mut_coord] = mut_to
                    else:
                        # bulge 
                        seq_compl[mut_coord] = mut_to
        # check that mutations did not occur on the very end or the very beginning of the seq/seq_compl
        if part1[len(part1) - 1] == '.':
            seq = seq[:len(part1) - 1]
            part1 = part1[:len(part1) - 1]
        elif part1[0] == '.':
            seq = seq[1:len(part1)]
            part1 = part1[1:len(part1)]
        if part2[len(part2) - 1] == '.':
            seq_compl = seq_compl[:len(part2) - 1]
            part2 = part2[:len(part2) - 1]
        elif part2[0] == '.':
            seq_compl = seq_compl[1:len(part2)]
            part2 = part2[1:len(part2)]
        seq = ''.join(seq)
        seq_compl = ''.join(seq_compl)
        part1 = ''.join(part1)
        part2 = ''.join(part2)
        return (seq, seq_compl, part1, part2)


def FindEnergyFromStructure(seq, seq_compl, part1, part2):
    energy = 0
    seq_length = len(seq)
    seq_compl_length = len(seq_compl)
    i = seq_compl_length - 1
    j = 0
    while i > 0 and j < seq_length - 1:
        if part1[j: j + 2] == '((' and part2[i - 1: i + 1] == '))':  # stem

            energy_add = stacking_matrix[Dic_bp.get(seq_compl[i - 1] + seq[j + 1], 6)][
                Dic_bp.get(seq[j] + seq_compl[i], 6)]
            i -= 1
            j += 1
        elif part1[j: j + 2] == '(.' and part2[i - 1: i + 1] == '))':
            if part1[j + 2] == '(':  # Bulge 1x0

                energy_add = bulge_list[1] + stacking_matrix[Dic_bp.get(seq_compl[i - 1] + seq[j + 2], 6)][
                    Dic_bp.get(seq[j] + seq_compl[i], 6)]
                i -= 1
                j += 2

            elif part1[j + 2] == '.':
                if part1[j + 3] == '(':  # Bulge 2x0

                    energy_add = bulge_list[2] + (
                        TerminalAU if Dic_bp.get(seq_compl[i - 1] + seq[j + 3], 6) > 2 or Dic_bp.get(
                            seq[j] + seq_compl[i], 6) > 2 else 0)
                    i -= 1
                    j += 3
                elif part1[j + 3] == '.':
                    # Bulge 3x0 
                    energy_add = bulge_list[3] + (
                        TerminalAU if Dic_bp.get(seq_compl[i - 1] + seq[j + 4], 6) > 2 or Dic_bp.get(
                            seq[j] + seq_compl[i], 6) > 2 else 0)
                    i -= 1
                    j += 4
        elif part1[j: j + 2] == '((' and part2[i - 1: i + 1] == '.)':
            if part2[i - 2] == ')':  # Bulge 0x1

                energy_add = bulge_list[1] + stacking_matrix[Dic_bp.get(seq_compl[i - 2] + seq[j + 1], 6)][
                    Dic_bp.get(seq[j] + seq_compl[i], 6)]
                i -= 2
                j += 1
            elif part2[i - 2] == '.':
                if part2[i - 3] == ')':  # Bulge 0x2

                    energy_add = bulge_list[2] + (
                        TerminalAU if Dic_bp.get(seq_compl[i - 3] + seq[j + 1], 6) > 2 or Dic_bp.get(
                            seq[j] + seq_compl[i], 6) > 2 else 0)
                    i -= 3
                    j += 1

                elif part2[i - 3] == '.':
                    # Bulge 0x3 

                    energy_add = bulge_list[3] + (
                        TerminalAU if Dic_bp.get(seq_compl[i - 4] + seq[j + 1], 6) > 2 or Dic_bp.get(
                            seq[j] + seq_compl[i], 6) > 2 else 0)
                    i -= 4
                    j += 1

        elif part1[j: j + 2] == '(.' and part2[i - 1: i + 1] == '.)':
            if part1[j + 2] == '(':
                if part2[i - 2] == ')':  # Loop 1x1

                    energy_add = \
                    intl11_matrix[Dic_bp.get(seq_compl[i - 2] + seq[j + 2], 7)][Dic_bp.get(seq[j] + seq_compl[i], 7)][
                        Dic_nt.get(seq_compl[i - 1], 5)][Dic_nt.get(seq[j + 1], 5)]
                    i -= 2
                    j += 2

                elif part2[i - 2] == '.':
                    if part2[i - 3] == ')':  # Loop 1x2

                        energy_add = intl12_matrix[Dic_bp.get(seq[j] + seq_compl[i], 7)][
                            Dic_bp.get(seq_compl[i - 3] + seq[j + 2], 7)][Dic_nt.get(seq[j + 1], 5)][
                            Dic_nt.get(seq_compl[i - 2], 5)][Dic_nt.get(seq_compl[i - 1], 5)]
                        i -= 3
                        j += 2
                    elif part2[i - 3] == '.':  # Loop 1x3

                        energy_add = internal_loop[4] + ninio + \
                                     intl13_matrix[Dic_bp.get(seq[j] + seq_compl[i], 7)][Dic_nt.get(seq[j + 1], 5)][
                                         Dic_nt.get(seq_compl[i - 1], 5)] + \
                                     intl13_matrix[Dic_bp.get(seq_compl[i - 4] + seq[j + 2], 7)][
                                         Dic_nt.get(seq_compl[i - 3], 5)][Dic_nt.get(seq[j + 1], 5)]
                        i -= 4
                        j += 2
            elif part1[j + 2] == '.':
                if part2[i - 2] == ')':
                    if part1[j + 3] == '(':  # Loop 2x1

                        energy_add = intl12_matrix[Dic_bp.get(seq_compl[i - 2] + seq[j + 3], 7)][
                            Dic_bp.get(seq[j] + seq_compl[i], 7)][Dic_nt.get(seq_compl[i - 1], 5)][
                            Dic_nt.get(seq[j + 1], 5)][Dic_nt.get(seq[j + 2], 5)]
                        i -= 2
                        j += 3
                    elif part1[j + 3] == '.':  # Loop 3x1
                        energy_add = internal_loop[4] + ninio + \
                                     intl13_matrix[Dic_bp.get(seq[j] + seq_compl[i], 7)][Dic_nt.get(seq[j + 1], 5)][
                                         Dic_nt.get(seq_compl[i - 1], 5)] + \
                                     intl13_matrix[Dic_bp.get(seq_compl[i - 2] + seq[j + 4], 7)][
                                         Dic_nt.get(seq_compl[i - 1], 5)][Dic_nt.get(seq[j + 3], 5)]
                        i -= 2
                        j += 4
                elif part2[i - 2] == '.':
                    if part2[i - 3] == ')':
                        if part1[j + 3] == '(':  # Loop 2x2
                            energy_add = intl22_matrix[Dic_bp.get(seq_compl[i - 3] + seq[j + 3], 7)][
                                Dic_bp.get(seq[j] + seq_compl[i], 7)][Dic_nt.get(seq_compl[i - 2], 5)][
                                Dic_nt.get(seq_compl[i - 1], 5)][Dic_nt.get(seq[j + 1], 5)][Dic_nt.get(seq[j + 2], 5)]
                            i -= 3
                            j += 3
                        elif part1[j + 3] == '.':  # Loop 3x2
                            energy_add = internal_loop[5] + ninio + \
                                         intl23_matrix[Dic_bp.get(seq[j] + seq_compl[i], 7)][Dic_nt.get(seq[j + 1], 5)][
                                             Dic_nt.get(seq_compl[i - 1], 5)] + \
                                         intl23_matrix[Dic_bp.get(seq_compl[i - 3] + seq[j + 4], 7)][
                                             Dic_nt.get(seq_compl[i - 2], 5)][Dic_nt.get(seq[j + 3], 5)]
                            i -= 3
                            j += 4
                    elif part2[i - 3] == '.':
                        if part1[j + 3] == '(':  # Loop2x3
                            energy_add = internal_loop[5] + ninio + \
                                         intl23_matrix[Dic_bp.get(seq[j] + seq_compl[i], 7)][Dic_nt.get(seq[j + 1], 5)][
                                             Dic_nt.get(seq_compl[i - 1], 5)] + \
                                         intl23_matrix[Dic_bp.get(seq_compl[i - 4] + seq[j + 3], 7)][
                                             Dic_nt.get(seq_compl[i - 3], 5)][Dic_nt.get(seq[j + 2], 5)]
                            i -= 4
                            j += 3

                        elif part1[j + 3] == '.':  # Loop 3x3
                            energy_add = internal_loop[6] + ninio + \
                                         intl33_matrix[Dic_bp.get(seq[j] + seq_compl[i], 7)][Dic_nt.get(seq[j + 1], 5)][
                                             Dic_nt.get(seq_compl[i - 1], 5)] + \
                                         intl33_matrix[Dic_bp.get(seq_compl[i - 4] + seq[j + 4], 7)][
                                             Dic_nt.get(seq_compl[i - 3], 5)][Dic_nt.get(seq[j + 3], 5)]

                            i -= 4
                            j += 4

        energy += energy_add
    energy = energy / 100
    structure = part1 + part2
    return (energy, structure, seq, seq_compl)