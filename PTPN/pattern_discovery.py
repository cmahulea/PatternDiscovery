######################################
#           Pattern discovery        #
#     Author :  Le Moigne            #
# theo.le-moigne@ens-paris-saclay.fr #
#              30/11/2023            #
#         modified 04/06/2024        #
######################################

"""
Implementation of Continuous Pattern Discovery from "Automated Generation of Models of Activities
of Daily Living" by Jérémie Saives and Gregory Faraut (https://hal.science/hal-00999505) and
"Automatic Discovery and Deviation Detection of Temporal Inhabitant Habits" by Kevin Viard
(unpublished).

Use a database of n sequences of event (max length: m) from an alphabet sigma (max length: l) and a minimal support supp_min
"""

import copy
from datetime import timedelta
import math


def eval_dates(timed_database, pattern):
    """
    Evaluate dates of a pattern
    timed_database: list of n sequences (lists containing events from alphabet sigma and time)
    pattern: a list of event to find in sequences
    return ld: the dates of the given pattern (start, end for each occurrence)
    """
    k = len(pattern)
    ld = []
    for sequence in timed_database:
        for n in range(len(sequence)-k+1):
            if pattern == map(lambda x: x[0], sequence[n:n+k]):
                ld.append([sequence[n][1], sequence[n+k][1]])
    return ld


def average_pattern_duration(timed_database, pattern):
    """
    Evaluate average duration of a pattern
    timed_database: list of n sequences (lists containing events from alphabet sigma and time)
    pattern: a list of event to find in sequences
    return: float
    """
    ld = eval_dates(timed_database, pattern)
    durations = map(lambda x: x[1]-x[0], ld)
    return sum(durations)/len(durations)


def average_event_duration(timed_database, event):
    """
    Evaluate average duration of an event
    timed_database: list of n sequences (lists containing events from alphabet sigma and time)
    return: float
    """
    durations = []
    for id_sequence in range(len(timed_database)):
        sequence = timed_database[id_sequence]
        for id_event in range(len(sequence)):
            if sequence[id_event][0] == event:
                if id_event != 0:
                    durations.append(sequence[id_event][1] - sequence[id_event-1][1])
                elif id_sequence != 0:
                    durations.append(sequence[id_event][1] - timed_database[id_sequence-1][-1][1])
                else:
                    durations.append(sequence[id_event][1])
    return sum(durations)/len(durations)


def interval_event_duration(timed_database, event, initial_time=None):
    """
    Evaluate the interval duration of an event
    timed_database: list of n sequences (lists containing events from alphabet sigma and time)
    return: interval, list of min and max duration
    """
    durations = []
    for id_sequence in range(len(timed_database)):
        sequence = timed_database[id_sequence]
        for id_event in range(len(sequence)):
            if sequence[id_event][0] == event:
                if id_event != 0:
                    durations.append(sequence[id_event][1] - sequence[id_event-1][1])
                elif id_sequence != 0:  # b careful, only possible if it's folowing sequences
                    durations.append(sequence[id_event][1] - timed_database[id_sequence-1][-1][1])
                elif initial_time is not None:  # if initial time is not defined, it's not possible to compute the duration of the first event
                    durations.append(sequence[id_event][1] - initial_time)
    return [min(durations), max(durations)]


def lognormal_event_duration(timed_database, event, initial_time=None):
    """
    Evaluate log-normal parameters of the duration of an event
    timed_database: list of n sequences (lists containing events from alphabet sigma and time)
    return: dict with "mu" and "sigma" parameters
    """
    following_event = False  # if True, compute duration by difference of following event time, else, the duration is already given
    durations = []
    for id_sequence in range(len(timed_database)):
        sequence = timed_database[id_sequence]
        for id_event in range(len(sequence)):
            if sequence[id_event][0] == event:
                if following_event:
                    if id_event != 0 and sequence[id_event][1] > sequence[id_event-1][1]:
                        durations.append(sequence[id_event][1] - sequence[id_event-1][1])
                    elif id_sequence != 0 and sequence[id_event][1] > timed_database[id_sequence-1][-1][1]:  # be careful, only possible for folowing sequences
                        durations.append(sequence[id_event][1] - timed_database[id_sequence-1][-1][1])
                    elif initial_time is not None:  # if initial time is not defined, it's not possible to compute the duration of the first event
                        durations.append(sequence[id_event][1] - initial_time)
                else:
                    durations.append(sequence[id_event][1])
    durations = [d.total_seconds() if type(d) is timedelta else d for d in durations]
    mu = sum([math.log(d) for d in durations])/len(durations)
    sigma = math.sqrt(sum([(math.log(d) - mu)**2 for d in durations])/(len(durations)-1))  # Bessel's correction
    return {"mu": mu, "sigma": sigma}


def normal_event_duration(timed_database, event, initial_time=None):
    """
    Evaluate normal parameters of the duration of an event
    timed_database: list of n sequences (lists containing events from alphabet sigma and time)
    return: dict with "mu" and "sigma" parameters
    """
    following_event = False  # if True, compute duration by difference of following events time, else, the duration is already given
    durations = []
    for id_sequence in range(len(timed_database)):
        sequence = timed_database[id_sequence]
        for id_event in range(len(sequence)):
            if sequence[id_event][0] == event:
                if following_event:
                    if id_event != 0 and sequence[id_event][1] > sequence[id_event-1][1]:
                        durations.append(sequence[id_event][1] - sequence[id_event-1][1])
                    elif id_sequence != 0 and sequence[id_event][1] > timed_database[id_sequence-1][-1][1]:  # be careful, only possible for folowing sequences
                        durations.append(sequence[id_event][1] - timed_database[id_sequence-1][-1][1])
                    elif initial_time is not None:  # if initial time is not defined, it's not possible to compute the duration of the first event
                        durations.append(sequence[id_event][1] - initial_time)
                else:
                    durations.append(sequence[id_event][1])
    durations = [d.total_seconds() if type(d) is timedelta else d for d in durations]
    mu = sum(durations)/len(durations)
    sigma = math.sqrt(sum([(d - mu)**2 for d in durations])/(len(durations)-1))  # Bessel's correction
    return {"mu": mu, "sigma": sigma}


def exponential_event_duration(timed_database, event, initial_time=None):
    """
    Evaluate exopnential parameter of the duration of an event
    timed_database: list of n sequences (lists containing events from alphabet sigma and time)
    return: dict with "lambd" parameter
    """
    following_event = False  # if True, compute duration by difference of following events time, else, the duration is already given
    durations = []
    for id_sequence in range(len(timed_database)):
        sequence = timed_database[id_sequence]
        for id_event in range(len(sequence)):
            if sequence[id_event][0] == event:
                if following_event:
                    if id_event != 0 and sequence[id_event][1] > sequence[id_event-1][1]:
                        durations.append(sequence[id_event][1] - sequence[id_event-1][1])
                    elif id_sequence != 0 and sequence[id_event][1] > timed_database[id_sequence-1][-1][1]:  # be careful, only possibles for folowing sequences
                        durations.append(sequence[id_event][1] - timed_database[id_sequence-1][-1][1])
                    elif initial_time is not None:  # if initial time is not defined, it's not possible to compute the duration of the first event
                        durations.append(sequence[id_event][1] - initial_time)
                else:
                    durations.append(sequence[id_event][1])
    durations = [d.total_seconds() if type(d) is timedelta else d for d in durations]
    lambd = len(durations)/sum(durations)
    return {"lambd": lambd}


def gamma_event_duration(timed_database, event, initial_time=None):
    """
    Evaluate gamme parameters of the duration of an event
    timed_database: list of n sequences (lists containing events from alphabet sigma and time)
    return: dict with "k" and "theta" parameter
    """
    following_event = False  # if True, compute duration by difference of following events time, else, the duration is already given
    durations = []
    for id_sequence in range(len(timed_database)):
        sequence = timed_database[id_sequence]
        for id_event in range(len(sequence)):
            if sequence[id_event][0] == event:
                if following_event:
                    if id_event != 0 and sequence[id_event][1] > sequence[id_event-1][1]:
                        durations.append(sequence[id_event][1] - sequence[id_event-1][1])
                    elif id_sequence != 0 and sequence[id_event][1] > timed_database[id_sequence-1][-1][1]:  # be careful, only possible for folowing sequences
                        durations.append(sequence[id_event][1] - timed_database[id_sequence-1][-1][1])
                    elif initial_time is not None:  # if initial time is not defined, it's not possible to compute the duration of the first event
                        durations.append(sequence[id_event][1] - initial_time)
                else:
                    durations.append(sequence[id_event][1])
    durations = [d.total_seconds() if type(d) is timedelta else d for d in durations]
    mu = sum(durations)/len(durations)
    sigma = math.sqrt(sum([(d - mu)**2 for d in durations])/(len(durations)-1))  # Bessel's correction
    k = (mu/sigma)**2
    theta = mu/k
    gamma_k = 1
    for i in range(1, 1000):
        gamma_k *= (1+1/i)**k/(1+k/i)
    gamma_k *= 1/k
    return {"k": k, "theta": theta, "gamma_k": gamma_k}


def check_inclusion(pattern, sequence):
    """
    Check if a pattern is include in a sequence (or an other pattern)
    pattern: the included or not pattern
    sequence: the sequence containing or not the pattern
    return: boolean
    complexity: O(m²/4)
    """
    # maybe optimisable
    for i in range(len(sequence)+1-len(pattern)):
        if pattern == sequence[i:i+len(pattern)]:
            return True
    return False


def generate_alphabet(database):
    """
    Generate the alphabet Sigma of event in a database
    database: list of n sequences (lists containing events from alphabet sigma)
    return sigma: the alphabet in database
    complexity: O(n*m*l)
    """
    sigma = []
    for sequence in database:
        for event in sequence:
            if event not in sigma:
                sigma.append(event)
    return sigma


def init_all_supports(database, sigma):
    """
    Initialize the list of pattern with pattern of elementary length (length 2) for all supports
    database: list of n sequences (lists containing events from alphabet sigma)
    sigma: an alphabet
    return list_elem: dict of list of patterns of elementary lenght indexed by support
           total: dict of number of sequence with occurence of an event
    complexity: O(n*l²*m)
    """
    patterns = {(e1, e2): 0 for e1 in sigma for e2 in sigma}
    total = {e: 0 for e in sigma}
    for sequence in database:
        patterns_found = {pattern: False for pattern in patterns}
        for i in range(len(sequence) - 1):
            patterns_found[(sequence[i], sequence[i+1])] = True
        for pattern, found in patterns_found.items():
            if found:
                patterns[pattern] += 1
                total[pattern[0]] += 1
    list_disc = {support: [] for support in range(1, len(database)+1)}
    for pattern, support in patterns.items():
        if support > 0:
            list_disc[support].append(list(pattern))
    return list_disc, total


def pattern_discovery(database, sigma=None):
    """
    Pattern Discovery which list pattern **of size 2** for each support found in a sequences
    database: list of n sequences (lists containing events from an alphabet sigma)
    return list_disc: dict of list of patterns of elementary lenght indexed by support
    """
    if not sigma:
        sigma = generate_alphabet(database)
    list_disc, _ = init_all_supports(database, sigma)                      # discovered patterns
    return list_disc


def pattern_discovery_95_wilson_confidence_intervals(database, sigma=None):
    """
    Pattern Discovery which list pattern **of size 2** for each support found in a sequences
    database: list of n sequences (lists containing events from an alphabet sigma)
    return list_disc: list of tuple (pattern, k(support), n, pi, ci95(95% Wilson confidence interval))
    """
    if not sigma:
        sigma = generate_alphabet(database)
    count, total = init_all_supports(database, sigma)                      # discovered patterns
    list_disc = []  # not the same as in pattern_discovery
    z = 1.96
    for support, patterns in count.items():
        for pattern in patterns:
            k = support
            n = total[pattern[0]]
            pi = k/n
            pi_tilde = (pi + z**2/(2*n))/(1 + z**2/n)
            half = z/(1 + z**2/n)*(pi*(1-pi)/n + z**2/(4*n**2))**0.5
            ci95 = [pi_tilde - half, pi_tilde + half]
            list_disc.append([pattern, k, n, pi, ci95])
    return list_disc
