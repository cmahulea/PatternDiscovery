
import os
import csv
from datetime import datetime
import json
import unicodedata
from collections import defaultdict
from PTPN import ProbabilisticTimePetriNet
from pattern_discovery import *

# Configuración
out_folder = "../examples/lca"
pnml_file = out_folder+"/base_model.pnml"
wilson_file = out_folder+"/wilson_intervals.json"
os.makedirs(out_folder, exist_ok=True)

csv_file = "database/lca/lca_eng.csv"  # Cambia por tu fichero real
patient_col = 0
state_col = 3
date_col = 4
time_col = 5

minimal_support = 1

def normalize(s):
    s = s.replace("/", "").replace(" ", "_")
    return ''.join(c for c in unicodedata.normalize('NFD', s) if unicodedata.category(c) != 'Mn')

# 3. Calcular intervalos de Wilson
def wilson_confidence_intervals(sup_count ,hypotetical_n):
    intervals = {}  # not the same as in pattern_discovery
    z = 1.96
    n = hypotetical_n
    for support, patterns in sup_count.items():
        for pattern in patterns:
            # print(f"Calculating Wilson interval for pattern {pattern} with support {support} and n={n}")
            k = support
            pi = k/n
            pi_tilde = (pi + z**2/(2*n))/(1 + z**2/n)

            sqrt_arg = max(0, pi*(1-pi)/n + z**2/(4*n**2)) # Avoid negative square root
            half = z/(1 + z**2/n)*(sqrt_arg)**0.5

            ci95 = [pi_tilde - half, pi_tilde + half]
            
            from_state, to_state = pattern
            if from_state not in intervals:
                intervals[from_state] = {}
            intervals[from_state][to_state] = ci95
    return intervals

def main():
    # 1. Read and process the CSV file
    print("Reading and processing the CSV file...")
    database = []
    timed_database = []
    current_seq = []
    with open(csv_file) as f:
        reader = csv.reader(f, delimiter=";")
        reader.__next__()  # skip the header line

        for row in reader:
            if row[patient_col]:
                print(f"Processing patient {row[patient_col]}...")
                current_seq = []
                timed_seq = []
                database.append(current_seq)
                timed_database.append(timed_seq)

            # Timestamp
            if not row[time_col].strip():
                time_str = row[date_col].strip()
                if time_str == "PRIVADA" or time_str == "":
                    continue  # Skip private entries
                time = datetime.strptime(time_str, "%d/%m/%Y")
            else:
                time_str = row[date_col] + ' ' + row[time_col]
                time = datetime.strptime(time_str, "%d/%m/%Y %H:%M:%S")
            # print(f"Time string: -{time_str}-")

            state = normalize(row[state_col])
            current_seq.append(state)
            timed_seq.append((state, time.timestamp()))

    nb_sequence = len(database)
    alphabet = generate_alphabet(database)

    ## Generate PTPN
    print("Generating PTPN model...")
    transitions = {}
    PTPN = ProbabilisticTimePetriNet("base_model")

    for transition_name in alphabet:
        interval = lognormal_event_duration(timed_database, transition_name)
        transition = PTPN.add_transition(transition_name, interval, None)
        transitions[transition_name] = transition
        print(transition_name, interval)

    ### Pattern discovery
    patterns = pattern_discovery(database, alphabet)
    distributions = {}
    for support in range(nb_sequence, minimal_support-1, -1):
        for pattern in patterns[support]:
            for i in range(len(pattern)-1):
                if transitions[pattern[i]] not in distributions:
                    distributions[transitions[pattern[i]]] = {transitions[pattern[i+1]]: support}
                elif transitions[pattern[i+1]] not in distributions[transitions[pattern[i]]]:
                    distributions[transitions[pattern[i]]][transitions[pattern[i+1]]] = support
    place_num = 0
    for transition in distributions:
        outcomes = []
        supports = []
        for post_transition in distributions[transition]:
            for potential_place in post_transition.pre:
                place = potential_place
                break
            else:
                place = PTPN.add_place("P" + str(place_num), 0)  # no marking
                place_num += 1
                PTPN.add_edge(place, post_transition, 1)  # no weight for the moment
                
            outcome = {place: 1}  # no weight for the moment
            outcomes.append(outcome)
            supports.append(distributions[transition][post_transition])
        probabilities = [support/sum(supports) for support in supports]
        PTPN.add_distribution(transition, outcomes, probabilities)

    start_place = PTPN.add_place("start", 1)  # 1 mark in starting place
    end_place = PTPN.add_place("end", 0)  # no marking

    for transition in PTPN.transitions:
        if not transition.pre:
            PTPN.add_edge(start_place, transition, 1)  # no weight for the moment
        if not transition.post:
            outcome = {end_place: 1}
            PTPN.add_distribution(transition, [outcome], [1])  # no weight for the moment, probability 1
    if start_place.post_transitions == []:
        print("Warning: start place has no outgoing transitions, removing it from the model.")
        # it is a loop: not playable without information about initial marking
        PTPN.remove_place(start_place)
    if end_place.pre_transitions == []:
        print("Warning: end place has no incoming transitions, removing it from the model.")
        # it is a loop: not playable without information about initial marking
        PTPN.remove_place(end_place)
        
    restart_transition = PTPN.add_transition("repeat", None, None)  # no duration, no weight
    PTPN.add_edge(end_place, restart_transition, 1)
    PTPN.add_distribution(restart_transition, [{start_place: 1}], [1])


    PTPN.export_pnml(pnml_file)


    print("Generating Wilson confidence intervals...")
    sup_count, _ = init_all_supports(database, alphabet)
    sample_sizes = [3, 10, 30, 100]
    intervals = {}
    for n_prime in sample_sizes:
        intervals[str(n_prime)] = wilson_confidence_intervals(sup_count, n_prime)
    # Guardar como JSON
    with open(wilson_file, "w") as f:
        json.dump(intervals, f, indent=2)
    print(f"Intervalos de Wilson exportados a {wilson_file}")

    print("Listo.")

main()