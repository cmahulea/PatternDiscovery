import os,sys
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../..")

import json
import random
import subprocess
from src.net.PTPN import PTPN
from src.ptpnbound import run_ptpnbound
import matplotlib.pyplot as plt

# Configuración
folder = "examples/lca/"
pnml_file = folder+"base_model.pnml"
pnml_file = folder+"given_example.pnml"
wilson_file = folder+"wilson_intervals.json"
results_file = folder+"mc_results.json"
results_txt = folder+"mc_results.txt"

MC_SAMPLES = 500  # Número de réplicas Monte Carlo

transition_name = "Discharge"  # Relevant transition for throughput

def load_intervals(wilson_file, sample_size):
    with open(wilson_file) as f:
        intervals = json.load(f)
    return intervals[str(sample_size)]

def sample_intervals(taget_places_for_src):
    probabilities = {}
    total_prob = 0.0
    for to_place, bounds in taget_places_for_src.items():
        # Asigna un valor aleatorio dentro del intervalo
        prob = random.uniform(bounds[0], bounds[1])
        probabilities[to_place] = prob
        total_prob += prob
    # Normaliza las probabilidades para que sumen 1
    for to_place in probabilities:
        probabilities[to_place] /= total_prob
    return probabilities

def find_next_transition(ptpn, place_name):
    transitions = []
    for arc in ptpn.get_arcs():
        if arc.get_source().get_name() == place_name:
            transitions.append(arc.get_target().get_name())
            # print(f"Found transition {arc.get_target().get_name()} after place '{place_name}'")
    
    if len(transitions) != 1:
        print(f"Warning: Expected 1 place before transition '{transition_name}', found {len(transitions)}")
        for place in transitions:
            print(f"\tTransition: {place.get_name()}")
    return transitions[0]

def assign_probabilities(ptpn, intervals):
    # Agrupa arcos salientes por origen
    for src_transition, to_places in intervals.items():
        probabilities = sample_intervals(to_places)
        nAssigned = 0
        totalAssigned = 0.0
        for arc in ptpn.get_arcs():
            if arc.get_source().get_name() == src_transition:
                transition_after_place = find_next_transition(ptpn, arc.get_target().get_name())
                if transition_after_place in probabilities.keys():
                    # Asigna la probabilidad al arco
                    arc._Arc__probability = probabilities[transition_after_place]
                    nAssigned += 1
                    totalAssigned += arc.get_prob()
        
        if nAssigned != probabilities.__len__():
            print(f"Warning: {nAssigned} arcs assigned, expected {probabilities.__len__()} for transition {src_transition}")
            input("Press Enter to continue...")

def extract_throughput(dot_file, transition_name):
    with open(dot_file) as f:
        for line in f:
            # Busca la línea del nodo de la transición
            if line.strip().startswith(f'{transition_name} ') and 'xlabel=' in line:
                next_line = next(f)
                # Busca la línea con el throughput
                number_str = next_line.split(',')[1].strip().strip(']')
                return float(number_str)
    return float('nan')


def plot_throughputs(metrics):
    sample_sizes = list(metrics.keys())
    means = [metrics[size]["mean"] for size in sample_sizes]
    variances = [metrics[size]["var"] for size in sample_sizes]

    fig = plt.figure(figsize=(10, 5))
    plt.plot(sample_sizes, means, marker='o', label='Mean Throughput')
    plt.fill_between(sample_sizes,[m - v**0.5 for m, v in zip(means, variances)], [m + v**0.5 for m, v in zip(means, variances)], alpha=0.2)
    plt.title('Monte Carlo Throughput Metrics')
    plt.xlabel('Sample Size')
    plt.ylabel('Throughput')
    plt.xticks(sample_sizes)
    plt.savefig(folder+"througput.png")

def plot_CI(metrics):
    sample_sizes = list(metrics.keys())
    ci_lows = [metrics[size]["CI_low"] for size in sample_sizes]
    ci_highs = [metrics[size]["CI_high"] for size in sample_sizes]

    fig = plt.figure(figsize=(10, 5))
    plt.fill_between(sample_sizes, ci_lows, ci_highs, color='blue', alpha=0.5, label='95% Confidence Interval')
    plt.title('Monte Carlo Throughput Confidence Intervals')
    plt.xlabel('Sample Size')
    plt.ylabel('Throughput')
    plt.xticks(sample_sizes)
    plt.legend()
    plt.savefig(folder+"confidence_intervals.png")
def main():
    txt = "Monte Carlo Efficiency Bounds Results\n=====================================\n"
    metrics = {}
    for sample_size in [3, 10, 30, 100]:
        metrics[sample_size] = {}
        print(f"Running Monte Carlo simulation with sample size: {sample_size}")
            
        intervals = load_intervals(wilson_file, sample_size)
        throughputs = []
        for b in range(MC_SAMPLES):
            ptpn = PTPN("net_"+str(b))

            ptpn.import_pnml(pnml_file)
            print("PTPN", b, "loaded")
            # ptpn.print_net()
            # input("Press Enter to continue...")

            assign_probabilities(ptpn, intervals)
            temp_pnml = "tmp"
            ptpn.export_pnml(temp_pnml + ".pnml")
            # Llama a ptpnbound.py
            
            run_ptpnbound(temp_pnml, transition_name, None, None, ["out_tmp", "dot"], False)

            throughput = extract_throughput("out_tmp.dot", transition_name)
            throughputs.append(throughput)

            print(f"Replica {b+1}/{MC_SAMPLES}: throughput={throughput}")

            # Clear tmp files
            os.remove(temp_pnml + ".pnml")
            os.remove("out_tmp.dot")
            os.remove("out_tmp.dot.pdf")

        # Compute metrics
        sorted_throughputs = sorted(throughputs)
        mean_throughput = sum(throughputs) / len(throughputs)
        metrics[sample_size]["mean"] = mean_throughput
        metrics[sample_size]["var"] = sum((x - mean_throughput) ** 2 for x in throughputs) / len(throughputs)
        metrics[sample_size]["CI_low"] = sorted_throughputs[round(0.025 * MC_SAMPLES)]
        metrics[sample_size]["CI_high"] = sorted_throughputs[min(round(0.975 * MC_SAMPLES), len(sorted_throughputs) - 1)]
        metrics[sample_size]["CI_width"]  = metrics[sample_size]["CI_high"] - metrics[sample_size]["CI_low"]

        print("Monte Carlo Throughput Metrics: ")
        txt += f"Sample Size: {sample_size}\n"
        for name, value in metrics[sample_size].items():
            print(f"{name}: {value}")
            txt += f"\t{name}: {value}\n"
    
    # Plot metrics 
    plot_throughputs(metrics)
    plot_CI(metrics)

    # Save metrics
    with open(results_file, 'w') as f:
        json.dump(metrics, f, indent=4)
    with open(results_txt, 'w') as f:
        f.write(txt)

if __name__ == "__main__":
    main()