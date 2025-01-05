import argparse
import numpy as np
import pandas as pd

def F1_score(a, b):
    
    return a * b / (a + b + 1e-12)

def product(a, b):
    
    return np.sqrt(a * b)

def weighted_jaccard_distance(feature_matrix):

    distance_matrix = np.zeros((feature_matrix.shape[0], feature_matrix.shape[0]))
    for i in range(feature_matrix.shape[0]):
        for j in range(i+1, feature_matrix.shape[0]):
            q = np.vstack((feature_matrix[i], feature_matrix[j]))
            value = 1 - np.sum(np.min(q, axis=0)) / (np.sum(np.max(q, axis=0)) + 1e-12)
            distance_matrix[i][j] = value
            distance_matrix[j][i] = value
    
    return distance_matrix

def euclidean_distance(feature_matrix):
    distance_matrix = np.zeros((feature_matrix.shape[0], feature_matrix.shape[0]))
    for i in range(feature_matrix.shape[0]):
        for j in range(i+1, feature_matrix.shape[0]):
            value = np.sqrt(np.sum((feature_matrix[i] - feature_matrix[j])**2))
            distance_matrix[i][j] = value
            distance_matrix[j][i] = value
    min_val = np.min(distance_matrix)
    max_val = np.max(distance_matrix)
    distance_matrix = (distance_matrix - min_val) / (max_val - min_val)

    return distance_matrix

def correlation_distance(feature_matrix):
    distance_matrix = np.zeros((feature_matrix.shape[0], feature_matrix.shape[0]))
    for i in range(feature_matrix.shape[0]):
        for j in range(i+1, feature_matrix.shape[0]):
            correlation = np.corrcoef(feature_matrix[i], feature_matrix[j])[0, 1]
            value = 1 - correlation
            distance_matrix[i][j] = value
            distance_matrix[j][i] = value
    min_val = np.min(distance_matrix)
    max_val = np.max(distance_matrix)
    distance_matrix = (distance_matrix - min_val) / (max_val - min_val)
    
    return distance_matrix

def manhattan_distance(feature_matrix):
    distance_matrix = np.zeros((feature_matrix.shape[0], feature_matrix.shape[0]))
    for i in range(feature_matrix.shape[0]):
        for j in range(i+1, feature_matrix.shape[0]):
            value = np.sum(np.abs(feature_matrix[i] - feature_matrix[j]))
            distance_matrix[i][j] = value
            distance_matrix[j][i] = value
    min_val = np.min(distance_matrix)
    max_val = np.max(distance_matrix)
    distance_matrix = (distance_matrix - min_val) / (max_val - min_val)
    
    return distance_matrix

def cosine_distance(feature_matrix):
    distance_matrix = np.zeros((feature_matrix.shape[0], feature_matrix.shape[0]))
    for i in range(feature_matrix.shape[0]):
        for j in range(i+1, feature_matrix.shape[0]):
            dot_product = np.dot(feature_matrix[i], feature_matrix[j])
            norm_i = np.linalg.norm(feature_matrix[i])
            norm_j = np.linalg.norm(feature_matrix[j])
            value = 1 - (dot_product / (norm_i * norm_j + 1e-12))
            distance_matrix[i][j] = value
            distance_matrix[j][i] = value
    min_val = np.min(distance_matrix)
    max_val = np.max(distance_matrix)
    distance_matrix = (distance_matrix - min_val) / (max_val - min_val)
    
    return distance_matrix


def generate_interaction_abundance(species_interaction_list, species_abundance_f):
    
    data = pd.read_csv(species_abundance_f, index_col=0, sep="\t")
    species_list = data.index.values.tolist()
    sample_list = data.columns.values.tolist()
    data_matrix = np.array(data)
    data_matrix = np.log(data_matrix + 1)
    data_matrix[data_matrix!=0] = 1
    interaction_abundance_matrix = np.zeros((len(species_interaction_list), len(sample_list)))
    for i in range(interaction_abundance_matrix.shape[0]):
        species_interaction = species_interaction_list[i]
        species_1 = species_interaction.split("&")[0]
        species_2 = species_interaction.split("&")[1]

        if species_1 in species_list and species_2 in species_list:
            interaction_abundance_matrix[i] = product(data_matrix[species_list.index(species_1)], data_matrix[species_list.index(species_2)])

    interaction_abundance_matrix = interaction_abundance_matrix / (np.sum(data_matrix, axis=0) + 1e-7) ** 2 * 100 

    return interaction_abundance_matrix, sample_list


def cal_interaction_redundancy(interaction_matrix, interaction_abundance_matrix, distance_measure):
          
    distance_matrix = eval(distance_measure)(interaction_matrix)
    
    interaction_redundancy_value = np.zeros((interaction_abundance_matrix.shape[1]))
    td_value = np.zeros((interaction_abundance_matrix.shape[1]))
    for i in range(interaction_abundance_matrix.shape[0]):
        for j in range(interaction_abundance_matrix.shape[0]):
            if i != j:
                interaction_redundancy_value += (interaction_abundance_matrix[i] * interaction_abundance_matrix[j] * (1 - distance_matrix[i, j]))
                td_value += (interaction_abundance_matrix[i] * interaction_abundance_matrix[j])

    return interaction_redundancy_value, td_value


def main():

    parser = argparse.ArgumentParser()
    # required parameters
    parser.add_argument("--ICN_ref", default=None, type=str, required=True, help="The reference interaction content network file (OTU interaction x KEGG interaction).")
    parser.add_argument("--abundance_file", default=None, type=str, required=True, help="The taxon abundance data (OTU x sample).")

    # other paramters
    parser.add_argument("--output_file", default="IR_output.txt", type=str, help="The output file (default: IR_output.txt).")
    parser.add_argument("--distance_measure", default="weighted_jaccard_distance", type=str, help="The distance measures used for calculating interaction redundancy, including weighted_jaccard_distance, euclidean_distance, correlation_distance, and manhattan_distance (default: weighted_jaccard_distance).")

    args = parser.parse_args()
    interaction_matrix = pd.read_csv(args.ICN_ref, index_col=0, sep=",")
    species_interaction_list = interaction_matrix.index.values.tolist()
    interaction_matrix = np.array(interaction_matrix)

    interaction_abundance_matrix, sample_list = generate_interaction_abundance(species_interaction_list, args.abundance_file)

    interaction_redundancy_value, td_value = cal_interaction_redundancy(interaction_matrix, interaction_abundance_matrix, args.distance_measure)
    
    f_out = open(args.output_file, "w")
    f_out.write("Sample" + "\t" + "Interaction Redundancy" + "\t" + "Interaction Diversity" + "\n")
    for i in range(len(sample_list)):
        f_out.write(sample_list[i] + "\t" + str(interaction_redundancy_value[i]) + "\t" + str(td_value[i]) + "\n")

    f_out.close()
    
    


if __name__ == "__main__":
    main()

