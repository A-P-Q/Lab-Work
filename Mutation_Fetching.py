import requests
import argparse
import sys
import csv
import os
import pandas as pd
from io import StringIO

#add prediction (Polyphen) for those that doesn't have clisigtypes 

#read the mutation table we made
def read_mut_table(table_path):
    return

def query_uniprot(query):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={query}+AND+reviewed:true&format=json&fields=accession,id&size=500"
    response = requests.get(url)
    if not response.ok:
        print(f"Error fetching data: {response.status_code}")
        return [], []

    data = response.json()
    results = data.get('results', [])
    filtered_accessions = []
    filtered_entry_names = []

    if not results:
        print(f"No results found for the query: {query}")
        return [], []

    # Check if the query is an accession or entry name
    is_accession = results[0].get('primaryAccession', '').lower() == query.lower()
    is_entry_name = results[0].get('uniProtkbId', '').lower() == query.lower()

    if not is_accession and not is_entry_name:
        # If the query is not an accession or entry name, return all results without filtering
        for result in results:
            accession = result.get('primaryAccession', '')
            entry_name = result.get('uniProtkbId', '')
            filtered_accessions.append(accession)
            filtered_entry_names.append(entry_name)

        print(f"All results retrieved without filtering.")
        return filtered_accessions, filtered_entry_names

    # Otherwise, proceed with filtering
    for result in results:
        accession = result.get('primaryAccession', '')
        entry_name = result.get('uniProtkbId', '')

        if query.lower() in accession.lower():
            filtered_accessions.append(accession)
            filtered_entry_names.append(entry_name)
        if query.lower() in entry_name.lower() and entry_name not in filtered_entry_names:
            filtered_accessions.append(accession)
            filtered_entry_names.append(entry_name)

    if filtered_accessions:
        print(f"Filtered Accessions: {', '.join(filtered_accessions)}")
    else:
        print(f"No accessions found containing the query string.")

    if filtered_entry_names:
        print(f"Filtered Entry Names: {', '.join(filtered_entry_names)}")
    else:
        print(f"No entry names found containing the query string.")

    return filtered_accessions, filtered_entry_names


# Function to fetch protein variations from UniProt API
def get_protein_variations(accession):
    base_url = 'https://www.ebi.ac.uk/proteins/api/variation'
    full_url = f'{base_url}/{accession}'
    response = requests.get(full_url)
    if response.status_code == 200:
        return response.json()
    else:
        print(f'Failed to retrieve variation data for accession {accession}. HTTP Status Code: {response.status_code}')
        return None

# Function to fetch and save clinically significant missense variants as CSV
def fetch_and_save_variants(accession, output_dir):
    protein_variations = get_protein_variations(accession)
    variant_list = []

    if protein_variations:
        features = protein_variations.get('features', [])
        if not features:
            print(f"No features found for accession {accession}. Skipping CSV generation.")
            return False

        for variant in features:
            begin = variant.get('begin', '')
            end = variant.get('end', '')
            wildType = variant.get('wildType', '')
            mutatedType = variant.get('alternativeSequence', '')
            # Check if clinical significances exist, if not, set to 'N/A'
            clinical_significance_dicts = variant.get('clinicalSignificances', [])
            types = ','.join([d.get('type', '') for d in clinical_significance_dicts]) if clinical_significance_dicts else 'N/A'
            # Include the variant regardless of clinical significance
            variant_entry = [f"{wildType}{begin}{mutatedType}", begin, end, wildType, mutatedType, types, 'proteinsAPI', 0.0]
            variant_list.append(variant_entry)

        if variant_list:
            csv_filename = os.path.join(output_dir, f'{accession}_all_variants.csv')
            with open(csv_filename, 'w', newline='') as f:
                csv_writer = csv.writer(f)
                csv_writer.writerow(['Variant', 'Begin', 'End', 'Wild Type', 'Mutated Type', 'Clinical Significance Types', 'Origin', 'Coefficient'])
                csv_writer.writerows(variant_list)

            print(f"Variant fetching completed for accession {accession}.")
            return True
        else:
            print(f"No variants found for accession {accession}. Skipping CSV generation.")
    else:
        print(f"Failed to retrieve protein variations for accession {accession}.")
    return False
def fetch_and_save_variants(accession, output_dir):
    protein_variations = get_protein_variations(accession)
    variant_list = []

    if protein_variations:
        features = protein_variations.get('features', [])
        if not features:
            print(f"No features found for accession {accession}. Skipping CSV generation.")
            return False

        for variant in features:
            begin = variant.get('begin', '')
            end = variant.get('end', '')
            wildType = variant.get('wildType', '')
            mutatedType = variant.get('alternativeSequence', '')
            # Check if clinical significances exist, if not, set to 'N/A'
            clinical_significance_dicts = variant.get('clinicalSignificances', [])
            types = ','.join([d.get('type', '') for d in clinical_significance_dicts]) if clinical_significance_dicts else 'N/A'
            # Include the variant regardless of clinical significance
            variant_entry = [f"{wildType}{begin}{mutatedType}", begin, end, wildType, mutatedType, types, 'proteinsAPI', 0.0]
            variant_list.append(variant_entry)

        if variant_list:
            csv_filename = os.path.join(output_dir, f'{accession}_all_variants.csv')
            with open(csv_filename, 'w', newline='') as f:
                csv_writer = csv.writer(f)
                csv_writer.writerow(['Variant', 'Begin', 'End', 'Wild Type', 'Mutated Type', 'Clinical Significance Types', 'Origin', 'Coefficient'])
                csv_writer.writerows(variant_list)

            print(f"Variant fetching completed for accession {accession}.")
            return True
        else:
            print(f"No variants found for accession {accession}. Skipping CSV generation.")
    else:
        print(f"Failed to retrieve protein variations for accession {accession}.")
    return False

# Function to fetch protein mutations from UniProt API
def get_protein_mutations(accession):
    base_url = 'https://www.ebi.ac.uk/proteins/api/mutagenesis'
    full_url = f'{base_url}/{accession}'
    response = requests.get(full_url)
    if response.status_code == 200:
        return response.json()
    else:
        print(f'Failed to retrieve mutation data for accession {accession}. HTTP Status Code: {response.status_code}')
        return None

# Function to fetch and save mutations as CSV
def fetch_and_save_mutations(accession, output_dir):
    protein_mutations = get_protein_mutations(accession)
    mutation_list = []

    if protein_mutations:
        seq = protein_mutations.get('sequence', '')
        features = protein_mutations.get('features', [])

        if not features:
            print(f"No mutations found for accession {accession}. Skipping CSV generation.")
            return False

        for variant in features:
            begin = int(variant['begin'])
            end = int(variant['end'])
            alternativeSequence = variant.get('alternativeSequence', '')
            description = variant.get('description', '')
            evidences = variant.get('evidences', [])
            origin = 'mutagenesis'
            coefficient = 0.0

            if 0 <= begin - 1 < len(seq):
                wild_type = seq[begin - 1]
            else:
                wild_type = ''

            variant_list = [f"{wild_type}{begin}{alternativeSequence}", begin, end, wild_type, alternativeSequence, description, evidences, origin, coefficient]
            mutation_list.append(variant_list)

        csv_filename = os.path.join(output_dir, f'{accession}_mutagenesis.csv')
        with open(csv_filename, 'w', newline='') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(['Variant', 'Begin', 'End', 'Wild Type', 'Mutated Type', 'Clinical Significance Types', 'evidences', 'Origin', 'Coefficient'])
            csv_writer.writerows(mutation_list)
        print(f"Mutation fetching completed.")
        return True
    else:
        print(f"No mutation data found for accession {accession}.")
    return False

# Function to extract variant information from a string
def extract_info(variant):
    digit_start = next(i for i, c in enumerate(variant) if c.isdigit())
    digit_end = next(i for i, c in enumerate(variant[digit_start:], start=digit_start) if not c.isdigit())

    position = int(variant[digit_start:digit_end])
    wild_type = variant[:digit_start]
    mutated_type = variant[digit_end:]

    return position, position, wild_type, mutated_type

def calculate_average_coefficient(df):
    # Group by the 'Begin' column and calculate the mean of the 'Coefficient' column
    df['Coefficient'] = pd.to_numeric(df['Coefficient'], errors='coerce')
    average_coefficients = df.groupby('Begin', as_index=False)['Coefficient'].mean()
    average_coefficients.rename(columns={'Coefficient': 'Average_Coefficient'}, inplace=True)

    return average_coefficients

# Function to reformat AlphaMissense data and save as TSV
def reformatAlphaMissense(accession_file, output_dir):
    s_GitHub = os.environ.get('HOME')
    s_MutDataDir = os.path.join(s_GitHub, "Documents/GitHub/Cx40/data/alphaMissense")

    fileNm = os.path.join(s_MutDataDir, f"{accession_file}.tsv")
    outNm = os.path.join(output_dir, f"{accession_file}_alphaMissense.csv")

    if not os.path.exists(fileNm):
        print(f"File {fileNm} does not exist. Skipping AlphaMissense processing.")
        return False

    with open(fileNm, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        with open(outNm, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter=',')
            writer.writerow(['Variant', 'Begin', 'End', 'Wild Type', 'Mutated Type', 'Clinical Significance Types', 'Origin', 'Coefficient'])

            for row in reader:
                variant = row[1]
                begin, end, wild_type, mutated_type = extract_info(variant)
                clinical_significance = row[3]
                coefficient = row[2]
                origin = "AlphaMissense"

                writer.writerow([variant, begin, end, wild_type, mutated_type, clinical_significance, origin, coefficient])

            print(f"Formatting completed. Data saved in {outNm}")

            alpha_missense_df = pd.read_csv(outNm)
            unique_begins_count = alpha_missense_df['Begin'].nunique()
            print(f"Unique 'Begin' values after reformatting: {unique_begins_count}")
            average_coefficients_df = calculate_average_coefficient(alpha_missense_df)
            average_coefficients_file = os.path.join(output_dir, f"{accession_file}_average_coefficients.csv")
            average_coefficients_df.to_csv(average_coefficients_file, index=False)
            print(f"Average coefficients saved in: {average_coefficients_file}")

            return True
        
def replace_coefficients(merged_df, alpha_missense_file):
    # Check if the alpha_missense_file exists before proceeding
    if os.path.exists(alpha_missense_file):
        alpha_missense_df = pd.read_csv(alpha_missense_file)
        alpha_coeff_dict = alpha_missense_df.set_index(['Begin', 'End', 'Wild Type', 'Mutated Type'])['Coefficient'].to_dict()

        for index, row in merged_df.iterrows():
            key = (row['Begin'], row['End'], row['Wild Type'], row['Mutated Type'])
            if key in alpha_coeff_dict:
                merged_df.at[index, 'Coefficient'] = alpha_coeff_dict[key]
    else:
        print(f"File {alpha_missense_file} does not exist. Skipping coefficient replacement.")

    return merged_df
        
# Function to convert TSV file to CSV file
def convert_tsv_to_csv(file_path, output_dir):
    csv_file_path = os.path.join(output_dir, os.path.basename(file_path).replace('.tsv', '.csv'))

    with open(file_path, 'r', newline='') as tsv_file:
        with open(csv_file_path, 'w', newline='') as csv_file:
            tsv_reader = csv.reader(tsv_file, delimiter='\t')
            csv_writer = csv.writer(csv_file)

            for row in tsv_reader:
                csv_writer.writerow(row)

    print(f"Converted TSV to CSV. Data saved in {csv_file_path}")

def combine_csv_files(accession, output_dir, files_generated):
    # Only generate the combined table if at least one file was created
    if not files_generated:
        print(f"No files generated for accession {accession}. Skipping combined table generation.")
        return

    file_paths = [
        os.path.join(output_dir, f'{accession}_all_variants.csv'),
        os.path.join(output_dir, f'{accession}_mutagenesis.csv'),
        os.path.join(output_dir, f'{accession}_alphaMissense.csv')
    ]
    
    dataframes = []
    for file_path in file_paths:
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            dataframes.append(df)

    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.insert(0, 'Accession', accession)

        combined_file = os.path.join(output_dir, f'{accession}_combined.csv')
        combined_df.to_csv(combined_file, index=False)
        print(f"Combined CSV file created at: {combined_file}")
        return combined_file # so that we can use it in load_and_merge_tables
    else:
        print(f"No valid files to combine for accession {accession}.")

def load_and_merge_tables(master_file_path, uniProtId):
    master_df = pd.read_csv(master_file_path)
    master_df.columns = master_df.columns.str.lower()  # Standardize column names to lowercase
    uniProtId = uniProtId.lower()

    s_GitHub = os.environ.get('HOME')
    s_MutDataDir = os.path.join(s_GitHub, "Documents/GPCRDB/")
    
    # Adjusted to use UniProtKB ID (uniProtId) instead of accession
    BW_fileNm = os.path.join(s_MutDataDir, f"cxbw_{uniProtId}.csv")
    
    if os.path.exists(BW_fileNm):
        BW_df = pd.read_csv(BW_fileNm)
        merged_df = master_df.merge(BW_df[['sequence_number', 'ballesteros', 'alignment_index']], 
                                    left_on='begin', 
                                    right_on='sequence_number', 
                                    how='left')
        merged_df.rename(columns={'ballesteros': 'CxBW'}, inplace=True)
    else:
        print(f"No cxbw file found for {uniProtId}. Adding placeholder columns.")
        merged_df = master_df
        for column in ['sequence_number', 'ballesteros', 'alignment_index']:
            merged_df[column] = pd.NA

    merged_df.columns = merged_df.columns.str.title()  # Ensure consistent column naming
    return merged_df


def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def parse_arguments():
    default_output_dir = os.path.join(os.environ.get('HOME', os.getcwd()), "Documents/CxVariations/")
    ensure_directory_exists(default_output_dir)

    parser = argparse.ArgumentParser(description="Fetch and save protein variants based on input query.")
    parser.add_argument("--input_query", help="Value to search for (e.g., gene name, entry name, protein name, or direct accession).")
    parser.add_argument("--output_dir", help="Directory to save output files", default=default_output_dir)
    args = parser.parse_args()

    return args.input_query, args.output_dir

if __name__ == "__main__":
    input_query, output_dir = parse_arguments()
    accessions, entryNames = query_uniprot(input_query)

    for accession, entryName in zip(accessions, entryNames):
        print(f"Processing UniProt accession: {accession} with entry name: {entryName}")
        files_generated = False
        if fetch_and_save_variants(accession, output_dir):
            files_generated = True
        if fetch_and_save_mutations(accession, output_dir):
            files_generated = True
        if reformatAlphaMissense(accession, output_dir):
            files_generated = True

        if files_generated:
            combined_file = combine_csv_files(accession, output_dir, files_generated)
            if combined_file:
                merged_df = load_and_merge_tables(combined_file, entryName) 
                alpha_missense_path = os.path.join(output_dir, f'{accession}_alphaMissense.csv')
                updated_combined_df = replace_coefficients(merged_df, alpha_missense_path)
                merged_file = os.path.join(output_dir, f'{accession}_combined.csv')
                updated_combined_df.to_csv(merged_file, index=False)
                print(f"Updated merged file with actual coefficients saved at: {merged_file}")
        else:
            print("Could not resolve the input query to this UniProt accession number.")  