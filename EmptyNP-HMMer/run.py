#!/usr/bin/env python
import os, argparse, subprocess, logging, re, shutil
from alive_progress import alive_bar; import time

def process_fasta_files(directory):
    """
    Processes each FASTA file in the input folder using awk and seqkit.
    """
    for filename in os.listdir(directory):
        # Ensure filename ends with .fasta or .fa but not with _clean.fasta
        if (filename.endswith(".fasta") or filename.endswith(".fa")) and not filename.endswith("_clean.fasta"):
            base_name = os.path.splitext(filename)[0]  # Extract base name without extension
            input_path = os.path.join(directory, filename)
            clean_path = os.path.join(directory, f"{base_name}_clean.txt")
            clean_fasta_path = os.path.join(directory, f"{base_name}_clean.fasta")

            # Skip processing if the cleaned file already exists
            if os.path.exists(clean_path) and os.path.exists(clean_fasta_path):
                logging.info(f"Skipped processing as cleaned file exists: {clean_path} and {clean_fasta_path}")
                continue

            # Run awk command
            logging.info(f"Processing: {filename}")

        # Using subprocess.Popen for commands that involve redirections
            with open(clean_path, "w") as outfile:
                subprocess.run(["awk", "{print $1}", input_path], stdout=outfile, check=True)

            seqkit_command = f"seqkit seq {clean_path} -o {clean_fasta_path}"
            subprocess.run(seqkit_command, shell=True, check=True)
            print(f"Processed: {filename}")

def list_files(directory):
    """
    List all files in a given directory.

    :param directory: The directory path to list files from.
    """
    if directory is None:
        logging.error("No database directory provided.")
        return None
    
    if not os.path.exists(directory):
        logging.error(f"The directory {directory} does not exist.")
        return []

    database_path = os.path.join(directory, 'database')
    database_output = open(database_path, "w")

    for filename in os.listdir(directory):
        if filename.endswith(".hmm"):
            input_path = os.path.join(directory, filename)
            base_name = filename.rsplit('.hmm', 1)[0]  # Remove the .hmm extension
            new_hmm_file = os.path.join(directory,f"{base_name}_new")  # New file name

            # Run hmmconvert command
            logging.info(f"Processing HMM file: {input_path}")
            subprocess.run(['hmmconvert', input_path], stdout=database_output)

    database_output.close()
    subprocess.run(f"hmmpress -f {database_path}", shell=True, check=True)
    return database_path

def run_hmm(directory, cpu_count, database_path, output, E, dome, visualization):
    """
    Processes each cleaned FASTA file in the input folder using HMMER to search against a database.
    """
    output_dir = output
    #if os.path.exists(output_dir):  
    #    shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

    database_alignment = 'Database_alignment'
    for filename in os.listdir(directory):
        if filename.endswith("_clean.fasta"):
            base_name = filename.rsplit('_clean.fasta', 1)[0]
            input_path = os.path.join(directory, filename)
            base_output_dir = os.path.join(output_dir, base_name)
            os.makedirs(base_output_dir, exist_ok=True)  # Create a subdirectory for each file's output

            output_file_path = os.path.join(base_output_dir, f"{base_name}.hmm_results")
            logging.info(f"Running HMMER for {filename}")
            #running_hmm_command = f'hmmsearch --noali --domE {dome} --domtblout, {output_file_path} -E {E} --cpu {str(cpu_count), database_path, input_path}'
            #print(f'Running HMM command: {running_hmm_command}')

            with alive_bar(enrich_print=False) as bar:
                hmm_command = ["hmmsearch", "--noali", "--domE", dome, "--domtblout", output_file_path, "-E", E, "--cpu", str(cpu_count), database_path, input_path]
                # Run hmmsearch command and log the result
                try:
                    print('\n')
                    print(f'Processing HMMsearch using {dome} (E-value) to report domains and {E} (E-value) to report sequences')
                    result = subprocess.run(hmm_command, check=True, capture_output=True, text=True)
                    logging.info(result.stdout)
                except subprocess.CalledProcessError as e:
                    logging.error(f"Error running HMMER: {e}")
                    logging.error(e.stderr)
                bar()

    afa_files_dict = {} #k = Afa file name : v = path to Afa file
    for afa_files in os.listdir(database_alignment):
        afa_file_name = afa_files.split('.afa')
        afa_file_key = afa_file_name[0]
        afa_files_dict[afa_file_key] = os.path.join(database_alignment, afa_files)

    #Output/Species_folder/hmm_results
    for species_folder_in_output in os.listdir(output_dir):
        species_folder_path = os.path.join(output_dir, species_folder_in_output)

        if not os.path.isdir(species_folder_path):
            continue

        # Non conflict
        list_file_dict = {}
        bed_file_dict = {}    

        # Conflict
        conflict_list_dict = {}
        conflict_bed_dict = {}

        # Table
        table_count_dict = {}

        for hmm_result_file in os.listdir(species_folder_path):
            if hmm_result_file.endswith('.hmm_results'):
                species_name = hmm_result_file.strip().split('.hmm_results')[0]
                hmm_result_file_path = os.path.join(species_folder_path, hmm_result_file)
                with open(hmm_result_file_path, 'r') as result_file:
                    for lines in result_file:
                        if lines.startswith('#') or lines.startswith('-'):
                            continue
                        else:
                            remove_empty_spaces = re.sub(r'\s+', '\t', lines)
                            columns = remove_empty_spaces.strip().split('\t')
                            gene_id, domain_id, ali_from, ali_to = columns[0], columns[3], columns[17], columns[18]
                            
                            # Create list_file_dict, Keys = domain_id : Values = gene_id
                            if domain_id not in list_file_dict:
                                list_file_dict[domain_id] = []
                            list_file_dict[domain_id].append(gene_id)

                            # Create bed_file_dict, Keys = domain_id : Values = gene_id \t ali_from \t ali_to
                            if domain_id not in bed_file_dict:
                                bed_file_dict[domain_id] = []
                            bed_file_dict[domain_id].append(f'{gene_id}\t{int(ali_from) - 1}\t{ali_to}')

                            # Create conflict_list_dict, keys = gene_id : values = domain_id
                            if gene_id not in conflict_list_dict:
                                conflict_list_dict[gene_id] = []
                            conflict_list_dict[gene_id].append(domain_id)

                            # Create conflict_bed_dict, Keys = gene_id : values = ali_from \t ali_to
                            if gene_id not in conflict_bed_dict:
                                conflict_bed_dict[gene_id] = []
                            conflict_bed_dict[gene_id].append(f'{int(ali_from) - 1}\t{ali_to}') 

                            # Create Table, Keys = species_name : Values = counts
                            if species_name not in table_count_dict:
                                table_count_dict[species_name] = {}
                            if domain_id not in table_count_dict[species_name]:
                                table_count_dict[species_name][domain_id] = 1
                            else:
                                table_count_dict[species_name][domain_id] += 1

                # Write '{species}_{domain}.list' files
                for k, v in list_file_dict.items():
                    output_list_file_name = f'{species_name}_{k}.list'
                    with open(os.path.join(species_folder_path, output_list_file_name), 'w') as writing_list_output_file:
                        for genes in v:
                            writing_list_output_file.write(f'{genes}\n')

                # Write '{species}_{domain}_domain.bed' files
                for k, v in bed_file_dict.items():
                    output_bed_file_name = f'{species_name}_{k}_domain.bed'
                    with open(os.path.join(species_folder_path, output_bed_file_name), 'w') as writing_bed_output_file:
                        for bed_format in v:
                            writing_bed_output_file.write(f'{bed_format}\n')

                # Write '{species}_conflict.list' files
                output_conflict_file_name = f'{species_name}_conflict.list'
                with open(os.path.join(species_folder_path, output_conflict_file_name), 'w') as writing_conflict_output_file:
                    for k, v in conflict_list_dict.items():
                        remove_duplicates = list(set(v))
                        if len(remove_duplicates) > 1:
                            conflicting_domain = '\t'.join(remove_duplicates)
                            #writing_conflict_output_file.write(f'{k}\t{conflicting_domain}\n')
                            writing_conflict_output_file.write(f'{k}\n')
                                      
###################################################################################################################
        for conflict_list in os.listdir(output_dir):
            if conflict_list.endswith('_conflict.list'):
                species_name = conflict_list.strip().split('_conflict.list')[0]
                conflict_file_path = os.path.join(output_dir, conflict_list)
                output_conflict_bed_name = f'{species_name}_conflict_domain.bed'
                with open(conflict_file_path, 'r') as result_file, open(os.path.join(output_dir, output_conflict_bed_name), 'w') as writing_conflict_bed_file:
                    test_dict = {}
                    for lines in result_file:
                        gene_id_in_conflict_file = lines.strip().split('\t')[0]
                        test_dict[gene_id_in_conflict_file] = []
                        for k, v in conflict_bed_dict.items():
                            #print(f'{k}\t{v}\n')
                            if k == gene_id_in_conflict_file:
                                writing_conflict_bed_file.write(f'{v}')
###################################################################################################################
                                
        for k, v in table_count_dict.items():
            output_table_file_name = f'{k}_counts.tab'
            with open(os.path.join(species_folder_path, output_table_file_name), 'w') as writing_table_file:
                # Write header
                writing_table_file.write(k + '\t'.join([''] + list(v.keys())) + '\n')
        
                # Write species name and counts
                writing_table_file.write(f'{k}\t')
                for count in v.values():
                    writing_table_file.write(f'{count}\t')
                writing_table_file.write('\n')

        # Create Fasta files using '.list' against 'clean.fasta'
        for list_file in os.listdir(species_folder_path):
            if list_file.endswith('.list'):
                list_file_path = os.path.join(species_folder_path, list_file)
                name = list_file.replace('.list', '')
                output_file_name = f'{list_file.split(".list")[0]}.fasta'
                output_file_name_path = os.path.join(species_folder_path, output_file_name)
                for clean_files in os.listdir(directory):
                    if clean_files.endswith('_clean.fasta') and clean_files.startswith(name.split('_')[0]):
                        clean_fasta_path = os.path.join(directory, clean_files)
                        running_seqkit = f'seqkit grep -f {list_file_path} {clean_fasta_path} -o {output_file_name_path}'
                        with alive_bar(enrich_print=False) as bar:
                            try:
                                print('\n')
                                print(f'Processing {list_file} and {clean_files}')
                                seqkit = subprocess.run(running_seqkit, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                                logging.info(seqkit.stdout)
                            except subprocess.CalledProcessError as e:
                                logging.error(f'Error running {running_seqkit} as {e}')
                                logging.erorr(e.stderr)
                            finally:
                                bar()

            elif list_file.endswith('_domain.bed'):
                list_file_path = os.path.join(species_folder_path, list_file)
                name = list_file.replace('_domain.bed', '')
                output_file_name = f'{list_file.split("_domain.bed")[0]}_domain.fasta'
                output_file_name_path = os.path.join(species_folder_path, output_file_name)  
                for clean_files in os.listdir(directory):
                    if clean_files.endswith('_clean.fasta') and clean_files.startswith(name.split('_')[0]):
                        clean_fasta_path = os.path.join(directory, clean_files)
                        running_seqkit_subseq = f'seqkit subseq --update-faidx --bed {list_file_path} {clean_fasta_path} -o {output_file_name_path}'
                        with alive_bar(enrich_print=False) as bar:
                            try:
                                print('\n')
                                print(f'Processing {list_file} and {clean_files}')
                                seqkit_subseq = subprocess.run(running_seqkit_subseq, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                                logging.info(seqkit_subseq.stdout)
                            except subprocess.CalledProcessError as e:
                                logging.error(f'Error running {running_seqkit_subseq} as {e}')
                                logging.error(e.stderr)
                            finally:
                                bar()

        for fasta_file in os.listdir(species_folder_path):
            if fasta_file.endswith('_domain.fasta'):
                fasta_file_path = os.path.join(species_folder_path, fasta_file)
                #fasta_name = fasta_file.replace('_domain.fasta', '_muscled_domain.fasta')
                fasta_name = f'{fasta_file.split("_domain.fasta")[0]}_muscled_domain.fa'
                output_muscled_file_name = os.path.join(species_folder_path, fasta_name)
                muscle_command = f'muscle -in {fasta_file_path} -out {output_muscled_file_name}'
                with alive_bar(enrich_print=False) as bar:
                    try:
                        print('\n')
                        print(f'Processing {fasta_file} for alignment using MUSCLE')
                        muscle = subprocess.run(muscle_command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                        logging.info(muscle.stdout)
                    except subprocess.CalledProcessError as e:
                        logging.error(f'Error running {muscle_command} as {e}')
                        logging.error(e.stderr)
                    finally:
                        bar()

        for fasta in os.listdir(species_folder_path):
            if fasta.endswith('_muscled_domain.fa'):
                for k, v in afa_files_dict.items():
                    if k in fasta:
                        #fasta_file_new_name = fasta.replace("_muscled_domain.fasta", "_muscled_combined_domain.afa")
                        fasta_file_new_name = f'{fasta.split("_muscled_domain.fa")[0]}_muscled_combined_domain.afa'
                        running_profile = f'muscle -profile -in1 {v} -in2 {os.path.join(species_folder_path, fasta)} -out {os.path.join(species_folder_path, fasta_file_new_name)}'
                        with alive_bar(enrich_print=False) as bar:
                            try:
                                print('\n')
                                print(f'Processing {fasta} to create Combined FASTA file using MUSCLE profile')
                                muscle_profile = subprocess.run(running_profile, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                                logging.info(muscle_profile.stdout)
                            except subprocess.CalledProcessError as e:
                                logging.error(f'Error running {running_profile} as {e}')
                                logging.error(e.stderr)
                            finally:
                                bar()

        for muscled_file in os.listdir(species_folder_path):
            if muscled_file.endswith('_muscled_domain.fa'):
                muscled_file_path = os.path.join(species_folder_path, muscled_file)
                #trimal_name = muscled_file.replace('_muscled_domain.fasta', '_muscled_trimal_domain.fasta')
                trimal_name = f'{muscled_file.split("_muscled_domain.fa")[0]}_muscled_trimal_domain.fa'
                output_trimal_file_name = os.path.join(species_folder_path, trimal_name)
                running_trimal = f'trimal -in {muscled_file_path} -out {output_trimal_file_name} -gt 0.50 -cons 60'
                try:
                    trimal = subprocess.run(running_trimal, shell = True, check = True)
                    logging.info(trimal.stdout)
                except subprocess.CalledProcessError as e:
                    logging.error(e.stderr)
            elif muscled_file.endswith('_muscled_combined_domain.afa'):
                combined_muscled_file_path = os.path.join(species_folder_path, muscled_file)
                trimal_combined_name = f'{combined_muscled_file_path.split("_muscled_combined_domain.afa")[0]}_muscled_combined_trimal_domain.afa'
                running_trimal = f'trimal -in {combined_muscled_file_path} -out {trimal_combined_name} -gt 0.50 -cons 60'
                subprocess.run(running_trimal, shell = True, check = True)
                try:
                    trimal = subprocess.run(running_trimal, shell = True, check = True)
                    logging.info(trimal.stdout)
                except subprocess.CalledProcessError as e:
                    logging.error(e.stderr)

        for seqkit_fasta_files in os.listdir(species_folder_path):
            if seqkit_fasta_files.endswith('.fasta') or seqkit_fasta_files.endswith('.fa'):
                seqkit_file_path = os.path.join(species_folder_path, seqkit_fasta_files)
                running_sed = f"sed -i 's/:.//g' {seqkit_file_path}"
                try:
                    sed = subprocess.run(running_sed, shell = True, check = True)
                    logging.info(sed.stdout)
                except subprocess.CalledProcessError as e:
                    logging.error(e.stderr)
#
    #species_dict = {}
    #for species_file in os.listdir(output_dir):
    #    if species_file.endswith('_domain.fasta')
                #

        if visualization:
            for muscle_or_trimal_file in os.listdir(species_folder_path):
                # Handle .afa extensions
                if muscle_or_trimal_file.endswith('_muscled_combined_trimal_domain.afa') or muscle_or_trimal_file.endswith('_muscled_combined_domain.afa'):
                    muscle_path = os.path.join(species_folder_path, muscle_or_trimal_file)
                    weblogo_muscle_name = muscle_or_trimal_file.replace('.afa', '_weblogo.pdf')
                    output_weblogo_muscle = os.path.join(species_folder_path, weblogo_muscle_name)
                    running_weblogo_for_muscle = f'weblogo -f {muscle_path} -D fasta -o {output_weblogo_muscle} -F pdf -n 80 --resolution 400'
                    with alive_bar(enrich_print=False) as bar:
                        try:
                            print('\n')
                            print(f'Processing {weblogo_muscle_name}')
                            weblogo = subprocess.run(running_weblogo_for_muscle, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                            logging.info(weblogo.stdout)
                        except subprocess.CalledProcessError as e:
                            logging.error(f'Error running {running_weblogo_for_muscle} as {e}')
                            logging.error(e.stderr)
                        finally:
                            bar()

                    pymsaviz_path = os.path.join(species_folder_path, muscle_or_trimal_file)
                    pymsaviz_name = muscle_or_trimal_file.replace('.afa', '_pymsaviz.pdf')
                    output_pymsaviz = os.path.join(species_folder_path, pymsaviz_name)
                    running_weblogo_for_trimal = f'pymsaviz -i {pymsaviz_path} -o {output_pymsaviz} --wrap_length 80 --color_scheme Taylor --show_consensus --show_count'
                    with alive_bar(enrich_print=False) as bar:
                        try:
                            print('\n')
                            print(f'Processing {pymsaviz_name}')
                            pymsaviz = subprocess.run(running_weblogo_for_trimal, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                            logging.info(pymsaviz.stdout)
                        except subprocess.CalledProcessError as e:
                            logging.error(f'Error running {running_weblogo_for_trimal} as {e}')
                            logging.error(e.stderr)
                        finally:
                            bar()
                # Handle .fa extensions
                if muscle_or_trimal_file.endswith('_muscled_domain.fa') or muscle_or_trimal_file.endswith('_muscled_trimal_domain.fa'):
                    muscle_path = os.path.join(species_folder_path, muscle_or_trimal_file)
                    weblogo_muscle_name = muscle_or_trimal_file.replace('.fa', '_weblogo.pdf')
                    output_weblogo_muscle = os.path.join(species_folder_path, weblogo_muscle_name)
                    running_weblogo_for_muscle = f'weblogo -f {muscle_path} -D fasta -o {output_weblogo_muscle} -F pdf -n 80 --resolution 400'
                    with alive_bar(enrich_print=False) as bar:
                        try:
                            print('\n')
                            print(f'Processing {weblogo_muscle_name}')
                            weblogo = subprocess.run(running_weblogo_for_muscle, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                            logging.info(weblogo.stdout)
                        except subprocess.CalledProcessError as e:
                            logging.error(f'Error running {running_weblogo_for_muscle} as {e}')
                            logging.error(e.stderr)
                        finally:
                            bar()

                    pymsaviz_path = os.path.join(species_folder_path, muscle_or_trimal_file)
                    pymsaviz_name = muscle_or_trimal_file.replace('.fa', '_pymsaviz.pdf')
                    output_pymsaviz = os.path.join(species_folder_path, pymsaviz_name)
                    running_weblogo_for_trimal = f'pymsaviz -i {pymsaviz_path} -o {output_pymsaviz} --wrap_length 80 --color_scheme Taylor --show_consensus --show_count'
                    with alive_bar(enrich_print=False) as bar:
                        try:
                            print('\n')
                            print(f'Processing {pymsaviz_name}')
                            pymsaviz = subprocess.run(running_weblogo_for_trimal, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                            logging.info(pymsaviz.stdout)
                        except subprocess.CalledProcessError as e:
                            logging.error(f'Error running {running_weblogo_for_trimal} as {e}')
                            logging.error(e.stderr)
                        finally:
                            bar()
def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="List files from specified directories and set number of CPUs.")
    # Add arguments for input folder, database folder, and number of CPUs
    parser.add_argument("--input", help="Path to the input folder", required=True)
    parser.add_argument("--db", help="Path to the database folder", required=True)
    parser.add_argument("--output", help="Name for your output folder, Default folder name is 'output'.", default='output', required=False)
    parser.add_argument("--E", help="HMMsearch option. Report sequence E-value threshold value is defaulted as 1e-5", default='1e-5', required=False)
    parser.add_argument("--domE", help="HMMsearch option. Report domains E-value threshold value is 1e-10", default='1e-10', required=False)
    parser.add_argument("--CPU", help="Number of CPUs", default=2, type=int, required=False)
    parser.add_argument("--logging", help="Enable logging", action='store_true', required=False)
    parser.add_argument("--visualization", help="Call Weblogo and Pymsaviz", action='store_true', required=False)
    
    # Parse the arguments
    args = parser.parse_args()

    # Configure logging based on user input
    if args.logging:
        logging.basicConfig(level=logging.INFO,filename='hmm_run.log', format='%(asctime)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(level=logging.CRITICAL)  # Print CRITICAL logs only

    logging.info("Starting processing of FASTA files.")
    process_fasta_files(args.input)

    logging.info("Creating HMM database.")
    database_path = list_files(args.db)
    if database_path:
        logging.info("Running HMM searches.")
        run_hmm(args.input, args.CPU, database_path, args.output, args.E, args.domE, args.visualization)
    else:
        logging.error("Database creation failed or was skipped.")

if __name__ == "__main__":
    main()
