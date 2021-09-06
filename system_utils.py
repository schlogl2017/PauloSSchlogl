#!usr/bin/env python
# -*- coding: utf-8 -*-


import os
from collections import defaultdict, Counter
import collections
import pandas as pd
import numpy as np


def get_files(dir_name):
    # create a list of file and sub directories
    # names in the given directory
    files = os.listdir(dir_name)
    all_files = []
    # Iterate over all the entries
    for entry in files:
        # Create full path
        full_path = os.path.join(dir_name, entry)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(full_path):
            all_files += get_files(full_path)
        else:
            all_files.append(full_path)
    return all_files


def get_names(path):
    return [(file.split('/')[3][:-3]) for file in get_files(path)]


def get_dir_names(filenames):
    dirn = set()
    subd = set()
    subsub = set()
    for filename in filenames:
        dir_splt = filename.split('/')
        dirn.add(dir_splt[0])
        subd.add(dir_splt[1])
        subsub.add(dir_splt[2])
    return next(iter(dirn)), next(iter(subd)), next(iter(subsub))


def get_dir_name(filename):
    subd, subsub = set(), set()
    names = filename.split('/')
    for i in range(len(names)):
        subd.add(names[1])
        subsub.add(names[2])
    return next(iter(subd)), next(iter(subsub))


def make_me_a_folder(folder_name):
    os.getcwd()
    try:
        os.makedirs(folder_name)
    except OSError:
        print("Directory already exists")
    return os.path.exists(folder_name)


def get_full_name(dir_root, sud_dir, ssub):
    return os.path.join(dir_root, sud_dir, ssub)


def get_fasta_files(dir_name):
    # tuple with the file extention to check
    ext = tuple(['fa.gz','.fa', '.fasta', '.fa.gz', '.fna', '.fna.gz'])
    infiles = []
    # lokking for file in the path
    for path, subdirs, files in os.walk(dir_name):
        for name in files:
            # getting the filenames
            input_files = os.path.join(path, name)
            if input_files.endswith(ext):
                infiles.append(input_files)
    return infiles


def get_files2(dir_name, subdir):
    # create a list of file and sub directories
    # names in the given directory
    files = os.listdir(dir_name)
    all_files = []
    # Iterate over all the entries
    for entry in files:
        # Create full path
        full_path = os.path.join(dir_name, entry, subdir)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(full_path):
            all_files += get_files(full_path)
        else:
            all_files.append(full_path)
    return all_files


def get_list_paths(spc_names, path, ssub):
    pwd = defaultdict()
    for name in spc_names:
        pwd[name] = get_full_name(path, name, ssub)
    return pwd


def get_fasta_files_paths(dict_paths):
    dict_fasta_paths = defaultdict(list)
    for name, path in dict_paths.items():
        dict_fasta_paths[name] = dict_fasta_paths.get(name, []) + glob.glob(path)
    return dict_fasta_paths


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z


def get_fasta_files_dict(dir_name):
    ext = tuple(['fa.gz','.fa', '.fasta', '.fa.gz', '.fna', '.fna.gz'])
    infiles = defaultdict(list)
    for path, subdirs, files in os.walk(dir_name):
        for name in files:
            infiles[name] = infiles.get(name, [])
            input_files = os.path.join(path, name)
            if input_files.endswith(ext):
                infiles[name].append(input_files)
    return infiles


def fasta_file_paths(fasta_dirs):
    fasta_paths = defaultdict(list)
    exts = tuple(['fa.gz', '.fa', '.fasta', '.fa.gz', '.fna', '.fna.gz'])
    for name, pwd in fasta_dirs.items():
        fasta_paths[name] = fasta_paths.get(name, [])
        for path, _, files in os.walk(pwd):
            for n in files:
                input_files = os.path.join(path, n)
                if input_files.endswith(exts):
                    fasta_paths[name].append(input_files)
    return fasta_paths


def get_names_list(text_file):
    """
    Function to get the names of species or genus in a text file.
    
    Inputs:
    
        text_file - a string representing a text file where the 
                    species or genus are listed one by line.
    
    Outputs:
    
        names_list - a list/array-like that countain the name of the
                     species or genus.
                
        > get_names_list('teste_genus.txt')
        
        ['Acidiphilium']
    """
    # initialize the list
    names_list = []
    # open the text file
    with open(text_file, 'r') as fh:
        # iterates through the file handle
        for line in fh:
            # clean the spaces
            name = line.strip()
            # add the name to the list
            names_list.append(name)
    return names_list


def find_csv_filenames(dir_name, spc_name, sub_dir, ssub, num_k):
    name = os.path.join(dir_name, spc_name, sub_dir)
    path_files = [os.path.join(name, f'{ssub}{str(i)}') for i in range(1, num_k)]
    files_dict = defaultdict(list)
    files_dict[spc_name] = files_dict.get(spc_name, [])
    for path in path_files:
        filenames = ''.join(os.listdir(path))
        files_dict[spc_name].append(os.path.join(path, filenames))
    return files_dict


def get_csvs_to_df_concatenation(spc_name, filenames_dict):
    """Function tha receives a species name and a dictionary with
    all csvs files to be concatenated.
    """
    dfs = []
    for filename in filenames_dict[spc_name]:
        dfs.append(pd.read_csv(filename, dtype={'kmers': str, 'counts': np.int32}))
    dfs = [df.set_index("kmers", drop=True) for df in dfs]
    concat = pd.concat(dfs, axis=0, copy=False).reset_index()
    return concat, len(dfs)


def csv_filenames(dir_name, spc_name, sub_dir, ssub):
    name = os.path.join(dir_name, spc_name, sub_dir)
    path_files = [os.path.join(name, ssub)]
    files_dict = defaultdict(list)
    files_dict[spc_name] = files_dict.get(spc_name, [])
    for path in path_files:
        filenames = ''.join(os.listdir(path))
        files_dict[spc_name].append(os.path.join(path, filenames))
    return files_dict


def path_generator(dir_name, sub_dir, genus, sub_sub_dir, type_file):
    """
    Function to generate one by one the complete paths to the 
    fasta files.
    The directory need to be designed like this:
        base directory/sub_directory/genus/sub_sub_directory/files.type
        Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000016725.1_chr.fna.gz
        
    Inputs:
        
        dir_name - string representing the base directory.
        sub_dir - string representing the sub directory.
        genus - string representing the bacterial genus.
        sub_sub_dir - string representing the sub sub directory, where the files
                      are located.
        type_file - represents the final file extention (csv, gz, etc)
    
    Outputs:
    
        Yields the complete path to the source files.
        
        >list(path_generator('Data', 'Genomes_splitted', 'Acidiphilium', 'Chromosomes', 'gz'))
        
        ['Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000016725.1_chr.fna.gz',
        'Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000202835.1_chr.fna.gz']
        
    """
    # generates the paths to the directory
    pwd = os.path.join(dir_name, sub_dir, genus, sub_sub_dir)
    # iterates through the directory
    for filename in os.listdir(pwd):
        # checks if the files have the correct type
        # and if so yields the paths
        if filename.endswith(type_file):
            yield os.path.join(pwd, filename)


def data_generator(file_path_dict):
    """
    Yields a tuples as names/genus and the full path to the source
    files.
    
    Inputs:
    
        file_path_dict - a dictionary-like object mapping the names/genus to
                         the source of files to be analized.
    
    Outputs:
    
        Yields a tuple with first element as the name or genus and the second
        element the full path to the source file.
    
    for data in data_generator(file_dict):
        print(data)
    
    ('Acidiphilium', 'Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000016725.1_chr.fna.gz')
    ('Acidiphilium', 'Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000202835.1_chr.fna.gz')
    """
    # iterates through the dictionary items
    for name, filenames in file_path_dict.items():
        # gets the name
        spc = name
        # and then iterates through the filenames
        # that is a list with the full path
        # then yields the name and the full path
        for filename in filenames:
            yield spc, filename

def file_dict(text_file, dir_name, sub_dir, sub_sub_dir, type_file):
    """
    Function to generate a dictionary-like object with the name/genus to
    the complete path to the fasta files.
    The directory need to be designed like this:
        base directory/sub_directory/genus/sub_sub_directory/files.type
        Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000016725.1_chr.fna.gz
        
    Inputs:
    
        text_file - a string representing a text file where the 
                    species or genus are listed one by line.        
        dir_name - string representing the base directory.
        sub_dir - string representing the sub directory.
        sub_sub_dir - string representing the sub sub directory, where the files
                      are located.
        type_file - represents the final file extention (csv, gz, etc)
    
    Outputs:
    
        paths_to_files - a dictionary-like object where keys are names and
                         values are the full path to the files

    >file_dict('teste_genus.txt', 'Data', 'Genomes_splitted', 'Chromosomes', 'gz')
    defaultdict(list,
            {'Acidiphilium': 
            ['Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000016725.1_chr.fna.gz',
              'Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000202835.1_chr.fna.gz']})
    """
    # initialze the dictionary with default values as lists
    paths_to_files = defaultdict(list)
    # get the names to the keys
    names = get_names_list(text_file)
    # iterates through the names list
    for name in names:
        # add the keys to the dictionary
        # add the full paths to the values items
        # as a list and return the dictionary
        paths_to_files[name] = paths_to_files.get(name, []) + list(path_generator(dir_name,
                                                                            sub_dir,
                                                                            name,
                                                                            sub_sub_dir,
                                                                            type_file))
    return paths_to_files


