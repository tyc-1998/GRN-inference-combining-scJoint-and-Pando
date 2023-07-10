import numpy as np
import scipy.sparse
import h5py
import sys
import os

os.chdir("./data")

def h5_to_npz(file_name):
    h5 = h5py.File(file_name, 'r')
    h5_data = h5['matrix/data']
    print('H5 dataset shape:', h5_data.shape)

    sparse_data = scipy.sparse.csr_matrix(np.array(h5_data).transpose())
    scipy.sparse.save_npz(file_name.replace('h5', 'npz'), sparse_data)


def label_parsing(label_file):
    total_rna_labels = []
    with open(label_file) as fp:
        lines = fp.readlines()
        lines = lines[1:]
    
        for line in lines:
            line = line.split(',')
            label = line[1].replace('\"', '').replace('\n', '')
            total_rna_labels.append(label)
    
        # label to idx
        label_idx_mapping = {}
        unique_labels = np.unique(total_rna_labels)
        for i, name in enumerate(unique_labels):
            label_idx_mapping[name] = i
        print(label_idx_mapping)
    with open("label_to_idx.txt", "w") as fp:
        for key in sorted(label_idx_mapping):
            fp.write(key + " " + str(label_idx_mapping[key]) + '\n')
    
        # gen rna label files
    with open(label_file) as fp:
        lines = fp.readlines()
        lines = lines[1:]
    
        with open(label_file.replace('csv', 'txt'), 'w') as rna_label_f:
            for line in lines:
                line = line.split(',')
                label = line[1].replace('\"', '').replace('\n', '')
                rna_label_f.write(str(label_idx_mapping[label]) + '\n')

h5_to_npz("exprs_rna.h5")
h5_to_npz("exprs_atac.h5")
label_parsing("celltype_rna.csv")