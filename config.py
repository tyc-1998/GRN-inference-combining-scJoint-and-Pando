import torch
import os

class Config(object):
    def __init__(self):
        self.use_cuda = False
        self.threads = 1

        if not self.use_cuda:
            self.device = torch.device('cpu')
        else:
            self.device = torch.device('cuda:0')
        
        # DB info
        self.number_of_class = 9 #celltype numbers.
        self.input_size = 1880 #rna and atac common input gene numbers.
        self.rna_paths = ['./data/exprs_rna.npz']
        self.rna_labels = ['./data/celltype_rna.txt']             
        self.atac_paths = ['./data/exprs_atac.npz']
        self.atac_labels = [] #Optional. If atac_labels are provided, accuracy after knn would be provided.
        self.rna_protein_paths = []
        self.atac_protein_paths = []
        
        # Training config            
        self.batch_size = 256
        self.lr_stage1 = 0.01
        self.lr_stage3 = 0.01
        self.lr_decay_epoch = 20
        self.epochs_stage1 = 20
        self.epochs_stage3 = 20
        self.p = 0.8
        self.embedding_size = 64
        self.momentum = 0.9
        self.center_weight = 1
        self.with_crossentorpy = True
        self.seed = 1
        self.checkpoint = ''