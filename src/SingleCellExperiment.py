#SingleCellExperiment
import os
import pandas as pd

class SingleCellExperiment:
    barcodes_file = ""
    seurat_file = ""
    metacell_file = ""
    ecDNA_file = ""
    df=pd.DataFrame()

    def __init__(self,barcodes_file, seurat_file=None, metacell_file=None, ecDNA_file=None):
        self.add_barcodes(barcodes_file)
        if seurat_file != None:
            self.add_seurat(seurat_file)
        if metacell_file != None:
            self.add_metacell(metacell_file)
        if ecDNA_file != None:
            self.add_ecDNA(ecDNA_file)

    def __repr__(self):
        return repr(self.df)
    
    def add_barcodes(self,barcodes_file):
        '''
        Makes a new df as an index of barcodes. 
        Currently erases previous df.
        '''
        with open(barcodes_file,'r') as f:
            barcodes = list(map(lambda x: x.strip(),f.readlines()))
        self.df=pd.DataFrame(index=barcodes)

    def add_seurat(self,seurat_file):
        '''
        Add a seurat metadata table
        '''
        seurat = pd.read_csv(seurat_file,sep='\t')
        seurat=seurat[['nCount_RNA','nFeature_RNA','percent.mt','nCount_ATAC','nFeature_ATAC','seurat_clusters']]
        seurat=seurat.rename(columns={'seurat_clusters':'seurat_cluster'})
        self.df['qc_pass_seurat'] = self.df.index.isin(seurat.index)
        self.df=self.df.join(seurat,how='left')
    
    def add_ecDNA(self,ecDNA_file):
        '''
        Add an ecDNA metadata table. See 
        2021-10-07_better-background/permutation-testing.ipynb
        to generate this file.
        '''
        ecDNA_annot = pd.read_csv(ecDNA_file,sep='\t',index_col=0)
        self.df=self.df.join(ecDNA_annot,how='left')

    def add_metacell(self,metacell_file):
        '''
        Add metacell annotations. See 2021-09-07_metacell/metacell.ipynb.
        '''
        mc_annot = pd.read_csv(metacell_file,sep='\t',index_col=0)
        mc_annot.columns = ['metacell_cluster']
        self.df = self.df.join(mc_annot,how='left')

class rcmb56pdx(SingleCellExperiment):

    root="/mnt/c/Users/ochapman/Documents/Mesirov/scRNA+ATAC/"
    barcodes_file = root+"RCMB56-pdx/RCMB56-pdx/outs/filtered_feature_bc_matrix/barcodes.tsv"
    seurat_file = root+"2021-08-29_seurat/rcmb56-pdx_seurat_metadata.tsv"
    metacell_file = ''
    ecDNA_file = root+"2021-10-07_better-background/rcmb56-pdx/ecDNA_status_gbg.tsv"
    cs_file = root+'2021-08-27_Kazachkova/Kazachkova_rotation_sp21/calls_5000_1e6.csv'
    df=pd.DataFrame()

    def __init__(self):
        super().__init__(
            self.barcodes_file,
            seurat_file=self.seurat_file,
            ecDNA_file=self.ecDNA_file
        )
        self.add_copyscat_dm(self.cs_file)
    
    def add_copyscat_dm(self,cs_file):
        df = pd.read_csv(cs_file,index_col=0)
        df = df.rename(columns={'chr1_83000000.chr1_96000000':'cs_ecDNA1'})
        df['cs_ecDNA2&3'] = df[['chr7_0.chr7_8000000','chr7_86000000.chr7_96000000','chr7_98000000.chr7_107000000',
                     'chr7_125000000.chr7_131000000','chr7_41000000.chr7_55000000','chr7_32000000.chr7_35000000']].any(axis='columns',bool_only=True)
        df = df[['cs_ecDNA1','cs_ecDNA2&3']]
        self.df = self.df.join(df,how='left')

class rcmb56ht(SingleCellExperiment):

    root="/mnt/c/Users/ochapman/Documents/Mesirov/scRNA+ATAC/"
    barcodes_file = root+"RCMB56-ht/cellranger-2.0.0/outs/filtered_feature_bc_matrix/barcodes.tsv"
    seurat_file = root+"2021-08-29_seurat/rcmb56-ht_seurat_metadata.tsv"
    metacell_file = root+"2021-09-07_metacell/rcmb56-ht_metacell_ids.tsv"
    ecDNA_file = root+"2021-10-07_better-background/rcmb56-ht/ecDNA_status_gbg.tsv"
    cs_file = root+"2021-08-27_Kazachkova/copyscat_dm.csv"
    df=pd.DataFrame()

    def __init__(self):
        super().__init__(
            self.barcodes_file,
            seurat_file=self.seurat_file,
            ecDNA_file=self.ecDNA_file,
            metacell_file=self.metacell_file
        )
        self.add_copyscat_dm(self.cs_file)

    def add_copyscat_dm(self,cs_file):
        df=pd.read_csv(cs_file,index_col=0)
        df=df[['chr1_83000000.chr1_97000000','chr7_0.chr7_6000000']]
        df=df.rename(columns={'chr1_83000000.chr1_97000000':'cs_ecDNA1','chr7_0.chr7_6000000':'cs_ecDNA2&3'})
        self.df = self.df.join(df,how='left')

class rcmb56ht_alt_bg(rcmb56ht):
    root="/mnt/c/Users/ochapman/Documents/Mesirov/scRNA+ATAC/"
    ecDNA_file = root+"2021-10-07_better-background/rcmb56-ht/ecDNA_status_ebg.tsv"
    def __init__(self):
        super().__init__()

class rcmb56pdx_alt_bg(rcmb56pdx):
    root="/mnt/c/Users/ochapman/Documents/Mesirov/scRNA+ATAC/"
    ecDNA_file = root+"2021-10-07_better-background/rcmb56-pdx/ecDNA_status_ebg.tsv"
    def __init__(self):
        super().__init__()