import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
from anndata import AnnData
from xml.dom.minidom import parse
import xml.dom.minidom

def gem2adata(gem_data_raw,bin,geneID='geneID',x='x',y='y',MIDCounts='MIDCounts'):
    gem_data = gem_data_raw.copy()
    gem_data['x'] = gem_data[x]
    gem_data['y'] = gem_data[y]
    gem_data['geneID'] = gem_data[geneID]
    gem_data['MIDCounts'] = gem_data[MIDCounts]
    
    gem_data['x'] = round(gem_data['x'] / bin, 0).astype(int)
    gem_data['y'] = round(gem_data['y'] / bin, 0).astype(int)
    gem_data = gem_data.groupby(['x', 'y', 'geneID']).agg(MIDCounts=('MIDCounts', 'sum')).reset_index()
    gem_data = gem_data[['geneID', 'x', 'y', 'MIDCounts']]
    gem_data['cell'] = gem_data['x'].map(str) + '_' + gem_data['y'].map(str)

    uniq_cell, uniq_gene = gem_data.cell.unique(), gem_data.geneID.unique()
    uniq_cell, uniq_gene = list(uniq_cell), list(uniq_gene)
    cell_dict = dict(zip(uniq_cell, range(0, len(uniq_cell))))
    gene_dict = dict(zip(uniq_gene, range(0, len(uniq_gene))))
    gem_data["csr_x_ind"] = gem_data["cell"].map(cell_dict)
    gem_data["csr_y_ind"] = gem_data["geneID"].map(gene_dict)

    matrix = csr_matrix((gem_data['MIDCounts'], (gem_data["csr_x_ind"], gem_data["csr_y_ind"])),
                        shape=((len(uniq_cell), len(uniq_gene))))

    var = pd.DataFrame({"gene_short_name": uniq_gene})
    var.set_index("gene_short_name", inplace=True)
    obs = pd.DataFrame({"cell": gem_data['cell'].unique().tolist()})
    obs["coor_x"] = pd.Series(obs['cell']).str.split('_', expand=True).iloc[:, 0].astype(float)
    obs['coor_y'] = pd.Series(obs['cell']).str.split('_', expand=True).iloc[:, 1].astype(float)
    obsm = {"spatial": obs.loc[:, ['coor_x', "coor_y"]].values}
    obs.set_index("cell", inplace=True)
    obs['nCount_RNA'] = np.sum(matrix,axis=1)
    tmp_X = matrix.copy()
    tmp_X[tmp_X>0] = 1
    obs['nFeature_RNA'] = np.sum(tmp_X,axis=1)
    adata = AnnData(matrix, dtype=matrix.dtype, obs=obs.copy(), var=var.copy(), obsm=obsm.copy())
    return adata

def get_xml_matrix(path):
    DOMTree = xml.dom.minidom.parse(path)
    collection = DOMTree.documentElement
    patchs = collection.getElementsByTagName("t2_patch")
    
    title_list = []
    matrix_list = []
    for patch in patchs:
    
        patch_list = []
        if patch.hasAttribute("title"):
            tit = patch.getAttribute("title")
            title_list.append(tit)
        if patch.hasAttribute("transform"):
            transform = patch.getAttribute("transform")
            transform = transform[7:-1]
            matrix_list.append(transform)
    
    #json_list=str(json_list).replace("'","\"")
    df = pd.DataFrame({'title':title_list,'matrix':matrix_list})
    return df


def trakEM_transform(ix, iy, matrix, matrix_bin,transform_bin):
    affine = np.zeros((3, 3))
    matrix = np.array(matrix.split(',')).astype(float).reshape((3, 2))
    affine[0:2, :] = matrix.T
    affine[2] = [0, 0, 1]
    
    gem_data = pd.DataFrame({'x': ix, 'y': iy})
    gem_data['i'] = 1
    scale = np.matrix([[transform_bin / matrix_bin, 0, 0],
                       [0, transform_bin / matrix_bin, 0],
                       [0, 0, 1]])

    affine = (scale.I).dot(affine).dot(scale)
    new_coor = (affine.dot([gem_data['x'], gem_data['y'], gem_data['i']])).T
    rx = new_coor[:, 0]
    ry = new_coor[:, 1]
    return rx, ry

