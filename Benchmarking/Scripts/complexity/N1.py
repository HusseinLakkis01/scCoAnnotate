import numpy as np
import itertools
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import distance
from sys import argv
from sklearn.preprocessing import LabelEncoder
import pandas as pd
import os
import seaborn as sns 
import matplotlib as plt

def precompute_fx(X: np.ndarray, y: np.ndarray):
    """Precompute some useful things to support complexity measures.
    Parameters
    ----------
    X : :obj:`np.ndarray`, optional
            Attributes from fitted data.
    y : :obj:`np.ndarray`, optional
            Target attribute from fitted data.
    Returns
    -------
    :obj:`dict`
        With following precomputed items:
        - ``ovo_comb`` (:obj:`list`): List of all class OVO combination, 
            i.e., [(0,1), (0,2) ...].
        - ``cls_index`` (:obj:`list`):  The list of boolean vectors
            indicating the example of each class. 
        - ``cls_n_ex`` (:obj:`np.ndarray`): The number of examples in
            each class. The array indexes represent the classes.
    """
    prepcomp_vals = {}
    classes, class_freqs = np.unique(y, return_counts=True)
    cls_index = [np.equal(y, i) for i in classes]
    #cls_n_ex = np.array([np.sum(aux) for aux in cls_index])
    cls_n_ex = list(class_freqs)
    ovo_comb = list(itertools.combinations(range(classes.shape[0]), 2))
    prepcomp_vals["ovo_comb"] = ovo_comb
    prepcomp_vals["cls_index"] = cls_index
    prepcomp_vals["cls_n_ex"] = cls_n_ex
    return prepcomp_vals

# Fraction of Borderline Points (N1)
def ft_N1(X: np.ndarray, y: np.ndarray, metric: str = "euclidean") -> np.ndarray:
    # 0-1 scaler to compute distances 
    scaler = MinMaxScaler(feature_range=(0, 1)).fit(X)
    X_ = scaler.transform(X)
    # compute the distance matrix and the minimum spanning tree.
    dist_m = np.triu(distance.cdist(X_, X_, metric), k=1)
    mst = minimum_spanning_tree(dist_m)
    node_i, node_j = np.where(mst.toarray() > 0)
    # which edges have nodes with different class
    which_have_diff_cls = y[node_i] != y[node_j]
    # number of different vertices connected
    aux = np.unique(np.concatenate([node_i[which_have_diff_cls],node_j[which_have_diff_cls]])).shape[0]
    return aux/X.shape[0]

# Fraction of Borderline Points (N1)
def get_borderlinecells(X: np.ndarray, y: np.ndarray, metric: str = "euclidean") -> np.ndarray:
    # 0-1 scaler to compute distances 
    scaler = MinMaxScaler(feature_range=(0, 1)).fit(X)
    X_ = scaler.transform(X)
    # compute the distance matrix and the minimum spanning tree.
    dist_m = np.triu(distance.cdist(X_, X_, metric), k=1)
    mst = minimum_spanning_tree(dist_m)
    node_i, node_j = np.where(mst.toarray() > 0)
    # which edges have nodes with different class
    which_have_diff_cls = y[node_i] != y[node_j]
    # number of different vertices connected
    aux = np.unique(np.concatenate([node_i[which_have_diff_cls],node_j[which_have_diff_cls]]))
    cell_status = np.zeros(y.shape[0]).astype(str)
    cell_status[aux] = "Border cell"
    cell_status[~np.isin(np.arange(cell_status.size),aux)] = "Non-border cell"
    return cell_status

def loop_over_classes(data, y):
    '''
    Function that loops over the possible permutations of cell-types (classes) to calculate pairwise N1 between them
    input: data: Pandas Dataframe with genes as columns and cells as rows
           y: pandas Dataframe with a column 'label' containing cell types
    '''
    # get unique classes
    classes = np.unique(y['label'])
    # encode it for multiclass problems
    le = LabelEncoder()
    # Generate all possible combinations
    combinations = list(itertools.combinations(classes, 2))
    results= []
    # loop over all combinations and calculate pairwise 
    for comb in combinations:
        subset_ = y[y['label'].isin(comb)]['label'].astype(str)
        data_subset = data.loc[subset_.index]
        subset_ = le.fit_transform(subset_.to_numpy())    
        X = data_subset.to_numpy()
        precomp_fx = precompute_fx(X, subset_)
        cls_index = precomp_fx['cls_index']
        cls_n_ex = precomp_fx['cls_n_ex']
        ovo_comb = precomp_fx['ovo_comb']
        N1 = ft_N1(X,subset_)
        results.append(N1)
    return results, combinations
        

def create_conf_matrix(combinations, vals, classes):
    N1_results = {combinations[i]: vals[i] for i in range(len(combinations))}
    df=pd.DataFrame(np.zeros([len(classes),len(classes)]))
    df.columns = list(classes)
    df.index = list(classes)
    for comb in N1_results:
        print(comb)
        df.loc[comb[0],comb[1]] = N1_results[comb]
        df.loc[comb[1],comb[0]] = N1_results[comb]
    return df


#def plot_heatmap(df):
    #df = df * 100
    #ax = sns.heatmap(df)
    #plt.savefig("pairwiseN1_heatmap.png",bbox_inches='tight')


def main():

    data = pd.read_csv(argv[1],index_col = 0)
    labels =  pd.read_csv(argv[2],index_col = 0)
    classes = np.unique(labels['label'].astype(str))
    results, combinations = loop_over_classes(data, labels)
    df = create_conf_matrix(combinations, results, classes)
    os.chdir(argv[3])
    df.to_csv("pairwise_N1.csv")
    #plot_heatmap(df)
    cells = get_borderlinecells(data.to_numpy(), labels['label'].to_numpy())
    cells = pd.DataFrame(cells)
    cells.to_csv("BorderlineCells.csv")
    X = data.to_numpy()
    # creating initial dataframe
    le = LabelEncoder()
    y = le.fit_transform(labels['label'].to_numpy())
    N1 = ft_N1(X,y)
    results ={'N1':[N1]}
    results = pd.DataFrame.from_dict(results)
    results.to_csv("N1.csv")
    os.remove("VariableFeatures.csv")


    
if __name__ == "__main__":
    main()





