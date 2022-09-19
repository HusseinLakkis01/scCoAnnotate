import numpy as np
import itertools
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import LeaveOneOut
from scipy.sparse.csgraph import minimum_spanning_tree
from sklearn.metrics import accuracy_score
from scipy.spatial import distance
from sys import argv
from sklearn.preprocessing import LabelEncoder
import pandas as pd
import os

def precompute_fx(X: np.ndarray, y: np.ndarray):

    """Precompute some useful things to support complexity measures.
    https://github.com/jcelias98/Complexity-Measures for the code
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

# Volume of Overlapping Region (F2)

def _minmax(X: np.ndarray, class1: np.ndarray, class2: np.ndarray) -> np.ndarray:
    """ This function computes the minimum of the maximum values per class
    for all features.
    """
    max_cls = np.zeros((2, X.shape[1]))
    max_cls[0, :] = np.max(X[class1], axis=0)
    max_cls[1, :] = np.max(X[class2], axis=0)
    aux = np.min(max_cls, axis=0)
    
    return aux

def _minmin(X: np.ndarray, class1: np.ndarray, class2: np.ndarray) -> np.ndarray:
    """ This function computes the minimum of the minimum values per class
    for all features.
    """
    min_cls = np.zeros((2, X.shape[1]))
    min_cls[0, :] = np.min(X[class1], axis=0)
    min_cls[1, :] = np.min(X[class2], axis=0)
    aux = np.min(min_cls, axis=0)
    
    return aux

def _maxmin(X: np.ndarray, class1: np.ndarray, class2: np.ndarray) -> np.ndarray:
    """ This function computes the maximum of the minimum values per class
    for all features.
    """
    min_cls = np.zeros((2, X.shape[1]))
    min_cls[0, :] = np.min(X[class1], axis=0)
    min_cls[1, :] = np.min(X[class2], axis=0)
    aux = np.max(min_cls, axis=0)
    
    return aux

def _maxmax(X: np.ndarray, class1: np.ndarray, class2: np.ndarray) -> np.ndarray:
    """ This function computes the maximum of the maximum values per class
    for all features. 
    """
    max_cls = np.zeros((2, X.shape[1]))
    max_cls[0, :] = np.max(X[class1], axis=0)
    max_cls[1, :] = np.max(X[class2], axis=0)
    aux = np.max(max_cls, axis=0)
    return aux

def ft_F2(X: np.ndarray, ovo_comb: np.ndarray, cls_index: np.ndarray) -> float:
    f2_list = []
    
    for idx1, idx2 in ovo_comb:
        y_class1 = cls_index[idx1]
        y_class2 = cls_index[idx2]
        zero_ = np.zeros(np.shape(X)[1])
        overlap_ = np.maximum(zero_, _minmax(X, y_class1, y_class2)-_maxmin(X, y_class1, y_class2))
        range_ = _maxmax(X, y_class1, y_class2)-_minmin(X, y_class1, y_class2)
        ratio = (overlap_)/(range_)
        f2_list.append(np.prod(ratio))
        
    return np.mean(f2_list)


def main():
    data = pd.read_csv(argv[1],index_col = 0)
    labels =  pd.read_csv(argv[2],index_col = 0)['label'].astype(str)
    X = data.to_numpy()
    # creating initial dataframe
    le = LabelEncoder()
    y = le.fit_transform(labels.to_numpy())
    precomp_fx = precompute_fx(X, y)
    cls_index = precomp_fx['cls_index']
    cls_n_ex = precomp_fx['cls_n_ex']
    ovo_comb = precomp_fx['ovo_comb']
    F2 = ft_F2(X, ovo_comb, cls_index)
    results ={'F2':[F2]}
    results = pd.DataFrame.from_dict(results)
    os.chdir(argv[3])
    results.to_csv("F2.csv")
    os.remove("pca.csv")


if __name__ == "__main__":
    main()