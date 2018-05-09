import numpy as np

def get_some_trait():
    mat = np.array([[1,0,0],[0,0,1],[0,1,0]])
    det = np.linalg.det(mat)
    print("the det of mat is:",det,"\n")

    inv_mat = np.linalg.inv(mat)
    print("the inv_mat is:")
    print(inv_mat,"\n")

    eig1,eig2 = np.linalg.eig(mat)
    print("the eig of mat is:",eig1,"\n")

    print("the feature vector of mat is:")
    print(eig2)

get_some_trait()
