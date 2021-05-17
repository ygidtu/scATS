from math import e, log
import numpy as np

def entropy(mtx):

    def __entropy__(labels):
        n_labels = len(labels)

        if n_labels <= 1:
            return 0

        probs = labels / n_labels
        ent = 0.

        for i in probs:
            if i > 0:
                ent -= i * log(i, e)

        return ent

    return np.array([__entropy__(x) for x in mtx])
