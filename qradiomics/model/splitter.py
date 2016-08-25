from sklearn.cross_validation import StratifiedKFold
from sklearn.cross_validation import KFold

def split(X, y, eval_size):
    # eval_size = 0.10
    is_labels = True

    if is_labels:
        kf = StratifiedKFold(y, round(1. / eval_size))
    else:
        kf = KFold(y, round(1. / eval_size))
    train_indices, valid_indices = next(iter(kf))

# X_train, y_train = X[train_indices], y[train_indices]
# X_valid, y_valid = X[valid_indices], y[valid_indices]
