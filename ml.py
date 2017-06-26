import matplotlib.pyplot as plt
import numpy as np
import time
import sklearn
import sys
import argparse

from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import train_test_split
from sklearn.externals import joblib

class_freq_n = 5

def data_preprocess(data):
  unique, counts = np.unique(data[:,data.shape[1]-3], return_counts=True)
  class_freq = dict(zip(unique, counts))
  for key in class_freq:
    if(class_freq[key] < class_freq_n):
      data = data[data[:,data.shape[1]-3] != key,:]
  return data

def trainSVM(train_data,var):

    x_coord_train = train_data.shape[0]
    y_coord_train = train_data.shape[1]
    train_data = data_preprocess(train_data[0:var,:])
    X_ = train_data[:,0:y_coord_train - 3]
    Y_ = train_data[:,y_coord_train-3:y_coord_train-2]
    unique, counts = np.unique(Y_, return_counts=True)
    class_freq = dict(zip(unique, counts))
    print("Class frequencies:")
    print(class_freq)

    #X_train, X_test, Y_train, Y_test = train_test_split(X_, Y_, test_size=0.5, random_state=42, stratify=Y_)
    X_train, X_test, Y_train, Y_test = train_test_split(X_, Y_, test_size=0.5, random_state=42)
    tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4], 'C': [1, 10, 100, 1000]}, {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]
    svm = SVC()
    clf = GridSearchCV(svm,tuned_parameters,iid=False)
    clf.fit(X_train,Y_train.ravel())

    joblib.dump(clf,'SVM.pkl')
    print(clf.score(X_test,Y_test))

def testSVM(test_data,var):

    x_coord_test = test_data.shape[0]
    y_coord_test = test_data.shape[1]
    test_data = data_preprocess(test_data[0:var,:])
    X_t = test_data[:,0:y_coord_test - 3]
    Y_t = test_data[:,y_coord_test-3:y_coord_test-2]
    unique, counts = np.unique(Y_t, return_counts=True)
    class_freq = dict(zip(unique, counts))
    print(class_freq)

    #_ , X_test, _ , Y_test = train_test_split(X_t, Y_t, test_size=0.99, random_state=42, stratify=Y_t)
    _ , X_test, _ , Y_test = train_test_split(X_t, Y_t, test_size=0.99, random_state=42)

    clf = joblib.load('SVM.pkl')
    print(clf.score(X_test,Y_test))

def predictSVM(predict_data):

    clf = joblib.load('SVM.pkl')
    y_c = predict_data.shape[1]
    X_predict = predict_data[:,0:y_c - 3]
    Y_predict = clf.predict(X_predict)
    predict_data[:,y_c-3:y_c-2] = Y_predict.reshape((Y_predict.shape[0],1))

    return predict_data

def trainRFC(train_data,var):

    x_coord_train = train_data.shape[0]
    y_coord_train = train_data.shape[1]
    train_data = data_preprocess(train_data[0:var,:])
    X_ = train_data[:,0:y_coord_train - 3]
    Y_ = train_data[:,y_coord_train-3:y_coord_train-2]
    unique, counts = np.unique(Y_, return_counts=True)
    class_freq = dict(zip(unique, counts))
    print("Class frequencies:")
    print(class_freq)

    #X_train, X_test, Y_train, Y_test = train_test_split(X_, Y_, test_size=0.5, random_state=42, stratify=Y_)
    X_train, X_test, Y_train, Y_test = train_test_split(X_, Y_, test_size=0.5, random_state=42)
    param_grid = {"n_estimators"      : [100,200,300],
           "criterion"         : ["gini"],
           "max_features"      : [3, 4, 5],
           "max_depth"         : [10, 12, 14],
           "min_samples_split" : [4, 5, 6]}
    rfc = RandomForestClassifier(oob_score = True)
    clf = GridSearchCV(rfc,param_grid,iid=False)
    clf.fit(X_train,Y_train.ravel())

    joblib.dump(clf,'RFC.pkl')
    print(clf.score(X_test,Y_test))

def testRFC(test_data,var):

    x_coord_test = test_data.shape[0]
    y_coord_test = test_data.shape[1]
    test_data = data_preprocess(test_data[0:var,:])
    X_t = test_data[:,0:y_coord_test - 3]
    Y_t = test_data[:,y_coord_test-3:y_coord_test-2]
    unique, counts = np.unique(Y_t, return_counts=True)
    class_freq = dict(zip(unique, counts))
    print(class_freq)

    #_ , X_test, _ , Y_test = train_test_split(X_t, Y_t, test_size=0.99, random_state=42, stratify=Y_t)
    _ , X_test, _ , Y_test = train_test_split(X_t, Y_t, test_size=0.99, random_state=42)

    clf = joblib.load('RFC.pkl')
    print(clf.score(X_test,Y_test))

def predictRFC(predict_data):

    clf = joblib.load('RFC.pkl')
    y_c = predict_data.shape[1]
    X_predict = predict_data[:,0:y_c - 3]
    Y_predict = clf.predict(X_predict)
    predict_data[:,y_c-3:y_c-2] = Y_predict.reshape((Y_predict.shape[0],1))

    return predict_data

def trainGBC(train_data,var):

    x_coord_train = train_data.shape[0]
    y_coord_train = train_data.shape[1]
    train_data = data_preprocess(train_data[0:var,:])
    X_ = train_data[:,0:y_coord_train - 3]
    Y_ = train_data[:,y_coord_train-3:y_coord_train-2]
    unique, counts = np.unique(Y_, return_counts=True)
    class_freq = dict(zip(unique, counts))
    print("Class frequencies:")
    print(class_freq)

    #X_train, X_test, Y_train, Y_test = train_test_split(X_, Y_, test_size=0.5, random_state=42, stratify=Y_)
    X_train, X_test, Y_train, Y_test = train_test_split(X_, Y_, test_size=0.5, random_state=42)
    gb_grid_params = {"n_estimators": [200,300,400],
             'learning_rate': [0.05, 0.02, 0.01],
             'max_depth': [4, 6, 8],
             'min_samples_leaf': [20, 50,100]}
    gb_gs = GradientBoostingClassifier()
    clf = GridSearchCV(gb_gs,gb_grid_params,iid=False);
    clf.fit(X_train,Y_train.ravel())

    joblib.dump(clf,'GBC.pkl')
    print(clf.score(X_test,Y_test))

def testGBC(test_data,var):

    x_coord_test = test_data.shape[0]
    y_coord_test = test_data.shape[1]
    test_data = data_preprocess(test_data[0:var,:])
    X_t = test_data[:,0:y_coord_test - 3]
    Y_t = test_data[:,y_coord_test-3:y_coord_test-2]
    unique, counts = np.unique(Y_t, return_counts=True)
    class_freq = dict(zip(unique, counts))
    print(class_freq)

    #_ , X_test, _ , Y_test = train_test_split(X_t, Y_t, test_size=0.99, random_state=42, stratify=Y_t)
    _ , X_test, _ , Y_test = train_test_split(X_t, Y_t, test_size=0.99, random_state=42)

    clf = joblib.load('RFC.pkl')
    print(clf.score(X_test,Y_test))

def predictGBC(predict_data):

    clf = joblib.load('RFC.pkl')
    y_c = predict_data.shape[1]
    X_predict = predict_data[:,0:y_c - 3]
    Y_predict = clf.predict(X_predict)
    predict_data[:,y_c-3:y_c-2] = Y_predict.reshape((Y_predict.shape[0],1))

    return predict_data

def main(arguments):

  train_file = open("train.csv")
  train_data = np.genfromtxt(train_file, dtype=int, delimiter=',')
  test_file = open("test.csv")
  test_data = np.genfromtxt(test_file, dtype=int, delimiter=',')

  predict_file = open("predict.csv")
  predict_data = np.genfromtxt(predict_file, dtype=int, delimiter=',')
  predict_data = trainSVM(train_data,test_data,predict_data)

  with open("predict_SVM.csv",'w') as f:
        np.savetxt(f, predict_data.astype(int), fmt='%i', delimiter=',')

  # predict_file = open("predict.csv")
  # predict_data = np.genfromtxt(predict_file, dtype=int, delimiter=',')
  # predict_data = trainRFCV(train_data,test_data,predict_data)

  # with open("predict_RFC.csv",'w') as f:
  #       np.savetxt(f, predict_data.astype(int), fmt='%i', delimiter=',')

  # predict_file = open("predict.csv")
  # predict_data = np.genfromtxt(predict_file, dtype=int, delimiter=',')
  # predict_data = trainGBCCV(train_data,test_data,predict_data)

  # with open("predict_GBC.csv",'w') as f:
  #       np.savetxt(f, predict_data.astype(int), fmt='%i', delimiter=',')


if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))