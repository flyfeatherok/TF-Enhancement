# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 15:37:34 2017
#TF image Test
@author: xuan.xia@siat.ac.cn
"""

# -*- coding: utf-8 -*-

#import modules
import tensorflow as tf
import numpy as np
import os
import cv2
import sys
 
SignalType = sys.argv[1]
DataPath = sys.argv[2]
ModlePath = sys.argv[3]
if SignalType=='1':
    Dat='/FH_Test.dat'
    Img='/FH_Test.bmp'
if SignalType=='2':
    Dat='/LFM_Test.dat'
    Img='/LFM_Test.bmp'
if SignalType=='3':
    Dat='/SFM_Test.dat'
    Img='/SFM_Test.bmp'
TFDir = DataPath + Dat
IMGDir = DataPath + Img
Graph = ModlePath + '/model.ckpt.meta'
Modle = ModlePath + '/model.ckpt'
FFTn=256
Sizeh=np.uint32(FFTn)
Sizew=np.uint32(FFTn)
Size=Sizeh*Sizew
Gap=255*np.ones((FFTn,1),np.uint8)

TF = np.fromfile(TFDir,np.float32)
TF = np.reshape(TF, (-1, Size) )
IMG = cv2.imread(IMGDir)
IMG = cv2.cvtColor(IMG,cv2.COLOR_BGR2GRAY)
print("Data is ready")

def ImgCon(TF, IMG, Y, Gap):
    Max=max(max(TF))
    I=TF*255/Max
    I[I<0]=0
    I[I>255]=255
    I=np.reshape( I, (Sizeh,Sizew) )
    a=np.hstack((np.uint8(I),Gap))
    a=np.hstack((a,IMG))
    a = np.hstack((a,Gap))
    Y[Y<0]=0
    Y[Y>255]=255
    Y = np.reshape( np.uint8(Y), (Sizeh,Sizew) )
    Y = np.hstack((a,Y))
    return Y

with tf.Session() as sess:  
  
    saver = tf.train.import_meta_graph(Graph)
    saver.restore(sess, Modle)
    y_image = tf.get_collection('y_image')[0]
    graph = tf.get_default_graph()
    input_x = graph.get_operation_by_name('input_x').outputs[0]
    train_flag = graph.get_operation_by_name('train_flag').outputs[0]
  
    for i in range(TF.shape[0]): 
        Ix = np.reshape(TF[i], (1,Size))
        Ic = IMG[0:Sizeh,i:i+Sizew]
        y = sess.run( y_image, feed_dict={input_x: np.reshape(Ix, (1,Sizeh,Sizew,1)),train_flag:False} )  
        Y = ImgCon(Ix, Ic, y, Gap)
        cv2.imshow("WVD+GroundTurth+Estimation",Y)
        cv2.waitKey(10) 
    cv2.waitKey(0)
    
print("DONE!")