#!/usr/bin/python
import copy
import string
import subprocess
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from jdftx_constants import *

try:
    defaulttext = subprocess.check_output(['jdftx', '-t'], universal_newlines=True)
    defaulttext = '\n'.join([line for line in defaulttext.splitlines() if line!='' and line[0]!='#'])
except:
    try:
        defaulttext = subprocess.check_output(['cat', '/home/main/jdftx/help.txt'], universal_newlines=True)
    except:
        defaulttext = ''

def removeComments(text):
    i=text.find('#')
    while i!=-1:
        j=text.find('\n',i)
        if j==-1:
            text = text[:i]
        else:
            text = text[:i]+text[j:]
        i=text.find('#')
    return text

def findCommand(text,cmd,includeDefaultText=False):
    """
    Searches the (text) for (cmd) and returns argument of command i.e. the line without (cmd).
    Takes care of \ character.
    """
    if includeDefaultText:
        text=text+'\n'+defaulttext
    text = removeComments(text)

    i = text.find(cmd)

    while i!=-1 and (text[i+len(cmd)] not in string.whitespace or (i!=0 and text[i-1]!='\n')):
        i = text.find(cmd,i+len(cmd))
    if i!=-1:
        text = text.replace('\\\n','  ')
        j = text.find('\n',i)  #find the next newline after cmd
        return (text[i+len(cmd):j]).lstrip()
    else:
        raise IOError('No such command:%s'%cmd)

def removeCommand(text,cmd):
    """
    Removes the first occurrence of the cmd. Removes the whole line.
    Returns a tuple of (modified text,position of the removed command)
    """
    i = text.find(cmd)
    while i!=-1 and (text[i+len(cmd)] not in string.whitespace \
                or (i!=0 and text[i-1] not in string.whitespace)):
        i = text.find(cmd,i+len(cmd))
    if i == -1:
        raise IOError('No such command:%s'%cmd)
    j = text.find('\n',i)
    while text[j-1] == '\\':
        j=text.find('\n',j)
    #~ if j==-1:
        #~ return(text[0:i],i)
    return (text[0:i]+text[j+1:],i)

def removeAllCommands(text,cmd):
    """
    Removes all of the occurences of the command (cmd). Returns remaining text.
    """
    while True:
        try:
            text,i = removeCommand(text,cmd)
        except IOError:
            break
    return text

def isCartesian(text,is_input=True):
    """
    Checks the input string 'text' if the coord-type is set to be cartesian or not.
    If it is not an input file, set is_input to be False so that it doesn't include the default inputs.
    Return 	True if cartesian
            False if lattice
    """
    coordstype = findCommand(text,'coords-type',includeDefaultText=is_input)
    if coordstype =='lattice':
        return False
    elif coordstype == 'cartesian':
        return True
    else:
        raise IOError('Couldn\'t decide if it is cartesian or not: %s'%coordstype)

def includeFiles(text):
    while True:
        try:
            fname = findCommand(text,'include')
            f=open(fname,'r')
            includedtext = f.read()
            f.close()
            (text,i)=removeCommand(text,'include')
            text = text[:i] + includedtext + '\n' + text[i:]
        except IOError as err:
            if err.errno==None:
                break
            elif err.errno==2:
                print('The file to include is not there!')
                raise err
    return text

def readR(text):
    """Reads first occurance of
    R =
    [R00 R01 R02]
    [R10 R11 R12]
    [R20 R21 R22]
    block in an output text. Returns numpy array in shape of 3x3.
    """
    if text.find('R = \n') == -1:
        raise IOError('R = [] block was not found.')
    R = '\n'.join(text.split('R = \n',1)[1].splitlines()[:3])
    R = R.replace('[','').replace(']','')
    return np.array(R.split(),float).reshape((3,3))

def readLatticeVectors(text):
    R = np.array(findCommand(text,'lattice').split(),dtype='float64')
    R = R.reshape((3,3))
    while True:
        try:
            latticeReScale = np.array(findCommand(text,'latt-scale').split(),dtype='float64')
            latticeReScale = np.diag(latticeReScale)
            R = np.dot(R,latticeReScale)
            text = removeCommand(text,'latt-scale')[0]
        except IOError as e:
            if format(e) == 'No such command:latt-scale':
                break
            else:
                raise
    return R
