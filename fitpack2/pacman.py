#!/usr/bin/env python
import sys,os
import subprocess
import argparse

def checkdir(path):
    if not os.path.exists(path): 
        os.makedirs(path)

def setup():
    global repos
    repos={}
    repos['core']=[  'analysis'
                    ,'database'
                    ,'examples'
                    ,'qcdlib'
                    ,'nuclib'
                    ,'fitlib*'
                    ,'tools*'
                  ]
    
    repos['obslib']=[]    
    repos['grids'] =[]    
    repos['acc']   =None
    L=open('repos.txt').readlines()
    L=[_.strip() for _ in L if _.startswith('#')==False if _.strip()!='']
    for entry in L:
        key,val=entry.split(':')
        if   key=='obs': repos['obslib'].append(val)
        elif key=='grd': repos['grids'].append(val)
        elif key=='acc': repos['acc']=val
        else:
            sys.exitt('ERR: entry at repos.txt file not interpretable')

def print_msg(msg):
    msg='\n%s:'%msg
    print(msg)
    print('='*len(msg[:-1]))

def install_core_libs():
    msg='install core'
    print_msg(msg)
    account='QCDHUB'
    if repos['acc']!=None:
        account=repos['acc']

    for rep in repos['core']:
        if os.path.exists(rep.replace('*','') ):
            print('%10s already installed'%rep) 
            continue
        if '*' in rep:
            r=rep.replace('*','')
            link='git@github.com:%s/%s.git'%(account,r)
            try:
                subprocess.call(['git','clone',link])
            except:
                msg='ERR: %s not found at github account %s'
                print(msg%(r,account))
                continue
        else:
            link='git@github.com:QCDHUB/%s.git'%(rep)
            subprocess.call(['git','clone',link])

def install_obslib():
    msg='install obslib'
    print_msg(msg)
    checkdir('obslib')
    subprocess.call(['touch','obslib/__init__.py'])
    for rep in repos['obslib']:
        if os.path.exists('obslib/%s'%rep):
            print('%10s already installed'%rep) 
            continue
        link='git@github.com:QCDHUB/%s.git'%(rep)
        subprocess.call(['git','clone',link])
        subprocess.call(['mv',rep,'obslib'])

def install_grids():
    msg='install grids'
    print_msg(msg)
    checkdir('grids')
    subprocess.call(['touch','grids/__init__.py'])
    for rep in repos['grids']:
        if os.path.exists('grids/%s'%rep):
            print('%10s already installed'%rep) 
            continue
        link='git@github.com:QCDHUB/%s.git'%(rep)
        subprocess.call(['git','clone',link])
        subprocess.call(['mv',rep,'grids'])

def update_core_libs():
    msg='update core'
    print_msg(msg)
    cwd=os.getcwd()
    for rep in repos['core']:
        r=rep.replace('*','')
        if not os.path.exists(r): continue 
        print('updating %s'%r)
        os.chdir(r)
        try:
            subprocess.call(['git','pull'])
        except:
            print('ERR: git pull %s failed'%r)
        os.chdir(cwd)

def update_obsib():
    msg='update obslib'
    print_msg(msg)
    cwd=os.getcwd()
    for rep in repos['obslib']:
        print('updating %s'%rep)
        os.chdir('obslib/%s'%rep)
        try:
            subprocess.call(['git','pull'])
        except:
            print('ERR: git pull %s failed'%rep)
        os.chdir(cwd)

def update_grids():
    msg='update grids'
    print_msg(msg)
    cwd=os.getcwd()
    for rep in repos['grids']:
        print('updating %s'%rep)
        os.chdir('grids/%s'%rep)
        try:
            subprocess.call(['git','pull'])
        except:
            print('ERR: git pull %s failed'%rep)
        os.chdir(cwd)

if __name__== "__main__":

    setup()

    ap = argparse.ArgumentParser()
    ap.add_argument('option', help='install, update',type=str)
    args = ap.parse_args()

    if args.option=='install':

        install_core_libs()
        install_obslib()
        install_grids()

    if args.option=='update':

        update_core_libs()
        update_obsib()
        update_grids()
