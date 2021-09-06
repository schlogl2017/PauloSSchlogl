import ftputil
import string
import os
import sys

#Where to put the genomes (script will create sub directories
#using the same names as the NCBI use).  This directory must
#exist already:
base_path="C:\\genomes\\Bacteria\\"

#People on linux try something like:
#base_path="~/genomes/Bacteria/"

host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')
host.chdir('/genomes/Bacteria/')

dir_list = host.listdir(host.curdir)
for dir_name in dir_list :
    host.chdir('/genomes/Bacteria/')
    if host.path.isdir(dir_name):
        print dir_name
        host.chdir('/genomes/Bacteria/' + dir_name + '/')
        file_list = host.listdir(host.curdir)
        for file_name in file_list :
            if file_name[-4:]==".gbk" :
                print "File " + file_name
                if not os.path.isdir(os.path.join(base_path,dir_name)) :
                    print "Making directory " + os.path.join(base_path,dir_name)
                    os.chdir(base_path)
                    os.mkdir(os.path.join(base_path,dir_name))
                if os.path.isfile(os.path.join(base_path,dir_name,file_name)) :
                    print "Skiping file " \
                          + os.path.join(base_path,dir_name,file_name)
                elif host.path.isfile(file_name) :
                    print "Downloading file " \
                          + os.path.join(base_path,dir_name,file_name)
                    host.download(file_name, \
                          os.path.join(base_path,dir_name,file_name), 't')
                    #Download arguments: remote filename, local filename, mode
                else :
                    print "ERROR - Not a file " + dir_name + "/" + file_name
