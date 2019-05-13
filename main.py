import sys

from get_data import get_data_homology,get_data_genome

arg=sys.argv
arg=arg[1:]

if len(arg)!=5:
    print("No. of arguments more or less. Please check")
    sys.exit(1)

a=[]
d={}
a_h=[]
d_h={}
ld=[]
ldg=[]
dir_g="data"

ld,ldg,a,d=get_data_genome(arg,dir_g,a,d,ld,ldg)
#print(a[0][0:10],"\n",a[2][0:10],"\n",d,"\n",ld[0][0:10],"\n",ld[2][0:10])
dir_hom="data_homology"

a_h,d_h=get_data_homology(arg,dir_hom,a_h,d_h)
#print(a_h[0][0:10],"\n",d_h)

x=input()
